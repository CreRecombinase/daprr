#include <RcppGSL.h>
#include "classdef.hpp"
#include "dap_classdef.hpp"
#include "logistic.hpp"
#include <RcppEigen.h>


#include <algorithm>
#include <cmath>

using namespace Eigen;


using	loc_it = Rcpp::IntegerVector::iterator;
using namespace Rcpp;
bool verbose=true;

template<typename T>
std::string_view sexp_f(const T tr){

  if  (tr == 0){ return "NILSXP";}    /* nil = NULL */
  if  (tr == 1){ return "SYMSXP";}    /* symbols */
  if  (tr == 2){ return "LISTSXP";}    /* lists of dotted pairs */
  if  (tr == 3){ return "CLOSXP";}    /* closures */
  if  (tr == 4){ return "ENVSXP";}    /* environments */
  if  (tr == 5){ return "PROMSXP";}    /* promises: [un]evaluated closure arguments */
  if  (tr == 6){ return "LANGSXP";}    /* language constructs (special lists) */
  if  (tr == 7){ return "SPECIALSXP";}    /* special forms */
  if  (tr == 8){ return "BUILTINSXP";}    /* builtin non-special forms */
  if  (tr == 9){ return "CHARSXP";}    /* "scalar" string type (internal only)*/
  if  (tr == 10){ return "LGLSXP";}   /* logical vectors */
  if  (tr == 13){ return "INTSXP";}   /* integer vectors */
  if  (tr == 14){ return "REALSXP";}   /* real variables */
  if  (tr == 15){ return "CPLXSXP";}   /* complex variables */
  if  (tr == 16){ return "STRSXP";}   /* string vectors */
  if  (tr == 17){ return "DOTSXP";}   /* dot-dot-dot object */
  if  (tr == 18){ return "ANYSXP";}   /* make "any" args work */
  if  (tr == 19){ return "VECSXP";}   /* generic vectors */
  if  (tr == 20){ return "EXPRSXP";}   /* expressions vectors */
  if  (tr == 21){ return "BCODESXP";}   /* byte code */
  if  (tr == 22){ return "EXTPTRSXP";}   /* external pointer */
  if  (tr == 23){ return "WEAKREFSXP";}   /* weak reference */
  if  (tr == 24){ return "RAWSXP";}   /* raw bytes */
  if  (tr == 25){ return "S4SXP";}
  if  (tr == 30){ return "NEWSXP";}   /* fresh node creaed in new page */
  if  (tr == 31){ return "FREESXP";}   /* node released by GC */
  if  (tr == 99){ return "FUNSXP";}
  return "Don't Know";
};



class factor_col{
public:
  Rcpp::IntegerVector vec;
  Rcpp::StringVector names;
  factor_col(SEXP tr,const bool need_names=false){
    auto t = TYPEOF(tr);
    int k=0;
    if(t ==INTSXP){
      vec = tr;
      Rcpp::Nullable<Rcpp::StringVector> n_names =  vec.attr("levels");
      if(need_names){
	if(n_names.isNull()){
	  Rcpp::stop("Features must be named");
	}else{
	  names=n_names;
	  k=names.size();
	}
      }
    } else{
      if(t!=STRSXP){
	Rcpp::Rcerr<<"t is of type"<<sexp_f(t)<<std::endl;
	Rcpp::stop("column  must be int(factor) or character");
      }
      std::unordered_map<std::string,int> u_names;
      Rcpp::StringVector feat_v = tr;
      vec = Rcpp::IntegerVector(feat_v.size());
      int i=0;
      std::string buffer;
      for(auto fv: feat_v ){
	buffer = fv;
	auto u_name_i = u_names.find(buffer);
	if(u_name_i==u_names.end()){
	  auto mp = u_names.insert({buffer,++k});
	  u_name_i = mp.first;
	}
	vec(i++)=u_name_i->second;
      }
      u_names.size();
      std::vector<std::string> tnames(k);
      for(auto [mf,ti] :u_names ){
	tnames[ti-1]=mf;
      }
      names = Rcpp::wrap(tnames);
    }
  }
};

template<typename T>
class ProxyTrip {
  
  Rcpp::IntegerVector::iterator rvec;
  Rcpp::IntegerVector::iterator cvec;
  using trip_t =   Eigen::Triplet<T,typename Eigen::SparseMatrix<T>::StorageIndex>;
  trip_t ret;

public:
  ProxyTrip(Rcpp::IntegerVector::iterator rvec_, Rcpp::IntegerVector::iterator cvec_):rvec(rvec_),cvec(cvec_){
  }
  trip_t* operator->(){
    ret=trip_t((*rvec)-1,(*cvec)-1,1);

    return(&ret);
  }
  ProxyTrip& operator++(){     // prefix
    rvec++;
    cvec++;
    return(*this);
  }

  ProxyTrip& operator--(){  // prefix
    rvec--;
    cvec--;
    return *this;

  }

  ProxyTrip operator++(int){ // postfix

    ProxyTrip temp(rvec,cvec);
    // Use prefix operator to increment this digit
    ++(*this); // apply operator

    // return temporary result
    return temp; // return saved state
  }

  bool operator ==(const ProxyTrip<T> other)const {
    return( (rvec==other.rvec)&&(cvec==other.cvec));
  }
  bool operator !=(const ProxyTrip<T> other)const {
    return !((rvec==other.rvec)&&(cvec==other.cvec)) ;
  }

  ProxyTrip operator--(int){ // postfix

    ProxyTrip temp(rvec,cvec);
    // Use prefix operator to increment this digit
    --(*this); // apply operator

    // return temporary result
    return temp; // return saved state
  }
 
};


class SparseDF{
  factor_col anno_row_id;
  factor_col anno_col_id;
  const size_t p;
  const int k;
public:
  SparseDF(factor_col &&row_id, factor_col &&col_id,const size_t p_):anno_row_id(row_id),
						     anno_col_id(col_id),
								     p(p_),
								     k(anno_col_id.names.size()){
  }
  SparseDF(Rcpp::DataFrame anno_df,const size_t p_,const std::string row_name="SNP",const std::string col_name = "feature"):
  anno_row_id(anno_df[row_name]),
  anno_col_id(anno_df[col_name],true),
  p(p_),
  k(anno_col_id.names.size()){

  }

  RcppGSL::matrix<int> getMat(){
    RcppGSL::matrix<int> anno_mat(p,std::max(k,1));

    //std::fill(anno_mat.begin(),anno_mat.end(),0);
        std::memset(anno_mat->data,0,anno_mat->size1*anno_mat->size2*sizeof(int));
    const size_t nr = anno_row_id.vec.size();
    for(int i=0; i<nr; i++){
      anno_mat(anno_row_id.vec(i)-1,anno_col_id.vec(i)-1)=1;
    }
    return(anno_mat);
  }
  Rcpp::StringVector names(){
    return anno_col_id.names;
  }
};




//[[Rcpp::export]]
Rcpp::List make_matrix(const size_t p,Rcpp::DataFrame anno_df){


  SparseDF spdf( anno_df, p,"SNP","feature");

  using namespace Rcpp;
  return(List::create(_["annomat"]=spdf.getMat(),_["names"]=spdf.names()));
}







// [[Rcpp::export]]
Rcpp::List run_torus(Rcpp::IntegerVector locus_id, Rcpp::NumericVector z_hat,RcppGSL::matrix<int> anno_mat,Rcpp::StringVector names,const bool prior=false,const bool do_verbose=false,bool use_glmnet=true){


    if(!do_verbose){
        verbose=false;
    }else{
        verbose=true;
    }
    double EM_thresh = 0.05;
    double init_pi1 = 1e-3;
    int print_avg = 0;
    auto split = donut::make_splitter(locus_id.begin(),locus_id.end());
    donut::Result_obj res(locus_id.size(),anno_mat.ncol()+1);
    donut::controller con(split,res,anno_mat,EM_thresh,init_pi1,print_avg);
    con.load_data_R(z_hat);
    con.load_annotations_R(anno_mat,Rcpp::as<std::vector<std::string>>(names));
    try{
      auto result = con.estimate(use_glmnet);
      using namespace Rcpp;
      if(!prior){
        return Rcpp::List::create(_["est"]=result);
      }else{

        return (Rcpp::List::create(_["prior"]=res.prior,_["est"]=result));           
      }

    }catch(std::exception e){
      Rcpp::Rcerr<<e.what()<<std::endl;
      Rcpp::stop("Caught exception from torus!");
    }catch(int e){
      Rcpp::Rcerr<<e<<std::endl;
      Rcpp::stop("Caught exception from GSL!");

    }
}


// [[Rcpp::export]]
Rcpp::List torus_df(Rcpp::IntegerVector locus_id, Rcpp::NumericVector z_hat,Rcpp::DataFrame anno_df,const bool prior=false,const bool do_verbose=false, bool use_glmnet=true){

    const size_t p = locus_id.size();
    if(!do_verbose){
        verbose=false;
    }else{
        verbose=true;
    }
    // double EM_thresh = 0.05;
    // double init_pi1 = 1e-3;
    // int print_avg = 0;
    SparseDF spdf( anno_df, p,"SNP","feature");

    auto anno_mat = spdf.getMat();
    auto names = spdf.names();
    return run_torus(locus_id,z_hat,anno_mat,names,prior,do_verbose,use_glmnet);
}


RcppGSL::vector<double>	logit_donut(RcppGSL::matrix<int> X,RcppGSL::vector<double> y,double lambdaL1=0,double lambdaL2=0){


  const size_t npar = X.ncol()+1;

  RcppGSL::vector<double> beta(npar);

  auto kd = X.ncol();
  RcppGSL::vector<int>nlev = gsl_vector_int_calloc(kd);
  for(int i=0; i<kd; i++){
    gsl_vector_int_set(nlev, i,2);
  }

  Logistic logit(npar,y.size());
  logit.fit(beta,X,nlev,y,lambdaL1,lambdaL2);
  return(beta);
}




  //[[Rcpp::export]]
std::string dap_torus(std::vector<std::string> argv){

  int argc=argv.size();

  // creating the grid

  //olist.push_back(0.1);
  //phlist.push_back(0.05);
  std::stringstream buff_o;

  char data_file[256];
  char gmap_file[256];
  char smap_file[256];
  char annot_file[256];
  char prior_dir[256];
  char output_pip[256];
  char lik_file[256];
  int csize = -1;
  int gsize = -1;
  int nthread = 1;
  int print_avg = 0;

  int fastqtl_use_dtss = 1;



  memset(data_file,0,256);


  char init_file[256];


  int find_egene = 0;
  int est = 0;
  double alpha = 0.05;

  double init_pi1 = 1e-3;

  char ci_file[256];
  memset(lik_file,0,256);
  memset(ci_file,0,256);
  memset(data_file,0,256);
  memset(gmap_file,0,256);
  memset(smap_file,0,256);
  memset(annot_file,0,256);
  memset(init_file,0,256);
  memset(prior_dir,0,256);
  memset(output_pip,0,256);

  int force_logistic = 0;
  int prob_annot = 0;

  int data_format = 1;

  double EM_thresh = 0.05;
  double dist_bin_size = -1;
  double l1_lambda = 0;
  double l2_lambda = 0;

  for(int i=0;i<argc;i++){

    if(strcmp(argv[i].c_str(), "-d")==0 || strcmp(argv[i].c_str(), "-data")==0){
      strcpy(data_file,argv[++i].c_str());
      continue;
    }



    if(strcmp(argv[i].c_str(), "-gmap")==0){
      strcpy(gmap_file,argv[++i].c_str());
      continue;
    }


    if(strcmp(argv[i].c_str(), "-smap")==0 ){
      strcpy(smap_file,argv[++i].c_str());
      continue;
    }

    if(strcmp(argv[i].c_str(), "-lik")==0){
      strcpy(lik_file,argv[++i].c_str());
      continue;
    }

    if(strcmp(argv[i].c_str(), "-annot")==0 ){
      strcpy(annot_file,argv[++i].c_str());
      continue;
    }

    if(strcmp(argv[i].c_str(), "-t")==0 || strcmp(argv[i].c_str(), "-thresh")==0){
      EM_thresh = atof(argv[++i].c_str());
      continue;
    }

    if(strcmp(argv[i].c_str(), "-l1_lambda")==0){
      l1_lambda = atof(argv[++i].c_str());
      continue;
    }

    if(strcmp(argv[i].c_str(), "-l2_lambda")==0){
      l2_lambda = atof(argv[++i].c_str());
      continue;
    }


    if(strcmp(argv[i].c_str(), "--single_fuzzy_annot")==0){
      prob_annot = 1;
      continue;
    }


    if(strcmp(argv[i].c_str(), "-init_pi1")==0){
      init_pi1 = atof(argv[++i].c_str());
      continue;
    }


    if(strcmp(argv[i].c_str(), "--force_logistic")==0){
      force_logistic = 1;
      continue;
    }

    if(strcmp(argv[i].c_str(), "--load_bf")==0 || strcmp(argv[i].c_str(), "--bf")==0 ){
      data_format = 2;
      continue;
    }


    if(strcmp(argv[i].c_str(), "--load_zval")==0 || strcmp(argv[i].c_str(), "--zval")==0){
      data_format = 3;
      continue;
    }


    if(strcmp(argv[i].c_str(), "--load_fastqtl")==0 || strcmp(argv[i].c_str(), "--fastqtl")==0){
      data_format = 4;
      continue;
    }

    if(strcmp(argv[i].c_str(), "--load_matrixeqtl")==0 || strcmp(argv[i].c_str(), "--matrixeqtl")==0){
      data_format = 1;
      continue;
    }

    if(strcmp(argv[i].c_str(), "--no_dtss") == 0) {
      fastqtl_use_dtss = 0;
      continue;
    }





    if(strcmp(argv[i].c_str(), "-dist_bin_size") == 0){
      dist_bin_size = atof(argv[++i].c_str());
      continue;
    }


    if(strcmp(argv[i].c_str(), "-est")==0){
      est = 1;
      continue;
    }


    if(strcmp(argv[i].c_str(), "-egene")==0 || strcmp(argv[i].c_str(), "-qtl")==0 ){
      find_egene = 1;
      continue;
    }

    if(strcmp(argv[i].c_str(), "-dump_prior")==0){
      strcpy(prior_dir, argv[++i].c_str());
      continue;
    }

    if(strcmp(argv[i].c_str(), "-dump_pip")==0){
      strcpy(output_pip, argv[++i].c_str());
      continue;
    }


    if(strcmp(argv[i].c_str(), "-h")==0 || strcmp(argv[i].c_str(), "-help")==0 ){
      torus::show_banner();
      continue;
    }


    if(strcmp(argv[i].c_str(), "--print_avg")==0 ){
      print_avg = 1;
      continue;
    }


    if(strcmp(argv[i].c_str(), "-alpha")==0){
      alpha = atof(argv[++i].c_str());
      continue;
    }


    fprintf(stderr, "Error: undefined option %s\n", argv[i].c_str());
    torus::show_banner();
    Rcpp::stop("Error: undefined option "+argv[i]);

  }



  // checking mandatory arguments
  if(strlen(data_file)==0){
    fprintf(stderr,"Error: data file unspecified\n");

    torus::show_banner();
    Rcpp::stop("Error: data file unspecified\n");
  }


  // a global variable
  torus::controller con(buff_o);
  con.EM_thresh = EM_thresh;

  if(dist_bin_size > 0){
    con.dist_bin_size = dist_bin_size;
  }

  if(force_logistic){
    con.force_logistic = 1;
  }

  if(prob_annot){
    con.single_fuzzy_annot = 1;
  }

  if(l1_lambda!=0){
    con.l1_lambda = l1_lambda;
    con.force_logistic = 1;
  }

  if(l2_lambda!=0){
    con.l2_lambda = l2_lambda;
    con.force_logistic =1;
  }


  con.init_pi1 = init_pi1;
  con.print_avg = print_avg;

  switch(data_format){
  case 1:
    con.load_data(data_file);
    break;
  case 2:
    con.load_data_BF(data_file);
    break;
  case 3:
    con.load_data_zscore(data_file);
    break;
  case 4:
    con.fastqtl_use_dtss = fastqtl_use_dtss;
    con.load_data_fastqtl(data_file);
    gmap_file[0]=smap_file[0] = 0;
    break;
  default:
    con.load_data(data_file);
    break;
  }


  con.load_map(gmap_file, smap_file);
  con.load_annotation(annot_file);
  fprintf(stderr,"Initializing ... \n");

  if(est==0 && find_egene==0){
    est = 1;
  }


  if(est)
    con.estimate();
  if(find_egene){
    con.find_eGene(alpha);
  }
  if(strlen(prior_dir)>0){
    con.dump_prior(prior_dir);
  }
  if(strlen(output_pip)>0){
    fprintf(stderr,"#-dump_pip output_pip file: %s \n",output_pip);
    con.dump_pip(output_pip);
  }

  if(strlen(lik_file)>0){
    std::string filename = lik_file;
    std::ofstream ostrm(filename);
    ostrm << con.final_log10_lik<<std::endl;
  }

  return con.buf_o.str();
}
