#include <RcppGSL.h>
#include "classdef.hpp"
#include "logistic.hpp"
#include <RcppEigen.h>


#include <algorithm>
#include <cmath>

using namespace Eigen;


using	loc_it = Rcpp::IntegerVector::iterator;
using namespace Rcpp;
bool verbose=true;



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
	Rcpp::stop("column `feature` in anno_df  must be int(factor) or character");
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
  ProxyTrip(  Rcpp::IntegerVector::iterator rvec_, Rcpp::IntegerVector::iterator cvec_):rvec(rvec_),cvec(cvec_){
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

  SparseDF(Rcpp::DataFrame anno_df,const size_t p_,const std::string row_name="SNP",const std::string col_name = "feature"):
  anno_row_id(anno_df[row_name]),
  anno_col_id(anno_df[col_name],true),
  p(p_),
  k(anno_col_id.names.size()){

  }

  // Eigen::SparseMatrix<double> getspMat(){

  //   Eigen::SparseMatrix<double> ret(p,std::max(k,1));
  //   ProxyTrip<double> ptb(anno_row_id.vec.begin(),anno_col_id.vec.begin());
  //   ProxyTrip<double> pte(anno_row_id.vec.end(),anno_col_id.vec.end());
  //   ret.setFromTriplets(ptb,pte);
  //   return(ret);
  // }
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
  //  auto [mat,names] = gen_anno_mat(p,anno_df);

  using namespace Rcpp;
  return(List::create(_["annomat"]=spdf.getMat(),_["names"]=spdf.names()));
}



// [[Rcpp::export]]
Rcpp::List torus(Rcpp::IntegerVector locus_id, Rcpp::NumericVector z_hat,RcppGSL::matrix<int> anno_mat,Rcpp::StringVector names,const bool prior=false,const bool do_verbose=false,bool use_glmnet=true){


    if(!do_verbose){
        verbose=false;
    }else{
        verbose=true;
    }
    double EM_thresh = 0.05;
    double init_pi1 = 1e-3;
    int print_avg = 0;
    auto split = make_splitter(locus_id.begin(),locus_id.end());
    Result_obj res(locus_id.size(),anno_mat.ncol()+1);
    controller con(split,res);
    con.EM_thresh = EM_thresh;
    con.init_pi1 = init_pi1;
    con.print_avg = print_avg;
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
    auto [anno_mat,names] = gen_anno_mat(p,anno_df);

    return torus(locus_id,z_hat,anno_mat,names,prior,do_verbose,use_glmnet);
}




//[[Rcpp::export]]
RcppGSL::vector<double>	logit_torus(RcppGSL::matrix<int> X,RcppGSL::vector<double> y,double lambdaL1=0,double lambdaL2=0){


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
RcppGSL::vector<double>	logit_qbinom(RcppGSL::matrix<int> X,RcppGSL::vector<double> y){


  const size_t npar = X.ncol()+1;

  RcppGSL::vector<double> beta(npar-1);

  auto kd = X.ncol();
  RcppGSL::vector<int>nlev = gsl_vector_int_calloc(kd);
  for(int i=0; i<kd; i++){
    gsl_vector_int_set(nlev, i,2);
  }

  Logistic logit(npar,y.size());
  logit.fit_qbinom(beta,X,y);
  return(beta);
}




// Rcpp::NumericVector logistic_cat_pred(RcppGSL::vector<double> beta,  // Vector of parameters length = 1 + Sum_k(C_k - 1)
//                          RcppGSL::matrix<int> X,  //Matrix Nobs x K
//                          RcppGSL::vector<int> nlev,
//                          RcppGSL::vector<double> yhat
//                                         )
// {

//   for(int i = 0; i < X->size1; ++i) {
//     double Xbetai=beta->data[0];
//     int iParm=1;
//     for(int k = 0; k < X->size2; ++k) {
//       if(gsl_matrix_int_get(X,i,k)>0)
        
//         Xbetai+=beta->data[gsl_matrix_int_get(X,i,k)-1+iParm];
//       iParm+=nlev->data[k]-1;
//     }
//     yhat->data[i]=1/(1 + std::exp(-Xbetai));
//   }
//   return(Rcpp::wrap(yhat));
// }
