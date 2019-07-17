#include <RcppGSL.h>
#include "classdef.hpp"
#include "logistic.hpp"
#include <RcppEigen.h>
#include <gsl/gsl_multifit.h>
#include <gsl/span>

#include <algorithm>
#include <cmath>

using namespace Eigen;


extern "C" void F77_SUB(splognet)(double *alpha,
	       int *n,
	       int *p,
	       int *nc,
	       double *x,
	       int *ix,
	       int *jx,
	       double *y,
	       double *g,
	       int *jd,
	       double *vp,
	       double *cl,
	       int *ne,
	       int *nx,
	       int *nlam,
	       double *flmin,
	       double *ulam,
	       double *thr,
	       int *isd,
	       int *intr,
	       int *maxit,
	       int *kopt,
	       int *lmu,
	       double *a0,
	       double *ca,
	       int *ia,
	       int *nin,
	       double *dev0,
	       double *dev,
	       double *alm,
	       int *nlp,
	       int *jerr);

extern "C" void F77_SUB(luncomp)(int* p,int *nx,int* nc,double* ca,int* ia,int *nin,double* a);

extern "C" void F77_SUB(lsolns)(int* p,int* nx,int* nc,int* lmu,double* ca,int* ia,int* nin,double* b);

extern "C" void F77_SUB(lognet) (double *parm,
		     int *no,
		     int *p,
		     int *nc,
		     double *x,
		     double *y,
		     double *g,
		     int *jd,
		     double *vp,
		     double *cl,
		     int *ne,
		     int *nx,
		     int *nlam,
		     double *flmin,
		     double *ulam,
		     double *thr,
		     int *isd,
		     int *intr,
		     int *maxit,
		     int *kopt,
		     int *lmu,
		     double *a0,
		     double *ca,
		     int *ia,
		     int *nin,
		     double *dev0,
		     double *dev,
		     double *alm,
		     int *nlp,
		     int *jerr
		     );


class comp_coeff{
  int p;
  int nlam;


public:
  comp_coeff(const int nlam_, const int p_,const bool zero=true):
  p(p_),
  nlam(nlam_),
  a0(nlam),
  ca(p,nlam),
  nin(nlam),
  ia(p),
  alm(nlam),
  lmu(0){

  if(zero){
    a0.setZero();
    ca.setZero();
    nin.setZero();
    alm.setZero();
  }
  }

  void setZero(){
    a0.setZero();
    ca.setZero();
    nin.setZero();
    alm.setZero();
  }
  void unpack_coeffs(){
    if(lmu==0){
      Rcpp::stop("number of lambda values is 0!");
    }

    int	one=1;
    retcoeff.resize(nlam,lmu);
    F77_SUB(lsolns)(&p,&nlam,&one,&lmu,ca.data(),ia.data(),nin.data(),retcoeff.data());
  }

  Eigen::ArrayXd a0;
  Eigen::MatrixXd ca;
  Eigen::ArrayXi nin;
  Eigen::ArrayXi ia;
  Eigen::ArrayXd alm;
  Eigen::MatrixXd retcoeff;
  int lmu;
};

// c call lognet (alpha,n,p,one,x,y,o,use_f,vp,cl,ne,nx,nlam,flmin,ulam,thr,isd,
// c              intr,maxit,kopt,lmu,a0,ca,ia,nin,dev0,fdev,alm,nlp,jerr)
// c call lognet (parm,no,ni,nc,x,y,o,jd,vp,cl,ne,nx,nlam,flmin,ulam,thr,isd,
// c              intr,maxit,kopt,lmu,a0,ca,ia,nin,dev0,fdev,alm,nlp,jerr)

template<typename T>
class Lognet{
  const int n;
  const int p;
  std::vector<double> ulam;
  int nlam;
  comp_coeff coeff;
  double alpha ;
  const T& x_o;
  T x;
  Eigen::MatrixXd y;
  Eigen::MatrixXd o;
  Eigen::ArrayXd one_vec;
  Eigen::MatrixXd interval_mat;
  int ne;
  int nx;
  double flmin;
  double thr;
  int isd;
  int intr;
  int maxit;
  int kopt;

  double dev0;
  Eigen::ArrayXd fdev;
  Eigen::ArrayXd alm;
  Eigen::ArrayXd coeffs;
  int nlp;
  int jerr;
public:
  Lognet(const T &X_, const double alpha_=0,std::vector<double> lambda= {0},const double thresh=1e-07,const int maxiter=100000):
  n(X_.rows()),
  p(X_.cols()),
  ulam(std::move(lambda)),
  nlam(ulam.size()),
  coeff(nlam,p),
  alpha(alpha_),
  x_o(X_),
  x(X_),
  y(n,2),
  o(n,1),
  one_vec(Eigen::ArrayXd::Ones(p)),
  interval_mat(2,p),
  ne(p+1),
  nx(p),
  flmin(1.0),
  thr(thresh),
  isd(1),
  intr(1),
  maxit(maxiter),
  kopt(0),
  dev0(0),
  fdev(nlam),
  alm(nlam),
  nlp(0),
  jerr(0){

  for(int i=0; i<p; i++){
    interval_mat(0,i)=-9e35;
    interval_mat(1,i)=9e35;
  }

  //  interval_mat.colwise() = (Eigen::VectorXd(2) <<,9e35).finished();

}
  void fit(const gsl::span<double> yd);
  void read_coeffs(gsl::span<double> beta,int index=0);
  void predict(const gsl::span<double> beta,gsl::span<double> y)const;
};


// template<>
// inline Lognet<Eigen::MatrixXd>::Lognet(const Eigen::MatrixXd &X_,const double alpha_, std::vector<double> lambda,const double thresh,const int maxiter)

// template<>
// inline Lognet<Eigen::SparseMatrix<double>>::Lognet(const Eigen::SparseMatrix<double> &X_,const double alpha_, std::vector<double> lambda,const double thresh,const int maxiter):
//   n(X_.rows()),
//   p(X_.cols()),
//   ulam(std::move(lambda)),
//   nlam(ulam.size()),
//   coeff(nlam,p),
//   alpha(alpha_),
//   x_o(X_),
//   x(X_),
//   y(n,2),
//   o(n,1),
//   one_vec(Eigen::ArrayXd::Ones(p)),
//   interval_mat(2,p),
//   ne(p+1),
//   nx(p),
//   flmin(1.0),
//   thr(thresh),
//   isd(1),
//   intr(1),
//   maxit(maxiter),
//   kopt(0),
//   dev0(0),
//   fdev(nlam),
//   alm(nlam),
//   nlp(0),
//   jerr(0){

//   for(int i=0; i<p; i++){
//     interval_mat(0,i)=-9e35;
//     interval_mat(1,i)=9e35;
//   }


// }

template<>
inline void Lognet<Eigen::MatrixXd>::fit(const gsl::span<double> yd){


  y.col(0)=  Eigen::Map<Eigen::ArrayXd>(yd.data(),yd.size());
  y.col(1)=1-y.col(0).array();
  x=x_o;
  o.setZero();
  int jd = 0;
  int nc = 1;
  //  lmu=0;
  lognet_(&alpha,
	  &n,
	  &p,
	  &nc,
	  x.data(),
	  y.data(),
	  o.data(),
	  &jd,
	  vp.data(),
	  cl.data(),
	  &ne,
	  &nx,
	  &nlam,
	  &flmin,
	  ulam.data(),
	  &thr,
	  &isd,
	  &intr,
	  &maxit,
	  &kopt,
	  &coeff.lmu,
	  coeff.a0.data(),
	  coeff.ca.data(),
	  coeff.ia.data(),
	  coeff.nin.data(),
	  &dev0,
	  fdev.data(),
	  alm.data(),
	  &nlp,
	  &jerr);

  coeff.unpack_coeffs();

  // Eigen::Map<Eigen::ArrayXd > coeff_v(ca.data(),lmu*ni);



}

template<>
inline void Lognet<Eigen::SparseMatrix<double>>::fit(const gsl::span<double> yd){

  y.col(0)=  Eigen::Map<Eigen::ArrayXd>(yd.data(),yd.size());
  y.col(1)=1-y.col(0).array();

  x=x_o;
  x.makeCompressed();
  o.setZero();

  splognet_(&parm,
	    &no,
	    &ni,
	    &nc,
	    x.valuePtr(),
	    x.innerIndexPtr(),
	    x.outerIndexPtr(),
	    y.data(),
	    o.data(),
	    &jd,
	    vp.data(),
	    cl.data(),
	    &ne,
	    &nx,
	    &nlam,
	    &flmin,
	    ulam.data(),
	    &thr,
	    &isd,
	    &intr,
	    &maxit,
	    &kopt,
	    &coeff.lmu,
	    coeff.a0.data(),
	    coeff.ca.data(),
	    coeff.ia.data(),
	    coeff.nin.data(),
	    &dev0,
	    fdev.data(),
	    alm.data(),
	    &nlp,
	    &jerr);

    coeff.unpack_coeffs();

}

template<typename T>
inline void Lognet<T>::read_coeffs(gsl::span<double> beta,int index){
  beta[0]=coeff.a0[index];
  Eigen::VectorXd tcol=coeff.retcoeff.col(index);
  std::copy_n(tcol.data(),tcol.size(),beta.begin()+1);
}

template<typename T>
inline void Lognet<T>::predict(const gsl::span<double> beta,gsl::span<double> y)const {

  const Eigen::Map<Eigen::VectorXd> beta_v(beta.data(),beta.size());
  Eigen::Map<Eigen::ArrayXd> y_v(y.data(),y.size());
  y_v = beta_v.tail(beta_v.size()-1)*x_o;

  y_v =	1/(1+(-(beta_v[0]+y_v).exp()));

}





double prior_to_pip(gsl::span<double> pip,const gsl::span<double> prior,const bool fix_pi0=true){

  double locus_pi0=1;
  std::transform(prior.begin(),prior.end(),pip.begin(),[&locus_pi0](double prior){
							 locus_pi0*=(1-prior);
							 return(prior/(1-prior));
						       });
  if(locus_pi0 < 1e-100 && fix_pi0){
    locus_pi0 = 1e-100;
  }

  for(auto &p_v : pip)
    p_v *= locus_pi0;
  return locus_pi0;

}


void transform_pip(gsl::span<double> pip,const gsl::span<double> BF,const double log10lik){

  std::transform(pip.begin(),pip.end(),BF.begin(),pip.begin(),
		 [log10lik](double prior,double log10_BF){
		   return pow(10,(log10(prior) + log10_BF - log10lik));
		 });
}


double log10_lik(const double locus_pi0,const gsl::span<double> BF, const gsl::span<double> pip,std::optional<double> &BF_max){
  if(!BF_max){
    *BF_max =  std::max(*(std::max_element(BF.begin(),BF.end())),0.0);
  }
  double max_el = *BF_max;
  double sum = locus_pi0*pow(10,0-max_el);

  double log10lik =	std::inner_product(BF.begin(),
					   BF.end(),
					   pip.begin(),
					   sum,std::plus<>(),
					   [max_el](double vec,double wts){
					     return(wts*pow(10,vec-max_el));});
  log10lik=max_el+log10(log10lik);
  return log10lik;
}



double E_step(gsl::span<double> pip,const gsl::span<double> BF, const gsl::span<double> prior,std::optional<double> &BF_max){

  double locus_pi0=  prior_to_pip(pip,prior);

  double lik10 = log10_lik(locus_pi0,BF,pip,BF_max);
  transform_pip(pip,BF,lik10);

  return(lik10);
}

double fdr(double log10lik,const gsl::span<double> BF, const gsl::span<double> prior,gsl::span<double> pip,std::optional<double> &BF_max){


  double locus_pi0=prior_to_pip(prior,pip,false);

  transform_pip(pip,BF,log10lik);

  log10lik = log10_lik(locus_pi0,BF,prior,BF_max);

  if(!BF_max){
    *BF_max =  std::max(*(std::max_element(BF.begin(),BF.end())),0.0);
  }
  double max_el = *BF_max;
  log10lik=max_el+log10(log10lik);

  double fdr =  std::pow(10,log10(locus_pi0)-log10lik);

  return fdr;
}

using SplitView = std::vector< gsl::span<double> >;

class GroupedView{
public:
  SplitView r_view;
  gsl::span<double> d_view;
  GroupedView(const splitter &splt_,
	      double* BF_p, const size_t p):r_view(splt_.split_view(BF_p)),
					    d_view(BF_p,p){}

};




//SumStatRegion is a non-owning view of the summary statistics
class SumStatRegion{
public:
  const splitter& splt;
  const GroupedView BF;
private:
  const size_t num_reg;
  // SplitView r_prior;
  // gsl::span<double> prior;
  // SplitView r_pip;
  // gsl::span<double> pip;
  // const SplitView r_BF;
  // const gsl::span<double> BF;
  mutable std::vector<std::optional<double>> BF_max;

public:
  SumStatRegion(const splitter &splt_,
		double*	BF_p,const size_t p):splt(splt_),
					     BF(splt,BF_p,p),
					     num_reg(splt_.num_regions())
  {

  }
  double E_steps(SplitView &r_pip,const SplitView &r_prior)const{
    double loglik=0;
    for(int i=0; i<num_reg; i++){
      loglik+=E_step(r_pip[i],BF.r_view[i],r_prior[i],BF_max[i]);
    }
    return loglik;
  }
  const size_t size() const {
    return BF.d_view.size();
  }
};












using namespace Eigen;
// template<typename T>
// class Lognet{
//   const T &X;
//   Lognet<T> lnm;
// public:
//   Lognet(const T &X_,const double alpha=0,std::vector<double> lambda= {0},const double thresh=1e-07,const int maxiter=100000):
//     X(X_),
//     lnm(X,alpha,lambda,thresh,maxiter){
//   }
//   void fit(gsl::span<double> beta, const gsl::span<double> y){
//     lnm.fit(Eigen::Map<Eigen::ArrayXd>(y.data(),y.size()));
//     lnm.read_coeffs(beta,0);
//   }


//   void predict(const gsl::span<double> beta,gsl::span<double> y) const{

//     lnm.predict(beta,y);

//   }
//   size_t num_params()const {
//     return X.cols()+1;
//   }
// };



double estimate_torus_sp(gsl::span<double> beta,GroupedView pip_v,GroupedView prior_v,const splitter& splt,const gsl::span<double> BF,const SparseMatrix<double> &X,const double EM_thresh=0.05){

  SumStatRegion sumstats(splt,BF.data(),prior_v.d_view.size());
  auto &pip_r = pip_v.r_view;
  auto &prior_r = prior_v.r_view;
  auto &prior =  prior_v.d_view;
  Lognet<SparseMatrix<double>> logistic(X);
  logistic.predict(beta,prior);
  int iter_ct=0;
  int iter_max=20;
  double last_log10_lik = -9999999;
  double curr_log10_lik = 0;
  while(iter_ct<iter_max && fabs(curr_log10_lik-last_log10_lik)>EM_thresh){
    last_log10_lik=curr_log10_lik;
    curr_log10_lik=sumstats.E_steps(pip_r,prior_r);
     logistic.fit(prior);
     logistic.read_coeffs(beta);
     //     lnm.read_coeffs(beta,0);
     //    logistic.fit(beta,prior);
    logistic.predict(beta,prior);
  }
  return curr_log10_lik;
}



template<typename T>
double evaluate_likelihood(GroupedView &pip,GroupedView &prior, const gsl::span<double> beta,SumStatRegion sumstats,const Lognet<T> &logistic){
  logistic.predict(beta,prior.d_view);
  return sumstats.E_steps(pip.r_view,prior.r_view);
}


template<typename T>
double estimate_torus_null(const gsl::span<double> beta,const size_t index ,SumStatRegion sumstats,const Lognet<T> &logistic,GroupedView &null_prior,GroupedView &null_pip){

  std::vector<double> null_beta_v(beta.size());
  std::copy(beta.begin(),beta.end(),null_beta_v.begin());
  null_beta_v[index]=0;
  gsl::span<double> null_beta(null_beta_v.data(),null_beta_v.size());
  logistic.predict(null_beta,null_prior.d_view);
  double null_log10_lik = sumstats.E_steps(null_pip.r_view,null_prior.r_view);
  return null_log10_lik;
}



template<typename T>
double fine_optimize_beta(double& curr_log10_lik,const gsl::span<double> beta,const size_t index ,SumStatRegion sumstats,const Lognet<T> &logistic,GroupedView &null_prior,GroupedView &null_pip){

  const double gr = (std::sqrt(5)-1)/2.0;
  std::vector<double> null_beta_v(beta.size());
  std::copy(beta.begin(),beta.end(),null_beta_v.begin());
  gsl::span<double> null_beta(null_beta_v.data(),null_beta_v.size());
  double a = null_beta[index];
  double b = 0;
  if(a > b){
    b = a;
    a = 0;
  }
  double c = b - gr*(b-a);
  double d = a + gr*(b-a);
  double thresh = 1e-3;
  if(fabs(c-d)<thresh){
    thresh = fabs(c-d)/2.0;
  }

  double fc;
  double fd;



  while((d-c)>thresh){

    null_beta[index]=c;
    fc = evaluate_likelihood(null_pip,null_prior,null_beta,sumstats,logistic);

    null_beta[index]=d;
    fd = evaluate_likelihood(null_pip,null_prior,null_beta,sumstats,logistic);

    if(fc > fd){

      b = d;
      d = c;
      c = b - gr*(b-a);

    }else{

      a = c;
      c = d;
      d = a + gr*(b-a);

    }
  }

  curr_log10_lik = fc;

  return (b+a)/2;

}


//
// template<typename T>
// Rcpp::List estimate_torus_glmnet(const gsl::span<int> region_id,const gsl::span<double>	z_hat, const T X,


    // update_prior(const gsl::span<double> beta){
    // for(int i = 0; i < X->size1; ++i) {
    //   double Xbetai=beta->data[0];
    //   int iParm=1;
    //   for(int k = 0; k < X->size2; ++k) {
    // 	if(gsl_matrix_int_get(X,i,k)>0)
    // 	  Xbetai+=beta->data[gsl_matrix_int_get(X,i,k)-1+iParm];
    // 	iParm+=nlev->data[k]-1;
    //   }
    //   yhat->data[i]=1/(1 + gsl_sf_exp(-Xbetai));
    // }

// template<typename Matrix>
// void fit_torus_glmnet(gsl::span<double> beta,const Matrix X,SumStatRegion &ssr){

//   Logistic logit(beta.size(),ssr.size());
//   double last_log10_lik = -9999999;

//   int itcc=0;
//   int iter_max=10;
//   while(itcc<iter_max){
//     if(verbose){
//       Rcpp::Rcerr<<"iter no:"<<itcc++<<std::endl;
//     }

//     double curr_log10_lik = 0;
//     for(auto &lv : locVec){
//       lv.EM_update();
//       curr_log10_lik += lv.log10_lik;
//     }
//     if(ncoef==1){
//       if(verbose){
// 	Rcpp::Rcerr<<"simple"<<std::endl;
//       }
//       simple_regression();
//     }
//     // only categrical annotations
//     else if(kd!=0){
//       if(kd == 1 && !force_logistic){
// 	if(verbose){
// 	  Rcpp::Rcerr<<"single_ct"<<std::endl;
// 	}
// 	single_ct_regression();
//       }else{
// 	if(verbose){
// 	  Rcpp::Rcerr<<"cat_fit"<<std::endl;
// 	}
// 	if(use_glmnet){
// 	  logit.fit_glmnet(beta_vec,Xd,pip_vec,l1_lambda,l2_lambda);
// 	}else{
// 	  logit.fit(beta_vec,Xd,dlevel,pip_vec,l1_lambda,l2_lambda);
// 	}
// 	//	logistic_cat_fit(beta_vec, Xd, dlevel,pip_vec, l1_lambda,l2_lambda);
//       }
//       if(verbose){
// 	Rcpp::Rcerr<<"cat_pred"<<std::endl;
// 	Rcpp::Rcerr<<"beta:\n";
//       }
//       	gsl::span<double> beta_sp(beta_vec->data,beta_vec->size);
// 	if(std::accumulate(beta_sp.begin(),beta_sp.end(),0)<1e-8){
// 	  Rcpp::Rcerr<<"Beta is entirely 0..."<<std::endl;
// 	}
// 	for(auto be : beta_sp){
// 	  if(verbose){

// 	    Rcpp::Rcerr<<be<<"\n";
// 	  }
// 	  if(isnan(be)){
// 	    Rcpp::stop("NaN encountered in beta");
// 	  }
// 	}
// 	if(verbose){
// 	  Rcpp::Rcerr<<std::endl;
// 	}
//       logistic_cat_pred(beta_vec, Xd, dlevel,prior_vec);
//     }


//     if(fabs(curr_log10_lik-last_log10_lik)<EM_thresh){
//       final_log10_lik = curr_log10_lik;
//       break;
//     }

//     last_log10_lik = curr_log10_lik;
//   }


//   double tp = 0;
//   for(int i=0;i<p;i++){
//     tp += gsl_vector_get(prior_vec,i);
//   }

// }





using	loc_it = Rcpp::IntegerVector::iterator;
using namespace Rcpp;
bool verbose=true;

std::pair<RcppGSL::matrix<int>,Rcpp::StringVector> gen_anno_mat(const size_t p,Rcpp::DataFrame anno_df){

    Rcpp::IntegerVector anno_row_id = anno_df["SNP"];
    SEXP tr = anno_df["feature"];

    auto t = TYPEOF(tr);
    int k=0;
    const size_t nr = anno_row_id.size();

    Rcpp::IntegerVector anno_col_id;
    Rcpp::StringVector names;
    if(t ==INTSXP){
        anno_col_id = anno_df["feature"];
        Rcpp::Nullable<Rcpp::StringVector> n_names =  anno_col_id.attr("levels");
        if(n_names.isNull()){
            Rcpp::stop("Features must be named");
        }
        names=n_names;
        k=names.size();
    } else{
        if(t!=STRSXP){
            Rcpp::stop("column `feature` in anno_df  must be int(factor) or character");
        }
        std::unordered_map<std::string,int> u_names;
        Rcpp::StringVector feat_v = anno_df["feature"];
        anno_col_id = Rcpp::IntegerVector(nr);

        int i=0;
	std::string buffer;
        for(auto fv: feat_v ){
            buffer = fv;
            auto u_name_i = u_names.find(buffer);
            if(u_name_i==u_names.end()){
                auto mp = u_names.insert({buffer,++k});
                u_name_i = mp.first;
            }
            anno_col_id(i++)=u_name_i->second;
        }
        u_names.size();
        std::vector<std::string> tnames(k);
        for(auto [mf,ti] :u_names ){
            tnames[ti-1]=mf;
        }
        names = Rcpp::wrap(tnames);
    }


    RcppGSL::matrix<int> anno_mat(p,std::max(k,1));

    std::memset(anno_mat->data,0,anno_mat->size1*anno_mat->size2*sizeof(int));

    for(int i=0; i<nr; i++){
        anno_mat(anno_row_id(i)-1,anno_col_id(i)-1)=1;
    }
    if(nr==0){
      names=anno_df.attr("feature");
    }
    return(std::make_pair(anno_mat,names));
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

  Eigen::SparseMatrix<double> getspMat(){

    Eigen::SparseMatrix<double> ret(p,std::max(k,1));
    ProxyTrip<double> ptb(anno_row_id.vec.begin(),anno_col_id.vec.begin());
    ProxyTrip<double> pte(anno_row_id.vec.end(),anno_col_id.vec.end());
    ret.setFromTriplets(ptb,pte);

//     SEXP asdgCMatrix_( SEXP XX_ ){
//       // typedef Eigen::SparseMatrix<double> SpMat;
//       // typedef Eigen::Map<Eigen::MatrixXd> MapMatd; // Input: must be double
//       // MapMatd X(Rcpp::as<MapMatd>(XX_));
//       // SpMat Xsparse = X.sparseView();
//       // Output: sparse matrix
//       S4 Xout(wrap(Xsparse));                      // Output: as S4 object
//       NumericMatrix Xin(XX_);                      // Copy dimnames
//       Xout.slot("Dimnames") = clone(List(Xin.attr("dimnames")));
      return(ret);
  }

  RcppGSL::matrix<int> getMat(){
    RcppGSL::matrix<int> anno_mat(p,std::max(k,1));

    std::memset(anno_mat->data,0,anno_mat->size1*anno_mat->size2*sizeof(int));
    const size_t nr = anno_row_id.vec.size();
    for(int i=0; i<nr; i++){
        anno_mat(anno_row_id.vec(i)-1,anno_col_id.vec(i)-1)=1;
    }
    // if(nr==0){
    //   names=anno_df.attr("feature");
    // }
    return(anno_mat);
  }
  Rcpp::StringVector names(){
    return anno_col_id.names;
  }
};







// std::pair<Eigen::SparseMatrix<int>,Rcpp::StringVector> gen_anno_spmat(const size_t p,Rcpp::DataFrame anno_df){

//     Rcpp::IntegerVector anno_row_id = anno_df["SNP"];
//     SEXP tr = anno_df["feature"];

//     auto t = TYPEOF(tr);
//     int k=0;
//     const size_t nr = anno_row_id.size();

//     Rcpp::IntegerVector anno_col_id;
//     Rcpp::StringVector names;
//     if(t ==INTSXP){
//         anno_col_id = anno_df["feature"];
//         Rcpp::Nullable<Rcpp::StringVector> n_names =  anno_col_id.attr("levels");
//         if(n_names.isNull()){
//             Rcpp::stop("Features must be named");
//         }
//         names=n_names;
//         k=names.size();
//     } else{
//         if(t!=STRSXP){
//             Rcpp::stop("column `feature` in anno_df  must be int(factor) or character");
//         }
//         std::unordered_map<std::string,int> u_names;
//         Rcpp::StringVector feat_v = anno_df["feature"];
//         anno_col_id = Rcpp::IntegerVector(nr);

//         int i=0;
// 	std::string buffer;
//         for(auto fv: feat_v ){
//             buffer = fv;
//             auto u_name_i = u_names.find(buffer);
//             if(u_name_i==u_names.end()){
//                 auto mp = u_names.insert({buffer,++k});
//                 u_name_i = mp.first;
//             }
//             anno_col_id(i++)=u_name_i->second;
//         }
//         u_names.size();
//         std::vector<std::string> tnames(k);
//         for(auto [mf,ti] :u_names ){
//             tnames[ti-1]=mf;
//         }
//         names = Rcpp::wrap(tnames);
//     }


//

//
// };






//  kopt = optimization flag
//     kopt = 0 => Newton-Raphson (recommended)
//     kpot = 1 => modified Newton-Raphson (sometimes faster)
//     kpot = 2 => nonzero coefficients same for each class (nc > 1)


//x(no,ni) = predictor data matrix flat file (overwritten)
//ka = algorithm flag
//     ka=1 => covariance updating algorithm
//     ka=2 => naive algorithm
//  parm = penalty member index (0 <= parm <= 1)
//       = 0.0 => ridge
//       = 1.0 => lasso
//  no = number of observations
//  ni = number of predictor variables
//  y(no) = response vector (overwritten)
//  w(no)= observation weights (overwritten)
//  jd(jd(1)+1) = predictor variable deletion flag
//     jd(1) = 0  => use all variables
//     jd(1) != 0 => do not use variables jd(2)...jd(jd(1)+1)
//  vp(ni) = relative penalties for each predictor variable
//     vp(j) = 0 => jth variable unpenalized
//  cl(2,ni) = interval constraints on coefficient values (overwritten)
//     cl(1,j) = lower bound for jth coefficient value (<= 0.0)
//     cl(2,j) = upper bound for jth coefficient value (>= 0.0)
//  ne = maximum number of variables allowed to enter largest model
//       (stopping criterion)
//  nx = maximum number of variables allowed to enter all models
//       along path (memory allocation, nx > ne).
//  nlam = (maximum) number of lamda values
//  flmin = user control of lamda values (>=0)
//     flmin < 1.0 => minimum lamda = flmin*(largest lamda value)
//     flmin >= 1.0 => use supplied lamda values (see below)
//  ulam(nlam) = user supplied lamda values (ignored if flmin < 1.0)
//  thr = convergence threshold for each lamda solution.
//     iterations stop when the maximum reduction in the criterion value
//     as a result of each parameter update over a single pass
//     is less than thr times the null criterion value.
//     (suggested value, thr=1.0e-5)
//  isd = predictor variable standarization flag:
//     isd = 0 => regression on original predictor variables
//     isd = 1 => regression on standardized predictor variables
//     Note: output solutions always reference original
//           variables locations and scales.
//  intr = intercept flag
//     intr = 0/1 => don't/do include intercept in model
//  maxit = maximum allowed number of passes over the data for all lambda
//     values (suggested values, maxit = 100000)


//  lmu = actual number of lamda values (solutions)
//  a0(lmu) = intercept values for each solution
//  ca(nx,lmu) = compressed coefficient values for each solution
//  ia(nx) = pointers to compressed coefficients
//  nin(lmu) = number of compressed coefficients for each solution
//  alm(lmu) = lamda values corresponding to each solution
//  nlp = actual number of passes over the data for all lamda values
//  jerr = error flag:
//     jerr  = 0 => no error
//     jerr > 0 => fatal error - no output returned
//        jerr < 7777 => memory allocation error
//        jerr = 7777 => all used predictors have zero variance
//        jerr = 10000 => maxval(vp) <= 0.0
//     jerr < 0 => non fatal error - partial output:
//        Solutions for larger lamdas (1:(k-1)) returned.
//        jerr = -k => convergence for kth lamda value not reached
//           after maxit (see above) iterations.
//        jerr = -10000-k => number of non zero coefficients along path
//           exceeds nx (see above) at kth lamda value.







//[[Rcpp::export]]
Rcpp::List make_matrix(const size_t p,Rcpp::DataFrame anno_df){
  auto [mat,names] = gen_anno_mat(p,anno_df);
  using namespace Rcpp;
  return(List::create(_["annomat"]=mat,_["names"]=names));
}


//[[Rcpp::export]]
Rcpp::List make_matrix_df(const size_t p,Rcpp::DataFrame anno_df){


  SparseDF spdf( anno_df, p,"SNP","feature");
  auto [mat,names] = gen_anno_mat(p,anno_df);

  using namespace Rcpp;
  return(List::create(_["annomat"]=spdf.getMat(),_["names"]=spdf.names()));
}


//[[Rcpp::export]]
Rcpp::List make_spmatrix_df(const size_t p,Rcpp::DataFrame anno_df){


  SparseDF spdf( anno_df, p,"SNP","feature");
  auto spm = spdf.getspMat();

  using namespace Rcpp;
  return(List::create(
      _["annomat"]=spm,
      _["names"]=spdf.names())
           );
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


  //  auto  my_t = locus_id.sexp_type();

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

//[[Rcpp::export]]
Rcpp::NumericMatrix logit_glmnet(Eigen::MatrixXd X,Rcpp::NumericVector y,double alpha=0,Rcpp::NumericVector lambda=Rcpp::NumericVector::create(0)){

  const size_t npar = X.cols()+1;
  const size_t nlambda=lambda.size();
  Lognet<Eigen::MatrixXd> logit(X,alpha,Rcpp::as<std::vector<double> >(lambda));
  Rcpp::NumericMatrix beta_m(npar,nlambda);
  gsl::span<double> yd(&y(0),y.size());

  logit.fit(yd);
  for(int i=0; i<nlambda; i++){
    gsl::span<double> beta(&beta_m(0,i),npar);
    logit.read_coeffs(beta);
  }
  return(beta_m);
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
