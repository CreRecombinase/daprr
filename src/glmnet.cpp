#include <glmnet.hpp>
#include "classdef.hpp"
namespace elasticdonut {

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
  GroupedView(const donut::splitter &splt_,
	      double* BF_p, const size_t p):r_view(splt_.split_view(BF_p)),
					    d_view(BF_p,p){}

};




//SumStatRegion is a non-owning view of the summary statistics
class SumStatRegion{
public:
  const donut::splitter& splt;
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
  SumStatRegion(const donut::splitter &splt_,
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



  double estimate_torus_sp(gsl::span<double> beta,GroupedView pip_v,GroupedView prior_v,const donut::splitter& splt,const gsl::span<double> BF,const Eigen::SparseMatrix<double> &X,const double EM_thresh=0.05){

  SumStatRegion sumstats(splt,BF.data(),prior_v.d_view.size());
  auto &pip_r = pip_v.r_view;
  auto &prior_r = prior_v.r_view;
  auto &prior =  prior_v.d_view;
  Lognet<Eigen::SparseMatrix<double>> logistic(X);
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


// 
// Rcpp::NumericMatrix logit_glmnet(Eigen::MatrixXd X,std::vector<double> y,double alpha=0,Rcpp::NumericVector lambda=Rcpp::NumericVector::create(0)){
// 
//   const size_t npar = X.cols()+1;
//   const size_t nlambda=lambda.size();
//   Lognet<Eigen::MatrixXd> logit(X,alpha,Rcpp::as<std::vector<double> >(lambda));
//   Rcpp::NumericMatrix beta_m(npar,nlambda);
// 
//   double *yp = y.data();
//   const size_t ys = y.size();
//   Rcpp::Rcerr<<yp<<std::endl;
//   gsl::span<double> yd(yp,ys);
// 
//   logit.fit(yd);
//   for(int i=0; i<nlambda; i++){
//     gsl::span<double> beta(&beta_m(0,i),npar);
//     logit.read_coeffs(beta);
//   }
//   return(beta_m);
// }

}
