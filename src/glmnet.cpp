#include <glmnet.hpp>
#include "classdef.hpp"
#include "logistic.hpp"
#include <fmt/printf.h>
#include <fmt/format.h>

// [[Rcpp::depends(RcppEigen)]]

gsl_vector *view_span(gsl::span<double> data) {

  const size_t n = data.size();
  gsl_vector *ret = new gsl_vector;
  ret->size = n;
  ret->stride = 1;
  ret->data = data.data();
  ret->block = nullptr;
  ret->owner = 0;
  return (ret);
}


namespace elasticdonut {

double prior_to_pip(gsl::span<double> pip,const gsl::span<double> prior,const bool fix_pi0){

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



double E_step(const gsl::span<double> pip,const gsl::span<double> BF, const gsl::span<double> prior,std::optional<double> &BF_max){

  double locus_pi0=  prior_to_pip(pip,prior);

  double lik10 = log10_lik(locus_pi0,BF,pip,BF_max);
  transform_pip(pip,BF,lik10);

  return(lik10);
}

double fdr(double log10lik,const gsl::span<double> BF, const gsl::span<double> prior,gsl::span<double> pip,std::optional<double> &BF_max){


  double locus_pi0=prior_to_pip(prior,pip,false);

  transform_pip(pip,BF,log10lik);

  log10lik = log10_lik(locus_pi0,BF,prior,BF_max);

  if (!BF_max) {
    *BF_max = std::max(*(std::max_element(BF.begin(), BF.end())), 0.0);
  }
  double max_el = *BF_max;
  log10lik = max_el + log10(log10lik);

  double fdr = std::pow(10, log10(locus_pi0) - log10lik);

  return fdr;
}

  //SumStatRegion is a non-owning view of the summary statistics
  SumStatRegion::SumStatRegion(const GroupedView BF_):
    BF(BF_),
    p_vec(BF.p),
    p_view(BF.copy_view(p_vec)),
    BF_max(BF.nr,std::nullopt)
  {

  }
  double SumStatRegion::E_steps(const SplitView &r_prior)const{
    double loglik=0;
    auto r_pip = p_view.r_view;
    for(int i=0; i<BF.nr; i++){
      loglik+=E_step(r_pip[i],BF.r_view[i],r_prior[i],BF_max[i]);
    }
    return loglik;
  }
  size_t SumStatRegion::size() const { return BF.p; }
  gsl::span<double> SumStatRegion::pip() { return p_view.d_view; }

  template<typename T>
  double fit_torus(ParameterBuffer &buff,const GroupedView BF_v,const Eigen::Map<T> &X,const double EM_thresh=0.05,const double alpha=0,const std::vector<double> lambda={0}){

    SumStatRegion sumstats(BF_v);


  }



  template<typename T>
  double fit_dap_donut(ParameterBuffer &buff,const GroupedView BF_v,const Eigen::Map<T> &X,const double EM_thresh=0.05,const double alpha=0,const std::vector<double> lambda={0}){

    SumStatRegion sumstats(BF_v);

    auto &pip_r = buff.pip_v.r_view;
    auto &prior_r = buff.prior_v.r_view;
    auto &prior =  buff.prior_v.d_view;

    auto beta =	buff.beta;
    auto g_pip = view_span(pip);
    auto g_beta = view_span(beta);
    auto g_prior = view_span(prior);
    if (lambda.size() == 0) {
      Rcpp::stop("number of lambda values is 0 (before logistic)!");
    }
    using ST=typename T::Scalar;
    Dap_logit logistic(copy_matrix<ST>(X),alpha,lambda.front());
    //    logistic.predict(beta,prior);
    int iter_ct=0;
    int iter_max=20;
    double last_log10_lik = -9999999;
    double curr_log10_lik = 0;
    Rcpp::Rcerr<<fmt::sprintf("Iter\tloglik\tIntercept\t");
    for(auto &n : buff.names){
      Rcpp::Rcerr<<n<<"\t";
    }
    Rcpp::Rcerr<<std::endl;



    while (iter_ct < iter_max &&
           fabs(curr_log10_lik - last_log10_lik) > EM_thresh) {
      last_log10_lik = curr_log10_lik;
      curr_log10_lik = sumstats.E_steps(pip_r, prior_r);
      logistic.fit(pip);
      Rcpp::Rcerr<<iter_ct<<"\t"<<iter_ct<<"\t"<<curr_log10_lik/log10(exp(1));
      for (auto bv : beta) {
        Rcpp::Rcerr << "\t" << bv;
      }
      Rcpp::Rcerr << std::endl;
      logistic.read_coeffs(beta);
      logistic.predict(beta, prior);
      iter_ct++;
    }
    return curr_log10_lik;
  }

  // torus::controller con();
  // const size_t p = SNP.size();
  // SparseDF spdf(anno,p);

  // FileZscoreParser zsp(data_file);
  // FileAnnotationParser ap(con.snp_hash,annot_file);
  //  torus::controller con(buff_o);



  template<typename T>
  double evaluate_likelihood(GroupedView &pip,GroupedView &prior, const gsl::span<double> beta,SumStatRegion sumstats,const Lognet<T> &logistic){
    logistic.predict(beta,prior.d_view);
    return sumstats.E_steps(pip.r_view,prior.r_view);
  }


  template<typename T>
  double estimate_torus_null(const gsl::span<double> beta,const size_t index ,const SumStatRegion sumstats,const Lognet<T> &logistic,GroupedView &null_prior,GroupedView &null_pip){

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

      } else {

        a = c;
        c = d;
        d = a + gr * (b - a);
      }
    }

    curr_log10_lik = fc;

    return (b + a) / 2;
  }

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
  //       torus::simple_regression();
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
  // Rcpp::NumericMatrix logit_glmnet(Eigen::MatrixXd X,std::vector<double>
  // y,double alpha=0,Rcpp::NumericVector lambda=Rcpp::NumericVector::create(0)){
  //
  //   const size_t npar = X.cols()+1;
  //   const size_t nlambda=lambda.size();
  //   Lognet<Eigen::MatrixXd> logit(X,alpha,Rcpp::as<std::vector<double>
  //   >(lambda)); Rcpp::NumericMatrix beta_m(npar,nlambda);
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

} // namespace elasticdonut


//' Estimate torus model using elasticnet
//'
//' @param locus
//' @param z
//' @param X
//'
//' @return
//' @param
//'
//' @export
// [[Rcpp::export]]
Rcpp::List elastic_donut(Rcpp::IntegerVector locus,Rcpp::NumericVector z, Rcpp::NumericMatrix X,double EM_thresh=0.05,const double alpha=0,Rcpp::NumericVector lambda=Rcpp::NumericVector::create(0),const double prior_init=1e-3){

  std::vector<std::string> names = Rcpp::as<std::vector<std::string>>(Rcpp::colnames(X));
  const size_t p = locus.size();
  auto BF_d =	donut::make_BF(z);
  elasticdonut::ParameterData param(p,names,prior_init,lambda.size());
  auto annotation= Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(X);
  auto splt = elasticdonut::make_splitter(locus.begin(),locus.end());
  elasticdonut::ParameterBuffer buff(param,splt);
  elasticdonut::GroupedView BF_v(splt.split_view(&BF_d.front()),p);

  auto final_res=elasticdonut::fit_torus<Eigen::MatrixXd>(buff,BF_v,annotation,EM_thresh,alpha,Rcpp::as<std::vector<double>>(lambda));

  using namespace Rcpp;
  names.insert(names.begin(),"Intercept");
  rownames(param.beta)=Rcpp::wrap(names);
  return List::create(_["beta"]=param.beta,
		      _["prior"]=param.prior,
		      _["lik"]=Rcpp::wrap(final_res));
}




//' Estimate torus model using dap
//'
//' @param locus
//' @param z
//' @param X
//'
//' @return
//' @param
//'
//' @export
// [[Rcpp::export]]
Rcpp::List dap_donut(Rcpp::IntegerVector locus,Rcpp::NumericVector z, Rcpp::NumericMatrix X,double EM_thresh=0.05,const double alpha=0,Rcpp::NumericVector lambda=Rcpp::NumericVector::create(0),const double prior_init=1e-3){

  std::vector<std::string> names = Rcpp::as<std::vector<std::string>>(Rcpp::colnames(X));
  const size_t p = locus.size();
  auto BF_d =	donut::make_BF(z);
  elasticdonut::ParameterData param(p,names,prior_init,lambda.size());
  auto annotation= Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(X);
  auto splt = elasticdonut::make_splitter(locus.begin(),locus.end());
  elasticdonut::ParameterBuffer buff(param,splt);
  elasticdonut::GroupedView BF_v(splt.split_view(&BF_d.front()),p);

  auto final_res=elasticdonut::fit_dap_donut<Eigen::MatrixXd>(buff,BF_v,annotation,EM_thresh,alpha,Rcpp::as<std::vector<double>>(lambda));

  using namespace Rcpp;
  names.insert(names.begin(),"Intercept");
  rownames(param.beta)=Rcpp::wrap(names);
  return List::create(_["beta"]=param.beta,
		      _["prior"]=param.prior,
		      _["lik"]=Rcpp::wrap(final_res));
}



//' Predict using dap
//'
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector predict_dap(Rcpp::NumericVector beta,
                                Rcpp::IntegerMatrix X) {

  Rcpp::NumericVector y(X.rows());
  auto mX = Rcpp::as<Eigen::Map<Eigen::MatrixXi>>(X);
  auto gsl_X = copy_matrix(mX);
  gsl::span<double> beta_s(&beta[0], beta.size());
  gsl::span<double> y_s(&y[0], y.size());

  Dap_logit dl(gsl_X, 0, 0);
  dl.predict(beta_s, y_s);
  return (y);
}



//' Predict using dap
//'
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector predict_glmnet(Rcpp::NumericVector beta,
                                Rcpp::NumericMatrix X) {

  Rcpp::NumericVector y(X.rows());
  auto Xd = Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(X);
  gsl::span<double> beta_s(&beta[0], beta.size());
  gsl::span<double> y_s(&y[0], y.size());

  elasticdonut::Lognet<Eigen::MatrixXd> dl(Xd, 0, {0});
  dl.predict(beta_s, y_s);
  return (y);
}
