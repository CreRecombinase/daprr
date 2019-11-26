#include "elasticdonut.hpp"
#include <fmt/printf.h>
#include <fmt/format.h>

namespace elasticdonut {

  using SplitView = std::vector< Eigen::Map<Eigen::ArrayXd> >;



  double prior_to_pip(Eigen::Map<Eigen::ArrayXd> pip,const Eigen::Map<Eigen::ArrayXd> prior,const bool fix_pi0=true){

  double locus_pi0=1;
  std::transform(begin(prior),end(prior),begin(pip),[&locus_pi0](double prior){
							 locus_pi0*=(1-prior);
							 return(prior/(1-prior));
						       });
  if(locus_pi0 < 1e-100 && fix_pi0){
    locus_pi0 = 1e-100;
  }

  pip *= locus_pi0;
  // for(auto &p_v : pip)
  //   p_v *= locus_pi0;
  return locus_pi0;

  }


void transform_pip(Eigen::Map<Eigen::ArrayXd> pip,const Eigen::Map<Eigen::ArrayXd> BF,const double log10lik){

  std::transform(begin(pip), end(pip), cbegin(BF), begin(pip),
                 [log10lik](double prior, double log10_BF) {
                   return pow(10, (log10(prior) + log10_BF - log10lik));
                 });
}

double log10_lik(const double locus_pi0,const Eigen::Map<Eigen::ArrayXd> BF, const Eigen::Map<Eigen::ArrayXd> pip,const double* BF_max){
  if(BF_max==nullptr){
    BF_max =  &(*(std::max_element(cbegin(BF),cend(BF))));
  }
  double max_el = std::max(*BF_max,0.0);
  double sum = locus_pi0*pow(10,0-max_el);

  double log10lik =	std::inner_product(begin(BF),
					   end(BF),
					   begin(pip),
					   sum,std::plus<>(),
					   [max_el](double vec,double wts){
					     return(wts*pow(10,vec-max_el));});
  log10lik=max_el+log10(log10lik);
  return log10lik;
}



double E_step(const Eigen::Map<Eigen::ArrayXd> pip,const Eigen::Map<Eigen::ArrayXd> BF, const Eigen::Map<Eigen::ArrayXd> prior,double* BF_max){

  double locus_pi0=  prior_to_pip(pip,prior);

  double lik10 = log10_lik(locus_pi0,BF,pip,BF_max);
  transform_pip(pip,BF,lik10);

  return(lik10);
}

double fdr(double log10lik,const Eigen::Map<Eigen::ArrayXd> BF, const Eigen::Map<Eigen::ArrayXd> prior,Eigen::Map<Eigen::ArrayXd> pip,const double* BF_max){


  double locus_pi0=prior_to_pip(prior,pip,false);

  transform_pip(pip,BF,log10lik);

  log10lik = log10_lik(locus_pi0,BF,prior,BF_max);

  if (!BF_max) {
    BF_max = &(*(std::max_element(cbegin(BF), cend(BF))));
  }
  double max_el = std::max(*BF_max,0.0);
  log10lik = max_el + log10(log10lik);

  double fdr = std::pow(10, log10(locus_pi0) - log10lik);

  return fdr;
}



// double evaluate_likelihood(GroupedView &prior, const Eigen::Map<Eigen::ArrayXd> beta,
//                            SumStatRegion &sumstats, const Net &logistic) {
//   logistic.predict(beta, prior.d_view);
//   return sumstats.E_steps(prior.r_view);
// }


  SumStatRegion::SumStatRegion(const GroupedView BF_):
    BF(BF_),
    p_vec(BF.p,0.0),
    p_view(BF.copy_view(p_vec)),
    BF_max(BF.nr,nullptr)
  {

  }
  double SumStatRegion::E_steps(const SplitView &r_prior){
    double loglik=0;
    auto r_pip = p_view.r_view;
    for(int i=0; i<BF.nr; i++){
      loglik+=E_step(r_pip[i],BF.r_view[i],r_prior[i],BF_max[i]);
    }
    return loglik;
  }
  size_t SumStatRegion::size() const { return BF.p; }
  Eigen::Map<Eigen::ArrayXd> SumStatRegion::pip() { return p_view.d_view; }





std::vector<double> make_BF(Rcpp::NumericVector z_hat){
    std::vector<double> BFv(z_hat.size());
    BF bf;
    std::transform(std::begin(z_hat),std::end(z_hat),BFv.begin(),[&bf](double z){
							   return(bf.compute_log10_BF(z));});
    return BFv;
}

std::vector<double> make_BF(const Eigen::Map<Eigen::ArrayXd> z_hat){
  std::vector<double> BFv(z_hat.size());
  BF bf;
  std::transform(begin(z_hat),end(z_hat),begin(BFv),[&bf](double z){
							 return(bf.compute_log10_BF(z));});
  return BFv;
}


  BF::BF():wt(1/4.0),
	 l10(log(10)){
  kv={1,4,16,25};
  for(int i=0; i<4; i++){
    kva[i]= 0.5*log(1/(1+kv[i]));
    kvb[i]=0.5*(kv[i]/(1+kv[i]));
  }
}

double BF::compute_log10_BF(const double z_score) const {
  const double z2 = pow(z_score, 2.0);
  double max = 0;

  for(int i=0;i<4;i++){
    double val = (kva[i] + z2*kvb[i])/l10;
    max = val > max ? val : max;
    kvc[i]=val;
  }
  double sum=0;
  for(int i=0; i<4; i++){
    sum+= wt*pow(10,(kvc[i]-max));
  }
  return(max+log10(sum));
}










  double ElasticDonut::fine_optimize_beta(Eigen::Map<Eigen::ArrayXd> new_beta,const size_t index,GroupedView &null_prior,double &log10_lik){

    const double gr = (std::sqrt(5)-1)/2.0;
    double a = new_beta[index];
    double b = 0;
    if (a > b) {
      b = a;
      a = 0;
    }
    double c = b - gr*(b-a);
    double d = a + gr*(b-a);
    //    double thresh = 1e-3;
    const double thresh	 = std::fabs(c-d) < 1e-3 ? std::fabs(c-d)/2.0 :	1e-3;
    // if(fabs(c-d)<thresh){
    //   thresh = fabs(c-d)/2.0;
    // }
    double fc;
    double fd;
    while((d-c)>thresh){

      new_beta[index]=c;
      logistic.predict(new_beta, null_prior.d_view);
      fc = sumstats.E_steps(null_prior.r_view);
      // fc = evaluate_likelihood(null_pip,null_prior,null_beta,sumstats,logistic);

      new_beta[index]=d;
      logistic.predict(new_beta, null_prior.d_view);
      fd = sumstats.E_steps(null_prior.r_view);

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

    log10_lik = fc;

    return (b + a) / 2;
  }


ElasticDonut::ElasticDonut(GroupedView BF_v, Net &logn,
                             const double prior_init, double EM_thresh_)
      : sumstats(BF_v), logistic(logn), prior_vec(BF_v.p, prior_init),
        null_prior_vec(0), prior_view(BF_v.copy_view(prior_vec)),
        names(logistic.get_names()), beta(logistic.feature_num(), 0.0),
        sd(logistic.feature_num(), 0.0), curr_log10_lik(0.0),
        EM_thresh(EM_thresh_) {
    names.insert(begin(names), "Intercept");
  }

double ElasticDonut::fit(bool keep_prior){

    auto &prior_r = prior_view.r_view;
    auto &prior =  prior_view.d_view;
    int iter_ct=0;
    int iter_max=10;
    double last_log10_lik = -9999999;
    Rcpp::Rcerr<<fmt::sprintf("Iter\tloglik\t");
    for(auto &n : names){
      Rcpp::Rcerr<<n<<"\t";
    }
    Rcpp::Rcerr<<std::endl;
    while (iter_ct < iter_max &&
           fabs(curr_log10_lik - last_log10_lik) > EM_thresh) {
      last_log10_lik = curr_log10_lik;
      curr_log10_lik = sumstats.E_steps(prior_r);
      logistic.fit(sumstats.pip());
      Rcpp::Rcerr<<iter_ct<<"\t"<<curr_log10_lik/log10(exp(1));
      logistic.read_coeffs(toMapXd(beta));
      for (auto bv : beta) {
        Rcpp::Rcerr << "\t" << bv;
      }
      Rcpp::Rcerr << std::endl;
      logistic.predict(toMapXd(beta), prior);
      iter_ct++;
    }
    std::vector<double> null_beta=beta;
    null_prior_vec.resize(prior_vec.size());
    auto null_pview =	prior_view.copy_view(null_prior_vec);
    const size_t bs = beta.size();
    for (int i = 0; i < bs; i++) {
      double tbeta=null_beta[i];
      null_beta[i]=0;
      logistic.predict(toMapXd(null_beta),null_pview.d_view);
      double null_lik = sumstats.E_steps(null_pview.r_view);
      double tdiff = (curr_log10_lik - null_lik)/log10(exp(1));
      if (tdiff < 0 && i > 0) {
        double log10_lik = curr_log10_lik;
	null_beta[i]=tbeta;
        beta[i] = fine_optimize_beta(toMapXd(null_beta), i, null_pview, log10_lik);
        tdiff = (log10_lik - null_lik) / log10(exp(1));
      }
      if (tdiff < 1e-8) {
        Rcpp::Rcerr << "log_lik-null_log_lik<1e-8: " << tdiff
                    << " for term: " << names[i] << std::endl;
        tdiff = 1e-8;
      }
      sd[i] = std::fabs(beta[i]) / std::sqrt(2 * tdiff);
      null_beta[i] = tbeta;
    }

    return curr_log10_lik;
  }


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
  auto BF_d =	elasticdonut::make_BF(z);

  //  int *locb =&(*std::begin(locus));
  auto splt = make_splitter((locus));
  GroupedView BF_v(splt.split_view(&BF_d.front()),p);

  elasticdonut::Lognet ln(X,alpha,Rcpp::as<std::vector<double> >(lambda));
  elasticdonut::ElasticDonut ed(BF_v,ln,prior_init,EM_thresh);
  auto lik = ed.fit(true);
  // auto final_res=elasticdonut::fit_torus<Eigen::MatrixXd>(buff,BF_v,annotation,EM_thresh,alpha,Rcpp::as<std::vector<double>>(lambda));

  using namespace Rcpp;
  NumericVector beta = wrap(ed.beta);
  beta.attr("names")=wrap(ed.names);
  return List::create(_["beta"]=beta,
		      _["prior"]=wrap(ed.prior_vec),
		      _["sd"]=wrap(ed.sd),
		      _["lik"]=wrap(lik));
}




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
Rcpp::List elastic_donut_sp(Rcpp::IntegerVector locus,Rcpp::NumericVector z, Rcpp::S4 X,double EM_thresh=0.05,const double alpha=0,Rcpp::NumericVector lambda=Rcpp::NumericVector::create(0),const double prior_init=1e-3){

  // std::vector<std::string> names = Rcpp::as<std::vector<std::string>>(Rcpp::colnames(X));
  const size_t p = locus.size();
  auto BF_d =	elasticdonut::make_BF(z);

  auto splt = make_splitter(locus);
  GroupedView BF_v(splt.split_view(&BF_d.front()),p);

  elasticdonut::spLognet ln(X,alpha,Rcpp::as<std::vector<double> >(lambda));

  elasticdonut::ElasticDonut ed(BF_v,ln,prior_init,EM_thresh);
  auto lik = ed.fit(true);
  // auto final_res=elasticdonut::fit_torus<Eigen::MatrixXd>(buff,BF_v,annotation,EM_thresh,alpha,Rcpp::as<std::vector<double>>(lambda));

  using namespace Rcpp;
  NumericVector beta = wrap(ed.beta);
  beta.attr("names")=wrap(ed.names);
  return List::create(_["beta"]=beta,
		      _["prior"]=wrap(ed.prior_vec),
		      _["sd"]=wrap(ed.sd),
		      _["lik"]=wrap(lik));
}

//' @export
//[[Rcpp::export]]
Rcpp::NumericVector zhat2BF(Rcpp::NumericVector z_hat){
  Rcpp::NumericVector BFv(z_hat.size());
  elasticdonut::BF bf;
  std::transform(std::begin(z_hat),std::end(z_hat),std::begin(BFv),[&bf](double z){
							 return(bf.compute_log10_BF(z));});
  return BFv;
}



//' @export
//[[Rcpp::export]]
Rcpp::IntegerVector new_splitter(Rcpp::IntegerVector locus_id) {
  auto splt = make_splitter(locus_id);
  return splt.export_rle();
}



//' @export
//[[Rcpp::export]]
Rcpp::NumericVector Esteps(Rcpp::NumericVector BF, Rcpp::NumericVector prior,Rcpp::IntegerVector rle) {


  splitter splt(rle);
  double *BFb =	&(BF[0]);
  Eigen::Map<Eigen::ArrayXd> BF_s(BFb,BF.size());
  elasticdonut::SumStatRegion BF_v(splt.group_view(&(*std::begin(BF))));
  auto prior_v = BF_v.BF.copy_view(prior);
  // Rcpp::NumericVector pip(prior.size(),0.0);
  // auto pip_v = prior_v.copy_view(pip);
  auto ret = BF_v.E_steps(prior_v.r_view);

  using namespace Rcpp;
  auto pip_v = BF_v.p_vec;
  pip_v.attr("lik")=ret;
  return pip_v;
}
