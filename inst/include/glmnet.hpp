#pragma once

#include <gsl/span>
#include <fmt/printf.h>
#include <fmt/format.h>
#include <RcppEigen.h>
#include "torus.hpp"

namespace elasticdonut {

using SplitView = std::vector< gsl::span<double> >;

inline size_t total_size(const SplitView& r_view){
  return std::accumulate(r_view.begin(),r_view.end(),0,[](size_t s,auto &sp){
							 return s+sp.size();
						       });
}

class GroupedView{
public:
  SplitView r_view;
  const size_t p;
  const size_t nr;
  gsl::span<double> d_view;
  GroupedView(SplitView r_view_, size_t p_)
    : r_view(r_view_), p(p_),nr(r_view.size()), d_view(r_view.begin()->data(), p) {}
  GroupedView(SplitView r_view_)
      : r_view(r_view_), p(total_size(r_view)), nr(r_view.size()),
        d_view(r_view.begin()->data(), p) {}
  GroupedView copy_view(gsl::span<double>	o_data) const {
    SplitView tr_view(nr);
    const size_t offset=0;
    auto ob=o_data.begin();
    std::transform(r_view.cbegin(),r_view.cend(),tr_view.begin(),[&ob](const auto sp) mutable{
								   auto ret = gsl::span<double>(ob,sp.size());
								   ob=ret.end();
								   return(ret);
								 });
    return GroupedView(std::move(tr_view),p);
  }
};

class ParameterData {
public:
  const size_t p;
  std::vector<std::string> colnames;
  const size_t k;
  const size_t num_params;
  size_t n_lambda;
  double prior_init;
  Rcpp::NumericMatrix beta;
  Rcpp::NumericMatrix prior;
  ParameterData(size_t p_, std::vector<std::string> colnames_, double prior_init=1e-3,size_t n_lambda_ = 1)
    : p(p_),colnames(colnames_),
      k(colnames.size()),
      num_params(k+1),
      n_lambda(n_lambda_),
      beta(num_params, n_lambda),
      prior(p, n_lambda) {

    std::fill(prior.begin(),prior.end(),prior_init);

  }

  // ParameterData(Rcpp::NumericMatrix beta_, Rcpp::NumericMatrix pip_,
  //               Rcpp::NumericMatrix prior_, Rcpp::StringVector cnames)
  //     : p(pip_.nrow()), colnames(Rcpp::as<std::vector<std::string>>(cnames)),
  //       k(beta_.nrow()-1),n_lambda(prior_.ncol()),prior_init(*prior_.begin()), beta(Rcpp::clone(beta_)), pip(Rcpp::clone(pip_)), prior(Rcpp::clone(prior_)) {}
  void reset(){
    std::fill(prior.begin(),prior.end(),prior_init);
    std::fill(beta.begin(),beta.end(),0);
  }
  ParameterData(Rcpp::NumericVector beta_,
                Rcpp::NumericVector prior_, Rcpp::StringVector cnames)
    : p(prior_.size()), colnames(Rcpp::as<std::vector<std::string>>(cnames)),
      k(beta_.size()-1),num_params(k+1),n_lambda(1),prior_init(*prior_.begin()), beta(k+1,n_lambda,beta_.cbegin()), prior(p,n_lambda,prior_.cbegin()) {}

};


  class ParameterBuffer {
  public:
    const size_t p;
    const size_t k;
    const size_t num_params;
    GroupedView prior_v;
    gsl::span<double> beta;
    std::vector<std::string> names;
    ParameterBuffer(ParameterData &data, splitter &splt)
      : p(data.p), k(data.k),num_params(data.num_params),
	prior_v(splt.split_view(data.prior.begin()), p),
	beta(data.beta.begin(), num_params), names(data.colnames) {}
  };

  class SumStatRegion {
  public:
  private:
    mutable std::vector<std::optional<double>> BF_max;
    mutable std::vector<double> p_vec;
    GroupedView p_view;
    const GroupedView BF;
    //  double E_step(
  public:
    SumStatRegion(const GroupedView BF_);
    double E_steps(const SplitView &r_prior) const;
    size_t size() const;
    gsl::span<double> pip();
  };

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





// c call lognet (alpha,n,p,one,x,y,o,use_f,vp,cl,ne,nx,nlam,flmin,ulam,thr,isd,
// c              intr,maxit,kopt,lmu,a0,ca,ia,nin,dev0,fdev,alm,nlp,jerr)
// c call lognet (parm,no,ni,nc,x,y,o,jd,vp,cl,ne,nx,nlam,flmin,ulam,thr,isd,
// c              intr,maxit,kopt,lmu,a0,ca,ia,nin,dev0,fdev,alm,nlp,jerr)

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
//c call lsolns(p,pp1,nc,lmu,ca,ia,nin,b)
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
    if(nlam==0){
      Rcpp::stop("number of lambda must be greater than 0");
    }else{
      //      Rcpp::Rcerr<<"number of lambda values is:	"<<nlam<<std::endl;
    }

  if(zero){
    a0.setZero();
    ca.setZero();
    nin.setZero();
    alm.setZero();
  }
  }
  // void predict(Eigen::Map<Eigen::ArrayXd> b
  void setZero(){
    a0.setZero();
    ca.setZero();
    nin.setZero();
    alm.setZero();
  }
  void unpack_coeffs(){
    if(lmu==0){
      Rcpp::stop("number of (actual) lambda values is 0!");
    }

    int	one=1;
    int pp1=p+1;
    retcoeff.resize(p,lmu);
    F77_SUB(lsolns)(&p,&p,&one,&lmu,ca.data(),ia.data(),nin.data(),retcoeff.data());
  }

  Eigen::ArrayXd a0;
  Eigen::MatrixXd ca;
  Eigen::ArrayXi nin;
  Eigen::ArrayXi ia;
  Eigen::ArrayXd alm;
  Eigen::MatrixXd retcoeff;
  int lmu;
};



// template<typename T>
// class Logdap{
//   Logdap(const T &X_, const double alpha_=0,std::vector<double> lambda= {0},const double thresh=1e-07,const int maxiter=100000):





template<typename T>
class Lognet{
public:
  int n;
  int p;
  std::vector<double> ulam;
  int nlam;
  comp_coeff coeff;
  double alpha ;
  const Eigen::Map<T> x_o;
  T x;
  mutable Eigen::MatrixXd y;
  Eigen::MatrixXd o;
  Eigen::ArrayXd one_vec;
  Eigen::MatrixXd interval_mat;
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
  Lognet(const Eigen::Map<T> X_, const double alpha_=0,std::vector<double> lambda= {0},const double thresh=1e-07,const int maxiter=100000):
    n(X_.rows()),
    p(X_.cols()),
    ulam(lambda),
    nlam(ulam.size()),
    coeff(nlam,p),
    alpha(alpha_),
    x_o(X_),
    x(X_),
    y(n,2),
    o(n,1),
    one_vec(Eigen::ArrayXd::Ones(p)),
    interval_mat(2,p),
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
    if(ulam.size()==0){
      static_assert(std::is_same_v<typename T::Scalar,double>,"Matrix must be of scalar type");
      Rcpp::stop("lambda cannot be empty");
    }

    for(int i=0; i<p; i++){
      interval_mat(0,i)=-9e35;
      interval_mat(1,i)=9e35;
    }

    //  interval_mat.colwise() = (Eigen::VectorXd(2) <<,9e35).finished();

  }
  void fit(const gsl::span<double> yd);
  void read_coeffs(gsl::span<double> beta,int index=0);
  void predict(const gsl::span<double> beta,gsl::span<double> yhat)const;
};






template<>
inline void Lognet<Eigen::MatrixXd>::fit(const gsl::span<double> yd){

  y.col(0)=  Eigen::Map<Eigen::ArrayXd>(yd.data(),yd.size());
  y.col(1)=1-y.col(0).array();
  x=x_o;
  o.setZero();
  int jd = 0;
  int nc = 1;
  int pp1 = p+1;
  // if(coeff.lmu==0){
  //   Rcpp::stop("number of lambda values is 0 (before fit)!");
  // }else{
  //   Rcpp::Rcerr<<"number of lambda values is :"<<coeff.lmu<<std::endl;
  // }


  F77_SUB(lognet)(&alpha,
		  &n,
		  &p,
		  &nc,
		  x.data(),
		  y.data(),
		  o.data(),
		  &jd,
		  one_vec.data(),
		  interval_mat.data(),
		  &pp1,
		  &p,
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

  if(jerr!=0){
    if (jerr > 0) {
      Rcpp::Rcerr << "fatal error	in splognet: " << jerr << std::endl;
      if (jerr < 7777)
        Rcpp::stop("memory allocation error");

      if (jerr > 8000 && jerr < 9000)
        Rcpp::stop("null probability < 1e-5 for class :" +
                   std::to_string(jerr - 8000));

      if (jerr > 9000 && jerr < 10000)
        Rcpp::stop("null probability > 1.0 - 1e-5 for class :" +
                   std::to_string(jerr - 9000));

      if (jerr == 10000)
        Rcpp::stop("maxval(vp) <= 0.0");

      if (jerr == 90000)
        Rcpp::stop("Bounds adjustnemt non convergence");
    }
    if (jerr < 0) {
      Rcpp::Rcerr << "noon-fatal error in splognet: " << jerr << std::endl;

      if (abs(jerr) < 1000)
        Rcpp::stop("convergence not reached for class: " +
                   std::to_string(abs(jerr)));

      if (abs(jerr) > 10000)
        Rcpp::stop("number of non zero coefficients along path exceeds nx at "
                   "lambda index :" +
                   std::to_string(abs(jerr) - 10000));

      if (abs(jerr) > 20000)
        Rcpp::stop("max(p*(1-p)) < 1e-06 at lambda index: " +
                   std::to_string(abs(jerr) - 20000));
    }
  }




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
  int pp1 = p+1;
  int jd = 0;
  int one = 1;

  splognet_(&alpha,
	    &n,
	    &p,
	    &one,
	    x.valuePtr(),
	    x.innerIndexPtr(),
	    x.outerIndexPtr(),
	    y.data(),
	    o.data(),
	    &jd,
	    one_vec.data(),
	    interval_mat.data(),
	    &pp1,
	    &p,
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

    if(jerr!=0){
    if (jerr > 0) {
      Rcpp::Rcerr << "fatal error	in splognet: " << jerr << std::endl;
      if (jerr < 7777)
        Rcpp::stop("memory allocation error");

      if (jerr > 8000 && jerr < 9000)
        Rcpp::stop("null probability < 1e-5 for class :" +
                   std::to_string(jerr - 8000));

      if (jerr > 9000 && jerr < 10000)
        Rcpp::stop("null probability > 1.0 - 1e-5 for class :" +
                   std::to_string(jerr - 9000));

      if (jerr == 10000)
        Rcpp::stop("maxval(vp) <= 0.0");

      if (jerr == 90000)
        Rcpp::stop("Bounds adjustnemt non convergence");
    }
    if (jerr < 0) {
      Rcpp::Rcerr << "noon-fatal error in splognet: " << jerr << std::endl;

      if (abs(jerr) < 1000)
        Rcpp::stop("convergence not reached for class: " +
                   std::to_string(abs(jerr)));

      if (abs(jerr) > 10000)
        Rcpp::stop("number of non zero coefficients along path exceeds nx at "
                   "lambda index :" +
                   std::to_string(abs(jerr) - 10000));

      if (abs(jerr) > 20000)
        Rcpp::stop("max(p*(1-p)) < 1e-06 at lambda index: " +
                   std::to_string(abs(jerr) - 20000));
    }
  }

    coeff.unpack_coeffs();

}



template<typename T>
void Lognet<T>::read_coeffs(gsl::span<double> beta,int index){
  beta[0]=coeff.a0[index];
  Eigen::VectorXd tcol=coeff.retcoeff.col(index);
  std::copy_n(tcol.data(),tcol.size(),beta.begin()+1);
}

template <typename T>
inline void Lognet<T>::predict(const gsl::span<double> beta,
                               gsl::span<double> yhat) const {

  if (beta.size() != (p + 1)) {
    Rcpp::stop("size of beta is " + std::to_string(beta.size()) + " in predict");
  }
  if (yhat.size() != n) {
    Rcpp::stop("size of y is " + std::to_string(yhat.size()) + " and not " +
               std::to_string(n) + "in predict");
  }
  Eigen::Map<Eigen::VectorXd> beta_v(beta.subspan(1, p).data(), p);
  Eigen::Map<Eigen::VectorXd> y_v(yhat.data(), yhat.size());

  this->y.col(0) = (x_o * beta_v);

  this->y.col(0) = 1 / (1 + (-(beta[0] + this->y.col(0).array())).exp());
  y_v = this->y.col(0);
}

  class	ElasticDonut{
  public:
    SumStatRegion sumstats;
    Lognet<Eigen::MatrixXd> logistic;
    std::vector<double> prior_vec;
    GroupedView prior_view;
    std::vector< std::string> names;
    std::vector<double> beta;
    std::vector<double> diff;

    const double EM_thresh;

    ElasticDonut(GroupedView BF_v, const Rcpp::NumericMatrix X,const double prior_init=1e-3,
                 double EM_thresh_ = 0.05, const double alpha = 0,
                 const std::vector<double> lambda = {0})
      : sumstats(BF_v), logistic(Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(X), alpha, lambda),
	prior_vec(BF_v.p,prior_init),prior_view(BF_v.copy_view(prior_vec)),names(Rcpp::as<std::vector<std::string> >(Rcpp::colnames(X))),beta(logistic.p+1,0.0),diff(logistic.p+1,0.0),EM_thresh(EM_thresh_) {}
    double fit(){

      auto &prior_r = prior_view.r_view;
      auto &prior =  prior_view.d_view;
      int iter_ct=0;
      int iter_max=20;
    double last_log10_lik = -9999999;
    double curr_log10_lik = 0;
    Rcpp::Rcerr<<fmt::sprintf("Iter\tloglik\tIntercept\t");
    for(auto &n : names){
      Rcpp::Rcerr<<n<<"\t";
    }
    Rcpp::Rcerr<<std::endl;
    while (iter_ct < iter_max &&
           fabs(curr_log10_lik - last_log10_lik) > EM_thresh) {
      last_log10_lik = curr_log10_lik;
      curr_log10_lik = sumstats.E_steps(prior_r);
      logistic.fit(sumstats.pip());
      Rcpp::Rcerr<<iter_ct<<"\t"<<iter_ct<<"\t"<<curr_log10_lik/log10(exp(1));
      logistic.read_coeffs(beta);
      for (auto bv : beta) {
        Rcpp::Rcerr << "\t" << bv;
      }
      Rcpp::Rcerr << std::endl;
      logistic.predict(beta, prior);
    }
    return curr_log10_lik;
    }
  };

} // namespace elasticdonut
