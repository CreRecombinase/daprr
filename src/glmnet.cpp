#include "glmnet.hpp"

namespace elasticdonut {



  dgCMatrix::dgCMatrix(Rcpp::S4 mat_):mat(mat_),spX(Rcpp::as<Eigen::Map<Eigen::SparseMatrix<double>>>(mat)){
      using namespace Rcpp;
      CharacterVector cl = mat.slot("class");
      if(as<std::string>(cl)!="dgCMatrix"){
	Rcerr<<"mat_ is of class: "<<cl<<" not of class dgCMatrix"<<std::endl;
	stop("invalid type passed to dgCMatrix()");
      }
      IntegerVector tix=mat.slot("p");
      IntegerVector tjx=mat.slot("i");
      NumericVector tx = mat.slot("x");
      ix=as<Eigen::ArrayXi>(tix)+1;
      jx=as<Eigen::ArrayXi>(tjx)+1;
      xd=as<Eigen::ArrayXd>(tx);
      SEXP dnls = mat.slot("Dimnames");

      if (Rf_isNull(dnls)) {
        stop("X slot 'Dimnames' cannot be null!");
      }
      List dnl(dnls);
      StringVector cn(dnl[1]);
      colnames = as<std::vector<std::string>>(cn);
    }



Lognet::Lognet(const Rcpp::NumericMatrix X_, const double alpha_,std::vector<double> lambda,const double thresh,const int maxiter):
  x_o(Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(X_)),
  n(x_o.rows()),
  p(x_o.cols()),
  ulam(lambda),
  nlam(ulam.size()),
  coeff(nlam,p),
  names(Rcpp::as<std::vector<std::string>>(Rcpp::colnames(X_))),
  alpha(alpha_),
  x(Rcpp::as<Eigen::MatrixXd>(X_)),
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
  jerr(0),
  stop_on_nonfatal(false){
  if (ulam.size() == 0) {
    Rcpp::stop("lambda cannot be empty");
  }

  for (int i = 0; i < p; i++) {
    interval_mat(0, i) = -9e35;
    interval_mat(1, i) = 9e35;
  }
}

  spLognet::spLognet(Rcpp::S4 X_, const double alpha_,std::vector<double> lambda,const double thresh,const int maxiter):
    x_o(X_),
    n(x_o.MapMatrix().rows()),
    p(x_o.MapMatrix().cols()),
    ulam(lambda),
    nlam(ulam.size()),
    coeff(nlam,p),
    alpha(alpha_),
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
    jerr(0),
    stop_on_nonfatal(false){
    if (ulam.size() == 0) {
      Rcpp::stop("lambda cannot be empty");
    }
    using namespace Rcpp;


    for (int i = 0; i < p; i++) {
      interval_mat(0, i) = -9e35;
      interval_mat(1, i) = 9e35;
    }
  }

void Lognet::fit(const gsl::span<double> yd){

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
      Rcpp::Rcerr << "fatal error	in lognet: " << jerr << std::endl;
      if (jerr < 7777)
        Rcpp::stop("memory allocation error");

      if(jerr == 7777){
	Rcpp::stop("all predictors have zero variance!");
      }

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
      Rcpp::Rcerr << "non-fatal error in lognet: " << jerr << std::endl;
      if (stop_on_nonfatal) {
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
      } else {
        if (abs(jerr) < 1000)
          Rcpp::Rcerr << "convergence not reached for class: " +
                             std::to_string(abs(jerr))
                      << std::endl;
        if (abs(jerr) > 10000)
          Rcpp::Rcerr
              << "number of non zero coefficients along path exceeds nx at "
                 "lambda index :" +
                     std::to_string(abs(jerr) - 10000)
              << std::endl;
        if (abs(jerr) > 20000)
          Rcpp::Rcerr << "max(p*(1-p)) < 1e-06 at lambda index: " +
                             std::to_string(abs(jerr) - 20000)
                      << std::endl;
      }
    }
  }

  coeff.unpack_coeffs();

  // Eigen::Map<Eigen::ArrayXd > coeff_v(ca.data(),lmu*ni);
}

void spLognet::fit(const gsl::span<double> yd){

  y.col(0)=  Eigen::Map<Eigen::ArrayXd>(yd.data(),yd.size());
  y.col(1)=1-y.col(0).array();

  //  x=x_o;
  o.setZero();
  int pp1 = p+1;
  int jd = 0;
  int one = 1;

  splognet_(&alpha,
	    &n,
	    &p,
	    &one,
	    x_o.x().data(),
	    x_o.p().data(),
	    x_o.i().data(),
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

      if (jerr == 7777) {
        Rcpp::stop("all predictors have zero variance!");
      }

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
      Rcpp::Rcerr << "non-fatal error in splognet: " << jerr << std::endl;
      if (stop_on_nonfatal) {
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
      } else {
        if (abs(jerr) < 1000)
          Rcpp::Rcerr << "convergence not reached for class: " +
                             std::to_string(abs(jerr))
                      << std::endl;
        if (abs(jerr) > 10000)
          Rcpp::Rcerr
              << "number of non zero coefficients along path exceeds nx at "
                 "lambda index :" +
                     std::to_string(abs(jerr) - 10000)
              << std::endl;
        if (abs(jerr) > 20000)
          Rcpp::Rcerr << "max(p*(1-p)) < 1e-06 at lambda index: " +
                             std::to_string(abs(jerr) - 20000)
                      << std::endl;
      }
    }
    }

    coeff.unpack_coeffs();
}

void Lognet::read_coeffs(gsl::span<double> beta,int index){
  beta[0]=coeff.a0[index];
  Eigen::VectorXd tcol=coeff.retcoeff.col(index);
  std::copy_n(tcol.data(),tcol.size(),beta.begin()+1);
}

void spLognet::read_coeffs(gsl::span<double> beta,int index){
  beta[0]=coeff.a0[index];
  Eigen::VectorXd tcol=coeff.retcoeff.col(index);
  std::copy_n(tcol.data(),tcol.size(),beta.begin()+1);
}

void Lognet::predict(const gsl::span<double> beta,
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



void spLognet::predict(const gsl::span<double> beta,
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

  this->y.col(0) = (x_o.MapMatrix() * beta_v);

  this->y.col(0) = 1 / (1 + (-(beta[0] + this->y.col(0).array())).exp());
  y_v = this->y.col(0);
}

} // namespace elasticdonut
