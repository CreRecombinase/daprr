#pragma once

#include <gsl/span>
#include <fmt/printf.h>
#include <fmt/format.h>
#include <RcppEigen.h>
#include "dataview.hpp"

using SplitView = std::vector< gsl::span<double> >;

namespace elasticdonut {






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
	       const double *x,
	       const int *ix,
	       const int *jx,
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




  class Net {
  public:
    virtual void fit(const gsl::span<double> yd) =0;
    virtual std::vector<std::string> get_names() const = 0;
    virtual void read_coeffs(gsl::span<double> beta,int index=0) = 0;
    virtual void predict(const gsl::span<double> beta,gsl::span<double> yhat) const = 0;
    virtual size_t snp_num() const = 0;
    virtual size_t feature_num() const = 0;
  };





  class Lognet: public Net{
  public:
    const Eigen::Map<Eigen::MatrixXd> x_o;
    int n;
    int p;
    std::vector<double> ulam;
    int nlam;
    comp_coeff coeff;
    std::vector<std::string> names;
    double alpha;

    Eigen::MatrixXd x;
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
    const bool stop_on_nonfatal;
  public:
    Lognet(const Rcpp::NumericMatrix X_, const double alpha_=0,std::vector<double> lambda= {0},const double thresh=1e-07,const int maxiter=100000);
    void fit(const gsl::span<double> yd);
    void read_coeffs(gsl::span<double> beta,int index=0);
    void predict(const gsl::span<double> beta,gsl::span<double> yhat)const;
    size_t snp_num() const { return n; }
    size_t feature_num() const { return p+1; }
    std::vector<std::string> get_names() const { return names; }
  };


  class dgCMatrix{
    Rcpp::S4 mat;
    Eigen::ArrayXi ix;
    Eigen::ArrayXi jx;
    Eigen::ArrayXd xd;
    Eigen::Map<Eigen::SparseMatrix<double>> spX;
    std::vector<std::string> colnames;
  public:
    dgCMatrix(Rcpp::S4 mat_);
    const Eigen::ArrayXi& p() const{
      return ix;
    }
    const Eigen::ArrayXi& i() const{
      return jx;
    }
    const Eigen::ArrayXd& x() const{
      return xd;
    }
    const std::vector<std::string>& names() const{
      return colnames;
    }
    const Eigen::Map<Eigen::SparseMatrix<double> > MapMatrix() const{
      return spX;
    }

  };

  class spLognet: public Net{
  public:
    const dgCMatrix x_o;
    int n;
    int p;
    std::vector<double> ulam;
    int nlam;
    comp_coeff coeff;

    double alpha ;


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
    const bool stop_on_nonfatal;
  public:
    spLognet(Rcpp::S4 X_, const double alpha_ = 0,
             std::vector<double> lambda = {0}, const double thresh = 1e-07,
             const int maxiter = 100000);
    void fit(const gsl::span<double> yd);
    void read_coeffs(gsl::span<double> beta, int index = 0);
    void predict(const gsl::span<double> beta, gsl::span<double> yhat) const;
    size_t snp_num() const { return n; }
    size_t feature_num() const { return p+1; }
    std::vector<std::string> get_names() const { return x_o.names(); }
  };

} // namespace elasticdonut
