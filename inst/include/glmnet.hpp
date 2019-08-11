#pragma once

#include <gsl/span>
#include <RcppEigen.h>




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
    int pp1=p+1;
    retcoeff.resize(nlam,lmu);
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
  int n;
  int p;
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






template<>
inline void Lognet<Eigen::MatrixXd>::fit(const gsl::span<double> yd){


  y.col(0)=  Eigen::Map<Eigen::ArrayXd>(yd.data(),yd.size());
  y.col(1)=1-y.col(0).array();
  x=x_o;
  o.setZero();
  int jd = 0;
  int nc = 1;
  int pp1 = p+1;
  //  lmu=0;
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
