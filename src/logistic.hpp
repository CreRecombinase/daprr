#pragma once
#include <RcppGSL.h>
#include <RcppEigen.h>
#include "gsl/span"
/* Mixed interface */
void logistic_mixed_pred(gsl_vector *beta  // Vector of parameters length = 1 + Sum_k(C_k - 1) + Kc
			 ,gsl_matrix_int *X  //Matrix Nobs x K 
			 ,gsl_vector_int *nlev // Vector with number categories
			 ,gsl_matrix *Xc   // continuous covariates  Matrix Nobs x Kc 
			 ,gsl_vector *yhat //Vector of prob. predicted by the logistic
			 );
 
int logistic_mixed_fit(gsl_vector *beta  // Vector of parameters length = 1 + Sum_k(C_k - 1) + Kc
		       ,gsl_matrix_int *X  //Matrix Nobs x K 
		       ,gsl_vector_int *nlev // Vector with number categories
		       ,gsl_matrix *Xc   // continuous covariates  Matrix Nobs x Kc 
		       ,gsl_vector *y //Vector of prob. to predict
		       ,double lambdaL1 // Regularization L1 0.0 if not used
		       ,double lambdaL2); // Regularization L2 0.0 if not used

double fLogit_mixed(gsl_vector *beta
		    ,gsl_matrix_int *X
		    ,gsl_vector_int *nlev
		    ,gsl_matrix *Xc   // continuous covariates  Matrix Nobs x Kc 
		    ,gsl_vector *y
		    ,double lambdaL1
		    ,double lambdaL2);


/* Categorical only interface */
void logistic_cat_pred(gsl_vector *beta  // Vector of parameters length = 1 + Sum_k(C_k - 1) + Kc
		       ,gsl_matrix_int *X  //Matrix Nobs x K 
		       ,gsl_vector_int *nlev // Vector with number categories
		       ,gsl_vector *yhat //Vector of prob. predicted by the logistic
		       );
 
int logistic_cat_fit(gsl_vector *beta  // Vector of parameters length = 1 + Sum_k(C_k - 1) + Kc
		     ,gsl_matrix_int *X  //Matrix Nobs x K 
		     ,gsl_vector_int *nlev // Vector with number categories
		     ,gsl_vector *y //Vector of prob. to predict
		     ,double lambdaL1 // Regularization L1 0.0 if not used
		     ,double lambdaL2); // Regularization L2 0.0 if not used

double fLogit_cat(gsl_vector *beta
		  ,gsl_matrix_int *X
		  ,gsl_vector_int *nlev
		  ,gsl_vector *y
		  ,double lambdaL1
		  ,double lambdaL2);


/* Continuous only interface */
void logistic_cont_pred(gsl_vector *beta  // Vector of parameters length = 1 + Sum_k(C_k - 1) + Kc
			,gsl_matrix *Xc   // continuous covariates  Matrix Nobs x Kc 
			,gsl_vector *yhat //Vector of prob. predicted by the logistic
			);
 
int logistic_cont_fit(gsl_vector *beta  // Vector of parameters length = 1 + Sum_k(C_k - 1) + Kc
		      ,gsl_matrix *Xc   // continuous covariates  Matrix Nobs x Kc 
		      ,gsl_vector *y //Vector of prob. to predict
		      ,double lambdaL1 // Regularization L1 0.0 if not used
		      ,double lambdaL2); // Regularization L2 0.0 if not used

double fLogit_cont(gsl_vector *beta
		   ,gsl_matrix *Xc   // continuous covariates  Matrix Nobs x Kc 
		   ,gsl_vector *y
		   ,double lambdaL1
		   ,double lambdaL2);


// I need to bundle all the data that goes to the function to optimze together.

typedef struct{
  gsl_matrix_int *X;
  gsl_vector_int *nlev;
  gsl_vector *y;
  double lambdaL1;
  double lambdaL2;
}fix_parm_cat_T;




template<typename T,int Row=Eigen::Dynamic,int Col=Eigen::Dynamic>
inline gsl_matrix_int*	copy_matrix(const Eigen::Map<Eigen::Matrix<T,Row,Col>> X){

  const size_t n = X.rows();
  const size_t p = X.cols();

  gsl_matrix_int* retX = gsl_matrix_int_calloc(n,p);
  Eigen::Map<Eigen::Matrix<int,Row,Col,Eigen::RowMajor>> wrap_retX(retX->data,n,p);
  wrap_retX=X.template cast<int>();
  return(retX);
}


gsl_vector *view_span(gsl::span<double> data);


class Dap_logit{
  gsl_matrix_int *X;
  const size_t p;
  const size_t npar;
  const size_t n;
  double lambdaL1;
  double lambdaL2;
  gsl_vector_int *nlev;
  gsl_vector *beta;
  gsl_matrix *myH;
  gsl_vector *stBeta;
  gsl_vector *myG;
  gsl_vector *tau;
  // RcppGSL::matrix<double> ty;


public:
  Dap_logit(gsl_matrix_int *X_,double alpha,double lambda):
    X(X_),p(X->size2),npar(p+1),n(X->size1),lambdaL1(alpha*lambda),lambdaL2((1-alpha)*lambda),
    nlev(gsl_vector_int_calloc(p)),
    beta(gsl_vector_calloc(npar)),
    myH(gsl_matrix_alloc(npar,npar)),
    stBeta(gsl_vector_alloc(npar)),
    myG(gsl_vector_alloc(npar)),
    tau(gsl_vector_alloc(npar))
  {

    for(int j=0; j <p; j++){
      std::map<int, int> rcd;
      for (int i = 0; i < n; i++) {
        int val = gsl_matrix_int_get(X, i, j);
        rcd[val] = 1;
      }
      gsl_vector_int_set(nlev, j, rcd.size());
      //      Rcpp::Rcerr << "dlevel: " << j << ": " << rcd.size() << std::endl;
    }
  }

  void fit(gsl::span<double> data);
  void read_coeffs(gsl::span<double> sy);
  void predict(gsl::span<double> sbeta, gsl::span<double> syhat);
  ~Dap_logit(){
    gsl_matrix_int_free (X);
    gsl_vector_free (beta);
    gsl_vector_int_free (nlev);
    gsl_vector_free (tau);
    gsl_vector_free (stBeta);
    gsl_vector_free (myG);
    gsl_matrix_free (myH);
  }
};

class Logistic {

  // gsl_vector *beta;
  const int npar;
  const size_t p;
  // const size_t p;
  // const size_t nfeat;
  //  gsl_matrix_int *X;
  //  gsl_vector_int *nlev;
  //  gsl_vector *y;
  // double lambdaL1;
  // double lambdaL2;
  RcppGSL::matrix<double> ty;
  gsl_matrix *myH;
  gsl_vector *stBeta;

  gsl_vector *myG;
  gsl_vector *tau;

 public:
  Logistic(const size_t npar_,const size_t p);
  void fit_glmnet(RcppGSL::vector<double> beta ,RcppGSL::matrix<int> X,RcppGSL::vector<double> y,double lambdaL1,double lambdaL2);
  void fit_qbinom(RcppGSL::vector<double> beta ,RcppGSL::matrix<int> X,RcppGSL::vector<double> y);

  void fit(gsl_vector* beta,gsl_matrix_int *X,gsl_vector_int *nlev,gsl_vector *y,double lambdaL1,double lambdaL2);
  ~Logistic(){
    gsl_vector_free (tau);
    gsl_vector_free (stBeta);
    gsl_vector_free (myG);
    gsl_matrix_free (myH);
  }
};
