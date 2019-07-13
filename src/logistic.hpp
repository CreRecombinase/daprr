#pragma once
#include <RcppGSL.h>
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



class Logistic{

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

  // double mLogLik;
  // int iter;
  // double maxchange;

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
