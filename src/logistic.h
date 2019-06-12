#ifndef LOGISTIC_H_   /* Include guard */
#define LOGISTIC_H_
#ifndef EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS
#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS
#endif
#include <RcppEigen.h>
using namespace Eigen;

/* Categorical only interface */
void logistic_cat_pred(Map<VectorXd>beta  // Vector of parameters length = 1 + Sum_k(C_k - 1)
		       ,SparseMatrix<double>X  //Matrix Nobs x K

		       ,Map<VectorXd>yhat //Vector of prob. predicted by the logistic
		       );


int logistic_cat_fit(Map<VectorXd>beta
		     ,SparseMatrix<double> X
		     ,Map<VectorXd>y
		     ,double lambdaL1
		     ,double lambdaL2);


double fLogit_cat(Map<VectorXd>beta
		  ,SparseMatrix<double>X
		  ,Map<VectorXd>y
		  ,double lambdaL1
		  ,double lambdaL2);

#endif // LOGISTIC_H_
