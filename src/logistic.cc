#include <stdio.h>
#include <math.h>

#include "logistic.h"

// I need to bundle all the data that goes to the function to optimze together. 

using namespace Eigen;

/***************/
/* Categorical */
/***************/

// I need to bundle all the data that goes to the function to optimze together. 
typedef struct{
  SparseMatrix<double>&X;

  Map<VectorXi> y;
  double lambdaL1;
  double lambdaL2;
}fix_parm_cat_T;


double fLogit_cat(Map<VectorXd> beta,
                  SparseMatrix<double> & X,
                  Map<VectorXd> y,
                  double lambdaL1,
                  double lambdaL2)
{
  int n = y.size();
  //  int k = X.size()2;
  int npar = beta.size();
  double total = 0;
  double aux = 0;
  /*   omp_set_num_threads(ompthr); */
  /*   /\* Changed loop start at 1 instead of 0 to avoid regularization of beta 0*\/ */

  total =beta.tail(npar-1).array().square().sum();


  total = (-total*lambdaL2/2);
  aux = beta.tail(npar-1).array().abs().sum();

  total = total-aux*lambdaL1;
  Eigen::VectorXd Xb = X.transpose()*beta.tail(npar-1);
  total += y.dot(Xb)-Xb.exp().unaryExpr([](const double x){
    return(std::log1p(x));
					}).sum();

//   for(int i = 0; i < n; ++i) {
//     double Xbetai=beta[0];
//     int iParm=1;
//     for(int k = 0; k < X.cols(); ++k) {
//       if(X.coeff(i,k)>0)
// 	Xbetai+=beta[X.coeff(i,k)-1+iParm];
//       iParm+=nlev[k]-1;
//     }
//     total += y[i]*Xbetai-gsl_sf_log_1plusx(gsl_sf_exp(Xbetai));
//   }
  return -total;
} 


void logistic_cat_pred(Map<VectorXd>beta  // Vector of parameters length = 1 + Sum_k(C_k - 1)
		       ,SparseMatrix<double>X  //Matrix Nobs x K

		       ,Map<VectorXd>yhat //Vector of prob. predicted by the logistic
		       )
{

  yhat = 1/(1+(-X.transpose()*beta.tail(beta.size()-1)).exp().array());
  // for(int i = 0; i < X.rows(); ++i) {
  //   double Xbetai=beta[0];
  //   int iParm=1;
  //   for(int k = 0; k < X.cols(); ++k) {
  //     if(X(i,k)>0)
  // 	Xbetai+=beta[X(i,k)-1+iParm];
  //     iParm+=nlev[k]-1;
  //   }
  //   yhat[i]=1/(1 + gsl_sf_exp(-Xbetai));
  // }
}


/* The gradient of f, df = (df/dx, df/dy). */
void 
wgsl_cat_optim_df (const Map<VectorXd>beta, void *params,
		   Map<VectorXd>out)
{
  fix_parm_cat_T *p = (fix_parm_cat_T *)params;
  int n = p->y.size();
  int K = p->X.cols();
  int npar = beta.size();
  // Intitialize gradient out necessary?
  out.setZero();
  auto np1 =	beta.size()-1;
  /* Changed loop start at 1 instead of 0 to avoid regularization of beta 0 */
  out.tail(np1)=p->lambdaL2*beta.tail(np1);

  for(int i = 1; i < npar; ++i)
    out[i]+= p->lambdaL1*((beta[i]>0)-(beta[i]<0));
  
  double pn = (-(p->y.array()-1/(1+ (-p->X.transpose()*beta).array().exp()))).sum();
  out[0]+=pn;
 
  for (int k=0; k<(p->X.outerSize()); ++k)
  for (SparseMatrix<double>::InnerIterator it(p->X,k); it; ++it)
  {
    out[1+it.col()]+=pn;
  }
  
}


/* The Hessian of f */
void 
wgsl_cat_optim_hessian (const Map<VectorXd>beta, void *params,
			Map<MatrixXd> out)
{
  fix_parm_cat_T *p = (fix_parm_cat_T *)params;
  int n = p->y.size(); 
  int K = p->X.cols();
  int npar = beta.size();
  // Intitialize Hessian out necessary ???
  out.setZero();
  
  /* Changed loop start at 1 instead of 0 to avoid regularization of beta 0*/
  for(int i = 1; i < npar; ++i)
    out(i,i)=p->lambdaL2;
  // L1 penalty not working yet, as not differentiable, I may need to do coordinate descent (as in glm_net)



  Eigen::VectorXd Xb = p->X.transpose()*beta.tail(npar-1);
  for(int i = 0; i < n; ++i) {
    double pn=0;
    double aux=0;
    double Xbetai=beta[0];
    int iParm2=1;
    int iParm1=1;



    for(int k = 0; k < K; ++k) {
      if(gsl_matrix_int_get(p->X,i,k)>0)
	Xbetai+=beta[gsl_matrix_int_get(p->X,i,k)-1+iParm1];
      iParm1+=p->nlev[k]-1;  //-1?
    }
    //    total += y[i]*Xbetai-log(1+gsl_sf_exp(Xbetai));
    pn= 1/(1 + gsl_sf_exp(-Xbetai));
    // Add a protection for pn very close to 0 or 1?
    aux=pn*(1-pn);
    *gsl_matrix_ptr(out,0,0)+=aux;
    iParm2=1;
    for(int k2 = 0; k2 < K; ++k2) {
      if(gsl_matrix_int_get(p->X,i,k2)>0)
	*gsl_matrix_ptr(out,0,gsl_matrix_int_get(p->X,i,k2)-1+iParm2)+=aux;
      iParm2+=p->nlev[k2]-1;   //-1?
    }
    iParm1=1;
    for(int k1 = 0; k1 < K; ++k1) {
      if(gsl_matrix_int_get(p->X,i,k1)>0)
	*gsl_matrix_ptr(out,gsl_matrix_int_get(p->X,i,k1)-1+iParm1,0)+=aux;
      iParm2=1;
      for(int k2 = 0; k2 < K; ++k2) {
	if((gsl_matrix_int_get(p->X,i,k1)>0) && (gsl_matrix_int_get(p->X,i,k2)>0))
	  *gsl_matrix_ptr(out
			  ,gsl_matrix_int_get(p->X,i,k1)-1+iParm1
			  ,gsl_matrix_int_get(p->X,i,k2)-1+iParm2
			  )+=aux;
	iParm2+=p->nlev[k2]-1;  //-1?
      }
      iParm1+=p->nlev[k1]-1; //-1?
    }
  }
}


double wgsl_cat_optim_f(Map<VectorXd>v, void *params)
{
  double mLogLik=0;
  fix_parm_cat_T *p = (fix_parm_cat_T *)params;
  mLogLik = fLogit_cat(v,p->X,p->y,p->lambdaL1,p->lambdaL2);
  return mLogLik; 
}


/* Compute both f and df together. */
void 
wgsl_cat_optim_fdf (Map<VectorXd>x, void *params,
		    double *f, Map<VectorXd>df)
{
  *f = wgsl_cat_optim_f(x, params); 
  wgsl_cat_optim_df(x, params, df);
}


int logistic_cat_fit(Map<VectorXd>beta
		     ,SparseMatrix<double> X
		     ,Map<VectorXd>y
		     ,double lambdaL1
		     ,double lambdaL2)
{
  double mLogLik=0;
  fix_parm_cat_T p;
  int npar = beta.size();
  int iter=0;
  double maxchange=0;

  //Intializing fix parameters
  p.X=X;

  p.y=y;
  p.lambdaL1=lambdaL1;
  p.lambdaL2=lambdaL2;
  
  //Initial fit
  //#ifdef _RPR_DEBUG_
  mLogLik = wgsl_cat_optim_f(beta,&p);
  //fprintf(stderr,"#Initial -log(Lik(0))=%lf\n",mLogLik);
  //#endif //_RPR_DEBUG

  gsl_matrix *myH = gsl_matrix_alloc(npar,npar); /* Hessian matrix*/
  VectorXd stBeta(npar); /* Direction to move */

  VectorXd myG(npar); /* Gradient*/
  VectorXd tau(npar); /* tau for QR*/

  for(iter=0;iter<100;iter++){ 
    wgsl_cat_optim_hessian(beta,&p,myH); //Calculate Hessian
    wgsl_cat_optim_df(beta,&p,myG);      //Calculate Gradient
    gsl_linalg_QR_decomp(myH,tau);   //Calculate next beta
    gsl_linalg_QR_solve(myH,tau,myG,stBeta);
    gsl_vector_sub(beta,stBeta);
    
    //Monitor convergence
    maxchange=0;
    for(int i=0;i<npar; i++)
      if(maxchange<fabs(stBeta[i]))
	maxchange=fabs(stBeta[i]);

#ifdef _RPR_DEBUG_
    mLogLik = wgsl_cat_optim_f(beta,&p);
    //fprintf(stderr,"#iter %d, -log(Lik(0))=%lf,%lf\n",(int)iter,mLogLik,maxchange);
#endif //_RPR_DEBUG

    if(maxchange<1E-4)
      break;
  }

#ifdef _RPR_DEBUG_
  //for (int i = 0; i < npar; i++)
  //  fprintf(stderr,"#par_%d= %lf\n",i,beta[i]);
#endif //_RPR_DEBUG

  //Final fit
  mLogLik = wgsl_cat_optim_f(beta,&p);
  //fprintf(stderr,"#Final %d) -log(Lik(0))=%lf, maxchange %lf\n",iter,mLogLik,maxchange);

  // gsl_vector_free (tau);
  // gsl_vector_free (stBeta);
  // gsl_vector_free (myG);
  // gsl_matrix_free (myH);

  return 0;
}



