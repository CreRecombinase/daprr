#include <stdio.h>
#include <math.h>
#include <RcppGSL.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_linalg.h>

#include "logistic.hpp"
#include <RcppEigen.h>
// I need to bundle all the data that goes to the function to optimze together. 
typedef struct{
  gsl_matrix_int *X;
  gsl_vector_int *nlev;
  gsl_vector *y;
  gsl_matrix *Xc;   // continuous covariates  Matrix Nobs x Kc (NULL if not used)
  double lambdaL1;
  double lambdaL2;
}fix_parm_mixed_T;


double fLogit_mixed(gsl_vector *beta
		    ,gsl_matrix_int *X
		    ,gsl_vector_int *nlev
		    ,gsl_matrix *Xc
		    ,gsl_vector *y
		    ,double lambdaL1
		    ,double lambdaL2)
{
  int n = y->size; 
  //  int k = X->size2; 
  int npar = beta->size; 
  double total = 0;
  double aux = 0;
  /*   omp_set_num_threads(ompthr); */
  /*   /\* Changed loop start at 1 instead of 0 to avoid regularization of beta_0*\/ */
  /*   /\*#pragma omp parallel for reduction (+:total)*\/ */
  for(int i = 1; i < npar; ++i)
    total += beta->data[i]*beta->data[i];
  total = (-total*lambdaL2/2);
  /*   /\*#pragma omp parallel for reduction (+:aux)*\/ */
  for(int i = 1; i < npar; ++i)
    aux += (beta->data[i]>0 ? beta->data[i] : -beta->data[i]);
  total = total-aux*lambdaL1;
  /* #pragma omp parallel for schedule(static) shared(n,beta,X,nlev,y) reduction (+:total) */
  for(int i = 0; i < n; ++i) {
    double Xbetai=beta->data[0];
    int iParm=1;
    for(int k = 0; k < X->size2; ++k) {
      if(gsl_matrix_int_get(X,i,k)>0)
	Xbetai+=beta->data[gsl_matrix_int_get(X,i,k)-1+iParm];
      iParm+=nlev->data[k]-1;
    }
    for(int k = 0; k < (Xc->size2); ++k) 
      Xbetai+= gsl_matrix_get(Xc,i,k)*beta->data[iParm++];
    total += y->data[i]*Xbetai-gsl_sf_log_1plusx(gsl_sf_exp(Xbetai));
  }
  return -total;
} 


void logistic_mixed_pred(gsl_vector *beta  // Vector of parameters length = 1 + Sum_k(C_k - 1)
			 ,gsl_matrix_int *X  //Matrix Nobs x K 
			 ,gsl_vector_int *nlev // Vector with number categories
			 ,gsl_matrix *Xc   // continuous covariates  Matrix Nobs x Kc (NULL if not used)
			 ,gsl_vector *yhat //Vector of prob. predicted by the logistic
			 )
{
  for(int i = 0; i < X->size1; ++i) {
    double Xbetai=beta->data[0];
    int iParm=1;
    for(int k = 0; k < X->size2; ++k) {
      if(gsl_matrix_int_get(X,i,k)>0)
	Xbetai+=beta->data[gsl_matrix_int_get(X,i,k)-1+iParm];
      iParm+=nlev->data[k]-1;
    }
    // Adding the continuous
    for(int k = 0; k < (Xc->size2); ++k) 
      Xbetai+= gsl_matrix_get(Xc,i,k)*beta->data[iParm++];
    yhat->data[i]=1/(1 + gsl_sf_exp(-Xbetai));
  }
}


/* The gradient of f, df = (df/dx, df/dy). */
void 
wgsl_mixed_optim_df (const gsl_vector *beta, void *params, 
		     gsl_vector *out)
{
  fix_parm_mixed_T *p = (fix_parm_mixed_T *)params;
  int n = p->y->size; 
  int K = p->X->size2; 
  int Kc = p->Xc->size2; 
  int npar = beta->size; 
  // Intitialize gradient out necessary?
  for(int i = 0; i < npar; ++i) 
    out->data[i]= 0; 
  /* Changed loop start at 1 instead of 0 to avoid regularization of beta 0 */
  for(int i = 1; i < npar; ++i)
    out->data[i]= p->lambdaL2*beta->data[i]; 
  for(int i = 1; i < npar; ++i)
    out->data[i]+= p->lambdaL1*((beta->data[i]>0)-(beta->data[i]<0));
  
  for(int i = 0; i < n; ++i) {
    double pn=0;
    double Xbetai=beta->data[0];
    int iParm=1;
    for(int k = 0; k < K; ++k) {
      if(gsl_matrix_int_get(p->X,i,k)>0)
	Xbetai+=beta->data[gsl_matrix_int_get(p->X,i,k)-1+iParm];
      iParm+=p->nlev->data[k]-1;
    }
    // Adding the continuous
    for(int k = 0; k < Kc; ++k) 
      Xbetai+= gsl_matrix_get(p->Xc,i,k)*beta->data[iParm++];

    pn= -( p->y->data[i] - 1/(1 + gsl_sf_exp(-Xbetai)) );

    out->data[0]+= pn;
    iParm=1;
    for(int k = 0; k < K; ++k) {
      if(gsl_matrix_int_get(p->X,i,k)>0)
	out->data[gsl_matrix_int_get(p->X,i,k)-1+iParm]+=pn;
      iParm+=p->nlev->data[k]-1;
    }
    // Adding the continuous
    for(int k = 0; k < Kc; ++k) {
      out->data[iParm++] += gsl_matrix_get(p->Xc,i,k)*pn;
    }
  }

}


/* The Hessian of f */
void 
wgsl_mixed_optim_hessian (const gsl_vector *beta, void *params, 
			  gsl_matrix *out)
{
  fix_parm_mixed_T *p = (fix_parm_mixed_T *)params;
  int n = p->y->size; 
  int K = p->X->size2; 
  int Kc = p->Xc->size2; 
  int npar = beta->size; 
  gsl_vector *gn = gsl_vector_alloc(npar); /* gn */
  // Intitialize Hessian out necessary ???
  gsl_matrix_set_zero(out);
  /* Changed loop start at 1 instead of 0 to avoid regularization of beta 0*/
  for(int i = 1; i < npar; ++i)
    gsl_matrix_set(out,i,i,(p->lambdaL2));  // Double check this
  // L1 penalty not working yet, as not differentiable, I may need to do coordinate descent (as in glm_net)
  
  for(int i = 0; i < n; ++i) {
    double pn=0;
    double aux=0;
    double Xbetai=beta->data[0];
    int iParm1=1;
    for(int k = 0; k < K; ++k) {
      if(gsl_matrix_int_get(p->X,i,k)>0)
	Xbetai+=beta->data[gsl_matrix_int_get(p->X,i,k)-1+iParm1];
      iParm1+=p->nlev->data[k]-1;  //-1?
    }
    // Adding the continuous
    for(int k = 0; k < Kc; ++k) 
      Xbetai+= gsl_matrix_get(p->Xc,i,k)*beta->data[iParm1++];

    pn= 1/(1 + gsl_sf_exp(-Xbetai));
    // Add a protection for pn very close to 0 or 1?
    aux=pn*(1-pn);

    // Calculate sub-gradient vector gn 
    gsl_vector_set_zero(gn);
    gn->data[0]= 1;
    iParm1=1;
    for(int k = 0; k < K; ++k) {
      if(gsl_matrix_int_get(p->X,i,k)>0)
	gn->data[gsl_matrix_int_get(p->X,i,k)-1+iParm1]=1;
      iParm1+=p->nlev->data[k]-1;
    }
    // Adding the continuous
    for(int k = 0; k < Kc; ++k) {
      gn->data[iParm1++] = gsl_matrix_get(p->Xc,i,k);
    }

    for(int k1=0;k1<npar; ++k1)
      if(gn->data[k1]!=0)
	for(int k2=0;k2<npar; ++k2)
	  if(gn->data[k2]!=0)
	    *gsl_matrix_ptr(out,k1,k2) += (aux * gn->data[k1] * gn->data[k2]);
  }
  gsl_vector_free(gn);
}


double wgsl_mixed_optim_f(gsl_vector *v, void *params)
{
  double mLogLik=0;
  fix_parm_mixed_T *p = (fix_parm_mixed_T *)params;
  mLogLik = fLogit_mixed(v,p->X,p->nlev,p->Xc,p->y,p->lambdaL1,p->lambdaL2);
  return mLogLik; 
}


/* Compute both f and df together. */
void 
wgsl_mixed_optim_fdf (gsl_vector *x, void *params, 
		      double *f, gsl_vector *df) 
{
  *f = wgsl_mixed_optim_f(x, params); 
  wgsl_mixed_optim_df(x, params, df);
}


int logistic_mixed_fit(gsl_vector *beta
		       ,gsl_matrix_int *X
		       ,gsl_vector_int *nlev
		       ,gsl_matrix *Xc   // continuous covariates  Matrix Nobs x Kc (NULL if not used)
		       ,gsl_vector *y
		       ,double lambdaL1
		       ,double lambdaL2)
{

  double mLogLik=0;
  fix_parm_mixed_T p;
  int npar = beta->size; 
  int iter=0;
  double maxchange=0;

  //Intializing fix parameters
  p.X=X;
  p.Xc=Xc;
  p.nlev=nlev;
  p.y=y;
  p.lambdaL1=lambdaL1;
  p.lambdaL2=lambdaL2;
  
  //Initial fit

  mLogLik = wgsl_mixed_optim_f(beta,&p);
#ifdef _RPR_DEBUG_
  fprintf(stderr,"#Initial -log(Lik(0))=%lf\n",mLogLik);
#endif //_RPR_DEBUG

  gsl_matrix *myH = gsl_matrix_alloc(npar,npar); /* Hessian matrix*/
  gsl_vector *stBeta = gsl_vector_alloc(npar); /* Direction to move */

  gsl_vector *myG = gsl_vector_alloc(npar); /* Gradient*/
  gsl_vector *tau = gsl_vector_alloc(npar); /* tau for QR*/

  for(iter=0;iter<100;iter++){ 
    wgsl_mixed_optim_hessian(beta,&p,myH); //Calculate Hessian
    wgsl_mixed_optim_df(beta,&p,myG);      //Calculate Gradient
    gsl_linalg_QR_decomp(myH,tau);   //Calculate next beta
    gsl_linalg_QR_solve(myH,tau,myG,stBeta);
    gsl_vector_sub(beta,stBeta);
    
    //Monitor convergence
    maxchange=0;
    for(int i=0;i<npar; i++)
      if(maxchange<fabs(stBeta->data[i]))
	maxchange=fabs(stBeta->data[i]);

#ifdef _RPR_DEBUG_
    //    mLogLik = wgsl_mixed_optim_f(beta,&p);
    fprintf(stderr,"#iter %d, -log(Lik(0))=%lf,%lf\n",(int)iter,mLogLik,maxchange);
#endif //_RPR_DEBUG

    if(maxchange<1E-4)
      break;
  }

#ifdef _RPR_DEBUG_
  for (int i = 0; i < npar; i++)
    fprintf(stderr,"#par_%d= %lf\n",i,beta->data[i]);
#endif //_RPR_DEBUG

  //Final fit
  mLogLik = wgsl_mixed_optim_f(beta,&p);
  //fprintf(stderr,"#Final %d) -log(Lik(0))=%lf, maxchange %lf\n",iter,mLogLik,maxchange);

  gsl_vector_free (tau);
  gsl_vector_free (stBeta);
  gsl_vector_free (myG);
  gsl_matrix_free (myH);

  return 0;
}

/***************/
/* Categorical */
/***************/


double fLogit_cat(gsl_vector *beta
		  ,gsl_matrix_int *X
		  ,gsl_vector_int *nlev
		  ,gsl_vector *y
		  ,double lambdaL1
		  ,double lambdaL2)
{
  int n = y->size; 
  //  int k = X->size2; 
  int npar = beta->size; 
  double total = 0;
  double aux = 0;
  /*   omp_set_num_threads(ompthr); */
  /*   /\* Changed loop start at 1 instead of 0 to avoid regularization of beta 0*\/ */
  /*   /\*#pragma omp parallel for reduction (+:total)*\/ */
  for(int i = 1; i < npar; ++i)
    total += beta->data[i]*beta->data[i];
  total = (-total*lambdaL2/2);
  /*   /\*#pragma omp parallel for reduction (+:aux)*\/ */
  for(int i = 1; i < npar; ++i)
    aux += (beta->data[i]>0 ? beta->data[i] : -beta->data[i]);
  total = total-aux*lambdaL1;
  /* #pragma omp parallel for schedule(static) shared(n,beta,X,nlev,y) reduction (+:total) */
  for(int i = 0; i < n; ++i) {
    double Xbetai=beta->data[0];
    int iParm=1;
    for(int k = 0; k < X->size2; ++k) {
      if(gsl_matrix_int_get(X,i,k)>0)
	Xbetai+=beta->data[gsl_matrix_int_get(X,i,k)-1+iParm];
      iParm+=nlev->data[k]-1;
    }
    total += y->data[i]*Xbetai-gsl_sf_log_1plusx(gsl_sf_exp(Xbetai));
  }
  return -total;
} 


void logistic_cat_pred(gsl_vector *beta  // Vector of parameters length = 1 + Sum_k(C_k - 1)
		       ,gsl_matrix_int *X  //Matrix Nobs x K 
		       ,gsl_vector_int *nlev // Vector with number categories
		       ,gsl_vector *yhat //Vector of prob. predicted by the logistic
		       )
{
  for(int i = 0; i < X->size1; ++i) {
    double Xbetai=beta->data[0];
    int iParm=1;
    for(int k = 0; k < X->size2; ++k) {
      if(gsl_matrix_int_get(X,i,k)>0)
	Xbetai+=beta->data[gsl_matrix_int_get(X,i,k)-1+iParm];
      iParm+=nlev->data[k]-1;
    }
    yhat->data[i]=1/(1 + gsl_sf_exp(-Xbetai));
  }
}


/* The gradient of f, df = (df/dx, df/dy). */
void 
wgsl_cat_optim_df (const gsl_vector *beta, void *params, 
		   gsl_vector *out)
{
  fix_parm_cat_T *p = (fix_parm_cat_T *)params;
  int n = p->y->size; 
  int K = p->X->size2; 
  int npar = beta->size; 
  // Intitialize gradient out necessary?
  for(int i = 0; i < npar; ++i) 
    out->data[i]= 0; 
  /* Changed loop start at 1 instead of 0 to avoid regularization of beta 0 */
  for(int i = 1; i < npar; ++i)
    out->data[i]= p->lambdaL2*beta->data[i]; 
  for(int i = 1; i < npar; ++i)
    out->data[i]+= p->lambdaL1*((beta->data[i]>0)-(beta->data[i]<0));
  
  for(int i = 0; i < n; ++i) {
    double pn=0;
    double Xbetai=beta->data[0];
    int iParm=1;
    for(int k = 0; k < K; ++k) {
      if(gsl_matrix_int_get(p->X,i,k)>0)
	Xbetai+=beta->data[gsl_matrix_int_get(p->X,i,k)-1+iParm];
      iParm+=p->nlev->data[k]-1;
    }
    //    total += y->data[i]*Xbetai-log(1+gsl_sf_exp(Xbetai));
    pn= -( p->y->data[i] - 1/(1 + gsl_sf_exp(-Xbetai)) );

    out->data[0]+= pn;
    iParm=1;
    for(int k = 0; k < K; ++k) {
      if(gsl_matrix_int_get(p->X,i,k)>0)
	out->data[gsl_matrix_int_get(p->X,i,k)-1+iParm]+=pn;
      iParm+=p->nlev->data[k]-1;
    }
  }
}


/* The Hessian of f */
void 
wgsl_cat_optim_hessian (const gsl_vector *beta, void *params, 
			gsl_matrix *out)
{
  fix_parm_cat_T *p = (fix_parm_cat_T *)params;
  int n = p->y->size; 
  int K = p->X->size2; 
  int npar = beta->size; 
  // Intitialize Hessian out necessary ???
  gsl_matrix_set_zero(out);
  /* Changed loop start at 1 instead of 0 to avoid regularization of beta 0*/
  for(int i = 1; i < npar; ++i)
    gsl_matrix_set(out,i,i,(p->lambdaL2));  // Double check this
  // L1 penalty not working yet, as not differentiable, I may need to do coordinate descent (as in glm_net)

  
  for(int i = 0; i < n; ++i) {
    double pn=0;
    double aux=0;
    double Xbetai=beta->data[0];
    int iParm2=1;
    int iParm1=1;
    for(int k = 0; k < K; ++k) {
      if(gsl_matrix_int_get(p->X,i,k)>0)
	Xbetai+=beta->data[gsl_matrix_int_get(p->X,i,k)-1+iParm1];
      iParm1+=p->nlev->data[k]-1;  //-1?
    }
    //    total += y->data[i]*Xbetai-log(1+gsl_sf_exp(Xbetai));
    pn= 1/(1 + gsl_sf_exp(-Xbetai));
    // Add a protection for pn very close to 0 or 1?
    aux=pn*(1-pn);
    *gsl_matrix_ptr(out,0,0)+=aux;
    iParm2=1;
    for(int k2 = 0; k2 < K; ++k2) {
      if(gsl_matrix_int_get(p->X,i,k2)>0)
	*gsl_matrix_ptr(out,0,gsl_matrix_int_get(p->X,i,k2)-1+iParm2)+=aux;
      iParm2+=p->nlev->data[k2]-1;   //-1?
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
	iParm2+=p->nlev->data[k2]-1;  //-1?
      }
      iParm1+=p->nlev->data[k1]-1; //-1?
    }
  }
}


double wgsl_cat_optim_f(gsl_vector *v, void *params)
{
  double mLogLik=0;
  fix_parm_cat_T *p = (fix_parm_cat_T *)params;
  mLogLik = fLogit_cat(v,p->X,p->nlev,p->y,p->lambdaL1,p->lambdaL2);
  return mLogLik; 
}


/* Compute both f and df together. */
void 
wgsl_cat_optim_fdf (gsl_vector *x, void *params, 
		    double *f, gsl_vector *df) 
{
  *f = wgsl_cat_optim_f(x, params); 
  wgsl_cat_optim_df(x, params, df);
}


int logistic_cat_fit(gsl_vector *beta
		     ,gsl_matrix_int *X
		     ,gsl_vector_int *nlev
		     ,gsl_vector *y
		     ,double lambdaL1
		     ,double lambdaL2)
{
  double mLogLik=0;
  fix_parm_cat_T p;
  int npar = beta->size;
  int iter=0;
  double maxchange=0;

  //Intializing fix parameters
  p.X=X;
  p.nlev=nlev;
  p.y=y;
  p.lambdaL1=lambdaL1;
  p.lambdaL2=lambdaL2;

  //Initial fit

  mLogLik = wgsl_cat_optim_f(beta,&p);
#ifdef _RPR_DEBUG_
  fprintf(stderr,"#Initial -log(Lik(0))=%lf\n",mLogLik);
#endif //_RPR_DEBUG

  gsl_matrix *myH = gsl_matrix_alloc(npar,npar); /* Hessian matrix*/
  gsl_vector *stBeta = gsl_vector_alloc(npar); /* Direction to move */

  gsl_vector *myG = gsl_vector_alloc(npar); /* Gradient*/
  gsl_vector *tau = gsl_vector_alloc(npar); /* tau for QR*/

  for(iter=0;iter<100;iter++){
    wgsl_cat_optim_hessian(beta,&p,myH); //Calculate Hessian
    wgsl_cat_optim_df(beta,&p,myG);      //Calculate Gradient
    gsl_linalg_QR_decomp(myH,tau);   //Calculate next beta
    gsl_linalg_QR_solve(myH,tau,myG,stBeta);
    gsl_vector_sub(beta,stBeta);

    //Monitor convergence
    maxchange=0;
    for(int i=0;i<npar; i++)
      if(maxchange<fabs(stBeta->data[i]))
	maxchange=fabs(stBeta->data[i]);

// #ifdef _RPR_DEBUG_
//     mLogLik = wgsl_cat_optim_f(beta,&p);
//     fprintf(stderr,"#iter %d, -log(Lik(0))=%lf,%lf\n",(int)iter,mLogLik,maxchange);
// #endif //_RPR_DEBUG

    if(maxchange<1E-4)
      break;
  }

#ifdef _RPR_DEBUG_
  for (int i = 0; i < npar; i++)
   fprintf(stderr,"#par_%d= %lf\n",i,beta->data[i]);
#endif //_RPR_DEBUG

  //Final fit
  mLogLik = wgsl_cat_optim_f(beta,&p);
  //fprintf(stderr,"#Final %d) -log(Lik(0))=%lf, maxchange %lf\n",iter,mLogLik,maxchange);

  gsl_vector_free (tau);
  gsl_vector_free (stBeta);
  gsl_vector_free (myG);
  gsl_matrix_free (myH);

  return 0;
}


/***************/
/* Continuous  */
/***************/

// I need to bundle all the data that goes to the function to optimze together. 
typedef struct{
  gsl_matrix *Xc;   // continuous covariates  Matrix Nobs x Kc 
  gsl_vector *y;
  double lambdaL1;
  double lambdaL2;
}fix_parm_cont_T;


double fLogit_cont(gsl_vector *beta
		   ,gsl_matrix *Xc
		   ,gsl_vector *y
		   ,double lambdaL1
		   ,double lambdaL2)
{
  int n = y->size; 
  int npar = beta->size; 
  double total = 0;
  double aux = 0;
  /*   omp_set_num_threads(ompthr); */
  /*   /\* Changed loop start at 1 instead of 0 to avoid regularization of beta_0*\/ */
  /*   /\*#pragma omp parallel for reduction (+:total)*\/ */
  for(int i = 1; i < npar; ++i)
    total += beta->data[i]*beta->data[i];
  total = (-total*lambdaL2/2);
  /*   /\*#pragma omp parallel for reduction (+:aux)*\/ */
  for(int i = 1; i < npar; ++i)
    aux += (beta->data[i]>0 ? beta->data[i] : -beta->data[i]);
  total = total-aux*lambdaL1;
  /* #pragma omp parallel for schedule(static) shared(n,beta,X,nlev,y) reduction (+:total) */
  for(int i = 0; i < n; ++i) {
    double Xbetai=beta->data[0];
    int iParm=1;
    for(int k = 0; k < (Xc->size2); ++k) 
      Xbetai+= gsl_matrix_get(Xc,i,k)*beta->data[iParm++];
    total += y->data[i]*Xbetai-gsl_sf_log_1plusx(gsl_sf_exp(Xbetai));
  }
  return -total;
} 


void logistic_cont_pred(gsl_vector *beta  // Vector of parameters length = 1 + Sum_k(C_k - 1)
			,gsl_matrix *Xc   // continuous covariates  Matrix Nobs x Kc (NULL if not used)
			,gsl_vector *yhat //Vector of prob. predicted by the logistic
			)
{
  for(int i = 0; i < Xc->size1; ++i) {
    double Xbetai=beta->data[0];
    int iParm=1;
    for(int k = 0; k < (Xc->size2); ++k) 
      Xbetai+= gsl_matrix_get(Xc,i,k)*beta->data[iParm++];
    yhat->data[i]=1/(1 + gsl_sf_exp(-Xbetai));
  }
}


/* The gradient of f, df = (df/dx, df/dy). */
void 
wgsl_cont_optim_df (const gsl_vector *beta, void *params, 
		    gsl_vector *out)
{
  fix_parm_cont_T *p = (fix_parm_cont_T *)params;
  int n = p->y->size; 
  int Kc = p->Xc->size2; 
  int npar = beta->size; 
  // Intitialize gradient out necessary?
  for(int i = 0; i < npar; ++i) 
    out->data[i]= 0; 
  /* Changed loop start at 1 instead of 0 to avoid regularization of beta 0 */
  for(int i = 1; i < npar; ++i)
    out->data[i]= p->lambdaL2*beta->data[i]; 
  for(int i = 1; i < npar; ++i)
    out->data[i]+= p->lambdaL1*((beta->data[i]>0)-(beta->data[i]<0));
  
  for(int i = 0; i < n; ++i) {
    double pn=0;
    double Xbetai=beta->data[0];
    int iParm=1;
    for(int k = 0; k < Kc; ++k) 
      Xbetai+= gsl_matrix_get(p->Xc,i,k)*beta->data[iParm++];

    pn= -( p->y->data[i] - 1/(1 + gsl_sf_exp(-Xbetai)) );

    out->data[0]+= pn;
    iParm=1;
    // Adding the continuous
    for(int k = 0; k < Kc; ++k) {
      out->data[iParm++] += gsl_matrix_get(p->Xc,i,k)*pn;
    }
  }
}


/* The Hessian of f */
void 
wgsl_cont_optim_hessian (const gsl_vector *beta, void *params, 
			 gsl_matrix *out)
{
  fix_parm_cont_T *p = (fix_parm_cont_T *)params;
  int n = p->y->size; 
  int Kc = p->Xc->size2; 
  int npar = beta->size; 
  gsl_vector *gn = gsl_vector_alloc(npar); /* gn */
  // Intitialize Hessian out necessary ???
  gsl_matrix_set_zero(out);
  /* Changed loop start at 1 instead of 0 to avoid regularization of beta 0*/
  for(int i = 1; i < npar; ++i)
    gsl_matrix_set(out,i,i,(p->lambdaL2));  // Double check this
  // L1 penalty not working yet, as not differentiable, I may need to do coordinate descent (as in glm_net)
  
  for(int i = 0; i < n; ++i) {
    double pn=0;
    double aux=0;
    double Xbetai=beta->data[0];
    int iParm1=1;
    for(int k = 0; k < Kc; ++k) 
      Xbetai+= gsl_matrix_get(p->Xc,i,k)*beta->data[iParm1++];

    pn= 1/(1 + gsl_sf_exp(-Xbetai));
    // Add a protection for pn very close to 0 or 1?
    aux=pn*(1-pn);

    // Calculate sub-gradient vector gn 
    gsl_vector_set_zero(gn);
    gn->data[0]= 1;
    iParm1=1;
    for(int k = 0; k < Kc; ++k) {
      gn->data[iParm1++] = gsl_matrix_get(p->Xc,i,k);
    }

    for(int k1=0;k1<npar; ++k1)
      if(gn->data[k1]!=0)
	for(int k2=0;k2<npar; ++k2)
	  if(gn->data[k2]!=0)
	    *gsl_matrix_ptr(out,k1,k2) += (aux * gn->data[k1] * gn->data[k2]);
  }
  gsl_vector_free(gn);
}


double wgsl_cont_optim_f(gsl_vector *v, void *params)
{
  double mLogLik=0;
  fix_parm_cont_T *p = (fix_parm_cont_T *)params;
  mLogLik = fLogit_cont(v,p->Xc,p->y,p->lambdaL1,p->lambdaL2);
  return mLogLik; 
}


/* Compute both f and df together. */
void 
wgsl_cont_optim_fdf (gsl_vector *x, void *params, 
		     double *f, gsl_vector *df) 
{
  *f = wgsl_cont_optim_f(x, params); 
  wgsl_cont_optim_df(x, params, df);
}


int logistic_cont_fit(gsl_vector *beta
		      ,gsl_matrix *Xc   // continuous covariates  Matrix Nobs x Kc (NULL if not used)
		      ,gsl_vector *y
		      ,double lambdaL1
		      ,double lambdaL2)
{

  double mLogLik=0;
  fix_parm_cont_T p;
  int npar = beta->size; 
  int iter=0;
  double maxchange=0;

  //Intializing fix parameters
  p.Xc=Xc;
  p.y=y;
  p.lambdaL1=lambdaL1;
  p.lambdaL2=lambdaL2;
  
  //Initial fit

  mLogLik = wgsl_cont_optim_f(beta,&p);
  #ifdef _RPR_DEBUG_
  fprintf(stderr,"#Initial -log(Lik(0))=%lf\n",mLogLik);
  #endif //_RPR_DEBUG

  gsl_matrix *myH = gsl_matrix_alloc(npar,npar); /* Hessian matrix*/
  gsl_vector *stBeta = gsl_vector_alloc(npar); /* Direction to move */

  gsl_vector *myG = gsl_vector_alloc(npar); /* Gradient*/
  gsl_vector *tau = gsl_vector_alloc(npar); /* tau for QR*/

  for(iter=0;iter<100;iter++){ 
    wgsl_cont_optim_hessian(beta,&p,myH); //Calculate Hessian
    wgsl_cont_optim_df(beta,&p,myG);      //Calculate Gradient
    gsl_linalg_QR_decomp(myH,tau);   //Calculate next beta
    gsl_linalg_QR_solve(myH,tau,myG,stBeta);
    gsl_vector_sub(beta,stBeta);
    
    //Monitor convergence
    maxchange=0;
    for(int i=0;i<npar; i++)
      if(maxchange<fabs(stBeta->data[i]))
	maxchange=fabs(stBeta->data[i]);

#ifdef _RPR_DEBUG_
    mLogLik = wgsl_cont_optim_f(beta,&p);
    //fprintf(stderr,"#iter %d, -log(Lik(0))=%lf,%lf\n",(int)iter,mLogLik,maxchange);
#endif //_RPR_DEBUG

    if(maxchange<1E-4)
      break;
  }

#ifdef _RPR_DEBUG_
  for (int i = 0; i < npar; i++)
   fprintf(stderr,"#par_%d= %lf\n",i,beta->data[i]);
#endif //_RPR_DEBUG

  //Final fit
  mLogLik = wgsl_cont_optim_f(beta,&p);
  //fprintf(stderr,"#Final %d) -log(Lik(0))=%lf, maxchange %lf\n",iter,mLogLik,maxchange);

  gsl_vector_free (tau);
  gsl_vector_free (stBeta);
  gsl_vector_free (myG);
  gsl_matrix_free (myH);

  return 0;
}

Logistic::Logistic(const size_t npar_,const size_t p_):npar(npar_),
						       p(p_),
						       ty(p,2),
						       myH(gsl_matrix_alloc(npar,npar)),
						       stBeta(gsl_vector_alloc(npar)),
						       myG(gsl_vector_alloc(npar)),
						       tau(gsl_vector_alloc(npar))
{
}

void Logistic::fit_qbinom(RcppGSL::vector<double> beta ,RcppGSL::matrix<int> X,RcppGSL::vector<double> y){

  using namespace Rcpp;
  Environment pkg = Environment::namespace_env("stats");
  Function qb =  pkg["quasibinomial"];

  auto qbf = qb("logit");
    // Picking up Matrix() function from Matrix package
  Function f = pkg["glm.fit"];

  auto ctrlL = List::create(_["epsilon"]=NumericVector::create(1e-08),
			    _["maxit"]=NumericVector::create(25),
			    _["trace"]=LogicalVector::create(true));
  List ret_l = f( X,ty, Named("family", qbf),Named("control",ctrlL),Named("start",beta));

  Rcpp::NumericVector m = (ret_l["coefficients"]);
  for(int i=0; i < m.size(); i++){
    beta[i]=m(i);
  }
}

void Logistic::fit_glmnet(RcppGSL::vector<double> beta ,RcppGSL::matrix<int> X,RcppGSL::vector<double> y,double alpha,double lambda){

  for(int i=0; i<p; i++){
    ty(i,0)=1-y[i];
    ty(i,1)=y[i];
  }

  Rcpp::Environment pkg = Rcpp::Environment::namespace_env("glmnet");

    // Picking up Matrix() function from Matrix package
  Rcpp::Function f = pkg["glmnet"];
  Rcpp::List ret_l = f( X,ty, Rcpp::Named("family", "binomial"),Rcpp::Named("alpha",alpha),Rcpp::Named("lambda",lambda),Rcpp::Named("type.logistic","modified.Newton"));
  auto m = Rcpp::as<Eigen::Map<Eigen::SparseMatrix<double> > >(ret_l["beta"]);
  beta[0] = Rcpp::as<Rcpp::NumericVector>(ret_l["a0"])[0];
  for(int i=0; i < m.size(); i++){
    beta[i+1]=m.coeff(i,0);
  }
  // return m;
}





  //  Rcpp::NumericMatrix ty(
void Logistic::fit(gsl_vector *beta,gsl_matrix_int *X,gsl_vector_int *nlev,gsl_vector *y,double lambdaL1,double lambdaL2){



  double mLogLik=0;
  int iter=0;
  double maxchange=0;

  fix_parm_cat_T p;
  int npar = beta->size;

  //Intializing fix parameters
  p.X=X;
  p.nlev=nlev;
  p.y=y;
  p.lambdaL1=lambdaL1;
  p.lambdaL2=lambdaL2;

  //Initial fit
  mLogLik = wgsl_cat_optim_f(beta,&p);
#ifdef _RPR_DEBUG_
  fprintf(stderr,"#Initial -log(Lik(0))=%lf\n",mLogLik);
#endif //_RPR_DEBUG

  for (iter = 0; iter < 100; iter++) {
    wgsl_cat_optim_hessian(beta, &p, myH); // Calculate Hessian
    wgsl_cat_optim_df(beta, &p, myG);      // Calculate Gradient
    gsl_linalg_QR_decomp(myH, tau);        // Calculate next beta
    gsl_linalg_QR_solve(myH, tau, myG, stBeta);
    gsl_vector_sub(beta, stBeta);

    // Monitor convergence
    maxchange = 0;
    for (int i = 0; i < npar; i++)
      if (maxchange < fabs(stBeta->data[i]))
        maxchange = fabs(stBeta->data[i]);
    if (maxchange < 1E-4)
      break;
  }

#ifdef _RPR_DEBUG_
  for (int i = 0; i < npar; i++)
    fprintf(stderr, "#par_%d= %lf\n", i, beta->data[i]);
#endif //_RPR_DEBUG

  // Final fit
  mLogLik = wgsl_cat_optim_f(beta, &p);
  // fprintf(stderr,"#Final %d) -log(Lik(0))=%lf, maxchange
  // %lf\n",iter,mLogLik,maxchange);

  // return 0;
}
void Dap_logit::fit(gsl::span<double> data) {
  const size_t n = data.size();
  gsl_vector retv;
  gsl_vector *y = &retv;
  y->size = n;
  y->stride = 1;
  y->data = data.data();
  y->block = nullptr;
  y->owner = 0;

  double mLogLik=0;
  int iter=0;
  double maxchange=0;

  fix_parm_cat_T p;
  int npar = beta->size;

  //Intializing fix parameters
  p.X=X;
  p.nlev=nlev;
  p.y=y;
  p.lambdaL1=lambdaL1;
  p.lambdaL2=lambdaL2;

  //Initial fit
  mLogLik = wgsl_cat_optim_f(beta,&p);
#ifdef _RPR_DEBUG_
  fprintf(stderr,"#Initial -log(Lik(0))=%lf\n",mLogLik);
#endif //_RPR_DEBUG

  for (iter = 0; iter < 100; iter++) {
    wgsl_cat_optim_hessian(beta, &p, myH); // Calculate Hessian
    wgsl_cat_optim_df(beta, &p, myG);      // Calculate Gradient
    gsl_linalg_QR_decomp(myH, tau);        // Calculate next beta
    gsl_linalg_QR_solve(myH, tau, myG, stBeta);
    gsl_vector_sub(beta, stBeta);

    // Monitor convergence
    maxchange = 0;
    for (int i = 0; i < npar; i++)
      if (maxchange < fabs(stBeta->data[i]))
        maxchange = fabs(stBeta->data[i]);
    if (maxchange < 1E-4)
      break;
  }

#ifdef _RPR_DEBUG_
  for (int i = 0; i < npar; i++)
    fprintf(stderr, "#par_%d= %lf\n", i, beta->data[i]);
#endif //_RPR_DEBUG

  // Final fit
  mLogLik = wgsl_cat_optim_f(beta, &p);
  // fprintf(stderr,"#Final %d) -log(Lik(0))=%lf, maxchange
  // %lf\n",iter,mLogLik,maxchange);
}


void Dap_logit::read_coeffs(gsl::span<double> obeta) {
  gsl::span<double> mbeta(beta->data,beta->size);
  if (obeta.size() != mbeta.size()) {
    Rcpp::stop("obeta.size()!=mbeta.size() in Dap_logit::read_coeffs");
  }
  std::copy(mbeta.begin(), mbeta.end(), obeta.begin());
}

void Dap_logit::predict(gsl::span<double> beta, gsl::span<double> syhat){

  const size_t ny = syhat.size();
  const size_t nb = beta.size();
  gsl_vector sy,sb;
  gsl_vector *yhat = &sy;
  gsl_vector *obeta = &sb;

  yhat->size = ny;
  obeta->size = nb;

  yhat->stride = 1;
  obeta->stride = 1;

  yhat->data = syhat.data();
  obeta->data = beta.data();

  yhat->block = nullptr;
  obeta->block = nullptr;

  yhat->owner = 0;
  obeta->owner = 0;


  logistic_cat_pred(obeta,X,nlev,yhat);


}
