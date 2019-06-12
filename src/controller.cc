
#include "classdef.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h>
#include <math.h>


//#include <omp.h>
#include <boost/config/warning_disable.hpp>
#include "logistic.h"
#include <wordexp.h>
#include "parser.hpp"
#include <sys/types.h>
#include <errno.h>

class BF{
  const double wt;
  const double l10;
  std::array<double,4> kv;
  std::array<double,4> kva;
  std::array<double,4> kvb;
  mutable std::array<double,4> kvc;
public:
  BF():wt(1/4.0),l10(log(10)){
    kv={1,4,16,25};
    for(int i=0; i<4; i++){
      kva[i]= 0.5*log(1/(1+kv[i]));
      kvb[i]=0.5*(kv[i]/(1+kv[i]));
    }
  }
  double compute_log10_BF(const double z_score) const{
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
};










void controller::load_data_R(Rcpp::NumericVector z_hat){

  std::vector<SNP> snp_vec;


  p=static_cast<int>(z_hat.size());
  log10_BF.resize(p);



 auto t_z = Rcpp::as<Eigen::Map<Eigen::ArrayXd>>(z_hat);


  //  assert_size(anno_mat.col(0),"anno_mat col 1",p);
  BF bf;

  log10_BF=t_z.unaryExpr([&bf](double zh){
							      return(bf.compute_log10_BF(zh));
							    });

//  std::transform(z_hat.begin(),z_hat.end(),log10_BF.begin(),

  const auto num_reg = split.num_regions();
  auto BF_d = log10_BF.data();
  auto prior_v_d = prior_vec.data();
  auto pip_v_d = pip_vec.data();



  for(int l=0; l<num_reg;l++){
    locVec.emplace_back(l,split.split_range(BF_d,l),split.split_range(prior_v_d,l),split.split_range(pip_v_d,l));
  }

}

void controller::load_annotations_R(std::vector<std::string> names){

  kd = Xd.cols();

  dvar_name_vec=names;
  //  Xd =  anno_mat;
}







void controller::init_params(){
  

  prior_vec = Eigen::ArrayXd::Constant(init_pi1,prior_vec.size());

  
}


void controller::run_EM(){  

  double last_log10_lik = -9999999;
  init_params();

  while(1){
    
    
    double curr_log10_lik = 0;


    for(auto &lv : locVec){
   
      lv.EM_update();
      curr_log10_lik += lv.log10_lik;
               
    }
    if(ncoef==1){
      simple_regression();
    }
    // only categrical annotations
    else if(kd!=0){
      if(kd == 1 && !force_logistic){
        single_ct_regression();
      }else{
        logistic_cat_fit(beta_vec, Xd, pip_vec, l1_lambda,l2_lambda);
      }
      
      logistic_cat_pred(beta_vec, Xd, prior_vec);
    }

    if(fabs(curr_log10_lik-last_log10_lik)<EM_thresh){
      final_log10_lik = curr_log10_lik;
      break;
    }

    last_log10_lik = curr_log10_lik;
  }
  
  
  double tp = prior_vec.sum();
  // for(int i=0;i<p;i++){
  //   tp += gsl_vector_get(prior_vec,i);
  // }
  
}



void controller::simple_regression(){
  
  double sum = pip_vec.sum();

  double new_prior = sum/p;
  prior_vec.setConstant(new_prior);
  // for(int i=0;i<p;i++){
  //   gsl_vector_set(prior_vec,i,new_prior);
  // }
  beta_vec[0]=log(new_prior/(1-new_prior));
  //  gsl_vector_set(beta_vec,0, );
}














// void controller::single_probt_est(){
  
//   if(!(kd==1) && !(kd==0))
//     return;
  
//   // double c00 = 0;
//   // double c11 = 0;
//   // double c10 = 0;
//   // double c01 = 0;

//   auto c11=Xd.col(0).dot(pip_vec);
//   auto c10 = (Xd.col(0)+(-1)).dot(pip_vec);
//   auto  c01 = Xd.col(0).dot(1-pip_vec);
//   auto c00 = (1-pip_vec).dot(1-Xd.col(0));
//   // for(int i=0;i<p;i++){
//   //   double qv = double(gsl_matrix_int_get(Xd,i,0));
    
//   //   double pv = gsl_vector_get(pip_vec,i);
//   //   c11 += pv*qv;
//   //   c10 += pv*(1-qv);
//   //   c01 += (1-pv)*qv;
//   //   c00 += (1-pv)*(1-qv);
//   // }
  
//   double sd = sqrt( 1/c11 + 1/c10 + 1/c01 + 1/c00);
//   if(l2_lambda==0){
//     l2_lambda = sd;
//     // f//printf(stderr, "Applying adaptive L2 penalty = %.3f\n\n", l2_lambda);
//     //// f//printf(stderr, "Set L2 penalty = %.3f\n", sd);
//     force_logistic = 1;
//   }
//   //// f//printf (stderr, "CT sd = %.3f  %.1f \n",sd, c11);
// }



void controller::single_ct_regression(){
    using namespace Eigen;  

  std::array<double,2> sum_pip{0,0};
  std::array<double,2> sum{0,0};

  int levels = 2;
  // std::cerr<<"levels:"<<levels<<std::endl;

  double tot_pip_sum = pip_vec.sum();
  sum_pip[1]= (Xd.col(0).cast<double>().dot(pip_vec));
  sum_pip[0]=tot_pip_sum-sum_pip[1];
  sum[1]=Xd.sum();
  sum[0]=p-sum[1];

  prior_vec.setConstant(sum_pip[0]/sum[0]);
//  SparseMatrix<double> mat(rows,cols);
  double pvd = (sum_pip[1]/sum[1])-sum_pip[0]/sum[0];
  prior_vec+=(pvd*Xd.cast<double>());


  beta_vec[0]=log(sum_pip[0]/sum[0]);
  beta_vec[1] = sum_pip[1]/sum[1]-beta_vec[0];

}




// option 1: find egene



// option 2: parameter estimation

Rcpp::DataFrame controller::estimate(){

  
  if(!finish_em){
    run_EM();
    finish_em = 1;
  }


  saved_beta_vec= beta_vec;
  saved_prior_vec=prior_vec;


  




  // CI for the intercept term
  double est = beta_vec[0];
  beta_vec[0]=0;

  if(kd!=0){
    logistic_cat_pred(beta_vec, Xd,prior_vec);
  }
  
  if(kd ==0){
    prior_vec.setConstant(0.5);
  }

  double null_log10_lik = 0;
  for(int k=0;k<locVec.size();k++){
    locVec[k].EM_update();
    null_log10_lik += locVec[k].log10_lik;
  }
  double diff = (final_log10_lik - null_log10_lik)/log10(exp(1));
  if(diff<1e-8){
    diff = 1e-8;
  }
  
  double sd = fabs(est)/sqrt(2*diff);

  estvec.push_back(est);
  low_vec.push_back(est-1.96*sd);
  high_vec.push_back(est+1.96*sd);
  name_vec.push_back("Intercept");

  beta_vec= saved_beta_vec;
  prior_vec=saved_prior_vec;





  int index = 1; //start with 1, 0 always intercept
  
  
  
  for(int i=0;i<dvar_name_vec.size();i++){
    
    //    int level = gsl_vector_int_get(dlevel,i);
    string prefix = dvar_name_vec[i];
    //    name_vec.push_back(prefix);
    for(int j=1;j<2;j++){
      ostringstream stream;
      stream <<prefix<<"."<<j;
      string label = stream.str();
      double est = beta_vec[index];
      beta_vec[index]=0;

      if(kd!=0){
	logistic_cat_pred(beta_vec, Xd,prior_vec);
      }

      double null_log10_lik = 0;
     
      for(int k=0;k<locVec.size();k++){
	locVec[k].EM_update();
	null_log10_lik += locVec[k].log10_lik;
      }
  
      double diff = (final_log10_lik - null_log10_lik)/log10(exp(1));
      


      if(diff<0){
	double curr_log10_lik = final_log10_lik;
	est = fine_optimize_beta(index, est, null_log10_lik,  curr_log10_lik);
	diff = fabs(curr_log10_lik - null_log10_lik)/log10(exp(1));
      }
      
      if(diff<1e-8){
	diff = 1e-8;
      }
      double sd = fabs(est)/sqrt(2*diff);

      estvec.push_back(est);
      low_vec.push_back(est-1.96*sd);
      high_vec.push_back(est+1.96*sd);
      name_vec.push_back(prefix);
      index++;

      beta_vec=saved_beta_vec;
      prior_vec=saved_prior_vec;
    }
  }
  
  

  
  // restore all
  for(int k=0;k<locVec.size();k++){
    locVec[k].EM_update();
  }
  
  //  gsl_vector_free(saved_beta_vec);
  //  gsl_vector_free(saved_prior_vec);


  using namespace Rcpp;
  auto ret = Rcpp::DataFrame::create(_["term"]=Rcpp::wrap(name_vec),
				     _["estimate"]=Rcpp::wrap(estvec),
				     _["high"]=Rcpp::wrap(high_vec),
				     _["low"]=Rcpp::wrap(low_vec));
  return ret;

}





						
    

double controller::fine_optimize_beta(int index, double est, double null_log10_lik, double &curr_log10_lik){
  
  double gr = (sqrt(5)-1)/2.0;
  
  /*
  double a = -fabs(est);
  double b = fabs(est);
  */

  double a = est;
  double b = 0;
  if(a > b){
    b = a;
    a = 0;
  }
  double c = b - gr*(b-a);
  double d = a + gr*(b-a);

  double thresh = 1e-3;
  if(fabs(c-d)<thresh){
    thresh = fabs(c-d)/2.0;
  }

  double fc;
  double fd;


  while((d-c)>thresh){
    
    fc = eval_likelihood(c, index);
    fd = eval_likelihood(d, index);

    
    ////printf("%f %f %f %f   %f %f   %f\n",a, d
    if(fc > fd){

      b = d;
      d = c;
      c = b - gr*(b-a);
    }else{
      /*
      if(fd > null_log10_lik){
        curr_log10_lik = fd;
	return d;
      }
      */
      a = c;
      c = d;
      d = a + gr*(b-a);
    }
  }
  

  curr_log10_lik = fc;

  return (b+a)/2;

  
  
}



double controller::eval_likelihood(double x, int index){

  beta_vec[index]=x;
  if(kd!=0){
    logistic_cat_pred(beta_vec, Xd, prior_vec);
  }

  double log10_lik = 0;
  for(int k=0;k<locVec.size();k++){
    locVec[k].EM_update();
    log10_lik += locVec[k].log10_lik;
  }
  return log10_lik;
}


void Locus::EM_update(){
    // compute log10_lik

  double locus_pi0=1;
  prior_sp=prior_sp.array()/(1-prior_sp.array());
  locus_pi0=(1-prior_sp.array()).prod()*locus_pi0;


  if(locus_pi0 < 1e-100){
    locus_pi0 = 1e-100;
  }

  prior_sp*=locus_pi0;


  double max_el = BF_sp.maxCoeff();

  max_el = std::max(max_el,0.0);

  double sum = locus_pi0*pow(10,0-max_el);
  log10_lik = sum+BF_sp.matrix().dot(Eigen::pow(10,prior_sp-max_el).matrix());

  log10_lik=max_el+log10(log10_lik);

  prior_sp=Eigen::pow(10,(Eigen::log10(prior_sp)+BF_sp-log10_lik));

}






void Locus::compute_fdr(){

  double locus_pi0 = 1;

  prior_sp=prior_sp.array()/(1-prior_sp.array());
  locus_pi0=(1-prior_sp.array()).prod()*locus_pi0;

  prior_sp*=locus_pi0;



  double & ll = log10_lik;
  ////printf("%s  log10_lik = %f\n",id.c_str(), log10_lik);

  double max_el = BF_sp.maxCoeff();

  max_el = std::max(max_el,0.0);

  double sum = locus_pi0*pow(10,0-max_el);
  log10_lik = sum+BF_sp.matrix().dot(Eigen::pow(10,prior_sp-max_el).matrix());

  log10_lik=max_el+log10(log10_lik);



  fdr =  pow(10,log10(locus_pi0)-log10_lik);

}



// void Locus::compute_fdr(){

//   double locus_pi0 = 1;

//   // compute log10_lik
//   vector<double> BF_vec;
//   vector<double> p_vec;
//   for(int i=0;i<snpVec.size();i++){
//     double prior = gsl_vector_get(prior_vec, snpVec[i].index);
//     BF_vec.push_back(snpVec[i].log10_BF);
//     p_vec.push_back(prior/(1-prior));
//     locus_pi0 *= (1-prior);
//   }

//   for(int i=0;i<snpVec.size();i++){
//     p_vec[i] *= locus_pi0;
//   }


//   BF_vec.push_back(0);
//   p_vec.push_back(locus_pi0);

//   log10_lik = log10_weighted_sum(BF_vec, p_vec);
//   fdr =  pow(10,log10(locus_pi0)-log10_lik);

// }






double log10_weighted_sum(vector<double> &vec, vector<double> &wts){


  double max = vec[0];
  for(size_t i=0;i<vec.size();i++){
    if(vec[i]>max)
      max = vec[i];
  }
  double sum = 0;
  for(size_t i=0;i<vec.size();i++){
    sum += wts[i]*pow(10, (vec[i]-max));
  }

  return (max+log10(sum));
}

bool rank_by_fdr (const Locus & loc1 , const Locus & loc2){
  return (loc1.fdr < loc2.fdr);
}
