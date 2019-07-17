#include "classdef.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h>
#include <math.h>
#include <RcppEigen.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <boost/config/warning_disable.hpp>
#include "logistic.hpp"
#include <wordexp.h>
#include <sys/types.h>
#include <errno.h>



BF::BF():wt(1/4.0),
	 l10(log(10)){
  kv={1,4,16,25};
  for(int i=0; i<4; i++){
    kva[i]= 0.5*log(1/(1+kv[i]));
    kvb[i]=0.5*(kv[i]/(1+kv[i]));
  }
}

double BF::compute_log10_BF(const double z_score) const {
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



void controller::load_data_R(Rcpp::NumericVector z_hat){

  if(verbose){
    Rcpp::Rcerr<<"loading data"<<std::endl;
  }
  std::vector<SNP> snp_vec;


  p=static_cast<int>(z_hat.size());
  log10_BF.resize(p);

  //  assert_size(anno_mat.col(0),"anno_mat col 1",p);
  BF bf;


  std::transform(z_hat.begin(),z_hat.end(),log10_BF.begin(),[&bf](double zh){
							      return(bf.compute_log10_BF(zh));
							    });

  const auto num_reg = split.num_regions();
  auto BF_d = log10_BF.data();
  auto prior_v_d = prior_vec->data;
  auto pip_v_d = pip_vec->data;



  for(int l=0; l<num_reg;l++){
    locVec.emplace_back(l,split.split_range(BF_d,l),split.split_range(prior_v_d,l),split.split_range(pip_v_d,l));
  }
  if(verbose){
    Rcpp::Rcerr<<"data loaded"<<std::endl;
  }

}

void controller::load_annotations_R(RcppGSL::matrix<int> anno_mat,std::vector<std::string> names){

  if(verbose){
    Rcpp::Rcerr<<"loading annotation"<<std::endl;
  }
  kd = anno_mat.ncol();
  if(anno_mat.ncol()!=names.size()){
    Rcpp::stop("names and anno_mat do not match in size:"+std::to_string(kd)+" vs "+std::to_string(names.size()));
  }
  dvar_name_vec=names;
  Xd =  anno_mat;
  dlevel = gsl_vector_int_calloc(kd);

  for(int i=0;i<kd;i++){
    int nl = count_factor_level(i);
    gsl_vector_int_set(dlevel, i,nl);
  }

  if(verbose){
    Rcpp::Rcerr<<"annotation loaded"<<std::endl;
  }

}



int controller::count_factor_level(int col){
 
  std::unordered_map<int, int> rcd;
  for(int i=0;i<p;i++){
    int val = gsl_matrix_int_get(Xd,i,col);
    rcd[val] = 1;
  }

  return rcd.size();
}



void controller::init_params(){
  
  if(verbose){
    Rcpp::Rcerr<<"initializing"<<std::endl;
  }
  for(int i=0;i<p;i++){
    gsl_vector_set(prior_vec, i, init_pi1);
  }


  // for(int i=0;i<locVec.size();i++){
  //   locVec[i].pip_vec = pip_vec;
  //   locVec[i].prior_vec = prior_vec;
  // }

  if(verbose){
    Rcpp::Rcerr<<"initialized"<<std::endl;
  }
}


void controller::run_EM(const bool use_glmnet){


  double last_log10_lik = -9999999;
  init_params();
  int itcc=0;
  Logistic logit(beta_vec->size,p);
  while(1){
    if(verbose){
      Rcpp::Rcerr<<"iter no:"<<itcc++<<std::endl;
    }
    double curr_log10_lik = 0;
    for(auto &lv : locVec){
      lv.EM_update();
      curr_log10_lik += lv.log10_lik;
    }
    if(ncoef==1){
      if(verbose){
	Rcpp::Rcerr<<"simple"<<std::endl;
      }
      simple_regression();
    }
    // only categrical annotations
    else if(kd!=0){
      if(kd == 1 && !force_logistic){
	if(verbose){
	  Rcpp::Rcerr<<"single_ct"<<std::endl;
	}
	single_ct_regression();
      }else{
	if(verbose){
	  Rcpp::Rcerr<<"cat_fit"<<std::endl;
	}
	if(use_glmnet){
	  logit.fit_glmnet(beta_vec,Xd,pip_vec,l1_lambda,l2_lambda);
	}else{
	  logit.fit(beta_vec,Xd,dlevel,pip_vec,l1_lambda,l2_lambda);
	}
	//	logistic_cat_fit(beta_vec, Xd, dlevel,pip_vec, l1_lambda,l2_lambda);
      }
      if(verbose){
	Rcpp::Rcerr<<"cat_pred"<<std::endl;
	Rcpp::Rcerr<<"beta:\n";
      }
      gsl::span<double> beta_sp(beta_vec->data,beta_vec->size);
      if(std::accumulate(beta_sp.begin(),beta_sp.end(),0)<1e-8){
	Rcpp::Rcerr<<"Beta is entirely 0..."<<std::endl;
      }
      for(auto be : beta_sp){
	if(verbose){
	  Rcpp::Rcerr<<be<<"\n";
	}
	if(isnan(be)){
	  Rcpp::stop("NaN encountered in beta");
	}
      }
      if(verbose){
	Rcpp::Rcerr<<std::endl;
      }
      logistic_cat_pred(beta_vec, Xd, dlevel,prior_vec);
    }


    if(fabs(curr_log10_lik-last_log10_lik)<EM_thresh){
      final_log10_lik = curr_log10_lik;
      break;
    }

    last_log10_lik = curr_log10_lik;
  }


  double tp = 0;
  for(int i=0;i<p;i++){
    tp += gsl_vector_get(prior_vec,i);
  }

}



void controller::simple_regression(){
  
  double sum = 0;
  for(int i=0;i<p;i++){
    sum += gsl_vector_get(pip_vec,i);
  }
  double new_prior = sum/p;
  for(int i=0;i<p;i++){
    gsl_vector_set(prior_vec,i,new_prior);
  }
  gsl_vector_set(beta_vec,0, log(new_prior/(1-new_prior)));
}



void controller::single_probt_est(){
  
  if(!(kd==1) && !(kd==0))
    return;
  
  double c00 = 0;
  double c11 = 0;
  double c10 = 0;
  double c01 = 0;
  
  for(int i=0;i<p;i++){
    double qv = double(gsl_matrix_int_get(Xd,i,0));
    
    double pv = gsl_vector_get(pip_vec,i);
    c11 += pv*qv;
    c10 += pv*(1-qv);
    c01 += (1-pv)*qv;
    c00 += (1-pv)*(1-qv);
  }
  
  double sd = sqrt( 1/c11 + 1/c10 + 1/c01 + 1/c00);
  if(l2_lambda==0){
    l2_lambda = sd;
    force_logistic = 1;
  }

}



void controller::single_ct_regression(){
 
  std::unordered_map<int,double> sum_pip;
  std::unordered_map<int,double> sum;
  
  int levels = gsl_vector_int_get(dlevel,0);
  // std::cerr<<"levels:"<<levels<<std::endl;

  for(int i=0;i<levels;i++){
    sum_pip[i] = sum[i] = 0;
  }

  for(int i=0;i<p;i++){
    int cat = gsl_matrix_int_get(Xd,i,0);
    sum_pip[cat] += gsl_vector_get(pip_vec,i);
    sum[cat] += 1;
  }
  
  
  for(int i=0;i<p;i++){
    int cat = gsl_matrix_int_get(Xd,i,0);
    gsl_vector_set(prior_vec,i,sum_pip[cat]/sum[cat]);
  }
  

  double baseline=0;
  for(int i=0;i<levels;i++){
    double new_prior = sum_pip[i]/sum[i];
    gsl_vector_set(beta_vec,i, log(new_prior/(1-new_prior))-baseline);
    if(i==0){
      baseline = log(new_prior/(1-new_prior));
    }
 
  }
  
  
}




// option 2: parameter estimation

Rcpp::DataFrame controller::estimate(bool use_glmnet){

  if(!finish_em){
    run_EM(use_glmnet);
    finish_em = 1;
  }


  gsl_vector_memcpy(saved_beta_vec, beta_vec);  
  gsl_vector_memcpy(saved_prior_vec, prior_vec);

  // CI for the intercept term
  double est = gsl_vector_get(beta_vec, 0);
  gsl_vector_set(beta_vec,0, 0.0);

  if(kd!=0){
    logistic_cat_pred(beta_vec, Xd, dlevel,prior_vec);
  }
  
  if(kd ==0){
    for(int i=0;i<p;i++){
      gsl_vector_set(prior_vec,i,0.50);
    }
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

  gsl_vector_memcpy(beta_vec, saved_beta_vec);
  gsl_vector_memcpy(prior_vec, saved_prior_vec);





  int index = 1; //start with 1, 0 always intercept
  
  
  
  for(int i=0;i<dvar_name_vec.size();i++){
    
    int level = gsl_vector_int_get(dlevel,i);
    std::string prefix = dvar_name_vec[i];
    //    name_vec.push_back(prefix);
    for(int j=1;j<level;j++){
      std::ostringstream stream;
      stream <<prefix<<"."<<j;
      std::string label = stream.str();
      double est = gsl_vector_get(beta_vec, index);
      gsl_vector_set(beta_vec,index, 0.0);

      if(kd!=0){
	logistic_cat_pred(beta_vec, Xd, dlevel,prior_vec);
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

      gsl_vector_memcpy(beta_vec, saved_beta_vec);  
      gsl_vector_memcpy(prior_vec, saved_prior_vec);
    }
  }
  
  

  
  // restore all
  for(int k=0;k<locVec.size();k++){
    locVec[k].EM_update();
  }
  
  gsl_vector_free(saved_beta_vec);
  gsl_vector_free(saved_prior_vec);

  Rcpp::NumericVector lik(name_vec.size(),final_log10_lik);
  using namespace Rcpp;
  auto ret = Rcpp::DataFrame::create(_["term"]=Rcpp::wrap(name_vec),
				     _["estimate"]=Rcpp::wrap(estvec),
				     _["high"]=Rcpp::wrap(high_vec),
				     _["low"]=Rcpp::wrap(low_vec),
				     _["lik"]=final_log10_lik);
  return ret;

}





						
    

double controller::fine_optimize_beta(int index, double est, double null_log10_lik, double &curr_log10_lik){
  
  double gr = (sqrt(5)-1)/2.0;


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

    if(fc > fd){

      b = d;
      d = c;
      c = b - gr*(b-a);
    }else{
      a = c;
      c = d;
      d = a + gr*(b-a);
    }
  }


  curr_log10_lik = fc;

  return (b+a)/2;


  
}



double controller::eval_likelihood(double x, int index){

  gsl_vector_set(beta_vec,index, x);
  if(kd!=0){
    logistic_cat_pred(beta_vec, Xd, dlevel,prior_vec);
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



  std::transform(prior_sp.begin(),prior_sp.end(),prior_sp.begin(),[&locus_pi0](double prior){
								     locus_pi0*=(1-prior);
								     return(prior/(1-prior));
								  });



  if(locus_pi0 < 1e-100){
    locus_pi0 = 1e-100;
  }


  for(auto &p_v : prior_sp)
    p_v *= locus_pi0;

  double max_el = *(std::max_element(BF_sp.begin(),BF_sp.end()));

  max_el = std::max(max_el,0.0);

  double sum = locus_pi0*pow(10,0-max_el);
  log10_lik =	std::inner_product(BF_sp.begin(),
				   BF_sp.end(),
				   prior_sp.begin(),
				   sum,std::plus<>(),
				   [max_el](double vec,double wts){
				     return(wts*pow(10,vec-max_el));});
  log10_lik=max_el+log10(log10_lik);


  double & ll = log10_lik;

  std::transform(prior_sp.begin(),prior_sp.end(),BF_sp.begin(),pip_sp.begin(),
		 [ll](double prior,double log10_BF){
		   return pow(10,(log10(prior) + log10_BF - ll));

		 });

}



void Locus::compute_fdr(){

  double locus_pi0 = 1;
  std::transform(prior_sp.begin(),prior_sp.end(),prior_sp.begin(),[&locus_pi0](double prior){
								     locus_pi0*=(1-prior);
								     return(prior/(1-prior));
								  });

  for(auto &p_v : prior_sp)
    p_v *= locus_pi0;


  double & ll = log10_lik;
  ////printf("%s  log10_lik = %f\n",id.c_str(), log10_lik);
  std::transform(prior_sp.begin(),prior_sp.end()-1,BF_sp.begin(),pip_sp.begin(),
		 [ll](double prior,double log10_BF){
		   return pow(10,(log10(prior) + log10_BF - ll));
		 });


   double max_el = *(std::max_element(BF_sp.begin(),BF_sp.end()));
  max_el = std::max(max_el,0.0);
  // for(size_t i=0;i<vec.size();i++){
  //   if(vec[i]>max)
  //     max = vec[i];
  // }
  double sum = locus_pi0*pow(10,0-max_el);
  log10_lik =	std::inner_product(BF_sp.begin(),
				   BF_sp.end(),
				   prior_sp.begin(),
				   sum,std::plus<>(),
				   [max_el](double wts,double vec){
				     return(wts*pow(10,vec-max_el));});
  log10_lik=max_el+log10(log10_lik);

  fdr =  pow(10,log10(locus_pi0)-log10_lik);

}






double log10_weighted_sum(std::vector<double> &vec, std::vector<double> &wts){


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
