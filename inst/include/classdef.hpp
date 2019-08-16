#pragma once
#include <RcppGSL.h>
#include <stdlib.h>
#include <string>
#include <string.h>
#include <vector>
#include <map>
#include <unordered_map>
#include "gsl/span"
#include "torus.hpp"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

extern bool verbose;

namespace elasticdonut {

  double prior_to_pip(gsl::span<double> pip,const gsl::span<double> prior,const bool fix_pi0=true);

  void transform_pip(gsl::span<double> pip,const gsl::span<double> BF,const double log10lik);
  double log10_lik(const double locus_pi0,const gsl::span<double> BF, const gsl::span<double> pip,std::optional<double> &BF_max);
  double E_step(gsl::span<double> pip,const gsl::span<double> BF, const gsl::span<double> prior,std::optional<double> &BF_max);
  double fdr(double log10lik,const gsl::span<double> BF, const gsl::span<double> prior,gsl::span<double> pip,std::optional<double> &BF_max);
}


namespace donut {

class BF{
  const double wt;
  const double l10;
  std::array<double,4> kv;
  std::array<double,4> kva;
  std::array<double,4> kvb;
  mutable std::array<double,4> kvc;

public:
  BF();
  double operator()(const double z_score) const{
    return(this->compute_log10_BF(z_score));
  }
  double compute_log10_BF(const double z_score) const;
};

  inline  std::vector<double> make_BF(Rcpp::NumericVector z_hat){
    std::vector<double> BFv(z_hat.size());
    BF bf;
    std::transform(z_hat.begin(),z_hat.end(),BFv.begin(),[&bf](double z){
							   return(bf.compute_log10_BF(z));});
    return BFv;
  }

  inline  std::vector<double> make_BF(const gsl::span<double> z_hat){
    std::vector<double> BFv(z_hat.size());
    BF bf;
    std::transform(z_hat.begin(),z_hat.end(),BFv.begin(),[&bf](double z){
							   return(bf.compute_log10_BF(z));});
    return BFv;
  }

class SNP {

 public:

  int id;
  int index;
  double log10_BF;
  SNP():id(-1),index(-1),log10_BF(-1){}
  SNP(int snp_id, double snp_log10_BF, int snp_index):id(snp_id),index(snp_index),log10_BF(snp_log10_BF){
  }
};


class Locus {
  
 public:
  
  // possible gene information/annotation
  int id;
  
  // gsl_vector *prior_vec;
  // gsl_vector *pip_vec;
  double log10_lik; // log10 of marginal likelihood
  double fdr;

  gsl::span<const double> BF_sp;
  gsl::span<double> prior_sp;
  gsl::span<double> pip_sp;


  Locus(int locus_id,  gsl::span<const double> BFv,gsl::span<double> pv,gsl::span<double> pipv):id(locus_id),
												BF_sp(BFv),
												prior_sp(pv),
												pip_sp(pipv){}
  Locus(){};

  
  //  void EM_update();
  void EM_update();
  void compute_fdr();


};






class Result_obj{

public:
  RcppGSL::vector<double> beta;
  RcppGSL::vector<double> pip;
  RcppGSL::vector<double> prior;
  RcppGSL::vector<double> o_prior;
  RcppGSL::vector<double> o_beta;


  Result_obj(const size_t	p,const size_t k):
    beta(k),
    pip(p),
    prior(p),
    o_prior(p),
    o_beta(k){}
  Result_obj(Rcpp::List obj):beta(obj["beta"]),
			     pip(obj["pip"]),
			     prior(obj["prior"]),
			     o_prior(obj["o_prior"]),
			     o_beta(obj["o_beta"]){}
};








class controller {
  elasticdonut::splitter split;
  
public:

  RcppGSL::vector<double> saved_beta_vec;// = gsl_vector_calloc(ncoef);
  RcppGSL::vector<double> saved_prior_vec;// = gsl_vector_calloc(p);
  RcppGSL::vector<double> prior_vec;
  RcppGSL::vector<double> pip_vec;
  RcppGSL::vector<double> beta_vec;


  controller(const elasticdonut::splitter &split_,Result_obj &obj,RcppGSL::matrix<int> anno_mat,double EM_thresh_=0.05,double init_pi1_=1e-3,int print_avg_=0):
    split(split_),
    saved_beta_vec(obj.o_beta),
    saved_prior_vec(obj.o_prior),
    prior_vec(obj.prior),
    pip_vec(obj.pip),
    beta_vec(obj.beta),
    EM_thresh(EM_thresh_),
    init_pi1(init_pi1_),
    print_avg(print_avg_),
    dlevel(anno_mat.ncol()),
    Xd(anno_mat)

  {
    p=kd= 0;
    force_logistic = 0;
    finish_em = 0;
    single_fuzzy_annot = 0;
    l1_lambda=l2_lambda=0;
    init_pi1 = 1e-3;
    print_avg = 0;
  }


  // storage
  std::vector<Locus> locVec;
  
  
  int p; // number of loc-SNP pairs
  
  //  int kc; // number of continuous covariate
  int kd; // number of discrete covariate

  std::vector<double> log10_BF;
  std::vector<double> estvec;
  std::vector<double> low_vec;
  std::vector<double> high_vec;
  std::vector<std::string> name_vec;

  std::vector<std::string> dvar_name_vec;

  double EM_thresh;
  RcppGSL::vector<int> dlevel; // (kd+1) entry levels of each factor
  RcppGSL::matrix<int> Xd; // p x kd

  int ncoef;

  double final_log10_lik;

  double init_pi1;
  int nthread;
  int finish_em;

  int print_avg;
  int force_logistic;
  int single_fuzzy_annot;

  double l1_lambda, l2_lambda; //shrinkage for enrich param est


  void load_data_R(Rcpp::NumericVector z_hat);
  void load_annotations_R(RcppGSL::matrix<int> anno_mat ,std::vector<std::string> names);

  int count_factor_level(int col);

  void simple_regression();
  void single_ct_regression();
  void single_probt_est();
  void init_params();
  void run_EM(const bool use_glmnet);
  std::vector<int> get_region_id();

  // template<typename T>
  // void serialize(T & archive){

  //   archive(saved_beta_vec,// = gsl_vector_calloc(ncoef);
  // 	    saved_prior_vec,prior_vec,pip_vec,beta_vec,)


  Rcpp::DataFrame estimate(bool use_glmnet);

 private:
  double eval_likelihood(double x, int index);
  double fine_optimize_beta(int index, double est, double null_log10_lik, double &curr_log10_lik);
  

  };




double log10_weighted_sum(std::vector<double> & val_vec, std::vector<double> & wts_vec);
double compute_log10_BF(double beta, double se_beta);

bool   rank_by_fdr (const Locus & loc1 , const Locus & loc2);
  
int classify_dist_bin(int snp_pos, int tss, double bin_size = -1);
double map_bin_2_dist(int bin, double bin_size=-1);


}
