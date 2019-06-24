#ifndef __CLASSDEF_H_
#define __CLASSDEF_H_

#include <RcppGSL.h>
#include <stdlib.h>
#include <string>
#include <string.h>
#include <vector>
#include <map>
#include <unordered_map>
#include "gsl/span"


#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

extern bool verbose;





class SNP {

 public:

  int id;
  int index;
  double log10_BF;
  SNP():id(-1),index(-1),log10_BF(-1){}
  SNP(int snp_id, double snp_log10_BF, int snp_index):id(snp_id),index(snp_index),log10_BF(snp_log10_BF){
    // id = snp_id;
    // log10_BF = snp_log10_BF;
    // index = snp_index;
  }
};






class Locus {
  
 public:
  
  // possible gene information/annotation
  int id;
  
  std::vector<SNP> snpVec;
  
  gsl_vector *prior_vec;
  gsl_vector *pip_vec;


  
  double log10_lik; // log10 of marginal likelihood
  
  double fdr;

  gsl::span<const double> BF_sp;
  gsl::span<double> prior_sp;
  gsl::span<double> pip_sp;


  Locus(int locus_id,  gsl::span<const double> BFv,gsl::span<double> pv,gsl::span<double> pipv):id(locus_id),
												BF_sp(BFv),
												prior_sp(pv),
												pip_sp(pipv){}


    // BF_vec(snpVec.size()+1),
    // min_index(std::numeric_limits<int>::max()),
    // max_index(std::numeric_limits<int>::min()){


    // std::transform(snpVec.begin(),snpVec.end(),BF_vec.begin(),[&min_index,&max_index](const SNP& snp){
    // 								min_index = snp.index < min_index ? snp.index : min_index;
    // 								max_index = snp.index > max_index ? snp.index : max_index;
    // 								return snp.log10_BF;
    // 							      });

  
  
  // Locus(int locus_id,  std::vector<SNP>  snpVec_):id(locus_id),
  // 					     snpVec(snpVec_),
  // 					     prior_vec(nullptr)
  // 					     pip_vec(nullptr),p_vec(snpVec.size()+1,0.0),
  //   BF_vec(snpVec.size()+1),
  //   min_index(std::numeric_limits<int>::max()),
  //   max_index(std::numeric_limits<int>::min()){


  //   std::transform(snpVec.begin(),snpVec.end(),BF_vec.begin(),[&min_index,&max_index](const SNP& snp){
  // 								min_index = snp.index < min_index ? snp.index : min_index;
  // 								max_index = snp.index > max_index ? snp.index : max_index;
  // 								return snp.log10_BF;
  // 							      });







  // }

  Locus(){};

  
  //  void EM_update();
  void EM_update();
  void compute_fdr();


};






class splitter{
  std::vector<size_t> ret_v;
public:
  splitter(std::vector<size_t> &&ret_v_):ret_v(ret_v_){}
  gsl::span<double> split_range(double* data_pt,int idx)const {
    int beg=ret_v[idx];
    int end_s=ret_v[idx+1];
    int	ret_s =	ret_v[idx+1]-ret_v[idx];
    return(gsl::span<double>(data_pt+beg,ret_s));
  }
  const int num_regions()const {
    return(static_cast<int>(ret_v.size()-1));
  }
};


template<typename T>
splitter make_splitter(T v_begin, T v_end){
  std::vector<size_t> ret_v;
    auto t_b = v_begin;
    auto t_i = t_b;
    ret_v.push_back(0);
    while(t_i != v_end){
      t_i = std::upper_bound(t_i,v_end,*t_i);
      ret_v.push_back(std::distance(t_b,t_i));
    }
    return(splitter(std::move(ret_v)));
}

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


  gsl_vector *saved_beta_vec(){
    return o_beta;// = gsl_vector_calloc(ncoef);
  }
  gsl_vector *saved_prior_vec(){
    return o_prior;// = gsl_vector_calloc(p);
  }
  gsl_vector *prior_vec(){
    return prior;
  }
  gsl_vector *pip_vec(){
    return pip;
  }
  gsl_vector *beta_vec(){
    return beta;
  }
};



class controller {
  const splitter &split;
  
public:


  controller(const splitter &split_,Result_obj &obj):split(split_),
						     saved_beta_vec(obj.saved_beta_vec()),
						     saved_prior_vec(obj.saved_prior_vec()),
						     prior_vec(obj.prior_vec()),
						     pip_vec(obj.pip_vec()),
						     beta_vec(obj.beta_vec())
  {
    p=kd= 0;
    force_logistic = 0;
    finish_em = 0;
    single_fuzzy_annot = 0;
    l1_lambda=l2_lambda=0;
    init_pi1 = 1e-3;
    print_avg = 0;
  }


  gsl_vector *saved_beta_vec;// = gsl_vector_calloc(ncoef);
  gsl_vector *saved_prior_vec;// = gsl_vector_calloc(p);
  gsl_vector *prior_vec;
  gsl_vector *pip_vec;
  gsl_vector *beta_vec;


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


  
  //  gsl_vector_int *dist_bin;
  double EM_thresh;
  
  gsl_vector_int *dlevel; // (kd+1) entry levels of each factor


  gsl_matrix *Xc;  // p x kc
  gsl_matrix_int *Xd; // p x kd


  

  
  int ncoef;

  double final_log10_lik;
  

  double init_pi1;
  int nthread; 


  int finish_em;


  int print_avg;


  // void load_data(char *filename);    // load data with MatrixeQTL format -- default
  // void load_data_zscore(char *filename); // load data with z-score/t-score
  // void load_data_BF(char *filename); // load data with pre-computed log10 Bayes factors
  // void load_data_fastqtl(char *filename); //load data with fastQTL format
  // void torus_d(Rcpp::IntegerVector locus_id, Rcpp::NumericVector z_hat, Rcpp::IntegerMatrix anno_mat);
  void load_data_R(Rcpp::NumericVector z_hat);
  void load_annotations_R(RcppGSL::matrix<int> anno_mat ,std::vector<std::string> names);

  // void load_map(char *gene_map_file, char *snp_map_file);
  // void load_annotation(char *annot_file);
  int count_factor_level(int col);
  
  void simple_regression();
  void single_ct_regression();
  void single_probt_est();
  int force_logistic;
  int single_fuzzy_annot;
  
  double l1_lambda, l2_lambda; //shrinkage for enrich param est


  void init_params();
    
  void run_EM();



  
  void find_eGene(double thresh=0.05);
  Rcpp::DataFrame estimate();
  //  void dump_prior(char *path);
  //  double dump_locus_prior(const Locus& loc,fs::path file);
  //  void dump_pip(char *file);


 private:
  double eval_likelihood(double x, int index);
  double fine_optimize_beta(int index, double est, double null_log10_lik, double &curr_log10_lik);
  



};




double log10_weighted_sum(std::vector<double> & val_vec, std::vector<double> & wts_vec);
double compute_log10_BF(double beta, double se_beta);

bool   rank_by_fdr (const Locus & loc1 , const Locus & loc2);
  
int classify_dist_bin(int snp_pos, int tss, double bin_size = -1);
double map_bin_2_dist(int bin, double bin_size=-1);
#endif
