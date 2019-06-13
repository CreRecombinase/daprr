#ifndef __CLASSDEF_H_
#define __CLASSDEF_H_
#ifndef EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS
#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS
#endif
#include <RcppEigen.h>
#include <stdlib.h>
#include <string>
#include <string.h>
#include <vector>
#include <map>
#include <unordered_map>




using namespace std;




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
  
  vector<SNP> snpVec;



  
  double log10_lik; // log10 of marginal likelihood
  
  double fdr;

  Eigen::Map<Eigen::ArrayXd> BF_sp;
  Eigen::Map<Eigen::ArrayXd> prior_sp;
  Eigen::ArrayXd prior_v;
  Eigen::Map<Eigen::ArrayXd> pip_sp;


  Locus(int locus_id,  Eigen::Map<Eigen::ArrayXd> BFv,Eigen::Map<Eigen::ArrayXd> pv,Eigen::Map<Eigen::ArrayXd> pipv):id(locus_id),
														     BF_sp(BFv),
														     prior_sp(pv),
														     prior_v(prior_sp),
														     pip_sp(pipv){}








  // }


  
  //  void EM_update();
  void EM_update();
  void compute_fdr();


};






class splitter{
  std::vector<size_t> ret_v;
public:
  splitter(std::vector<size_t> &&ret_v_):ret_v(ret_v_){}
  Eigen::Map<Eigen::ArrayXd> split_range(double* data_pt,int idx)const {
    int beg=ret_v[idx];
    int	ret_s =	ret_v[idx+1]-ret_v[idx];
    return(Eigen::Map<Eigen::ArrayXd>(data_pt+beg,ret_s));
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
    Rcpp::NumericVector r_beta;  
  Rcpp::NumericVector r_pip;
  Rcpp::NumericVector r_prior;
  Rcpp::NumericVector r_o_prior;
  Rcpp::NumericVector r_o_beta;
  
  Eigen::Map<Eigen::VectorXd> beta;
  Eigen::Map<Eigen::VectorXd> pip;
  Eigen::Map<Eigen::VectorXd> prior;
  Eigen::Map<Eigen::VectorXd> o_prior;
  Eigen::Map<Eigen::VectorXd> o_beta;


  Result_obj(const size_t	p,const size_t k):
        r_beta(k),  
    r_pip(p),
    r_prior(p),
    r_o_prior(p),
    r_o_beta(k),
    beta(Rcpp::as<Eigen::Map<Eigen::VectorXd> >(r_beta)),
    pip(Rcpp::as<Eigen::Map<Eigen::VectorXd> >(r_pip)),
    prior(Rcpp::as<Eigen::Map<Eigen::VectorXd> >(r_prior)),
    o_prior(Rcpp::as<Eigen::Map<Eigen::VectorXd> >(r_o_prior)),
    o_beta(Rcpp::as<Eigen::Map<Eigen::VectorXd> >(r_o_beta)){}


  Eigen::Map<Eigen::VectorXd>saved_beta_vec(){
    return o_beta;// = gsl_vector_calloc(ncoef);
  }
  Eigen::Map<Eigen::VectorXd>saved_prior_vec(){
    return o_prior;// = gsl_vector_calloc(p);
  }
  Eigen::Map<Eigen::VectorXd>prior_vec(){
    return prior;
  }
  Eigen::Map<Eigen::VectorXd>pip_vec(){
    return pip;
  }
  Eigen::Map<Eigen::VectorXd>beta_vec(){
    return beta;
  }
};



class controller {
  const splitter &split;
public:


  controller(const splitter &split_,Result_obj &obj,Eigen::Map<Eigen::SparseMatrix<double>> Xd_):split(split_),
						     saved_beta_vec(Rcpp::as<Eigen::Map<Eigen::VectorXd> >(obj.r_o_beta)),
						     saved_prior_vec(Rcpp::as<Eigen::Map<Eigen::VectorXd> >(obj.r_o_prior)),
						     prior_vec(Rcpp::as<Eigen::Map<Eigen::VectorXd> >(obj.r_prior)),
						     pip_vec(Rcpp::as<Eigen::Map<Eigen::VectorXd> >(obj.r_pip)),
						     beta_vec(Rcpp::as<Eigen::Map<Eigen::VectorXd> >(obj.r_beta)),
						     Xd(Xd_)
  {
    p=kd= 0;
    force_logistic = 0;
    finish_em = 0;
    single_fuzzy_annot = 0;
    l1_lambda=l2_lambda=0;
    init_pi1 = 1e-3;
    print_avg = 0;
  }


  Eigen::Map<Eigen::VectorXd>saved_beta_vec;// = gsl_vector_calloc(ncoef);
  Eigen::Map<Eigen::VectorXd>saved_prior_vec;// = gsl_vector_calloc(p);
  Eigen::Map<Eigen::VectorXd>prior_vec;
  Eigen::Map<Eigen::VectorXd>pip_vec;
  Eigen::Map<Eigen::VectorXd>beta_vec;


  // storage
  vector<Locus> locVec;
  
  
  int p; // number of loc-SNP pairs
  
  //  int kc; // number of continuous covariate
  int kd; // number of discrete covariate

  Eigen::VectorXd log10_BF;

  std::vector<double> estvec;
  std::vector<double> low_vec;
  std::vector<double> high_vec;
  std::vector<std::string> name_vec;

  vector<string> dvar_name_vec;
  Eigen::SparseMatrix<double> Xd;

  
  //  gsl_vector_int *dist_bin;
  double EM_thresh;
  



  

  
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
  void load_annotations_R(std::vector<std::string> names);

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

  Rcpp::DataFrame estimate();
  //  void dump_prior(char *path);
  //  double dump_locus_prior(const Locus& loc,fs::path file);
  //  void dump_pip(char *file);


 private:
  double eval_likelihood(double x, int index);
  double fine_optimize_beta(int index, double est, double null_log10_lik, double &curr_log10_lik);
  



};




double log10_weighted_sum(vector<double> & val_vec, vector<double> & wts_vec);
double compute_log10_BF(double beta, double se_beta);

bool   rank_by_fdr (const Locus & loc1 , const Locus & loc2);
  
int classify_dist_bin(int snp_pos, int tss, double bin_size = -1);
double map_bin_2_dist(int bin, double bin_size=-1);
#endif
