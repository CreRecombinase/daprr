#pragma once
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

template<class Archive>
void save(Archive & archive,
	  RcppGSL::matrix<int> const & m){

  Rcpp::IntegerMatrix rm = Rcpp::wrap(m);
  archive(rm);

}


template<class Archive>
void save(Archive & archive,
	  RcppGSL::vector<double> const & m){

  Rcpp::NumericVector rm = Rcpp::wrap(m);
  archive(rm);
}




template<class Archive>
void save(Archive & archive,
	  RcppGSL::matrix<double> const & m){

  Rcpp::NumericMatrix rm = Rcpp::wrap(m);
  archive(rm);

}


template<class Archive>
void save(Archive & archive,
	  RcppGSL::vector<int> const & m){

  Rcpp::IntegerVector rm = Rcpp::wrap(m);
  archive(rm);
}



///////


template<class Archive>
void load(Archive & archive,
	  RcppGSL::matrix<int> & m){

  Rcpp::IntegerMatrix rm;
  archive(rm);
  m = Rcpp::as<RcppGSL::matrix<int> >(rm);

}


template<class Archive>
void load(Archive & archive,
	  RcppGSL::vector<double> & m){

  Rcpp::NumericVector rm;
  archive(rm);
  m = Rcpp::as<RcppGSL::vector<double> >(rm);

}




template<class Archive>
void load(Archive & archive,
	  RcppGSL::matrix<double> & m){
  Rcpp::NumericMatrix rm;
  archive(rm);
  m = Rcpp::as<RcppGSL::matrix<double> >(rm);


}


template<class Archive>
void load(Archive & archive,
	  RcppGSL::vector<int> & m){

  Rcpp::IntegerVector rm;
  archive(rm);
  m = Rcpp::as<RcppGSL::vector<int> >(rm);
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
  double compute_log10_BF(const double z_score) const;
};



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

class splitter{
  std::vector<size_t> ret_v;
public:
  splitter(Rcpp::IntegerVector ser):ret_v(Rcpp::as<std::vector<size_t>>(ser)){
  }
  splitter(std::vector<size_t> &&ret_v_):ret_v(ret_v_){}
  gsl::span<double> split_range(double* data_pt,int idx)const {
    int beg=ret_v[idx];
    int	ret_s =	ret_v[idx+1]-ret_v[idx];
    return(gsl::span<double>(data_pt+beg,ret_s));
  }
  const int num_regions()const {
    return(static_cast<int>(ret_v.size()-1));
  }
  std::vector<gsl::span<double> > split_view(double* data_pt) const{
    const int tot_r= num_regions();
    std::vector< gsl::span<double> > groupv(tot_r);
    for(int i=0; i<tot_r; i++){
      groupv[i]=split_range(data_pt,i);
    }
    return groupv;
  }
  Rcpp::IntegerVector export_object() const{
    return Rcpp::wrap(ret_v);
  }
  // template<typename T>
  // void serialize (Archive& archive){
  //   archive(ret_v);
  // }
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
  Result_obj(Rcpp::List obj):beta(obj["beta"]),
			     pip(obj["pip"]),
			     prior(obj["prior"]),
			     o_prior(obj["o_prior"]),
			     o_beta(obj["o_beta"]){}



  // gsl_vector *saved_beta_vec(){
  //   return o_beta;// = gsl_vector_calloc(ncoef);
  // }
  // gsl_vector *saved_prior_vec(){
  //   return o_prior;// = gsl_vector_calloc(p);
  // }
  // gsl_vector *prior_vec(){
  //   return prior;
  // }
  // gsl_vector *pip_vec(){
  //   return pip;
  // }
  // gsl_vector *beta_vec(){
  //   return beta;
  // }
  // template<class Archive>
  // void serialize(Archive & archive){
  //   archive(beta,pip,prior,o_prior,o_beta);
  // }





  // Rcpp::List export_object() const{
  //   using namespace Rcpp;

  //   return List::create(_["beta"]=beta,
  // 		 _["pip"]=pip,
  // 		 _["prior"]=prior,
  // 		 _["o_prior"]=o_prior,
  // 		 _["o_beta"]=o_beta);
  // }

};








class controller {
  splitter split;
  
public:

  RcppGSL::vector<double> saved_beta_vec;// = gsl_vector_calloc(ncoef);
  RcppGSL::vector<double> saved_prior_vec;// = gsl_vector_calloc(p);
  RcppGSL::vector<double> prior_vec;
  RcppGSL::vector<double> pip_vec;
  RcppGSL::vector<double> beta_vec;


  controller(const splitter &split_,Result_obj &obj,RcppGSL::matrix<int> anno_mat,double EM_thresh_=0.05,double init_pi1_=1e-3,int print_avg_=0):
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
