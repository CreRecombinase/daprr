#pragma once
#include "glmnet.hpp"



namespace elasticdonut {
using SplitView = std::vector< gsl::span<double> >;



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

  std::vector<double> make_BF(Rcpp::NumericVector z_hat);
  std::vector<double> make_BF(const gsl::span<double> z_hat);

// SumStatRegion is a non-owning view of the summary statistics
class SumStatRegion {

  public:
  const GroupedView BF;
  Rcpp::NumericVector p_vec;
  private:
  GroupedView p_view;
  mutable std::vector<double*> BF_max;


    //  double E_step(
  public:
    SumStatRegion(const GroupedView BF_);
    double E_steps(const SplitView &r_prior);
    size_t size() const;
    gsl::span<double> pip();
  };




  class	ElasticDonut{
  public:
    SumStatRegion sumstats;
    Net& logistic;
    std::vector<double> prior_vec;
    std::vector<double> null_prior_vec;
    GroupedView prior_view;
    std::vector< std::string> names;
    std::vector<double> beta;
    std::vector<double> sd;
    double curr_log10_lik;

    const double EM_thresh;
    ElasticDonut(GroupedView BF_v, Net &logn,
                 const double prior_init = 1e-3, double EM_thresh_ = 0.05);
    double fit(bool keep_prior=true);
    double fine_optimize_beta(gsl::span<double> new_beta,const size_t index,GroupedView &null_prior,double &log10_lik);
  };
  double E_step(const gsl::span<double> pip, const gsl::span<double> BF,
		const gsl::span<double> prior, double* BF_max);


    } // namespace elasticdonut
