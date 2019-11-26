#pragma once

#include <RcppEigen.h>

using SplitView = std::vector< Eigen::Map<Eigen::ArrayXd> >;

class GroupedView{
public:
  SplitView r_view;
  const size_t p;
  const size_t nr;
  Eigen::Map<Eigen::ArrayXd> d_view;
  GroupedView(SplitView r_view_, size_t p_);
  GroupedView(SplitView r_view_);
  GroupedView copy_view(Eigen::Map<Eigen::ArrayXd>	o_data) const;
  GroupedView copy_view(std::vector<double> &o_data) const;

  GroupedView copy_view(Rcpp::NumericVector	o_data) const;
};


size_t total_size(const SplitView& r_view);

class splitter{
  std::vector<size_t> ret_v;
public:
  splitter(Rcpp::IntegerVector ret_v_):ret_v(Rcpp::as<std::vector<size_t> >(ret_v_)){}
  splitter(std::vector<size_t> &&ret_v_):ret_v(ret_v_){}
  Eigen::Map<Eigen::ArrayXd> split_range(double* data_pt,int idx)const {
    int beg=ret_v[idx];
    int	ret_s =	ret_v[idx+1]-ret_v[idx];
    return(Eigen::Map<Eigen::ArrayXd>(data_pt+beg,ret_s));
  }
  const int num_regions()const {
    return(static_cast<int>(ret_v.size()-1));
  }
  std::vector<Eigen::Map<Eigen::ArrayXd> > split_view(double* data_pt) const{
    const int tot_r= num_regions();
    std::vector< Eigen::Map<Eigen::ArrayXd> > groupv;
    groupv.reserve(tot_r);

    for (int i = 0; i < tot_r; i++) {
      groupv.emplace_back(split_range(data_pt, i));
    }
    return groupv;
  }
  GroupedView group_view(double* data_pt) const{
    const int tot_r= num_regions();
    std::vector< Eigen::Map<Eigen::ArrayXd> > groupv;
    groupv.reserve(tot_r);
    size_t p=0;
    for (int i = 0; i < tot_r; i++) {
      groupv.emplace_back(split_range(data_pt, i));
      p+=groupv.back().size();
    }
    return GroupedView(groupv,p);

  }

  std::vector<Eigen::Map<Eigen::ArrayXd>> split_view(Eigen::Map<Eigen::ArrayXd> data) const {
    return(this->split_view(data.data()));
  }
  Rcpp::IntegerVector export_rle() const{
    return Rcpp::wrap(ret_v);
  }




};


splitter make_splitter(Rcpp::IntegerVector v);
