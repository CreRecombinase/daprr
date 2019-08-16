#pragma once
#include "gsl/span"
#include "Rcpp.h"

namespace elasticdonut {
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

  std::vector<gsl::span<double>> split_view(gsl::span<double> data) const {
    return(this->split_view(data.data()));
  }
  Rcpp::IntegerVector export_object() const{
    return Rcpp::wrap(ret_v);
  }

};

template<typename T>
inline elasticdonut::splitter make_splitter(T v_begin, T v_end){
  std::vector<size_t> ret_v;
    auto t_b = v_begin;
    auto t_i = t_b;
    ret_v.push_back(0);
    while (t_i != v_end) {
      t_i = std::upper_bound(t_i, v_end, *t_i);
      ret_v.push_back(std::distance(t_b, t_i));
    }
    return (splitter(std::move(ret_v)));
}
}
