#include "dataview.hpp"
#include <numeric>






size_t total_size(const SplitView& r_view){
  return std::accumulate(r_view.begin(),r_view.end(),0,[](size_t s,auto &sp){
							 return s+sp.size();
						       });
}


GroupedView::GroupedView(SplitView r_view_, size_t p_)
  : r_view(r_view_), p(p_),nr(r_view.size()), d_view(r_view.begin()->data(), p) {}
GroupedView::GroupedView(SplitView r_view_)
  : r_view(r_view_), p(total_size(r_view)), nr(r_view.size()),
        d_view(r_view.begin()->data(), p) {}

GroupedView GroupedView::copy_view(Rcpp::NumericVector o_data) const {
  gsl::span<double> sp(&(*o_data.begin()), o_data.size());
  return (copy_view(sp));
}

GroupedView GroupedView::copy_view(gsl::span<double>	o_data) const {
  SplitView tr_view;
  tr_view.reserve(this->nr);
  size_t offset=0;
  // auto ob=o_data.data();
  std::transform(r_view.cbegin(), r_view.cend(), std::back_inserter(tr_view),
                 [&offset, &o_data](const auto sp) mutable {
                   auto ret = o_data.subspan(offset, sp.size());
                   //                     auto ret = gsl::span<double>(ob,
                   //                     sp.size());
                   offset += ret.size();
                   //                     ob = ret.
                   return (ret);
                 });
  return GroupedView(std::move(tr_view), this->p);
}

// splitter make_splitter( v){


// }

splitter make_splitter(gsl::span<int> v){
  std::vector<size_t> ret_v;
  auto v_begin = v.begin();
  auto v_end = v.end();

    auto t_b = v_begin;
    auto t_i = t_b;
    ret_v.push_back(0);
    while (t_i != v_end) {
      t_i = std::upper_bound(t_i, v_end, *t_i);
      ret_v.push_back(std::distance(t_b, t_i));
    }
    return (splitter(std::move(ret_v)));
}
