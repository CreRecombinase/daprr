/*
 * This file uses the Catch unit testing library, alongside
 * testthat's simple bindings, to test a C++ function.
 *
 * For your own packages, ensure that your test files are
 * placed within the `src/` folder, and that you include
 * `LinkingTo: testthat` within your DESCRIPTION file.
 */

// All test files should include the <testthat.h>
// header file.
#include "elasticdonut.hpp"
#include <testthat.h>
#include "dap_classdef.hpp"

#include "dap_logistic.hpp"
#include <RcppEigen.h>
#include "glmnet.hpp"
#include "io.hpp"
#include <random>
// Normally this would be a function from your package's
// compiled library -- you might instead just include a header
// file providing the definition, and let R CMD INSTALL
// handle building and linking.


// Initialize a unit test context. This is similar to how you
// might begin an R test file with 'context()', expect the
// associated context should be wrapped in braced.

Rcpp::Function system_file("system.file");
const std::string gwas_file = Rcpp::as<std::string>(system_file("gwas_z_t.txt.gz",Rcpp::Named("package")="daprcpp"));
const std::string annotation_file = Rcpp::as<std::string>(system_file("gwas_anno_t.txt.gz",Rcpp::Named("package")="daprcpp"));


template<typename T,int Row=Eigen::Dynamic,int Col=Eigen::Dynamic>
inline gsl_matrix_int*	copy_matrix(const Eigen::Map<Eigen::Matrix<T,Row,Col>> X){

  const size_t n = X.rows();
  const size_t p = X.cols();

  gsl_matrix_int* retX = gsl_matrix_int_calloc(n,p);
  Eigen::Map<Eigen::Matrix<int,Row,Col,Eigen::RowMajor>> wrap_retX(retX->data,n,p);
  wrap_retX=X.template cast<int>();
  return(retX);
}




template<class T>
typename std::enable_if<!std::numeric_limits<T>::is_integer, bool>::type
    almost_equal(T x, T y, int ulp=5)
{
  bool ret = std::abs(x-y) <= std::numeric_limits<T>::epsilon() * std::abs(x+y) * ulp
        // unless the result is subnormal
        || std::abs(x-y) < std::numeric_limits<T>::min();
  if(!ret){
    //    Rcpp::Rcerr<<x<<"!="<<y<<std::endl;
  }
    // the machine epsilon has to be scaled to the magnitude of the values used
    // and multiplied by the desired precision in ULPs (units in the last place)
  return ret;
}




template<typename T>
gsl_vector *view_span(T data) {

  const size_t n = data.size();
  gsl_vector *ret = new gsl_vector;
  ret->size = n;
  ret->stride = 1;
  ret->data = data.data();
  ret->block = nullptr;
  ret->owner = 0;
  return (ret);
}

context("E step") {

torus::controller con(gwas_file,annotation_file);

FileZscoreParser zsp(gwas_file.c_str());


 test_that("single step"){
   auto tLoc = con.locVec.begin();
   auto pip = tLoc->get_pip();
   auto prior = tLoc->get_prior();
   auto BFd = tLoc->get_BF();
   tLoc->EM_update();
   const auto log10_lik = tLoc->log10_lik;
   double* BFm = nullptr;
   auto comp_log10_lik = elasticdonut::E_step(toMapXd(pip),toMapXd(BFd),toMapXd(prior),BFm);

   expect_true(almost_equal(comp_log10_lik,log10_lik));



   int i=0;
   bool all_eq = std::all_of(tLoc->snpVec.begin(), tLoc->snpVec.end(),
                             [&i, &pip, &tLoc](auto &snp) mutable {
                               double opip =
                                   gsl_vector_get(tLoc->pip_vec, snp.index);
                               double cpip = pip[i++];
                               return (almost_equal(opip, cpip));
                             });
   expect_true(all_eq);
 }

 test_that("each locus"){

   auto locus = con.get_region_id();
   int *locb =&(*locus.begin());
   auto splt = make_splitter(Rcpp::wrap(locus));
   auto tLoc = con.locVec.begin();

   auto BFd = con.get_BF();
   const size_t p = BFd.size();

   // elasticdonut::ParameterData param(Rcpp::wrap(con.beta_vec),Rcpp::wrap(con.pip_vec),Rcpp::wrap(con.prior_vec),Rcpp::wrap(con.dvar_name_vec));
   // elasticdonut::ParameterBuffer buff(param,splt);
   GroupedView BF_v(splt.split_view(&BFd.front()),p);
   std::vector<double> prior(p,0.0);
   std::vector<double> pip(p,0.0);
   auto	prior_view = BF_v.copy_view(prior);
   auto	pip_view = BF_v.copy_view(pip);
   elasticdonut::SumStatRegion sumstats(BF_v);
   double result_a = std::accumulate(con.locVec.begin(), con.locVec.end(), 0.0,
                                     [](const double &t_lik, const torus::Locus &a) {
				       torus::Locus newa = a;
                                       newa.EM_update();
				       double tloc =newa.log10_lik;
				       double new_sum =	(t_lik + newa.log10_lik);

                                       return new_sum;
                                     });

   auto prior_c = con.prior_vec;
   auto pip_c = con.pip_vec;

   double result_b = sumstats.E_steps(prior_view.r_view);
   expect_true(almost_equal(result_a,result_b));
   auto	pip_d =	sumstats.pip();
   bool all_eq_pip = toMapXd(pip).isApprox(pip_d);
   expect_true(all_eq_pip);
 }
}

// // context("M step") {


// //   test_that("Logistic"){

// //     const size_t n =100;
// //     const size_t p = 3;

// //     Eigen::VectorXd y= (Eigen::ArrayXd::Random(n)+1)/2;
// //     gsl::span<double> s_y(y.data(),y.size());

// //     Eigen::MatrixXi X(n,p);
// //     Eigen::Map<Eigen::MatrixXi> mX(X.data(),n,p);

// //     std::random_device rd;  //Will be used to obtain a seed for the random number engine
// //     std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
// //     std::uniform_int_distribution<> dis(0, 1);

// //     X = X.unaryExpr([&dis, &gen](const int t) { return (dis(gen)); });
// //     auto gsl_X=copy_matrix(mX);
// //     auto nlev = gsl_vector_int_calloc(p);

// //     for(int j=0; j <p; j++){
// //       std::map<int, int> rcd;
// //       for (int i = 0; i < n; i++) {
// //         int val = gsl_matrix_int_get(gsl_X, i, j);
// //         rcd[val] = 1;
// //       }
// //       gsl_vector_int_set(nlev, j, rcd.size());
// //     }

// //     std::vector<double> beta_a(p + 1);
// //     std::vector<double> beta_b(p + 1);

// //     torus::logistic_mixed_fit(view_span(beta_a), gsl_X, nlev, nullptr,
// //                               view_span(s_y), 0, 0);
// //     Dap_logit logistic(gsl_X, 0, 0);
// //     logistic.fit(s_y);
// //     expect_false(almost_equal_span(beta_a, beta_b));
// //     logistic.read_coeffs(beta_b);

// //     expect_true(almost_equal_span(beta_a, beta_b));
// //   }
// //   test_that("Prediction"){

// //     const size_t n =100;
// //     const size_t p = 3;

// //     Eigen::VectorXd y= (Eigen::ArrayXd::Random(n)+1)/2;
// //     gsl::span<double> s_y(y.data(),y.size());

// //     Eigen::MatrixXi X(n,p);

// //     const Eigen::Map<Eigen::MatrixXi> mX(X.data(),n,p);

// //     std::random_device rd;  //Will be used to obtain a seed for the random number engine
// //     std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
// //     std::uniform_int_distribution<> dis(0, 1);

// //     X = X.unaryExpr([&dis, &gen](const int t) { return (dis(gen)); });
// //     Eigen::MatrixXd Xd = X.cast<double>();
// //     Eigen::Map<Eigen::MatrixXd> mXd(Xd.data(),n,p);
// //     auto gsl_X=copy_matrix(mX);
// //     std::vector<double> beta(p + 1);

// //     std::vector<double> y_a(n);
// //     std::vector<double> y_b(n);


// //     Dap_logit logistic(gsl_X, 0, 0);

// //     logistic.fit(s_y);
// //     logistic.read_coeffs(beta);
// //     elasticdonut::Lognet<Eigen::MatrixXd> ln(mXd);
// //     ln.predict(beta,y_a);
// //     logistic.predict(beta,y_b);

// //     expect_true(almost_equal_span(y_a, y_b));
// //   }



// }

// auto line_t =	zsp.getline();
// test_that("parsing zscore file works correctly") {
//   expect_true(std::get<0>(*line_t) == "1:701835:T:C");
// }

  






