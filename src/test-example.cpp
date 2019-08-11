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
#include <testthat.h>
#include "dap_classdef.hpp"
#include "classdef.hpp"
#include <RcppEigen.h>
#include "io.hpp"

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


template<class T>
typename std::enable_if<!std::numeric_limits<T>::is_integer, bool>::type
    almost_equal(T x, T y, int ulp=3)
{
    // the machine epsilon has to be scaled to the magnitude of the values used
    // and multiplied by the desired precision in ULPs (units in the last place)
    return std::abs(x-y) <= std::numeric_limits<T>::epsilon() * std::abs(x+y) * ulp
        // unless the result is subnormal
        || std::abs(x-y) < std::numeric_limits<T>::min();
}

context("E step") {

torus::controller con(gwas_file,annotation_file);

FileZscoreParser zsp(gwas_file.c_str());
FileAnnotationParser ap(con.snp_hash,annotation_file.c_str());


 test_that("single step"){
   auto tLoc = con.locVec.begin();
   auto pip = tLoc->get_pip();
   auto prior = tLoc->get_prior();
   auto BF = tLoc->get_BF();
   tLoc->EM_update();
   const auto log10_lik = tLoc->log10_lik;
   std::optional<double> BFm = std::nullopt;
   auto comp_log10_lik = elasticdonut::E_step(pip,BF,prior,BFm);

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

 }
}

  // auto line_t =	zsp.getline();
  // test_that("parsing zscore file works correctly") {
  //   expect_true(std::get<0>(*line_t) == "1:701835:T:C");
  // }

  








context("comparing donut and torus") {

  using namespace Rcpp;


  // IntegerVector locus_id  = wrap(read_vec<int>("locus_id.txt.gz"));
  // NumericVector z_hat  = wrap(read_vec<double>("z_hat.txt.gz"));

  // auto anno_df = DataFrame::create(_["SNP"]=wrap(read_vec<int>("row_id.txt.gz")),
  // 				   _["feature"]=wrap(read_vec<std::string>("col_id.txt.gz")));
  // SparseDF spdf( anno_df, p,"SNP","feature");
  // auto anno_mat = spdf.getMat();

  // const size_t p = locus_id.size();

  // torus::controller cont = cont_wrapper(gwas_file,annotation_file);
  // auto split = donut::make_splitter(locus_id.begin(),locus_id.end());

  // donut::Result_obj res(locus_id.size(),anno_mat.ncol()+1);
  // donut::controller do_cont(split,res);

  // cont.init_params();

  test_that("two plus two equals four") {
    expect_true(4 == 4);
  }

}
