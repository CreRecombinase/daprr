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

// Normally this would be a function from your package's
// compiled library -- you might instead just include a header
// file providing the definition, and let R CMD INSTALL
// handle building and linking.


// Initialize a unit test context. This is similar to how you
// might begin an R test file with 'context()', expect the
// associated context should be wrapped in braced.
/*
context("comparing donut and torus") {

  using namespace Rcpp;
  const std::string gwas_file = "gwas.txt.gz";
  const std::string annotation_file = "annotation.txt.gz";

  IntegerVector locus_id  = wrap(read_vec<int>("locus_id.txt.gz"));
  NumericVector z_hat  = wrap(read_vec<double>("z_hat.txt.gz"));

  auto anno_df = DataFrame::create(_["SNP"]=wrap(read_vec<int>("row_id.txt.gz")),
				   _["feature"]=wrap(read_vec<std::string>("col_id.txt.gz")));
  SparseDF spdf( anno_df, p,"SNP","feature");
  auto anno_mat = spdf.getMat();

  const size_t p = locus_id.size();

  torus::controller cont = cont_wrapper(gwas_file,annotation_file);
  auto split = donut::make_splitter(locus_id.begin(),locus_id.end());

  donut::Result_obj res(locus_id.size(),anno_mat.ncol()+1);
  donut::controller do_cont(split,res);

  cont.init_params();

// The format for specifying tests is similar to that of
  // testthat's R functions. Use 'test_that()' to define a
  // unit test, and use 'expect_true()' and 'expect_false()'
  // to test the desired conditions.
test_that("two plus two equals four") {
    expect_true(4 == 4);
  }

}
*/
