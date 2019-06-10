// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppGSL.h>
#include <Rcpp.h>

using namespace Rcpp;

// torus
Rcpp::List torus(Rcpp::IntegerVector locus_id, Rcpp::NumericVector z_hat, RcppGSL::matrix<int> anno_mat, Rcpp::StringVector names, const bool prior);
RcppExport SEXP _daprcpp_torus(SEXP locus_idSEXP, SEXP z_hatSEXP, SEXP anno_matSEXP, SEXP namesSEXP, SEXP priorSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type locus_id(locus_idSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type z_hat(z_hatSEXP);
    Rcpp::traits::input_parameter< RcppGSL::matrix<int> >::type anno_mat(anno_matSEXP);
    Rcpp::traits::input_parameter< Rcpp::StringVector >::type names(namesSEXP);
    Rcpp::traits::input_parameter< const bool >::type prior(priorSEXP);
    rcpp_result_gen = Rcpp::wrap(torus(locus_id, z_hat, anno_mat, names, prior));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_daprcpp_torus", (DL_FUNC) &_daprcpp_torus, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_daprcpp(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
