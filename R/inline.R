inlineCxxPlugin <- function(...) {
  plugin <-
    Rcpp::Rcpp.plugin.maker(
      include.before = "#include <gsl/span>",
      package        = "daprcpp"
    )
  settings <- plugin()
  settings$env$PKG_CPPFLAGS <- paste("-I../inst/include")
  settings$env$USE_CXX17 <- "yes"
  settings
}
