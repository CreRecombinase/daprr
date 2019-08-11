library(testthat)
library(daprcpptest)
context("utils")


test_that("works like cmd", {

  gf <- system.file("gwas_z_t.txt.gz", package = "daprcpp")
  af <- system.file("gwas_anno_t.txt.gz", package = "daprcpp")

  ret <- tidyr::unnest(daprcpp::run_torus_cmd(gf = gf, af = af)$df,cols=c(data))
  cmd <- daprcpptest::gen_torus_cmd(gf = gf, af = af)
  cmd <- cmd[!cmd%in%c("--load_zval","-est")]
  bt <- daprcpptest::dap_torus(cmd)
  btl <- tidyr::unnest(parse_torus_s(bt$s,bt$lik),cols=c(data))
  testthat::expect_equal(unclass(ret),unclass(btl),tolerance=1e-5)

})




test_that("daprcpp and torus give similar results", {
  
    snpnum <- 10000 
    p <- 3
    
    anno_r <- 100
    anno_df <- tibble::tibble(SNP = sample(snpnum, anno_r, replace = T),
                              feature = sample(paste0(letters,"_d")[1:p],anno_r, replace = T))
    
    
    gw_df <- tibble::tibble(SNP = 1:snpnum, region_id = cut(1:snpnum, breaks = 3L, labels = F),`z-hat`=rnorm(n = snpnum))
    
  
  
  
})
# 
# 
# test_that("torus_df and torus give comparable results", {
#   
# 
#   anno_df <- tibble::tibble(SNP = c(1, 3,4),feature = factor(c("aa","aa","ab")))
#   annomat_l <- make_matrix(p, anno_df)
#   y <- c(0.684902017234619, 2.44814199215588, 0.116413627382255, 1.4847548363253
# )
#   true_res_a <- structure(list(term = structure(c(3L, 1L, 2L), .Label = c("aa",
# "ab", "Intercept"), class = "factor"), estimate = c(-4.18133509474218,
# -4.40803025445066, -2.97845638309911), high = c(57946.1665030528,
# 41.5930822848575, 83.1053718780682), low = c(-57954.5291732423,
# -50.4091427937588, -89.0622846442664), lik = c(0.00463521179289106,
# 0.00463521179289106, 0.00463521179289106)), class = "data.frame", row.names = c(NA,
# -3L))
#   
#   res_a <- torus(locus_id = rep(1L, p),z_hat = y, anno_mat = annomat_l$annomat, names = annomat_l$names, use_glmnet = F)
#   res_b <- torus_df(locus_id = rep(1L, p),z_hat = y, anno_df = anno_df, use_glmnet = F)
#   expect_equal(res_a, res_b)
#   expect_equal(true_res_a, res_b$est, tolerance = 1e-2)
# })
# 
# 
# test_that("dap logistic works like glmnet in the absence of regularization", {
#   
#   n <- 30
#   p <- 4
#   y <- runif(n)
#   
#   X <- matrix(sample(2, n*p, replace = T)-1, n,p)
#   
#   lr <- logit_torus(X, y)
#   lrg <- logit_glmnet(X, y)
#   expect_equal(lr, lrg, tolerance = 1e-4)
#   
#   
# })



