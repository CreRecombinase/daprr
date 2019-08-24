library(testthat)
library(daprcpptest)
context("utils")


test_that("works like cmd", {

  gf <- system.file2("gwas_z_t.txt.gz", package = "daprcpptest")
  af <- system.file("gwas_anno_t.txt.gz", package = "daprcpptest")

  ret <- tidyr::unnest(daprcpp::run_torus_cmd(gf = gf, af = af)$df,cols=c(data))
  cmd <- daprcpptest:::gen_torus_cmd(gf = gf, af = af)
  cmd <- cmd[!cmd%in%c("--load_zval","-est")]
  bt <- daprcpptest:::dap_torus(cmd)
  btl <- tidyr::unnest(daprcpptest:::parse_torus_s(bt$s,bt$lik),cols=c(data))
  testthat::expect_equal(unclass(ret),unclass(btl),tolerance=1e-5)

})




test_that("daprcpp and torus give similar results", {

    snpnum <- 10000
    p <- 3
    anno_r <- 250
    set.seed(100)
    ln <- letters[1:p]
    names(ln) <- paste0(letters[1:p],"_d")
    ct <- paste0(letters,"_d")[1:p]
    fl <- as.list(rep(0,p))
    names(fl) <- ct
    anno_df <- tibble::tibble(SNP = sample(snpnum, anno_r, replace = T),
                              feature = sample(paste0(letters,"_d")[1:p],anno_r, replace = T))
    anno_df <- dplyr::distinct(anno_df,SNP,feature)
    w_anno_df <- dplyr::mutate(anno_df,value=1L)
    w_anno_df <- tidyr::spread(data=w_anno_df,key="feature",value="value",fill=0L)
  
    gw_df <- tibble::tibble(SNP = 1:snpnum, locus = cut(1:snpnum, breaks = 3L, labels = F),`z-val`=rnorm(n = snpnum))
    #mymay <- daprcpptest:::make_matrix(snpnum,anno_df)
    
    anno_mat <- data.matrix(dplyr::select(tidyr::replace_na(dplyr::left_join(dplyr::select(gw_df,SNP),w_anno_df),replace=fl),-SNP))
    colnames(anno_mat) <- ln[colnames(anno_mat)]
    
    #colnames(mymay$annomat) <- mymay$names

    library(Matrix)
    spX <- as(anno_mat,"sparseMatrix")
    y <- runif(snpnum)
    ym <- cbind(y,1-y)
    tr <- glmnet::cv.glmnet(spX,ym,family="binomial")
     (elm_ret <- daprcpptest:::elastic_donut(locus = gw_df$locus,z = gw_df$`z-val`,X = anno_mat)[-2])
    (selm_ret <- daprcpptest:::elastic_donut_sp(locus = gw_df$locus,z = gw_df$`z-val`,X = spX)[-2])
    # 
    # elm_ret <- daprcpptest:::dap_donut(locus = gw_df$locus,z = gw_df$`z-val`,X = mymay$annomat)
    # 
    tgf <- tempfile(fileext="tsv.gz")
    readr::write_delim(gw_df,path = tgf)
    taf <- tempfile(fileext="txt.gz")
    readr::write_delim(w_anno_df,path = taf)
    
    
    ret <- tidyr::unnest(daprcpp::run_torus_cmd(gf = tgf, af = taf)$df,cols=c(data))
    cmd <- daprcpptest:::gen_torus_cmd(gf = tgf, af = taf)
    cmd <- cmd[!cmd%in%c("--load_zval","-est")]
    bt <- daprcpptest:::dap_torus(cmd)
    btl <- tidyr::unnest(daprcpptest:::parse_torus_s(bt$s,bt$lik),cols=c(data))
    testthat::expect_equal(unclass(ret),unclass(btl),tolerance=1e-5)
  
  
})


test_that("we can give lambda to both",{

  
  snpnum <- 10000
  p <- 3
  anno_r <- 250
  set.seed(100)
  ln <- letters[1:p]
  names(ln) <- paste0(letters[1:p],"_d")
  ct <- paste0(letters,"_d")[1:p]
  fl <- as.list(rep(0,p))
  names(fl) <- ct
  anno_df <- tibble::tibble(SNP = sample(snpnum, anno_r, replace = T),
                            feature = sample(paste0(letters,"_d")[1:p],anno_r, replace = T))
  anno_df <- dplyr::distinct(anno_df,SNP,feature)
  w_anno_df <- dplyr::mutate(anno_df,value=1L)
  w_anno_df <- tidyr::spread(data=w_anno_df,key="feature",value="value",fill=0L)
  
  gw_df <- tibble::tibble(SNP = 1:snpnum, locus = cut(1:snpnum, breaks = 3L, labels = F),`z-val`=rnorm(n = snpnum))

  
  #mymay <- daprcpptest:::make_matrix(snpnum,anno_df)
  
  anno_mat <- data.matrix(dplyr::select(tidyr::replace_na(dplyr::left_join(dplyr::select(gw_df,SNP),w_anno_df),replace=fl),-SNP))
  colnames(anno_mat) <- ln[colnames(anno_mat)]
  
  #colnames(mymay$annomat) <- mymay$names
  
  library(Matrix)
  spX <- as(anno_mat,"sparseMatrix")
 
  
  
  
  (elm_ret <- daprcpptest:::elastic_donut(locus = gw_df$locus,z = gw_df$`z-val`,X = anno_mat,alpha = 0,lambda = 0.5)[-2])
  (selm_ret <- daprcpptest:::elastic_donut_sp(locus = gw_df$locus,z = gw_df$`z-val`,X = spX,alpha=0,lambda=0.5)[-2])
 
  tgf <- tempfile(fileext="tsv.gz")
  readr::write_delim(gw_df,path = tgf)
  taf <- tempfile(fileext="txt.gz")
  readr::write_delim(w_anno_df,path = taf)
  
  
  ret <- tidyr::unnest(daprcpp::run_torus_cmd(gf = tgf, af = taf,l2 = 0.5)$df,cols=c(data))
  cmd <- daprcpptest:::gen_torus_cmd(gf = tgf, af = taf)
  cmd <- cmd[!cmd%in%c("--load_zval","-est")]
  bt <- daprcpptest:::dap_torus(cmd)
  btl <- tidyr::unnest(daprcpptest:::parse_torus_s(bt$s,bt$lik),cols=c(data))
  testthat::expect_equal(unclass(ret),unclass(btl),tolerance=1e-5)
  
  
  
  
})


test_that("prediction works as it should",{
  
  n <- 100
  p <- 2
  
  Xi <- matrix(sample(0L:1L,n*p,replace=T),n,p)
  Xd <- matrix(as.numeric(Xi),n,p)
  colnames(Xd) <- letters[1:p]
  Xdf <- tibble::as_tibble(Xd) 
  Xdf <- dplyr::mutate(Xdf,y=runif(n))
  beta <- rnorm(p+1)
  tfit <- glm(y~offset(beta[1]*rep(1,n))+offset(beta[2]*a)+offset(beta[3]*b)+0,data=Xdf,family=quasibinomial(link="logit"))
  #tfit <- glm(y~offset(beta[1]*rep(1,n))+offset(beta[2]*a)+0,data=Xdf,family=quasibinomial(link="logit"))
  #tfit <- glm(y~offset(beta[1]*rep(1,n))+offset(beta[2]*a)+offset(beta[3]*b)+offset(beta[4]*c)+0,data=Xdf,family=quasibinomial(link="logit"))
  ty <- c(predict.glm(tfit,type="response"))
  
  yd <- daprcpptest:::predict_dap(beta,Xi)
  yg <- daprcpptest:::predict_glmnet(beta,Xd)
  testthat::expect_equal(ty,yg,check.attributes=F)  
  testthat::expect_equal(ty,yd,check.attributes=F)  
  
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



