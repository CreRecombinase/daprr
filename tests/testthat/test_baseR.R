library(testthat)
library(daprcpp)
context("utils")


# test_that("R quasibinomial logit regression works like dap",{
#   
#   n <- 30
#   p <- 4
#   y <- runif(n)
#   X <- matrix(sample(2,n*p,replace=T)-1,n,p)
#   sX <- as(X,"CsparseMatrix")
#   lr <- logit_torus(X,y)
#   names(lr) <- c("(Intercept)",paste0("X",1:p))
#   tglm <- glm(y~sX,family=quasibinomial(link="logit"))
#   lrR <- glm.fit(x = X,y = y,family=quasibinomial(link="logit"),control = glm.control(trace=TRUE))
#   $coefficients
#   
#   expect_equal(lr,lrR)
#   
# })


test_that("torus_df and torus give comparable results",{
  
  p <- 4
  anno_df <- tibble::tibble(SNP=c(1,3,4),feature=factor(c("aa","aa","ab")))  
  annomat_l <- make_matrix(p,anno_df)
  y <- c(0.684902017234619, 2.44814199215588, 0.116413627382255, 1.4847548363253
)
  true_res_a <- structure(list(term = structure(c(3L, 1L, 2L), .Label = c("aa", 
"ab", "Intercept"), class = "factor"), estimate = c(-4.18133509474218, 
-4.40803025445066, -2.97845638309911), high = c(57946.1665030528, 
41.5930822848575, 83.1053718780682), low = c(-57954.5291732423, 
-50.4091427937588, -89.0622846442664), lik = c(0.00463521179289106, 
0.00463521179289106, 0.00463521179289106)), class = "data.frame", row.names = c(NA, 
-3L))
  
  res_a <- torus(locus_id = rep(1L,p),z_hat = y,anno_mat = annomat_l$annomat,names = annomat_l$names)
  res_b <- torus_df(locus_id = rep(1L,p),z_hat = y,anno_df = anno_df)
  expect_equal(res_a,res_b)  
  expect_equal(true_res_a,res_b$est,tolerance=1e-2)
})


test_that("dap logistic works like glmnet in the absence of regularization",{
  
  n <- 30
  p <- 4
  y <- runif(n)
  
  X <- matrix(sample(2,n*p,replace=T)-1,n,p)
  
  lr <- logit_torus(X,y)
  lrg <- logit_glmnet(X,y)
  expect_equal(lr,lrg,tolerance=1e-4)
  
  
})



