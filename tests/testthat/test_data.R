context("preprocessing")


test_that("annotation matrix is generated correctly",{
  
  p <- 4
  anno_df <- tibble::tibble(SNP=as.integer(c(1,3,4)),feature=c("aa","aa","aa"))  
  annomat_l <- make_matrix(p,anno_df)
  expect_equal(annomat_l,list(annomat=matrix(c(1,0,1,1),4,1),
                              names=c("aa")))
  
  
  p <- 4
  anno_df <- tibble::tibble(SNP=c(1L,3L,4L),feature=c("aa","aa","aa"))  
  annomat_l <- make_matrix(p,anno_df)
  expect_equal(annomat_l,list(annomat=matrix(c(1,0,1,1),4,1),
                              names=c("aa")))
  
  p <- 4
  anno_df <- tibble::tibble(SNP=c(1L,3L,4L),feature=c("aa","aa","aa"))  
  annomat_l <- make_spmatrix_df(p,anno_df)
  spm <- sparseMatrix(i =anno_df$SNP,j = as.integer(factor(anno_df$feature)),x=rep(1.0,nrow(anno_df)) )
  
  expect_equal(annomat_l,list(annomat=spm,
                              names=c("aa")))
  
  
  
  anno_df <- tibble::tibble(SNP=c(1,3,4,1),feature=c("aa","aa","aa","ab"))  
  annomat_l <- make_matrix(p,anno_df)
  expect_equal(annomat_l,list(annomat=matrix(c(c(1,0,1,1),
                                               c(1,0,0,0)),4,2),
                              names=c("aa","ab")))
  
  
  anno_df <- tibble::tibble(SNP=c(1,3,4),feature=c("aa","aa","ab"))  
  annomat_l <- make_matrix(p,anno_df)
  expect_equal(annomat_l,list(annomat=matrix(c(1,0,1,0,
                                               0,0,0,1),4,2),
                              names=c("aa","ab")))
  
    anno_df <- tibble::tibble(SNP=c(1,3,4),feature=factor(c("aa","aa","ab")))  
  annomat_l <- make_matrix(p,anno_df)
  expect_equal(annomat_l,list(annomat=matrix(c(1,0,1,0,
                                               0,0,0,1),4,2),
                              names=c("aa","ab")))
  
    })

