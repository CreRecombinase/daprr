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
  
 
    })

