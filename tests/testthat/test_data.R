context("preprocessing")


test_that("annotation matrix is generated correctly",{
  
  p <- 4
  anno_df <- tibble::tibble(SNP=as.integer(c(1,3,4)),feature=c("aa","aa","aa"))  
  annomat <- spread_sparse(anno_df,row_total = p)
  testthat::expect_equal(as.matrix(annomat),matrix(c(1,0,1,1),4,1,dimnames=list(NULL,"aa")))
    })

