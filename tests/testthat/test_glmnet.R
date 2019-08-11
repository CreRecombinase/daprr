context("elastic-net")

# 
# test_that("glmnet really shrinks",{
#   
#   
#   library(glmnet)
#   library(daprcpp)
#   n <- 5000
#   p <- 3 
#   
#   X <- matrix(rnorm(n*p),n,p)
#   beta <- c(0.2,1.5,0,-4)
#   xb <- cbind(1,X)%*%beta
#   
#   oy <- 1/(1+exp(-xb))
#   y <- exp(xb)/(exp(xb)+1)
#   ym <- cbind(1-y,y)
#   
#   mg <- glmnet::glmnet(X,ym,family="binomial",alpha=1)
#   acoef <- rbind(mg$a0,as.matrix(mg$beta))
#   lam <- mg$lambda
#   
#   c2=logit_glmnet(X,y =y, alpha=1,lambda = lam)
#     
#   
# })
