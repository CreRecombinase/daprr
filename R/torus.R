torus_est <- function(zhat,region_id,annot){
  
  n <- 30
  p <- 3
  y <- runif(n)
  Ym <- cbind(1-y,y)    
  X <- matrix(runif(n*p),n,p)  
  
  
  result <- glmnet::glmnet(X,Ym,"binomial",alpha=1,lambda=0.0)
  
}
