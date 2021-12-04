


# Tensor-on-Tensor Regression where B has no format
# S has two different structures: diagonal (IndpS) and Kronecker-separable (KronS)
#' @export
ToT_unformat <- function(Yall,Xall,it = 100,corrs = "N",err = 1e-8){
  p <- length(dim(Yall))-1
  l <- length(dim(Xall))-1
  n <- dim(Yall)[p+1]
  hs <- dim(Xall)[-(l+1)]
  ms <- dim(Yall)[-(p+1)]
  mn <- prod(dim(Yall))
  sol <- svd(array(Xall,c(prod(hs),n)))
  B <- amprod(Yall,sol$u %*% diag(1/sol$d) %*% t(sol$v),p+1)
  B <- aperm(B,c(p+1,1:p))
  dim(B) <- c(hs,ms)
  res <- amprod(Yall,diag(n) - tcrossprod(sol$v),p+1)
  fit1 <- MLE_tensnorm(res,it = it,corrs = corrs,err =err)
  fit2 <- apply(res^2,1:p,mean)
  return(list(B=B,IndpS = fit2,KronS = fit1))
}
