
mat <- function (A, k) {
  # taken from tensr, no need to import the whole package
  Ak <- t(apply(A, k, "c"))
  if (nrow(Ak) != dim(A)[k]) {
    Ak <- t(Ak)
  }
  return(Ak)
}
amprod <- function (A, M, k) {
# taken from tensr, no need to import the whole package
  K <- length(dim(A))
  AM <- M %*% mat(A, k)
  AMA <- array(AM, dim = c(dim(M)[1], dim(A)[-k]))
  return(aperm(AMA, match(1:K, c(k, (1:K)[-k]))))
}
norm <- function(x){
  #norm of an array
  sqrt(sum(x^2))
}
normalize <- function(x){
  #returns a matrix with normalized columns
  apply(x,2,function(y)  return(y/norm(y)))
}
ginvS <- function(X,tol=1e-10){
  #generalized inverse of non-negative definite matrix X; faster than MASS::ginv
  ei <- eigen(X,symmetric = T)
  r <- sum(ei$values>tol)
  crossprod((ei$values[1:r])^(-1/2)* t(ei$vectors[,1:r]))
}
Styp <- function(x){
#returns useful quantities from a positive definite matrix
  s <- eigen(x,symmetric = T)
  return(list(
    sqr  = s$vectors %*% diag(s$values^0.5) %*% t(s$vectors),
    isqr = s$vectors %*% diag(s$values^-0.5) %*% t(s$vectors),
    inv  = s$vectors %*% diag(s$values^-1) %*% t(s$vectors),
    ldet = sum(log(s$values)),
    norm = norm(s$values),
    orig = x
  ))
}
tprod <- function (A, B, modes = 1:length(B)) {
  # tucker product
  # taken from package mcmcFunc, tho it should be attributed to tensr
  X <- A
  for (k in modes) {
    X <- amprod(X, B[[k]], k)
  }
  X
}
matmult_kdim <- function(A,S,nam){
  # This multiplies the matrix S to the mode of the tensor A named nam, but keeps the dimension name
  # S must be square, and its second dimension is the one contracted (as in the kth mode product)
  A <- mark(A,i=which(nam == names(dim(A))))
  Sdim <- dim(S)
  names(Sdim) <- c(nam,paste0(nam,"'"))
  A %e% as.tensor(S,dims = Sdim)
}
tarray <- function(A){
  #transforms a tensor to an array
  array(A,dim(A))
}
sqmode <- function(X,k,t = T){
  # faster alternative to tcrossprod(mat(X,k))
  ms <- dim(X)
  X <- aperm(X,c(k,(1:length(ms))[-k]))
  dim(X) <- c(ms[k],prod(ms[-k]))
  if(t)  return(tcrossprod(X))
  return(crossprod(X))
}
unorm_solver <- function(Xy,XX){
  #given Xy and XX the LS estimator is bhat = XX^-Xy
  #here we perform LS but with contraint ||bhat|| = 1
  XXeig <- eigen(XX,symmetric = T)
  z <- t(XXeig$vectors) %*% t(t(Xy))
  d <- XXeig$values
  objB <- Vectorize(function(l) log(sum((z/(l-d))^2)))
  # I need the root of objB, here I will get creative:
  # indeterminacy at lambda = eigenvalue for positive value, get the largest one: dmax
  eigEV <- objB(d+1e-3)
  dmax <- d[which(eigEV==max(eigEV))]
  if(length(dmax)!=1) dmax <- dmax[1]
  # get farther and farther from dmax by factors of 10 until you get to something negative: dmin
  dmin <- abs(dmax)
  while(objB(dmin)>0){
    dmin <- 10*dmin
  }
  lam <- uniroot(objB,c(dmax,dmin))$root
  XXeig$vectors %*%(z/(lam-d))
}
