#' Tensor-on-Tensor Regression where B has no format
#' and S has two different structures: diagonal (IndpS) and Kronecker-separable (KronS)
#'
#' Tensor-on-tensor regression
#' \eqn{Y_i = < X_i | B > + E_i}\cr
#' with kronecker-separable covariance
#' \eqn{Var(vec(E_i)) = \sigma^2 \otimes_k S_k}
#' The size of the ith tensor covariate \eqn{X_i} is \eqn{h_1 \times .. \times h_l},
#' The size of the ith tensor response \eqn{Y_i} is \eqn{m_1 \times .. \times m_p},
#' the size of the regression coefficient B is of size \eqn{h_1 \times .. \times h_l \times m_1 \times .. \times m_p},
#' and \eqn{i=1,...,n} .\cr
#' the matrix \eqn{M_k} is of size \eqn{m_k \times R},\cr
#' and matrix \eqn{S_k} is positive definite of size \eqn{m_k \times m_k} and is restricted to (1,1) element equal to 1.
#'
#' Convergence is achieved when the relative difference in snorm(B) + snorm(S) is less than err.
#' snorm(Y) here means norm(Y)/sqrt(size of Y) \cr
#' Here the number of dimensions of Y and X must match. If they don't, then one can create dummy dimensions of size 1
#' to have them match using \code{dim(Xall) <- c(1,dim(Xall))}
#'
#' @param Yall Array containing the n tensor responses along the last mode, so that it is of size \eqn{m_1 \times .. \times m_p \times n}.
#' The last dimension must match the last dimension of Xall.
#' @param Xall Array containing the n tensor covariates along the last mode, so that it is of size \eqn{h_1 \times .. \times h_l \times n}.
#' The last dimension must match the last dimension of Yall.
#' @param it maximum number of iterations.
#' @param err relative error used to assess convergence.
#' @param corrs Character vector of size p indicating the types of covariance matrices desired for S_1 ,.., S_p.
#' CHANGE THIS ONE!!
#' @return A list containing the following elements: \cr\cr
#' \code{B} -  the estimated coefficient B. \cr\cr
#' \code{sig2} - the estimate of \eqn{\sigma^2}. \cr\cr
#' \code{covs} - a list with the estimated matrices \eqn{S_1 ,.., S_p}. \cr\cr
#' \code{allconv} - a vector with all the convergence criteria.(????) \cr\cr
#' \code{allik} - a vector with all the loglikelihoods, should be monotone increasing.(????) \cr\cr
#' \code{it} - the number of iterations taken \cr\cr
#' @export
#' @examples
#' # Tensor-on-Tensor Regression on 6x7x8x9 responses and 3x4x5 covariates
#' set.seed(1234)
#' RmsT <- 6:9
#' RhsT <- 3:5
#' Rbt <- 100
#' Rnn <- 2
#' Rcorrs = c("AR(1)","EQC","MA(1)","N")
#' Rsig2t <- 20
#' dat <- diagdat_sim(msT=RmsT,hsT=RhsT,bt=Rbt,nn=Rnn,sig2t=Rsig2t,corrs=Rcorrs)
#' dim(dat$Xall) <- c(1,dim(dat$Xall))
#' fit <- ToT_unformat(Yall = dat$Yall,Xall = dat$Xall,corrs = Rcorrs)
#' par(mfrow = c(1,2))
#' hist(fit$B,main = "estimates in B (true in blue)")
#' abline(v= Rbt,col = "blue",lwd=3)
#' # plot(ts(fit$allik),main = "loglikelihood")
#' par(mfrow = c(1,1))
#' covars <- rbind("true" = c(Rsig2t,sapply(dat$SST,function(x)x[2,1])),
#'                 "est" = c(fit$sig2,sapply(fit$covs,function(x)x[2,1])))
#' colnames(covars) <- c("sig2",paste0("S-",Rcorrs))
#' covars
#' @author Carlos Llosa-Vite, \email{llosacarlos2@@gmail.com}
#' @references \url{https://arxiv.org/abs/2012.10249}
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
  # CHECK: In the comment for this function, it is written #MLE algorithm with covariances set to AR(1) autocorrelations 
  # Then how can it accomodate other type cvariates???
  # It seems that MLE_transform has all the cases. Still check!
  fit2 <- apply(res^2,1:p,mean)
  return(list(B=B,IndpS = fit2,KronS = fit1))
}
