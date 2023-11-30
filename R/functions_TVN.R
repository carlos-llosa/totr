
#' Maximum likelihood estimation for the tensor normal distribution
#'
#' Tensor normal distribution 
#' \eqn{Y_i }
#' is Gaussian with mean \eqn{M }
#' and with kronecker-separable covariance
#' \eqn{Var(vec(E_i)) = \sigma^2 \otimes_k S_k}
#' The size of the ith tensor response \eqn{Y_i} is \eqn{m_1 \times .. \times m_p},
#' and \eqn{i=1,...,n} .\cr
#' and matrix \eqn{S_k} is positive definite of size \eqn{m_k \times m_k} and is restricted to (1,1) element equal to 1.
#'
#' Convergence is achieved when the relative difference in loglikelihood is less than err.
#'
#' @param Yall Array containing the n tensor responses along the last mode, so that it is of size \eqn{m_1 \times .. \times m_p \times n}.
#' The last dimension must match the last dimension of Xall.
#' @param it maximum number of iterations.
#' @param err relative error used to assess convergence.
#' @param corrs Character vector of size p indicating the types of covariance matrices desired for \eqn{S_1 ,.., S_p}.
#' Options are "AR(1)", "MA(1)", "ARMA"/"ARMA(p,q)"/"ARMA(p, q)", "EQC"  for
#' AR(1), MA(1), ARMA(p, q) and equivariance correlation matrices, and
#' "N" for general covariance with element (1,1) equal to 1.
#' If corrs is of size 1, then \eqn{S_1 ,.., S_p} will all have the same correlation structure.
#' @param centered If true then the mean \eqn{M }\cr is assumed to be zero.
#' @param initS list containing initial values. init$sig2 contains the initial value of \eqn{\sigma^2}, init$covs is a list of p matrices
#' containing initial values for \eqn{S_1 ,.., S_p}.
#' If init = NULL then the initial values will be set to identity matrix. 
#' @param arma_param A list of size \code{length(dim(Yall))}, each of which contains the
#' ARMA parameter orders (p, q) for that corresponding mode.
#' p is the AR parameter order and q is the MA parameter order
#' If some other mode has some other kind of correlation structure
#' and you still want to specify the ARMA orders,
#' you can input a list of size p with other cases as NULL.
#' The default ARMA order is (1, 1).
#' @return A list containing the following elements: \cr\cr
#' \code{M} -  the estimated mean M.  \cr\cr
#' \code{covs} - a list with the estimated matrices \eqn{S_1 ,.., S_p}. \cr\cr
#' \code{sig2} - the estimate of \eqn{\sigma^2}. \cr\cr
#' \code{allik} - a vector with all the loglikelihoods, should be monotone increasing. \cr\cr
#' \code{it} - the number of iterations taken \cr\cr
#' @export
#' @examples
#' # Tensor normal distribution on 6x7x8x9 responses 
#' set.seed(1234)
#' RmsT <- 6:9
#' RhsT <- 3:5
#' Rbt <- 100
#' Rnn <- 2
#' Rcorrs = c("AR(1)","EQC","ARMA","N")
#' arma_params <- list(NULL, NULL, c(1, 2), NULL)
#' Rsig2t <- 20
#' dat <- diagdat_sim(msT=RmsT,hsT=RhsT,bt=Rbt,nn=Rnn,
#'                    sig2t=Rsig2t,corrs=Rcorrs, arma_param = arma_params)
#' fit <- MLE_tensnorm(Yall = dat$Yall, centered = F,
#'                    corrs = Rcorrs, arma_param = arma_params)
#' par(mfrow = c(1,2))
#' hist(fit$M,main = "estimates in M (true in blue) ")
#' abline(v= Rbt,col = "blue",lwd=3)
#' plot(ts(fit$allik),main = "convergence",ylab = "loglikelihood")
#' par(mfrow = c(1,1))
#' covars <- rbind("true" = c(Rsig2t,sapply(dat$SST,function(x)x[2,1])),
#'                 "est" = c(fit$sig2,sapply(fit$covs,function(x)x[2,1])))
#' colnames(covars) <- c("sig2",paste0("S-",Rcorrs))
#' covars
#' @author Carlos Llosa-Vite, \email{llosacarlos2@@gmail.com}
#' @author Subrata Pal, \email{SubrataluPal@@gmail.com}
#' @author Ranjan Maitra, \email{maitra@@iastate.edu}
#' @references Llosa-Vite, C., & Maitra, R. (2022). 
#'   \href{https://doi.org/10.1109/TPAMI.2022.3164836}{Reduced-Rank Tensor-on-Tensor Regression and Tensor-variate Analysis of Variance}
#'   \emph{IEEE Transactions on Pattern Analysis and Machine Intelligence}, 45(2), 2282 - 2296.  
#' @references Manceur, A.M., Dutilleul, P. (2013).
#'   \href{https://doi.org/10.1016/j.cam.2012.09.017}{Maximum likelihood estimation for the tensor normal distribution: Algorithm, minimum sample size, and empirical bias and dispersion}
#'   \emph{Journal of Computational and Applied Mathematics}, 239, 37-49.  
MLE_tensnorm <- function(Yall,it = 100,err = 1e-8, corrs = "N",centered = T,initS = NULL, arma_param = NULL){
  #MLE algorithm with covariances set to AR(1) autocorrelations
  p <- length(dim(Yall))-1
  n <- dim(Yall)[p+1]
  ms <- dim(Yall)[-(p+1)]
  mn <- prod(dim(Yall))
  # center + centering
  if(!centered){
    M <- apply(Yall,1:p,mean)
    Yall <- Yall - replicate(n,M)
  }
  #setting up correlations types
  if(length(corrs) == 1) corrs <- rep(corrs,p)
  Sgen <- as.list(1:p)
  if_arma <- rep(FALSE, p)
  if (is.null(arma_param)) arma_param <- as.list(rep(NA,p))
  
  for(k in 1:p){
    if (corrs[k] == "ARMA" || corrs[k] == "ARMA(p,q)" || corrs[k] == "ARMA(p, q)") {
      if_arma[k] <- TRUE
      Sgen[[k]] <- arma
      if(length(arma_param[[k]]) != 2) { 
        arma_param[[k]] <- c(1, 1)
      }
    } else if(corrs[k] == "AR(1)"){
      Sgen[[k]] <- ar1
    } else if(corrs[k] == "EQC"){
      Sgen[[k]] <- eqc
    } else if(corrs[k] == "MA(1)"){
      Sgen[[k]] <- ma1
    } else if(corrs[k] == "independent"){
      Sgen[[k]] <- indep
    } else if(corrs[k] == "N"){
      Sgen[[k]] <- ADJUST
    } else {
      cat(k,"TVN initialization: unknown corr",corrs[k],". Will assume unconstrained Sigma \n")
      Sgen[[k]] <- ADJUST
    }
  }
  if(is.null(initS)){
    Stypa <- lapply(ms,function(i)Styp(diag(i)))
    sig2 <- norm(Yall)^2/mn
  }else{
    Stypa <- lapply(initS$covs,Styp)
    sig2 <- initS$sig2
    Yall <- tprod(Yall,lapply(Stypa,function(x)x$isqr))
  }

  prev <- 1
  allconv <- NULL
  for(i in 1:it){
    #inverse square root of covariance matrices
    for(j in 1:p){
      #find S2 along with sig2
      S <-  Stypa[[j]]$sqr %*% sqmode(Yall,j) %*% Stypa[[j]]$sqr
      if (if_arma[j]) {
        SR  <- Styp(arma(n*prod(ms[-j]), sig2, S, arma_param[[j]][1], arma_param[[j]][2]))
      } else {
        SR  <- Styp(Sgen[[j]](n*prod(ms[-j]),sig2,S))
      } 
      sig2 <- sum(SR$inv*S)/mn
      #update Yall and the list
      Yall <- amprod(Yall,SR$isqr%*%Stypa[[j]]$sqr,j)
      Stypa[[j]] <- SR
    }
    #convergence
    conv <-  log(sig2) + sum(sapply(Stypa,function(x)x$ldet)/ms)
    allconv <- c(allconv,conv)
    if(abs((conv-prev)/prev) < err) break
    else prev <- conv
  }
  covs <- lapply(Stypa,function(x)x$orig)
  allik <- -mn*(1+log(2*pi)+allconv)/2
  toret <- list(covs = covs,sig2 = sig2,allik = allik,it = i)
  if(!centered) toret$M <- M
  return(toret)
}
