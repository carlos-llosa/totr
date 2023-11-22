
#' Tensor-on-Tensor Regression with Kronecker Separable Covariance and OP Format (outer product) Coefficient
#'
#' Tensor-on-tensor regression
#' \eqn{Y_i = < X_i | B > + E_i}\cr
#' with kronecker-separable covariance
#' \eqn{Var(vec(E_i)) = \sigma^2 \otimes_k S_k}
#' and OP-formatted
#' \eqn{B = o[[  M_1 ,..., M_p ]]}.\cr
#' The size of the ith tensor covariate \eqn{X_i} is \eqn{h_1 \times .. \times h_p},
#' The size of the ith tensor response eqn{Y_i} is \eqn{m_1 \times .. \times m_p},
#' the size of the regression coefficient B is of size \eqn{h_1 \times .. \times h_l \times m_1 \times .. \times m_p},
#' and \eqn{i=1,...,n} .\cr
#' The matrix \eqn{M_k} is of size \eqn{h_k \times m_k},
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
#' @param init list containing initial values. init$covssig2 contains the initial value of \eqn{\sigma^2}, init$covs$covs is a list of p matrices
#' containing initial values for \eqn{S_1 ,.., S_p}, init$TK$Ls is a list of l matrices containing the initial values for \eqn{L_1 ,..,L_l}, and
#' init$TK$Ms is a list of p matrices containing the initial values for \eqn{M_1 ,..,M_p}.
#' If init = NULL then the elements in init$covs will be initiated from the TVN model fitted on the unconstrained B residuals
#' and init$Ls, init$Ms will contain elements generated randomly from the uniform(0,1) distribution.
#' @param corrs Character vector of size p indicating the types of covariance matrices desired for \eqn{S_1 ,.., S_p}.
#' Options are "AR(1)", "MA(1)", "ARMA"/"ARMA(p,q)"/"ARMA(p, q)", "EQC"  for
#' AR(1), MA(1), ARMA(p, q) and equivariance correlation matrices, and
#' "N" for general covariance with element (1,1) equal to 1.
#' If corrs is of size 1, then \eqn{S_1 ,.., S_p} will all have the same correlation structure.
#' @param arma_param A list of size \code{length(dim(Yall))}, each of which contains the
#' ARMA parameter orders (p, q) for that corresponding mode.
#' p is the AR parameter order and q is the MA parameter order
#' If some other mode has some other kind of correlation structure
#' and you still want to specify the ARMA orders,
#' you can input a list of size p with other cases as NULL.
#' The default ARMA order is (1, 1).
#' @return A list containing the following elements: \cr\cr
#' \code{B} -  the estimated coefficient B. \cr\cr
#' \code{OP} - A list with the estimated \eqn{M_1 ,.., M_p}.\cr\cr
#' \code{sig2} - the estimate of \eqn{\sigma^2}. \cr\cr
#' \code{covs} - a list with the estimated matrices \eqn{S_1 ,.., S_p}. \cr\cr
#' \code{allconv} - a vector with all the convergence criteria. \cr\cr
#' \code{allik} - a vector with all the loglikelihoods, should be monotone increasing. \cr\cr
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
#' fit <- OP_normal(Yall = dat$Yall,Xall = dat$Xall,corrs = Rcorrs)
#' par(mfrow = c(1,2))
#' hist(fit$B,main = "estimates in B (true in blue)")
#' abline(v= Rbt,col = "blue",lwd=3)
#' plot(ts(fit$allik),main = "loglikelihood")
#' par(mfrow = c(1,1))
#' covars <- rbind("true" = c(Rsig2t,sapply(dat$SST,function(x)x[2,1])),
#'                 "est" = c(fit$sig2,sapply(fit$covs,function(x)x[2,1])))
#' colnames(covars) <- c("sig2",paste0("S-",Rcorrs))
#' covars
#' @author Carlos Llosa-Vite, \email{llosacarlos2@@gmail.com}
#' @references Llosa-Vite, C., & Maitra, R. (2022). 
#'   \href{https://doi.org/10.1109/TPAMI.2022.3164836}{Reduced-Rank Tensor-on-Tensor Regression and Tensor-variate Analysis of Variance}
#'   \emph{IEEE TPAMI}, 45(2), 2282 - 2296.  
#' 
OP_normal <- function(Yall,Xall,it = 100, err = 1e-7,init = NULL,corrs = "N", arma_param = NULL){
  #setting up dimensions
  p <- length(dim(Yall))-1
  l <- length(dim(Xall))-1
  if(dim(Xall)[l+1] != dim(Yall)[p+1]) stop("sample size of X and Y do not match")
  n <- dim(Yall)[p+1]

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
      cat("unknown corr",corrs[k],". Will assume unconstrained Sigma \n")
      Sgen[[k]] <- ADJUST
    }
  }

  #setting up dimensions (very important for tensor algebra, see package tensorA)
  mdims <- dim(Yall)
  names(mdims) <- c(paste0("m",1:p),"n")
  ms <- mdims[1:p]
  ldims <- dim(Xall)
  names(ldims) <- c(paste0("l",1:l),"n")
  hs <- ldims[1:l]
  mn <- prod(dim(Yall))

  #initial values for covariance component
  if(is.null(init$covs)){
    res_des <- diag(n) - tcrossprod(svd(array(Xall,c(prod(hs),n)))$v)
    init2 <- MLE_tensnorm(amprod(Yall,res_des,p+1),it = 4,corrs = corrs)
    sig2 <- init2$sig2
    Stypa <- lapply(init2$covs,Styp)
  }else{
    sig2 <- init$covs$sig2
    Stypa <- lapply(init$covs$covs,Styp)
  }

  #initial values for OP component
  if(is.null(init$OP)){
    Ms <- lapply(1:p,function(i){
      M <- matrix(runif(ms[i]*hs[i]),ms[i])
      M/norm(M)
    })
  }else{
    Ms <- init$OP
  }
  MSi <- lapply(1:p,function(i)Stypa[[i]]$isqr %*% Ms[[i]])
  Yall <- tprod(Yall,lapply(Stypa,function(x)x$isqr))

  #convergence
  prev <- 1
  allconv <- NULL
  allik <- NULL
  #iterate
  for(j in 1:it){

    for(k in 1:p){
      #create M and S
      Sqrprev <- Stypa[[k]]$sqr
      Xd <- mat(tprod(Xall,MSi,c(1:p)[-k]),k)
      Xd2 <- tcrossprod(Xd)
      Yd <- Sqrprev %*%mat(Yall,k)
      M <- Yd %*% t(Xd) %*% ginvS(Xd2)
      S <- tcrossprod(Yd) - M %*% Xd2 %*% t(M)
      #save them and sig2
      if(k != p) M <- M/norm(M)
      if (if_arma[k]) {
        Stypa[[k]]  <- Styp(arma(prod(mdims[-k]), sig2, S, arma_param[[k]][1], arma_param[[k]][2]))
      } else {
      Stypa[[k]] <- Styp(Sgen[[k]](prod(mdims[-k]), sig2, S))
      }
      MSi[[k]]   <- Stypa[[k]]$isqr %*% M
      Ms[[k]]    <- M
      Yall       <- amprod(Yall,Stypa[[k]]$isqr%*%Sqrprev,k)
      sig2 <- sum(S * Stypa[[k]]$inv)/mn
    }


    #convergence
    normB <- prod(sapply(Ms,norm))/sqrt(prod(c(hs,ms)))
    normS <- sqrt(sig2)*prod(sapply(Stypa,function(x)x$norm))/prod(ms)
    conv <- normB + normS
    allik <- c(allik,   log(sig2) + sum(sapply(Stypa,function(x)x$ldet)/ms)  )
    allconv <- c(allconv,conv)
    if(abs((conv-prev)/prev) < err){
      #    print(paste("converged at",i,"iterations with a relative tolerance of ", err));
      break
    } else prev <- conv
  }
  allik <- -mn*(1+log(2*pi)+allik)/2
  B <- aperm(Reduce("%o%",Ms),c(2*(1:p),2*(1:p)-1))
  covs <- lapply(Stypa,function(x)x$orig)
  toret <- list( B = B, OP = Ms, sig2=sig2 ,covs = covs, allconv = allconv,allik=allik,it = j)
  return(toret)
}
