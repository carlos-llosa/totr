
#' Balanced designs TANOVA Data generation
#' 
#' Generates a TANOVA balanced design, where each covariate \eqn{X_i} is an 
#' indicator tensor with 1 at one position and zeroes everywhere else. Here, 
#' null hypothesis of equality of groups is true for all factor combinations.
#'
#' @rdname diagdat_sim
#'
#' @param msT The size of each \eqn{Y_i}, i.e., \eqn{(m_1, ..., m_p)}.
#' @param hsT The size of each \eqn{X_i}, i.e., \eqn{(h_1, ..., h_l)}.
#' @param bt The unique value repeated in the entries of  the Beta tensor of size \eqn{(h_1, ..., h_l,m_1, ..., m_p)}
#' @param nn The number of repetitions per factor combination.
#' @param sig2t General variance \eqn{\sigma^2}
#' @param corrs Character vector of size p indicating the
#' types of covariance matrices desired for \eqn{S_1 ,.., S_p}.
#' Options are "AR(1)", "MA(1)", "ARMA"/"ARMA(p,q)"/"ARMA(p, q)", "EQC"  for
#' AR(1), MA(1), ARMA(p, q) and equivariance correlation matrices, and
#' "N" for general covariance with element (1,1) equal to 1.
#' If corrs is of size 1, then \eqn{S_1 ,.., S_p} will all have the
#' same correlation structure.
#' @param arma_param A list of size \code{length(dim(Yall))},
#' each of which contains the ARMA parameter orders (p, q)
#' for that corresponding mode.q
#' p is the AR parameter order and q is the MA parameter order.
#' If some other mode has some other kind of correlation structure
#' and you still want to specify the ARMA orders,
#' you can input a list of size p with other cases as NULL.
#' The default ARMA order is (1, 1).
#' @param covars The covariance matrices \eqn{[S_1 ,.., S_p]} explicitly specified.
#' It would be a list of size p, where each element is a covariance matrix of
#' size \eqn{m_i \times m_i}.
#' @returns A list containing the following elements: \cr\cr
#' \code{Yall} - Array containing the n tensor responses along the last mode, so that it is of size \eqn{m_1 \times .. \times m_p \times n}.
#'   The last dimension must match the last dimension of \code{Xall}. \cr\cr
#' \code{Xall} Array containing the n tensor covariates along the last mode, so that it is of size \eqn{h_1 \times .. \times h_l \times n}.
#'   The last dimension must match the last dimension of \code{Yall}. \cr\cr
#' \code{SST} - The list of matrices \eqn{[S_1 ,.., S_p]}, \eqn{S_i} is of size \eqn{m_i \times m_i}. \cr\cr
#' @export
#' @examples
#' # Tensor-on-Tensor Regression on 6x7x8x9 responses and 3x4x5 covariates
#' set.seed(1234)
#' RmsT <- 6:9
#' RhsT <- 3:5
#' Rbt <- 100
#' Rnn <- 2
#' Rcorrs <- c("AR(1)","EQC","MA(1)","N")
#' Rsig2t <- 20
#' dat <- diagdat_sim(msT=RmsT,hsT=RhsT,bt=Rbt,nn=Rnn,sig2t=Rsig2t,corrs=Rcorrs)
#' print(dim(dat$Xall))
#' print(dim(dat$Yall))
#' print(dim(dat$SST[[1]]))
#' @author Carlos Llosa-Vite, \email{llosacarlos2@@gmail.com}
#' @author Subrata Pal, \email{subrata@@iastate.edu}
#' @references Llosa-Vite, C., & Maitra, R. (2022). 
#'   \href{https://doi.org/10.1109/TPAMI.2022.3164836}{Reduced-Rank Tensor-on-Tensor Regression and Tensor-variate Analysis of Variance}
#'   \emph{IEEE TPAMI}, 45(2), 2282 - 2296.  
#' 
diagdat_sim <- function(msT, hsT, bt, nn, sig2t, corrs, 
                        arma_param = NULL, covars = NULL){

  #size of response
  msp <- prod(msT)
  #size f covariate
  hsp <- prod(hsT)
  #sample size
  n <- prod(hsT)*nn
  #True parameters
  # M and B
  Bt <- array(bt,c(hsT,msT))

  # covariances
  if(is.null(covars)){
    SST <- list()
    for(k in 1:length(msT)){
      if(corrs[k] == "AR(1)") {
        SST[[k]] <- covAR1(runif(1,-.5,.5), msT[k])
      } else if(corrs[k] == "MA(1)"){
        SST[[k]] <- covMA1(runif(1,-.5,.5), msT[k])
      } else if(corrs[k] == "EQC") {
        SST[[k]] <- covEQC(runif(1,-.1,.5), msT[k])
      } else if(corrs[k] == "ARMA"){
        if(is.null(arma_param[[k]]) | (length(arma_param[[k]]) != 2)) arma_param[[k]] <- c(1, 1)
        SST[[k]] <- covARMA(runif(arma_param[[k]][1],-.5,.5),
                            runif(arma_param[[k]][2],-.5,.5),
                            msT[k]) # why -ve half to half?
      } else {
        if(corrs[k] != "N")  cat("unknown corr",corrs[k],". Will assume unconstrained Sigma \n")
        S <- drop(rWishart(1,1+msT[k],diag(msT[k])))
        SST[[k]] <- S/S[1,1]
      }
    }
  } else {
    SST <- covars
  }

  #simulation
  Xuniq <- array(diag(hsp),c(hsT,hsp))
  Xall <- list()
  for(i in 1:nn) Xall[[i]] <- Xuniq
  Xall <- simplify2array(Xall)
  dim(Xall) <- c(dim(Xall)[1:length(hsT)],n)
  Yall <- array(t(array(Bt,c(hsp,msp))) %*% array(Xall,c(hsp,n)),c(msT,n))
  Yall <- Yall + tprod(array(rnorm(msp*n),c(msT,n)),lapply(SST,function(x)Styp(x)$sqr))*sqrt(sig2t)

  return(list(Yall = Yall,Xall = Xall,SST=SST))

}
