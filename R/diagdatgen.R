

#' Diagonal Data generation
#' 
#' @rdname diagdat_sim
#'
#' @param msT The size of each Y_i, i.e., (m_1, ..., m_p).
#' @param hsT The size of each X_i, i.e., (h_1, ..., h_l).
#' @param bt The value generating the Beta tensor of size hsT x msT
#' @param nn sample size/prod(hsT) (??)
#' @param sig2t General variance
#' @param corrs Character vector of size p indicating the
#' types of covariance matrices desired for S_1 ,.., S_p.
#' Options are "AR(1)", "MA(1)", "ARMA"/"ARMA(p,q)"/"ARMA(p, q)", "EQC"  for
#' AR(1), MA(1), ARMA(p, q) and equivariance correlation matrices, and
#' "N" for general covariance with element (1,1) equal to 1.
#' If corrs is of size 1, then S_1 ,.., S_p will all have the
#' same correlation structure.
#' @param arma_param A list of size length(dim(Yall)),
#' each of which contains the ARMA parameter orders (p, q)
#' for that corresponding mode.
#' p is the AR parameter order and q is the MA parameter order.
#' If some other mode has some other kind of correlation structure
#' and you still want to specify the ARMA orders,
#' you can input a list of size p with other cases as NULL.
#' The default ARMA order is (1, 1).
#' @param covars The covariance matrices [S_1 ,.., S_p] explicitly specified.
#' It would be a list of size p, where each element is a covariance matrix of
#' size $m_i x m_i$.
#' @returns A list containing the following elements: \cr\cr
#' * \code{Yall} - Array containing the n tensor responses along the last mode, so that it is of size m_1 x .. x m_p x n.
#'   The last dimension must match the last dimension of Xall. \cr\cr
#' * \code{Xall} Array containing the n tensor covariates along the last mode, so that it is of size h_1 x .. x h_l x n.
#'   The last dimension must match the last dimension of Yall. \cr\cr
#' * \code{SST} - The list of matrices [S_1 ,.., S_p], S_i is of size m_i x m_i. \cr\cr
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
#' dim(dat$Xall)
#' dim(dat$Yall)
#' dim(dat$SST[[1]])
#' @author Carlos Llosa-Vite, \email{llosacarlos2@@gmail.com}
#' @references \url{https://arxiv.org/abs/2012.10249}
diagdat_sim <- function(msT, hsT, bt, nn, sig2t, corrs, arma_param = NULL, covars = NULL){

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
