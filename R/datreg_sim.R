
#' Data generation when X and Beta are given with types of covariance matrices
#'
#' @rdname datreg_sim
#'
#' @param msT The size of each \eqn{Y_i}, i.e., \eqn{(m_1, ..., m_p)}.
#' @param Xall Array containing the n tensor covariates along the last mode,
#' so that it is of size \eqn{h_1 \times .. \times h_l \times n}.
#' @param Bt The coefficient array of size \eqn{h_1 \times .. \times h_l \times m_1 \times .. \times m_p}
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
#' for that corresponding mode.
#' p is the AR parameter order and q is the MA parameter order.
#' If some other mode has some other kind of correlation structure
#' and you still want to specify the ARMA orders,
#' you can input a list of size p with other cases as NULL.
#' The default ARMA order is (1, 1).
#' @param covars The covariance matrices \eqn{[S_1 ,.., S_p]} explicitly specified.
#' It would be a list of size p, where each element is a covariance matrix of
#' size \eqn{m_i \times m_i}.
#' @returns
#' A list containing the following elements: \cr\cr
#' \code{Yall}: Array containing the n tensor responses along the last mode,
#'    so that it is of size \eqn{m_1 \times .. \times m_p \times n}.
#'    The last dimension must match the last dimension of \code{Xall}. \cr\cr
#' \code{Xall}: Array containing the n tensor covariates along the last mode,
#'    so that it is of size \eqn{h_1 \times .. \times h_l \times n}. \cr\cr
#' \code{SST}: The list of matrices \eqn{[S_1 ,.., S_p]}, \eqn{S_i} is of size \eqn{m_i \times m_i}.
#'
#' @export
#' @examples
#' # Tensor-on-Tensor Regression on 6x7x8 responses and 2x3x4 covariates
#' set.seed(1234)
#' RmsT <- 6:8
#' Rnn <- 10
#' X <- array(1:prod(2:4, Rnn), dim = c(2:4, Rnn))
#' Bt <- array(1:prod(dim(X), RmsT), c(dim(X), RmsT))/1e5
#' Rcorrs <- c("ARMA", "ARMA", "N")
#' arma_params <- list(c(1,1), c(1, 2), NULL)
#' Rsig2t <- 10
#' dat <- datreg_sim(msT=RmsT, Xall=X, Bt=Bt, sig2t=Rsig2t,
#'                    corrs=Rcorrs, arma_param = arma_params)
#' dim(dat$Yall)
#' dim(dat$SST[[1]])
#' @author Carlos Llosa-Vite, \email{llosacarlos2@@gmail.com}
#' @author Subrata Pal, \email{SubrataluPal@@gmail.com}
#' @author Ranjan Maitra, \email{maitra@@iastate.edu}
#' @seealso \code{\link{diagdat_sim}}
#' @references Llosa-Vite, C., & Maitra, R. (2022). 
#'   \href{https://doi.org/10.1109/TPAMI.2022.3164836}{Reduced-Rank Tensor-on-Tensor Regression and Tensor-variate Analysis of Variance}
#'   \emph{IEEE Transactions on Pattern Analysis and Machine Intelligence}, 45(2), 2282 - 2296.  
#' 
datreg_sim <- function(msT, Xall, Bt, sig2t, corrs, arma_param = NULL, covars = NULL){

  #size of response
  msp <- prod(msT)
  #size of covariate
  hsT <- dim(Xall)
  #sample size
  n <- hsT[length(hsT)]
  hsT <- hsT[-length(hsT)]
  hsp <- prod(hsT)

  # covariances
  if(is.null(covars)){
    SST <- list()
    if (is.null(arma_param)) arma_param <- as.list(1:length(msT))

    for(k in 1:length(msT)){
      if(length(arma_param[[k]]) != 2){
        arma_param[[k]] <- NA
      }
      if(corrs[k] == "AR(1)") {
        SST[[k]] <- covAR1(runif(1,-.5,.5), msT[k])
      } else if(corrs[k] == "MA(1)"){
        SST[[k]] <- covMA1(runif(1,-.5,.5), msT[k])
      } else if(corrs[k] == "EQC") {
        SST[[k]] <- covEQC(runif(1,-.1,.5), msT[k])
      } else if(corrs[k] == "ARMA"){
        if(length(arma_param[[k]]) != 2) {
          arma_param[[k]] <- c(1, 1)
        }
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
  Yall <- array(t(array(Bt, c(hsp, msp))) %*% array(Xall, c(hsp, n)), c(msT, n))
  eps <- array(rnorm(msp * n), c(msT, n))
  sigmas <- lapply(SST, function(x){Styp(x)$sqr})
  Yall <- Yall + tprod(eps, sigmas) * sqrt(sig2t)

  return(list(Yall = Yall, Xall = Xall, SST = SST))
}
