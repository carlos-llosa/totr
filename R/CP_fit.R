
#' Tensor-on-Tensor Regression with Kronecker Separable Covariance and CP Format Coefficient
#'
#' Tensor-on-tensor regression
#' \eqn{Y_i = < X_i | B > + E_i}\cr
#' with kronecker-separable covariance
#' \eqn{Var(vec(E_i)) = \sigma^2 \otimes_k S_k}
#' and CP-formatted
#' \eqn{B = [[  L_1 ,.., L_l , M_1 ,..., M_p ]]}.\cr
#' The size of the ith tensor covariate \eqn{X_i} is \eqn{h_1 \times .. \times h_l},
#' The size of the ith tensor response \eqn{Y_i} is \eqn{m_1 \times .. \times m_p},
#' the size of the regression coefficient B is of size \eqn{h_1 \times .. \times h_l \times m_1 \times .. \times m_p},
#' and \eqn{i=1,...,n}.\cr
#' The matrix \eqn{L_k} is of size \eqn{h_k \times R}, where R is the CP rank,
#' the matrix \eqn{M_k} is of size \eqn{m_k \times R},
#' and matrix \eqn{S_k} is positive definite of size \eqn{m_k \times m_k} and is restricted to (1,1) element equal to 1.
#'
#' Convergence is achieved when the relative difference in snorm(B) + snorm(S) is less than err.
#' snorm(Y) here means norm(Y)/sqrt(size of Y).
#'
#' @param Yall Array containing the n tensor responses along the last mode, so that it is of size \eqn{m_1 \times .. \times m_p \times n}.
#' The last dimension must match the last dimension of Xall.
#' @param Xall Array containing the n tensor covariates along the last mode, so that it is of size \eqn{h_1 \times .. \times h_l \times n}.
#' The last dimension must match the last dimension of Yall.
#' @param R positive whole number indicating the CP rank.
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
#' @param retB Return the dense tensor B of low rank? 
#' @return A list containing the following elements: \cr\cr
#' \code{B} -  the estimated coefficient B. \cr\cr
#' \code{P} - A list with the estimated \eqn{L_1 ,.., L_l} in the list \code{Ls} and estimated \eqn{M_1 ,.., M_p} in the list \code{Ms}.\cr\cr
#' \code{sig2} - the estimate of \eqn{\sigma^2}. \cr\cr
#' \code{covs} - a list with the estimated matrices \eqn{S_1 ,.., S_p}. \cr\cr
#' \code{allconv} - a vector with all the convergence criteria. \cr\cr
#' \code{allik} - a vector with all the loglikelihoods, should be monotone increasing. \cr\cr
#' \code{it} - the number of iterations taken \cr\cr
#' @export
#' @importFrom stats optimize rWishart rnorm runif toeplitz uniroot
#' @importFrom utils tail
#' @examples
#' # Tensor-on-Tensor Regression on 6x7x8x9 responses and 3x4x5 covariates
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
#' fit <- CP_normal(Yall = dat$Yall,Xall = dat$Xall, R = 2,
#'                    corrs = Rcorrs, arma_param = arma_params)
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
#' @author Subrata Pal, \email{SubrataluPal@@gmail.com}
#' @author Ranjan Maitra, \email{maitra@@iastate.edu}
#' @references Llosa-Vite, C., & Maitra, R. (2022). 
#'   \href{https://doi.org/10.1109/TPAMI.2022.3164836}{Reduced-Rank Tensor-on-Tensor Regression and Tensor-variate Analysis of Variance}
#'   \emph{IEEE Transactions on Pattern Analysis and Machine Intelligence}, 45(2), 2282 - 2296.  
#' 
CP_normal <- function(Yall, Xall, R, it = 100, err = 1e-8, init = NULL,
                      corrs = "N", arma_param = NULL, retB = T) {

  #setting up dimensions
  p <- length(dim(Yall))-1
  l <- length(dim(Xall))-1
  if(dim(Xall)[l+1] != dim(Yall)[p+1]) stop("sample size of X and Y do not match")
  n <- dim(Yall)[p+1]
  hs <- dim(Xall)[-(l+1)]
  ms <- dim(Yall)[-(p+1)]
  mn <- prod(dim(Yall))

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
      cat(k,"unknown corr",corrs[k],". Will assume unconstrained Sigma \n")
      Sgen[[k]] <- ADJUST
    }
  }

  supd1 <- function(Gs,X){
    #this function takes Ls,Xs, and returns the superdiagonals (for each i) of Tucker(Xi,t(Ls))
    #the i=1,2,..,n are in the last mode of the array Xs
    # this is performed efficiently by avoiding matrix products
    res <- NULL
    for(r in 1:R){
      Ns <- Reduce("%o%",lapply(Gs,function(x)x[,r]))                 # problem if hst <- 1:3
      res <- rbind(res,apply(X,length(dim(X)),function(x)sum(x*Ns)))
    }
    return(res)
  }
  supd2 <- function(Gs,X,k){
    #this function is similar to supd1, except it also takes an argument k
    #k is the mode that should be fixed, so not superdiagonal everywhere
    res <- list()
    for(r in 1:R){
      Ns <- Reduce("%o%",lapply(Gs,function(x)x[,r]))
      resH <- apply(X,c(k,length(dim(X))),function(x)sum(x*Ns))
      res[[r]] <- resH
    }
    ph <- length(dim(resH))
    aperm(simplify2array(res),c((ph+1),1:ph))
  }
  #initial values for CP component
  if(is.null(init$CP)){
    Ls <- lapply(hs[-l],function(l){
      LR <- normalize(array(runif(l*R),c(l,R)))
    })
    Ls[[l]] <- array(runif(hs[l]*R),c(hs[l],R)) # the last L is not normalized
    Ms <- lapply(ms,function(l){
      LR <- normalize(array(runif(l*R),c(l,R)))
    })
  }else{
    Ls <- init$CP$Ls
    Ms <- init$CP$Ms
  }
  #initial values for covariance component
  if(is.null(init$covs)){
    res_des <- diag(n) - tcrossprod(svd(array(Xall,c(prod(hs),n)))$v)
    init2 <- MLE_tensnorm(amprod(Yall,res_des,p+1),it = 10,corrs = corrs)
    sig2 <- init2$sig2
    STTs <- lapply(init2$covs,Styp)
  }else{
    sig2 <- init$covs$sig2
    STTs <- lapply(init$covs$covs,Styp)
  }
  #setting lists based on initial values
  MSs <- lapply(1:p,function(i) STTs[[i]]$isqr %*% Ms[[i]])
  MMs <- lapply(MSs,crossprod)
  Yall <- tprod(Yall,lapply(STTs,function(x)x$isqr))
  allconv <- NULL
  allik <- NULL
  prev <- 1


  ## HERE BEGIN ITERATIONS



  for(j in 1:it){
    ##################### STAGE 1 ######################
    # Ms and Ss
    GallR <- supd1(Ls,Xall)
    GallR2 <- tcrossprod(GallR)

    if(p == 1){ # for vector responses
      Sqrprev <- STTs[[1]]$sqr
      #M and S
      Mr <- Sqrprev %*% Yall %*% t(GallR) %*% ginvS(GallR2)
      SR <- tcrossprod(Sqrprev %*% Yall) - Mr %*% GallR2 %*% t(Mr)
      if (if_arma[k]) {
        STTs[[1]] <- Styp(arma(n, sig2, SR, arma_param[[k]][1], arma_param[[k]][2]))
      } else {
        STTs[[1]] <- Styp(Sgen[[1]](n,sig2,SR))
      }
      #sig2 and update list, Yall
      Yall <- STTs[[1]]$isqr %*% Sqrprev %*% Yall
      sig2 <- sum(STTs[[1]]$inv*SR)/mn
      Ms[[1]] <- normalize(Mr)
      MSs[[1]]  <- STTs[[1]]$isqr %*% Ms[[1]]
      MMs[[1]]  <- crossprod(MSs[[1]])
      HR <- GallR2
    } else{
      for(k in 1:p){  #otherwise
        Sqrprev <- STTs[[k]]$sqr #previous matrix square root, since standardized Yall is stored
        #find M and S (Mr and S1r)
        HR <- GallR2 * Reduce("*",MMs[-k])
        Yz <- supd2(MSs[-k],Yall,k)
        Mr <- apply(apply(t(1:n),2,function(i)GallR[,i]*Yz[,,i]),1,sum)
        Mr <- Sqrprev %*% t(array(Mr,c(R,ms[k]))) %*% ginvS(HR)
        SR <- Sqrprev %*% sqmode(Yall,k) %*% Sqrprev - Mr %*% HR %*% t(Mr)
        if (if_arma[k]) {
          STTs[[k]]  <- Styp(arma(n*prod(ms[-k]), sig2, SR, arma_param[[k]][1], arma_param[[k]][2]))
        } else {
          STTs[[k]]  <- Styp(Sgen[[k]](n*prod(ms[-k]),sig2,SR))
        }
        #update lists, Yall and s2
        Yall <- amprod(Yall,STTs[[k]]$isqr %*% Sqrprev,k)
        sig2 <- sum(STTs[[k]]$inv*SR)/mn
        Ms[[k]] <- normalize(Mr)
        MSs[[k]]  <- STTs[[k]]$isqr %*% Ms[[k]]
        MMs[[k]]  <- crossprod(MSs[[k]])
      }
    }

    ##################### STAGE 2 ######################
    MAM <- Reduce("*",MMs)
    YnMs <- supd1(MSs,Yall)
    # Ls
    if(l == 1){ #for vector covariates
      #find L1
      gA <- ginvS(tcrossprod(Xall))
      gB <- ginvS(MAM)
      gC <- tcrossprod(Xall,YnMs)
      LR <-  gA %*% gC %*% gB
      #update list and sig2
      Ls[[1]] <- LR
      sig2 <- sum(STTs[[p]]$inv*(SR + Mr %*% HR %*% t(Mr))) - sum(diag(gA %*% gC %*% gB %*% t(gC)))
      sig2 <- sig2/mn
    } else { #for all higher orders
      for(k in 1:l){
        #find the superdiagonals
        Xwn <- supd2(Ls[-k],Xall,k)
        #find \sum H_\Sigma^{-1}H_i' in Xw2
        Xw2 <- sqmode(Xwn,3,t=F)
        dim(Xw2) <-  c(R,hs[k],R,hs[k])
        Xw2 <- apply(Xw2,c(2,4),function(x)x*MAM)
        dim(Xw2) <- c(R,R,hs[k],hs[k])
        Xw2 <- aperm(Xw2,c(1,3,2,4))
        dim(Xw2) <- c(hs[k]*R,hs[k]*R)
        #now perform the multiplication with the inverse
        LRv <- apply(apply(t(1:n),2,function(i)YnMs[,i]*Xwn[,,i]),1,sum) %*% ginvS(Xw2)
        LR <- matrix(LRv,ncol = R,byrow = T)
        if(k != l) LR <- normalize(LR) #don't normalize the last L
        #update list
        Ls[[k]] <- LR
      }
      #update sig2 after all Ls
      sig2 <- sum(STTs[[p]]$inv*(SR + Mr %*% HR %*% t(Mr))) - sum((t(LRv) %*% LRv)*Xw2)
      sig2 <- sig2/mn
    }
    ##################### CONVERGENCE ######################
    #the following criteria is Rnorm(B) + norm(S), where Rnorm is norm divided by dimension, and S is precision
    normB <-  sum( Reduce("*",lapply(Ls,crossprod))*Reduce("*",lapply(Ms,crossprod)) )/sqrt(prod(hs)*prod(ms))
    normS <-  prod(sapply(STTs,function(x)x$norm))*sqrt(sig2)/prod(ms)
    conv <-  normB + normS
    allconv <- c(allconv,conv)
    allik <- c(allik, log(sig2) +sum(sapply(STTs,function(x)x$ldet)/ms)       )
    if(abs((conv-prev)/prev)<err){
      break
    } else prev <- conv
  }
  #create tensor B, with all Ms and Ls unit column norms, and norms in lam
  if(retB) {
    B <- array(0,c(hs,ms))
    for(r in 1:R) B <- B + Reduce("%o%",lapply(Ls,function(x)x[,r])) %o% Reduce("%o%",lapply(Ms,function(x)x[,r]))
  }
  lam <- apply(Ls[[l]],2,norm)
  Ls[[l]] <- normalize(Ls[[l]])
  #set to descending order in norms (in lam)
  ord <- order(lam,decreasing = T)
  lam <- lam[ord]
  Ls <- lapply(Ls,function(x)x[,ord])
  Ms <- lapply(Ms,function(x)x[,ord])
  #save
  CP = list(lam=lam,Ls=Ls,Ms=Ms)
  covs = lapply(STTs,function(x)x$orig)
  allik <- -mn*(1+log(2*pi)+allik)/2
  toret <- list(CP = CP, sig2=sig2 ,covs = covs, allconv = allconv,it = j,allik=allik)
  if(retB) toret$B <- B
  return(toret)
}

