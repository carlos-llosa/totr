
#' Tensor-on-Tensor Regression with Kronecker Separable Covariance and Tensor Ring Format Coefficient
#'
#' Tensor-on-tensor regression
#' \eqn{Y_i = < X_i | B > + E_i}\cr
#' with kronecker-separable covariance
#' \eqn{Var(vec(E_i)) = \sigma^2 \otimes_k S_k}
#' and Tensor-Ring formatted
#' \eqn{B = tr( L_1 \times^1 .. \times^1 L_l \times^1 M_1 \times^1 .. \times^1 M_p )}.\cr
#' The size of the ith tensor covariate \eqn{X_i} is \eqn{h_1 \times .. \times h_l},
#' The size of the ith tensor response \eqn{Y_i} is \eqn{m_1 \times .. \times m_p},
#' the size of the regression coefficient B is of size \eqn{h_1 \times .. \times h_l \times m_1 \times .. \times m_p},
#' and \eqn{i=1,...,n} .\cr
#' The TR rank (g) is \eqn{(s_1 ,.., s_l,g_1 ,..,g_p)}, where \eqn{s_0 = g_p} and \eqn{g_0 = s_l},\cr
#' the tensor \eqn{L_k} is of size \eqn{s_k-1 \times h_k \times s_k}  ,
#' the tensor \eqn{M_k} is of size \eqn{g_k-1 \times m_k \times g_k},
#' and matrix \eqn{S_k} is positive definite of size \eqn{m_k \times m_k} and is restricted to (1,1) element equal to 1.
#'
#' Convergence is achieved when the relative difference in snorm(B) + snorm(S) is less than err.
#' snorm(Y) here means norm(Y)/sqrt(size of Y)
#'
#' @param Yall Array containing the n tensor responses along the last mode, so that it is of size \eqn{m_1 \times .. \times m_p \times n}.
#' The last dimension must match the last dimension of Xall.
#' @param Xall Array containing the n tensor covariates along the last mode, so that it is of size \eqn{h_1 \times .. \times h_l \times n}.
#' The last dimension must match the last dimension of Yall.
#' @param g vector of size p+l containing the Tensor Ring rank
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
#' \code{TR} - A list containing the estimated \eqn{L_1 ,.., L_l} in the list \code{Ls} and
#' the estimated \eqn{M_1 ,.., M_p} in the list \code{Ms}.\cr\cr
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
#' Rcorrs = c("AR(1)","EQC","ARMA","N")
#' Rsig2t <- 20
#' dat <- diagdat_sim(msT=RmsT,hsT=RhsT,bt=Rbt,nn=Rnn,sig2t=Rsig2t,corrs=Rcorrs)
#' # No arma_param supplied means ARMA(1, 1) by default
#' 
#' g <- c(2,2,2,2,2,2,2)
#' fit <- TR_normal(Yall = dat$Yall,Xall = dat$Xall, g = g,corrs = Rcorrs,it=30)
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
#' @references \url{https://arxiv.org/abs/2012.10249}
#' @import tensorA
TR_normal <- function(Yall,Xall,g,it = 100,init = NULL, err = 1e-7,corrs = "N", arma_param = NULL){
  #setting up dimensions
  p <- length(dim(Yall))-1
  l <- length(dim(Xall))-1
  if(dim(Xall)[l+1] != dim(Yall)[p+1]) stop("sample size of X and Y do not match")
  n <- dim(Yall)[p+1]

  #setting up correlations types
  if(length(corrs) == 1) corrs <- rep(corrs,p)
  Sgen <- as.list(1:p)
  if_arma <- rep(FALSE, p)
  if (is.null(arma_param)) arma_param <- as.list(1:p)

  for(k in 1:p){
    if (corrs[k] == "ARMA" || corrs[k] == "ARMA(p,q)" || corrs[k] == "ARMA(p, q)") {
      if_arma[k] <- TRUE
    } else {
      arma_param[[k]] <- NA
    }

    if(corrs[k] == "AR(1)"){
      Sgen[[k]] <- ar1
    } else if(corrs[k] == "EQC"){
      Sgen[[k]] <- eqc
    } else if(corrs[k] == "MA(1)"){
      Sgen[[k]] <- ma1
    } else if (if_arma[k]) {
      Sgen[[k]] <- arma
      if(length(arma_param[[k]]) != 2) { 
        arma_param[[k]] <- c(1, 1)
      }
    } else Sgen[[k]] <- ADJUST
  }


  #setting up dimensions (very important for tensor algebra, see package tensorA)
  mdims <- dim(Yall)
  names(mdims) <- c(paste0("m",1:p),"n")
  ms <- mdims[1:p]
  ldims <- dim(Xall)
  names(ldims) <- c(paste0("l",1:l),"n")
  hs <- ldims[1:l]
  names(g) <- c(paste0("s",1:l),paste0("g",1:p))
  mn <- prod(dim(Yall))
  Ldims <- lapply(1:l,function(i){
    ii <- ifelse(i==1,p+l,i-1)
    c(hs[i],g[ii],g[i])
  })
  Mdims <- lapply(1:p,function(i){
    c(ms[i],g[i+l-1],g[i+l])
  })
  Xall <- as.tensor(Xall,dims = ldims)
  Yall <- as.tensor(Yall,dims = mdims)
  #initial values for TR component
  if(is.null(init$TR)){
    Ms <- lapply(1:p,function(i){
      M <- array(runif(prod(Mdims[[i]])),Mdims[[i]])
      M <- M/norm(M)
      as.tensor(M,dims = Mdims[[i]])
    })
    Ls <- lapply(1:l,function(i){
      L <- array(runif(prod(Ldims[[i]])),Ldims[[i]])
      if(i != l) L <- L/norm(L)
      L <- as.tensor(L,dims = Ldims[[i]])
    })
  }else{
    Ls <- init$TR$Ls
    Ms <- init$TR$Ms
  }
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
  MZs <- lapply(1:p,function(i){
    as.tensor(amprod(Ms[[i]],Stypa[[i]]$isqr,1),dims = Mdims[[i]])
  })

  #convergence monitoring
  allconv <- NULL
  allik <- NULL
  prev <- 1
  normg <- function(A){
    # This function is useful for finding the norm of B in TC format
    tor <- Reduce("+",lapply(1:dim(A)[1],function(i)A[i,,] %x% A[i,,]))
    dim(tor) <- dim(A)[2:3]^2
    tor
  }

  Yall <- tprod(Yall,lapply(Stypa,function(x)x$isqr))
  XXL <- Xall
  for(kk in 1:l) XXL <- XXL %e% Ls[[kk]]

  #iterate till convergence
  for(j in 1:it){

    ######################### STAGE 1


    # needed for M1and M2
    W <- XXL
    #now loop for each (M_k,S_k)
    for(k in 1:p){
      #create M and S
      prevSqr <- Stypa[[k]]$sqr
      Yallz <- prevSqr %*% mat(Yall,k) #remove the dependence on the kt mode
      if(k != p){
        DD <- W %e% Reduce("%e%",MZs[(k+1):p]) #only do forward ones (back ones are done)
      } else {
        DD <- W #if last one do nothing
      }
      H <- untensor(DD,list(names(mdims[-k]),names(Mdims[[k]][-1]))) #this is like X
      HH <- tcrossprod(H) #this is like XX'
      MR <- Yallz %*% t(H) %*% ginvS(HH) #this is like YX'(XX')^{-1}
      S <- tcrossprod(Yallz) - MR %*% HH %*% t(MR)
      MR <- MR/norm(MR)
      #update lists,Yall and sig2
      if (if_arma[k]) {
        Stypa[[k]]  <- Styp(arma(prod(mdims[-k]), sig2, S, arma_param[[k]][1], arma_param[[k]][2]))
      } else {
      Stypa[[k]] <- Styp(Sgen[[k]](prod(mdims[-k]),sig2,S))
      }
      MZs[[k]] <- as.tensor(amprod(MR,Stypa[[k]]$isqr,1),dims = Mdims[[k]])
      Ms[[k]] <- as.tensor( MR , dims= Mdims[[k]])
      sig2 <- sum(Stypa[[k]]$inv*S)/mn
      Yall <- amprod(Yall,Stypa[[k]]$isqr%*%prevSqr,k) #update Yall such that it remains standardized (not sig2)
      W <- W %e% MZs[[k]]
    }

    ######################### STAGE 2

    V <- Reduce("%e%",MZs)

    YMs <- V %e% as.tensor(Yall,dims=mdims)
    V2 <- mul.tensor(V,names(ms),mark(V))
    XXL <- Xall
    for(k in 1:l){
      # create L
      XL <- XXL
      if(k != l ) for(kk in (k+1):l) XL <- XL %e% Ls[[kk]] #only do forward ones (back ones are done)
      XL2 <- mul.tensor(XL,names(mdims)[p+1],mark(XL))
      YG <- untensor(YMs %e% XL,list(names(Ldims[[k]]))) #this is like YX'
      GG <- untensor(XL2 %e% V2,list(names(Ldims[[k]]),paste0(names(Ldims[[k]]),"'"))) #this is like XX'
      LR <- to.tensor(c(YG %*% ginvS(GG)),Ldims[[k]]) #this is like YX'(XX')^{-1}
      #save the result and update lists
      if(k != l) LR <- LR/norm(LR)
      Ls[[k]] <- LR
      XXL <- XXL %e% Ls[[k]]
    }
    #########################  convergence
    # sum of relative norms of B and S

    normB <- sqrt(sum(Reduce("%*%",lapply(Ls,normg)) * t(Reduce("%*%",lapply(Ms,normg)))))/sqrt(prod(c(hs,ms)))
    normS <- prod(sapply(Stypa,function(x)x$norm))*sqrt(sig2)/prod(ms)
    conv <- normB + normS
    allconv <- c(allconv,conv)
    allik <- c(allik, log(sig2) +sum(sapply(Stypa,function(x)x$ldet)/ms)       )
    if(abs((conv-prev)/prev)<err){
      #  print(paste("converged at",i,"iterations with a relative tolerance of ", err))
      break
    } else prev <- conv
  }
  B <- Reduce("%e%",Ls) %e% Reduce("%e%",Ms)

  #turn into arrays and aperm
  B <- tarray(B)
  Ls <- lapply(Ls,function(x) aperm(tarray(x),c(2,1,3)))
  Ms <- lapply(Ms,function(x) aperm(tarray(x),c(2,1,3)))

  #return
  TR = list(Ls=Ls,Ms=Ms)
  covs = lapply(Stypa,function(x)x$orig)
  allik <- -mn*(1+log(2*pi)+allik)/2
  toret <- list( B = B, TR = TR, allik=allik,sig2=sig2 ,covs = covs, allconv = allconv,it = j)
  return(toret)
}



