
#' Tensor-on-Tensor Regression with Kronecker Separable Covariance and Tucker Format Coefficient
#'
#' Tensor-on-tensor regression
#' \eqn{Y_i = < X_i | B > + E_i}\cr
#' with kronecker-separable covariance
#' \eqn{Var(vec(E_i)) = \sigma^2 \otimes_k S_k}
#' and Tucker-formatted
#' \eqn{B = [[ V; L_1 ,.., L_l , M_1 ,..., M_p ]]}.\cr
#' The size of the ith tensor covariate \eqn{X_i} is \eqn{h_1 \times .. \times h_l},
#' The size of the ith tensor response \eqn{Y_i} is \eqn{m_1 \times .. \times m_p},
#' the size of the regression coefficient B is of size \eqn{h_1 \times .. \times h_l \times m_1 \times .. \times m_p},
#' and \eqn{i=1,...,n} .\cr
#' The size of V is the Tucker rank \eqn{(c_1 ,.., c_l , d_1 ,.., d_p)}, while
#' the matrix \eqn{L_k} is of size \eqn{h_k \times c_k},
#' the matrix \eqn{M_k} is of size \eqn{m_k \times d_k}, \eqn{M_k' S_k^-1 M_k} is an identity matrix,
#' and matrix \eqn{S_k} is positive definite of size \eqn{m_k \times m_k} and is restricted to (1,1) element equal to 1.\cr
#' If the Tucker rank is chosen such that \eqn{c_k = h_k}, then \eqn{L_k} will be an identity matrix,
#' and if \eqn{d_k = m_k}, then \eqn{M_k} will bean identity matrix.
#'
#' Convergence is achieved when the relative difference in snorm(B) + snorm(S) is less than err.
#' snorm(Y) here means norm(Y)/sqrt(size of Y)
#'
#' @param Yall Array containing the n tensor responses along the last mode, so that it is of size \eqn{m_1 \times .. \times m_p \times n}.
#' The last dimension must match the last dimension of Xall.
#' @param Xall Array containing the n tensor covariates along the last mode, so that it is of size \eqn{h_1 \times .. \times h_l \times n}.
#' The last dimension must match the last dimension of Yall.
#' @param pdims vector of size \eqn{p+l} containing the Tucker rank, or size of V.
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
#' \code{Tucker} - A list containing the estimated tensor V, along with the estimated \eqn{L_1 ,.., L_l}
#' in the list Ls and estimated \eqn{M_1 ,.., M_p} in the list Ms.\cr\cr
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
#' pdims <- c(2,2,2,2,2,2,2)
#' fit <- TK_normal(Yall = dat$Yall,Xall = dat$Xall, pdims = pdims,corrs = Rcorrs,it=30)
#' par(mfrow = c(1,2))
#' hist(fit$B,main = "estimates in B (true in blue) ")
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
#' @import tensorA
TK_normal <- function( Yall, Xall, pdims, it = 100, err = 1e-7,init = NULL,corrs = "N", arma_param = NULL, retB = T){
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
      cat(k,"unknown corr",corrs[k],". Will assume unconstrained Sigma \n")
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
  names(pdims) <- c(paste0("gl",1:l),paste0("gm",1:p))
  pdimsM <- pdims[(l+1):(l+p)]
  pdimsL <- pdims[1:l]
  mn <- prod(dim(Yall))

  Ldims <- lapply(1:l,function(i)c(hs[i],pdims[i]))
  Mdims <- lapply(1:p,function(i)c(ms[i],pdims[i+l]))
  #set up the cases where L is identity (when the rank matches the X dimension)
  idenL <- which(pdimsL != hs)
  names(ldims)[which(pdimsL == hs)] <- names(pdimsL)[which(pdimsL == hs)]
  Xall <- as.tensor(Xall,dims = ldims)
  #set up the cases where M is identity (when the rank matches the Y dimension)
  idenM <- which(pdimsM != ms)
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

  #initial values for Tucker component
  if(is.null(init$TK)){
    Ls <- lapply(1:l,function(i){
      if(i %in% idenL){
        return(as.tensor(normalize(matrix(runif(prod(Ldims[[i]])),hs[i])),dims = Ldims[[i]]))
      }else{ #set up Li as identity if rank matches dimension of X
        return(diag(hs[i]))
      }
    })
    Ms <- lapply(1:p,function(i){
      if(i %in% idenM){
        M <- matrix(runif(ms[i]^2),ms[i])
        M <- Stypa[[i]]$sqr %*% svd(Stypa[[i]]$isqr %*% M)$u[,1:pdimsM[i]]
        return(as.tensor(M,dims = Mdims[[i]]))
      }else{
        return(diag(ms[i]))
      }
    })
  }else{
    Ls <- init$TK$Ls
    Ms <- init$TK$Ms
  }

  MSIs <- lapply(1:p,function(i)  as.tensor(Stypa[[i]]$isqr %*% Ms[[i]],dims = Mdims[[i]]))
  XL <- Xall
  for(i in idenL) XL <- XL %e% Ls[[i]]
  Yall <- as.tensor(tprod(Yall,lapply(Stypa,function(x)x$isqr)),dims=mdims)
  allconv <- NULL
  allik <- NULL

  prev <- 1
  for(j1 in 1:it){
    ##################### Block 1 ######################
    W <- svd(untensor(XL,list(names(pdimsL),names(ldims)[l+1])))
    WW <- tcrossprod(W$u) # W' (WW')^- W
    rr <- sum(W$d>1e-10)
    Wi <- as.tensor(W$v[,1:rr] %*% (1/W$d[1:rr] * t(W$u[,1:rr])),dims = c(pdimsL,mdims[1+p])) # W' (WW')^-
    YMSI <- Yall
    ##################### M1 ... Mp ######################

    cout <- 1
    if(length(idenM) != 0 ){ # if at least one idenL needs to be estimated
      for(i in idenM){
        # S and M
        YMR <-  matmult_kdim(YMSI,WW,"n")
        if(i != tail(idenM,1)) for(ii in tail(idenM,length(idenM)-cout)) YMR <- YMR %e% MSIs[[ii]]
        Q  <- Stypa[[i]]$sqr %*% to.matrix.tensor(YMR,names(ms)[i])
        MR <- (Stypa[[i]]$sqr %*% svd(Stypa[[i]]$isqr%*%Q)$u)[,1:pdimsM[i]]
        Ms[[i]]   <- as.tensor(MR,dims = Mdims[[i]])
        MSIs[[i]] <- as.tensor(Stypa[[i]]$isqr %*% Ms[[i]],dims = Mdims[[i]])
        YMSI <- MSIs[[i]] %e% YMSI
        cout <- cout + 1
      }
    }

    ##################### V ######################
    V <- YMSI %e% Wi
    # in the above line V is not complete, because it has \Sigma_k^{-1/2} whenever M_k=I (it shouldnt)
    # however, this is the V that I want when estimating L1..Ll as well as S1 ... Sp
    # so I will fix V at the end, AFTER estimating finishing all the iterations



    ##################### L1 ... Ll ######################


    if(length(idenL) != 0 ){ # if at least one idenL needs to be estimated
      VYm <- V %e% YMSI
      XL <- Xall
      cout <- 1
      for(i in idenL){
        G <- XL
        if(i != tail(idenL,1)) for(ii in tail(idenL,length(idenL)-cout)) G <- G %e% Ls[[ii]]
        hdum <- untensor(V %e% G ,list(names(Ldims[[i]])))
        Lr <- c(G %e% VYm) %*% ginvS(sqmode(hdum,1))
        Ls[[i]] <- as.tensor(Lr,dims = Ldims[[i]])
        XL <- XL %e% Ls[[i]]
        cout <- cout +1
      }
    }

    ##################### S1 ... Sp ######################
    Err <- V %e% XL
    for(i in idenM)  Err <- MSIs[[i]] %e% Err
    Err <- Yall - Err
    #saveRDS(list(Yall = Yall,V=V,XL=XL,MSIs=MSIs,Err=Err),"/home/carlos/Desktop/writes/tot_norm/code/estTK2/dnwork.rds")
    for(j in 1:p){
      #find S2 along with sig2
      S <-  Stypa[[j]]$sqr %*% sqmode(Err,j) %*% Stypa[[j]]$sqr
      if (if_arma[j]) {
        SR  <- Styp(arma(n*prod(ms[-j]), sig2, S, arma_param[[j]][1], arma_param[[j]][2]))
      } else {
        SR <- Styp(Sgen[[j]](n*prod(ms[-j]),sig2,S))
      }
      sig2 <- sum(SR$inv*S)/mn
      #update Yall and the list
      Scor <- SR$isqr%*%Stypa[[j]]$sqr
      Yall <- amprod(Yall,Scor,j)
      Err <- amprod(Err,Scor,j)
      Stypa[[j]] <- SR
    }
    Yall <- as.tensor(Yall,dims = mdims)

    ################### CONVERGENCE #####################
    #the norm of M1 and M2 are not important, L1 and L2 only though the R in QR decomp.
    Vdd <- V
    for(i in Ls[idenL]) Vdd <- Vdd %e% qr.R(qr(i))
    normB <- norm( Vdd )/sqrt(prod(c(ms,hs)))
    normS <- sqrt(sig2)*prod(sapply(Stypa,function(x)x$norm))/prod(ms)
    conv <- normB + normS
    allik <- c(allik,   log(sig2) + sum(sapply(Stypa,function(x)x$ldet)/ms)  )

    allconv <- c(allconv,conv)
    if(abs((conv-prev)/prev)<err){
      #  print(paste("converged at",i,"iterations with a relative tolerance of ", err))
      break
    } else prev <- conv
  }
  #This loop fixes V, which before has \Sigma_k^{-1/2} multiplied whenever M_k is identity
  for( k in (1:p)[!(1:p %in% idenM)] ) # for k not in idenM
    V <- matmult_kdim(V,Stypa[[k]]$sqr,names(ms)[k])

  
  if(retB) {
    B <- V
    for(i in idenL) B <- B %e% Ls[[i]]
    for(i in idenM) B <- B %e% Ms[[i]]
    B <- reorder.tensor(B,order(names(B),names(c(hs,ms))))
  }

  Ls[idenL] <- lapply(Ls[idenL],qr) #do the QR trick to all except the ones that are identity
  #V <- aperm(V,c((p+1):(l+p),1:p))
  for(i in idenL) V <- matmult_kdim(V,qr.Q(Ls[[i]]),names(pdimsL)[i])
  Ls[idenL] <- lapply(Ls[idenL],qr.Q)

  #reorder the dimensions of V
  names(pdimsM)[(1:p)[!(1:p %in% idenM)]] <- names(ms)[(1:p)[!(1:p %in% idenM)]]
  V <- reorder.tensor(V,order(names(V),names(c(pdimsL,pdimsM))))

  #turn into arrays
  V <- tarray(V)
  Ls <- lapply(Ls,tarray)
  Ms <- lapply(Ms,tarray)

  #return
  Tucker = list(V=V,Ls=Ls,Ms=Ms)
  covs <- lapply(Stypa,function(x)x$orig)
  allik <- -mn*(1+log(2*pi)+allik)/2
  toret <- list(Tucker = Tucker, sig2=sig2 ,covs = covs, allconv = allconv,allik=allik,it = j1)
  if(retB) toret$B <- B
  return(toret)
}

