

MLE_tensnorm <- function(Yall,it = 100,err = 1e-4, corrs = "N",centered = T,initS = NULL, arma_param = NULL){
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
  toret <- list(covs = covs,sig2 = sig2,allconv = allconv)
  if(!centered) toret$M <- M
  return(toret)
}
