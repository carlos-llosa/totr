

MLE_tensnorm <- function(Yall,it = 100,err = 1e-4, corrs = "N",centered = T,initS = NULL){
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
  #setting up types of correlations
  if(length(corrs) == 1) corrs <- rep(corrs,p)
  Sgen <- as.list(1:p)
  for(k in 1:p){
    if(corrs[k] == "AR(1)"){
      Sgen[[k]] <- ar1
    } else if(corrs[k] == "EQC"){
      Sgen[[k]] <- eqc
    } else if(corrs[k] == "MA(1)"){
      Sgen[[k]] <- ma1
    } else Sgen[[k]] <- ADJUST
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
      SR <- Styp(Sgen[[j]](n*prod(ms[-j]),sig2,S))
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
