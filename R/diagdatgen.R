

#' @export
diagdat_sim <- function(msT,hsT,bt,nn,sig2t,corrs){

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
  SST <- list()
  for(k in 1:length(msT)){
    if(corrs[k] == "AR(1)") SST[[k]] <- covAR1(runif(1,-.5,.5),msT[k])
    else if(corrs[k] == "MA(1)") SST[[k]] <- covMA1(runif(1,-.5,.5),msT[k])
    else if(corrs[k] == "EQC") SST[[k]] <- covEQC(runif(1,-.1,.5),msT[k])
    else{
      if(corrs[k] != "N")  cat("unknown corr",corrs[k],". Will assume unconstrained Sigma \n")
      S <- drop(rWishart(1,1+msT[k],diag(msT[k])))
      SST[[k]] <- S/S[1,1]
    }
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
