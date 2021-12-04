

#ADJUST procedure as in Glanz and Carvalho 2018
ADJUST <- function(n,sig2,S){
  mm <- ncol(S)
  if(mm == 1) return(matrix(1,1,1))
  nu <- S[1,1]/(n*sig2)
  S <- S/S[1,1]
  S[2:mm,2:mm] <- nu*S[2:mm,2:mm] + (1-nu)*S[1,2:mm] %*% t(S[1,2:mm])
  return(S)
}

covMA1 <- function(rho,m){
  if(m == 1) return(matrix(1,1,1))
  #Generates the m x m MA(1)-correlation matrix with coef rho
  toeplitz(c(1,rho,rep(0,m-2)))
}
ma1 <- function(m,sig2,Sq){
  # For a squre matrix Sq
  # finds the MA(1) matrix Minimizer of  m(log(det(A))) + tr(Sq*A^-1)/sig2
  # where A is the correlation matrix of the process
  # the MA coefficient is constrained to be between -0.5 to 0.5 only
  mm <- ncol(Sq)
  if(mm == 1) return(matrix(1,1,1))
  minlik <- function(ma1){ # negative likelihood of an AR(1) process
    R <- chol(covMA1(ma1,mm))
    ldA <- 2*sum(log(diag(R))) #log(det(A)) from cholesky
    Ai <- chol2inv(R)#inverse from cholesky
    return(m*ldA+sum(Ai* Sq)/sig2)    # tr(AB) has closed form sum(A*B)
  }
  par1 <- optimize(minlik,c(-.5,.5))$minimum
  return(covMA1(par1,mm))
}

covEQC <- function(rho,m){
  if(m == 1) return(matrix(1,1,1))
  #Generates the m x m equi-correlation matrix with coef rho
  toeplitz(c(1,rep(rho,m-1)))
}
eqc <- function(m,sig2,Sq){
  # For a squre matrix Sq
  # finds the equi-correlation matrix (A) that minimizes  m(log(det(A))) + tr(Sq*A^-1)/sig2
  mm <- ncol(Sq)
  if(mm == 1) return(matrix(1,1,1))
  minlik <- function(rho){ # negative likelihood
    ldA <- (mm-1)*log(1-rho) + log(1-rho+rho*mm) #log(det(A)) has closed form
    rd <- rho/(1-rho+mm*rho)
    Ai <- toeplitz(c(1-rd,rep(-rd,mm-1))/(1-rho))
    return(m*ldA+sum(Ai* Sq)/sig2)    # tr(AB) has closed form sum(A*B)
  }
  par1 <- optimize(minlik,c(-1/(m-1)+1e-5,1-1e-5))$minimum
  return(covEQC(par1,mm))
}
covAR1 <- function(ar1,m){
  if(m == 1) return(matrix(1,1,1))
  #Generates the m x m AR(1)-correlation matrix with coef ar1
  toeplitz(ar1^(0:(m-1)))
}
ar1 <- function(m,sig2,Sq){
  # For a squre matrix Sq
  # finds the AR(1) matrix Minimizer of  m(log(det(A))) + tr(Sq*A^-1)/sig2
  # where A is the correlation matrix of the process
  mm <- ncol(Sq)
  if(mm == 1) return(matrix(1,1,1))
  if(mm>3){ #closed forms for inverse and determinant
    minlik <- function(ar1){ # negative likelihood of an AR(1) process
      ldA <- (mm-1)*log(1-ar1^2) #log(det(A)) has closed form
      Ai <- toeplitz(c(1,-ar1,rep(0,mm-2))/(1-ar1^2))#inverse of AR(1)
      Ai[replicate(2,2:(mm-1))] <- Ai[replicate(2,2:(mm-1))] + ar1^2/(1-ar1^2)
      return(m*ldA+sum(Ai* Sq)/sig2)    # tr(AB) has closed form sum(A*B)
    }
  }else{
    minlik <- function(ar1){ # negative likelihood of an AR(1) process
      R <- chol(covAR1(ar1,mm))
      ldA <- 2*sum(log(diag(R))) #log(det(A)) from cholesky
      Ai <- chol2inv(R)#inverse from cholesky
      return(m*ldA+sum(Ai* Sq)/sig2)    # tr(AB) has closed form sum(A*B)
    }
  }
  par1 <- optimize(minlik,c(-1+1e-5,1-1e-5))$minimum
  return(covAR1(par1,mm))
}
