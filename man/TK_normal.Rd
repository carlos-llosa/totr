% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TK_fit.R
\name{TK_normal}
\alias{TK_normal}
\title{Tensor-on-Tensor Regression with Kronecker Separable Covariance and Tucker Format Coefficient}
\usage{
TK_normal(
  Yall,
  Xall,
  pdims,
  it = 100,
  err = 1e-07,
  init = NULL,
  corrs = "N",
  arma_param = NULL,
  retB = T
)
}
\arguments{
\item{Yall}{Array containing the n tensor responses along the last mode, so that it is of size \eqn{m_1 \times .. \times m_p \times n}.
The last dimension must match the last dimension of Xall.}

\item{Xall}{Array containing the n tensor covariates along the last mode, so that it is of size \eqn{h_1 \times .. \times h_l \times n}.
The last dimension must match the last dimension of Yall.}

\item{pdims}{vector of size \eqn{p+l} containing the Tucker rank, or size of V.}

\item{it}{maximum number of iterations.}

\item{err}{relative error used to assess convergence.}

\item{init}{list containing initial values. init$covssig2 contains the initial value of \eqn{\sigma^2}, init$covs$covs is a list of p matrices
containing initial values for \eqn{S_1 ,.., S_p}, init$TK$Ls is a list of l matrices containing the initial values for \eqn{L_1 ,..,L_l}, and
init$TK$Ms is a list of p matrices containing the initial values for \eqn{M_1 ,..,M_p}.
If init = NULL then the elements in init$covs will be initiated from the TVN model fitted on the unconstrained B residuals
and init$Ls, init$Ms will contain elements generated randomly from the uniform(0,1) distribution.}

\item{corrs}{Character vector of size p indicating the types of covariance matrices desired for \eqn{S_1 ,.., S_p}.
Options are "AR(1)", "MA(1)", "ARMA"/"ARMA(p,q)"/"ARMA(p, q)", "EQC"  for
AR(1), MA(1), ARMA(p, q) and equivariance correlation matrices, and
"N" for general covariance with element (1,1) equal to 1.
If corrs is of size 1, then \eqn{S_1 ,.., S_p} will all have the same correlation structure.}

\item{arma_param}{A list of size \code{length(dim(Yall))}, each of which contains the
ARMA parameter orders (p, q) for that corresponding mode.
p is the AR parameter order and q is the MA parameter order
If some other mode has some other kind of correlation structure
and you still want to specify the ARMA orders,
you can input a list of size p with other cases as NULL.
The default ARMA order is (1, 1).}

\item{retB}{Return the dense tensor B of low rank?}
}
\value{
A list containing the following elements: \cr\cr
\code{B} -  the estimated coefficient B. \cr\cr
\code{Tucker} - A list containing the estimated tensor V, along with the estimated \eqn{L_1 ,.., L_l}
in the list Ls and estimated \eqn{M_1 ,.., M_p} in the list Ms.\cr\cr
\code{sig2} - the estimate of \eqn{\sigma^2}. \cr\cr
\code{covs} - a list with the estimated matrices \eqn{S_1 ,.., S_p}. \cr\cr
\code{allconv} - a vector with all the convergence criteria. \cr\cr
\code{allik} - a vector with all the loglikelihoods, should be monotone increasing. \cr\cr
\code{it} - the number of iterations taken \cr\cr
}
\description{
Tensor-on-tensor regression
\eqn{Y_i = < X_i | B > + E_i}\cr
with kronecker-separable covariance
\eqn{Var(vec(E_i)) = \sigma^2 \otimes_k S_k}
and Tucker-formatted
\eqn{B = [[ V; L_1 ,.., L_l , M_1 ,..., M_p ]]}.\cr
The size of the ith tensor covariate \eqn{X_i} is \eqn{h_1 \times .. \times h_l},
The size of the ith tensor response \eqn{Y_i} is \eqn{m_1 \times .. \times m_p},
the size of the regression coefficient B is of size \eqn{h_1 \times .. \times h_l \times m_1 \times .. \times m_p},
and \eqn{i=1,...,n} .\cr
The size of V is the Tucker rank \eqn{(c_1 ,.., c_l , d_1 ,.., d_p)}, while
the matrix \eqn{L_k} is of size \eqn{h_k \times c_k},
the matrix \eqn{M_k} is of size \eqn{m_k \times d_k}, \eqn{M_k' S_k^-1 M_k} is an identity matrix,
and matrix \eqn{S_k} is positive definite of size \eqn{m_k \times m_k} and is restricted to (1,1) element equal to 1.\cr
If the Tucker rank is chosen such that \eqn{c_k = h_k}, then \eqn{L_k} will be an identity matrix,
and if \eqn{d_k = m_k}, then \eqn{M_k} will bean identity matrix.
}
\details{
Convergence is achieved when the relative difference in snorm(B) + snorm(S) is less than err.
snorm(Y) here means norm(Y)/sqrt(size of Y)
}
\examples{
# Tensor-on-Tensor Regression on 6x7x8x9 responses and 3x4x5 covariates
set.seed(1234)
RmsT <- 6:9
RhsT <- 3:5
Rbt <- 100
Rnn <- 2
Rcorrs = c("AR(1)","EQC","ARMA","N")
Rsig2t <- 20
dat <- diagdat_sim(msT=RmsT,hsT=RhsT,bt=Rbt,nn=Rnn,sig2t=Rsig2t,corrs=Rcorrs)
# No arma_param supplied means ARMA(1, 1) by default

pdims <- c(2,2,2,2,2,2,2)
fit <- TK_normal(Yall = dat$Yall,Xall = dat$Xall, pdims = pdims,corrs = Rcorrs,it=30)
par(mfrow = c(1,2))
hist(fit$B,main = "estimates in B (true in blue) ")
abline(v= Rbt,col = "blue",lwd=3)
plot(ts(fit$allik),main = "loglikelihood")
par(mfrow = c(1,1))
covars <- rbind("true" = c(Rsig2t,sapply(dat$SST,function(x)x[2,1])),
                "est" = c(fit$sig2,sapply(fit$covs,function(x)x[2,1])))
colnames(covars) <- c("sig2",paste0("S-",Rcorrs))
covars
}
\references{
Llosa-Vite, C., & Maitra, R. (2022). 
  \href{https://doi.org/10.1109/TPAMI.2022.3164836}{Reduced-Rank Tensor-on-Tensor Regression and Tensor-variate Analysis of Variance}
  \emph{IEEE Transactions on Pattern Analysis and Machine Intelligence}, 45(2), 2282 - 2296.
}
\author{
Carlos Llosa-Vite, \email{llosacarlos2@gmail.com}

Subrata Pal, \email{SubrataluPal@gmail.com}

Ranjan Maitra, \email{maitra@iastate.edu}
}
