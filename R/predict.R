#' Predict Method for Tensor-on-Tensor Regression with Tucker Decomposition
#'
#' Computes the tensor product \eqn{< X_i | B >} without explicitly constructing the full B tensor.
#' Uses the Tucker decomposition \eqn{B = [[V; L_1,...,L_p, M_1,...,M_q]]} for efficient computation.
#'
#' @param X Tensor covariate of size \eqn{X_i} is \eqn{h_1 \times .. \times h_l}
#' @param Tucker List containing Tucker decomposition components:
#'   - V: Core tensor of size \eqn{(c_1 ,.., c_l , d_1 ,.., d_p)}
#'   - Ls: List of l matrices \eqn{L_k} of size \eqn{h_k \times c_k} 
#'   - Ms: List of p matrices \eqn{M_k} of size \eqn{m_k \times d_k}
#' @param retB Logical indicating whether to return full tensor B (for small problems)
#'
#' @return Predicted tensor \eqn{Y_i} of size \eqn{m_1 \times .. \times m_p}
#' @export
#' @examples
#' # After fitting TK_normal model:
#' # fit <- TK_normal(Yall = Yall, Xall = Xall, pdims = pdims, corrs = corrs, it=2)
#' # Prediction for new observation new_X
#' prediction_fast <- TK_predict(new_X, fit$Tucker)
#' prediction_full <- TK_predict_full(new_X, fit$B)
#' norm(prediction - prediction2)
#' @author Carlos Llosa-Vite, \email{llosacarlos2@@gmail.com}
#' @author Subrata Pal, \email{SubrataluPal@@gmail.com}
#' @author Ranjan Maitra, \email{maitra@@iastate.edu}
#' @references Llosa-Vite, C., & Maitra, R. (2022). 
#'   \href{https://doi.org/10.1109/TPAMI.2022.3164836}{Reduced-Rank Tensor-on-Tensor Regression and Tensor-variate Analysis of Variance}
#'   \emph{IEEE Transactions on Pattern Analysis and Machine Intelligence}, 45(2), 2282 - 2296.  
#' 
#' @import tensorA
TK_predict <- function(X, Tucker, retB = FALSE) {
  # Input validation
  if (!all(c("V", "Ls", "Ms") %in% names(Tucker))) {
    stop("Tucker must contain V, Ls, Ms components")
  }
  
  V <- Tucker$V
  Ls <- Tucker$Ls
  Ms <- Tucker$Ms
  
  l <- length(Ls)
  p <- length(Ms)
  
  # Dimension compatibility check
  if (length(dim(X)) != l) {
    stop("X must have", l, "modes matching Ls length")
  }
  if (!all(sapply(1:l, function(i) dim(X)[i] == nrow(Ls[[i]])))) {
    stop("X dimensions must match Ls matrix rows")
  }
  
    # Efficient Tucker contraction
    h_dims <- sapply(1:l, function(i) dim(X)[i])
    names(h_dims) <- paste0("l", 1:l)
    X_tensor <- tensorA::as.tensor(X, dims = h_dims)
    
    # Left contraction: X ×_1 L_1^T ×_2 ... ×_l L_l^T
    Z <- X_tensor
    for (i in 1:l) {
      Z <- Z %e% Ls[[i]]  # tensorA mode product
    }
    
    # Core contraction
    W <- Z %e% V
    
    Y_pred <- W
    for (j in 1:p) {
      Y_pred <- Y_pred %e% Ms[[j]]
    }
    
  return(as.array(Y_pred))
}

#' Dense Prediction (Direct Tensor Contraction)
TK_predict_full <- function(X, B) {
  if (!all(dim(X) == dim(B)[1:length(dim(X))])) {
    stop("X and B dimension mismatch")
  }
  X_tensor <- tensorA::as.tensor(X)
  B_tensor <- tensorA::as.tensor(B)
  return(as.array(X_tensor %e% B_tensor))  # Full tensor contraction
}
