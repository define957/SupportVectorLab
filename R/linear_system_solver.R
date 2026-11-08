#' Solve Linear Systems Use Cholesky Decomposition.
#'
#' @author Zhang Jiaqi.
#' @param a a symmetric square matrix containing the coefficients of the linear system.
#' @param b a vector or matrix giving the right-hand side(s).
#' @export
cholsolve <- function(a, b) {
  Lt <- chol(a)
  x <- LLtsolve(Lt, b)
  return(x)
}

# cholsolve <- function(
#     a,
#     b,
#     check_symmetric = FALSE,
#     check_solution = FALSE,
#     tol_solution = 1e-7,
#     singular = c("zero", "NA")
# ) {
#   singular <- match.arg(singular)
#   if (check_symmetric && !isSymmetric(a)) {
#     stop("cholsolve: 'a' must be a symmetric matrix")
#   }
#   b      <- as.matrix(b)
#   n      <- ncol(a)
#   Lt     <- chol(a, pivot = TRUE)
#   pivot  <- attr(Lt, "pivot")
#   rank   <- attr(Lt, "rank")
#   x      <- matrix(0, nrow = n, ncol = ncol(b))
#   if (rank > 0) {
#     bp <- b[pivot, , drop = FALSE]
#     y  <- matrix(0, n, ncol(b))
#     y[seq_len(rank), ] <- LLtsolve(Lt[seq_len(rank), seq_len(rank), drop = FALSE],
#                                    bp[seq_len(rank), , drop = FALSE])
#     x[pivot, ] <- y
#   }
#   if (check_solution) {
#     resid <- max(abs(a %*% x - b))
#     scale <- tol_solution * max(1, max(abs(a)) * max(abs(x)))
#     if (resid > scale) {
#       warning("cholsolve: residual ", format(resid, digits = 3),
#               " exceeds scaled tolerance ", format(scale, digits = 3),
#               " (rank ", rank, "/", n, "); system may be inconsistent")
#     }
#   }
#   if (singular == "NA" && rank < n) {
#     x[pivot[(rank + 1):n], ] <- NA
#   }
#   return(x)
# }

#' Solve Linear Systems.
#'
#' \code{LLtsolve} solves the equation {\eqn{L L' x = b}}, where \eqn{L'} is a upper
#' triangula square matrix
#'
#' @author Zhang Jiaqi.
#' @param Lt a upper triangula square matrix containing the coefficients of the linear system.
#' @param b a vector or matrix giving the right-hand side(s).
#' @export
LLtsolve <- function(Lt, b) {
  x <- backsolve(Lt, forwardsolve(Lt, b, upper.tri = T, transpose = T))
  return(x)
}

stack_qrsolve <- function(blocks, rhs, ridge_stack = 0) {
  p     <- ncol(blocks)
  S     <- blocks
  if (ridge_stack > 0) {
    S <- rbind(S, diag(sqrt(ridge_stack), p))
  }
  QRS   <- qr(S, LAPACK = TRUE)
  if (QRS$rank < p) {
    stop("stack_qrsolve: blocks are rank-deficient!")
  }
  R          <- qr.R(QRS)
  pivot      <- QRS$pivot
  x          <- matrix(0, p, ncol(rhs))
  x[pivot, ] <- LLtsolve(R, rhs[pivot, , drop = FALSE])
  x
}

## Factor the constant block [x; sqrt(ridge_stack) * I_p] once, so the IRLS
## loops can stack it as rbind(changing_rows, R) instead of
## rbind(changing_rows, x, sqrt(ridge_stack) * I_p).
##
## The QR is column-pivoted: R'R = P'(x'x + ridge I)P, so callers un-pivot it
## with  R_u <- R[, order(pivot)]  (R_u'R_u = x'x + ridge I, dense p x p) and
## stack R_u with the un-permuted data rows (no rhs / solution permutation).
## R_u is dense: only use it to build the stacked matrix, never in LLtsolve().
qr_aug <- function(x, ridge_stack) {
  p  <- ncol(x)
  qS <- qr(rbind(x, diag(sqrt(ridge_stack), p)), LAPACK = TRUE)
  list(R = qr.R(qS), pivot = qS$pivot)
}

svd_solve <- function(lhs, rhs, lambda) {
  lhs_svd <- svd(lhs)
  V       <- lhs_svd$v
  d2      <- lhs_svd$d^2
  V %*% ((1 / (d2 + lambda)) * (t(V) %*% rhs))
}
