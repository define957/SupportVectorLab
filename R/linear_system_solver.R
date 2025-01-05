#' Solve Linear Systems Use Cholesky Decomposition.
#'
#' @author Zhang Jiaqi.
#' @param a a symmetric square matrix containing the coefficients of the linear system.
#' @param b a vector or matrix giving the right-hand side(s).
#' @export
cholsolve <- function(a, b) {
  Lt <- chol(a, pivot = T)
  pivot <- attr(Lt, "pivot")
  oo <- order(pivot)
  x <- LLtsolve(Lt, b[pivot,])
  dim(x) <- c(ncol(a), ncol(b))
  x <- x[oo, ]
  return(x)
}

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
