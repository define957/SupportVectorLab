#' Stopping Criterion Based on Function Value Tolerance.
#'
#' \code{function_tolerance} tests whether the change of an objective (or loss)
#' function value between two successive iterations is small enough to stop.
#'
#' The criterion is a mixed absolute-relative tolerance:
#' \deqn{|f_{new} - f_{old}| \le \epsilon (1 + |f_{old}|)}
#' When \eqn{|f_{old}|} is large it behaves as a relative tolerance of
#' \eqn{\epsilon}; when \eqn{|f_{old}|} is close to zero it behaves as an
#' absolute tolerance of \eqn{\epsilon}.
#'
#' @author Zhang Jiaqi.
#' @param f_new a numeric value, the function value at the current iteration.
#' @param f_old a numeric value, the function value at the previous iteration.
#' @param eps a non-negative numeric value, the tolerance \eqn{\epsilon}.
#' @return a logical value, \code{TRUE} if the function value change satisfies
#'   the tolerance (i.e. convergence is declared).
#' @export
#' @examples
#' function_tolerance(1.0001, 1.0, 1e-3)  # TRUE, converged
#' function_tolerance(1.01, 1.0, 1e-3)    # FALSE, not converged
function_tolerance <- function(f_new, f_old, eps) {
  return(abs(f_new - f_old) <= eps*(1 + abs(f_old)))
}
