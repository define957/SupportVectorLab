#' Lagrangian Support Vector Machine Optimizer
#'
#' @author Zhang Jiaqi.
#' @param H,Hinv,q input matrices.
#' @param eps,max.steps error and maximum iterations.
#' @param u init solution.
#' @param beta parameter for speeding up convergence.
#' @param ... redundant parameters.
#' @return return results list
#' @export
lsvm_optimizer <- function(H, Hinv, q, eps, max.steps, u, beta, ...) {
  return(arma_lsvm_optimizer(H, Hinv, q, eps, max.steps, u, beta))
}
