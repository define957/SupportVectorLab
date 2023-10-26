#' Clipping Dual Coordinate Decent Optimizer
#'
#' @author Zhang Jiaqi.
#' @param H,q input matrices.
#' @param lb,ub lower bound and upper bound.
#' @param eps,max.steps error and maximum iterations.
#' @param u init solution.
#' @param ... redundant parameters.
#' @return return results list
#' @export
clip_dcd_optimizer <- function(H, q, lb, ub,
                               eps = 1e-5, max.steps = 200,
                               u = (lb + ub) / 2, ...) {

  # Clipping dual coordinate decent optimizer
  # solve quadratic programming like:
  #       min t(x) %*% H %*% x - t(q) %*% x
  #       s.t lb < x < ub
  H <- as.matrix(H)
  q <- as.matrix(q)
  lb <- as.matrix(lb)
  ub <- as.matrix(ub)
  clip_dcd_res <- cpp_clip_dcd_optimizer(H, q, lb, ub, eps, max.steps, u)
  return(clip_dcd_res)
}
