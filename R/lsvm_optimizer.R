lsvm_optimizer <- function(H_left, H_right, a, A, B, q, eps, max.steps, u, beta, ...) {
  return(arma_lsvm_optimizer(H_left, H_right, a, A, B, q, eps, max.steps, u, beta))
}
