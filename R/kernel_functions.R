r_linear_kernel <- function(x1, x2) {
  return(x1 %*% t(x2))
}

r_poly_kernel <- function(x1, x2, gamma, degree = 3, coef0 = 0) {
  K <- (gamma * x1 %*% t(x2) + coef0)^(degree)
  return(K)
}

r_rbf_kernel <- function(x1, x2, gamma= 1/ncol(x2), symmetric = FALSE) {
  x1 <- as.matrix(x1)
  x2 <- as.matrix(x2)
  n1 <- nrow(x1)
  norms1 <- rowSums(x1^2, 1)#as.matrix(apply(x1^2, 1, sum))
  e2 <- matrix(1, 1, n1)
  if (symmetric == FALSE) {
    n2 <- nrow(x2)
    norms2 <- rowSums(x2^2, 1)#as.matrix(apply(x2^2, 1, sum))
    e1 <- matrix(1, 1, n2)
  } else {
    norms2 <- norms1
    e1 <- e2
  }
  K <- exp(gamma*(-norms1 %*% e1 - t(e2) %*% t(norms2) + 2*(x1 %*% t(x2))))
  return(K)
}


#' Kernel Function
#'
#' @author Zhang Jiaqi.
#' @param x1,x2 input matrices.
#' @param kernel.type choose kernel type.
#' @param gamma parameter for \code{'rbf'} and \code{'poly'} kernel. Default \code{gamma = 1/ncol(X)}.
#' @param degree parameter for polynomial kernel, default: \code{degree = 3}.
#' @param coef0 parameter for polynomial kernel,  default: \code{coef0 = 0}.
#' @param symmetric if \code{x1 == x2} you can set \code{symmetric == TRUE}.
#' @param ... redundant parameters.
#' @export
kernel_function <- function(x1, x2,
                            kernel.type = c('linear', 'rbf', 'poly'),
                            gamma = 1/ncol(x1), degree = 3, coef0 = 0,
                            symmetric = FALSE, ...) {
  kernel.type <- match.arg(kernel.type)
  if (kernel.type == 'linear') {
    K <- r_linear_kernel(x1, x2)
  }else if (kernel.type == 'rbf') {
    K <- r_rbf_kernel(x1, x2, gamma, symmetric = symmetric)
  }else if (kernel.type == 'poly') {
    K <- r_poly_kernel(x1, x2, gamma, degree = 3, coef0 = 0)
  }
  return(K)
}


kernel_select_option <- function(X, kernel, solver, randx,
                                 gamma, degree, coef0, ...) {
  n <- nrow(X)
  if (kernel == "linear" & solver == "primal") {
    KernelX <- X
  } else if (kernel != "linear" & solver == "primal") {
    if (randx > 0 && randx < 1) {
      randX <- X[sample(n, floor(randx*n)),]
    } else {
      randX <- X
    }
    KernelX <- kernel_function(X, randX,
                               kernel.type = kernel,
                               gamma = gamma, degree = degree, coef0 = coef0)
    X <- randX
  } else if (solver == "dual") {
    KernelX <- kernel_function(X, X,
                               kernel.type = kernel,
                               gamma = gamma, degree = degree, coef0 = coef0,
                               symmetric = T)
  }
  K <- list("X" = X, "KernelX" = KernelX)
  return(K)
}
