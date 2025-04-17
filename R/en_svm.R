en_svm_dual_solver <- function(KernelX, y, C1, C2, eps, max.steps) {

  n <- nrow(KernelX)
  e <- matrix(1, n, 1)

  H0 <- calculate_svm_H(KernelX, y)
  diag(H0) <- diag(H0) + 1/C1
  IC1 <- diag(1/C1, n)
  dualH <- rbind(cbind(H0, IC1), cbind(IC1, IC1))

  q <- matrix(0, 2*n)
  q[1:n] <- 1
  lb <- matrix(-C2, 2*n)
  lb[1:n] <- 0
  ub <- matrix(Inf, 2*n)
  u0 <- lb
  u <- clip_dcd_optimizer(dualH, q, lb, ub, eps, max.steps, u0)$x

  coef <- y * u[1:n]
  BaseDualENSVMClassifier <- list(coef = as.matrix(coef))
  class(BaseDualENSVMClassifier) <- "BaseDualENSVMClassifier"
  return(BaseDualENSVMClassifier)
}

#' Hinge Support Vector Machine
#'
#' \code{en_svm} is an R implementation of EN-SVM
#'
#' @author Zhang Jiaqi.
#' @param X,y dataset and label.
#' @param C1,C2 plenty term.
#' @param kernel kernel function. The definitions of various kernel functions are as follows:
#' \describe{
#'     \item{linear:}{\eqn{u'v}{u'*v}}
#'     \item{poly:}{\eqn{(\gamma u'v + coef0)^{degree}}{(gamma*u'*v + coef0)^degree}}
#'     \item{rbf:}{\eqn{e^{(-\gamma |u-v|^2)}}{exp(-gamma*|u-v|^2)}}
#' }
#' @param gamma parameter for \code{'rbf'} and \code{'poly'} kernel. Default \code{gamma = 1/ncol(X)}.
#' @param degree parameter for polynomial kernel, default: \code{degree = 3}.
#' @param coef0 parameter for polynomial kernel,  default: \code{coef0 = 0}.
#' @param eps the precision of the optimization algorithm.
#' @param max.steps the number of iterations to solve the optimization problem.
#' @param fit_intercept if set \code{fit_intercept = TRUE},
#'                      the function will evaluates intercept.
#' @param ... unused parameters.
#' @return return \code{SVMClassifier} object.
#' @export
en_svm <- function(X, y, C1 = 1, C2 = 1,
                   kernel = c("linear", "rbf", "poly"),
                   gamma = 1 / ncol(X), degree = 3, coef0 = 0,
                   eps = 1e-5, max.steps = 4000,
                   fit_intercept = TRUE, ...) {
  X <- as.matrix(X)
  y <- as.matrix(y)
  class_set <- sort(unique(y))
  idx <- which(y == class_set[1])
  y[idx] <- 1
  y[-idx] <- -1
  y <- as.matrix(as.numeric(y))
  if (length(class_set) > 2) {
    stop("The number of class should less 2!")
  }
  kernel <- match.arg(kernel)
  solver <- "dual"
  reduce_set <- NULL
  if (fit_intercept == TRUE) {
    X <- cbind(X, 1)
  }
  reduce_flag <- is.null(reduce_set) == FALSE
  if (solver == "dual" && reduce_flag == TRUE) {
    reduce_flag <- FALSE
    reduce_set <- NULL
    cat("The dual solver does not support the reduce set; it has been set to NULL.\n")
  }
  kso <- kernel_select_option_(X, kernel, reduce_set, gamma, degree, coef0)
  KernelX <- kso$KernelX
  if (solver == "dual") {
    solver.res <- en_svm_dual_solver(KernelX, y, C1, C2,eps, max.steps)
  }
  SVMClassifier <- list("X" = X, "y" = y, "class_set" = class_set,
                        "C1" = C1, "C2" = C2, "kernel" = kernel,
                        "gamma" = gamma, "degree" = degree, "coef0" = coef0,
                        "solver" = solver, "coef" = solver.res$coef,
                        "fit_intercept" = fit_intercept)
  class(SVMClassifier) <- "SVMClassifier"
  return(SVMClassifier)
}

