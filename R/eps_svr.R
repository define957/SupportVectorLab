eps_svr_dual_solver <- function(KernelX, y, C = 1, epsilon = 0.2,
                                eps = 1e-5, max.steps = 80) {

  Q1 <- cbind(KernelX, -KernelX)
  Q2 <- cbind(-KernelX, KernelX)
  H <- rbind(Q1, Q2)
  e <- rep(1, nrow(y))
  q1 <- epsilon * e - y
  q2 <- epsilon * e + y

  q <- t(rbind(c(q1, q2)))
  nH <- nrow(H)
  lb <- matrix(0, nrow = nH)
  ub <- matrix(C, nrow = nH)

  alphas <- clip_dcd_optimizer(H, -q, lb, ub, eps, max.steps)$x
  coef <- c(alphas[1:(nH/2)] - alphas[-c(1:(nH/2))])
  BaseDualEpsSVMRegressor <- list(coef = as.matrix(coef))
  class(BaseDualEpsSVMRegressor) <- "BaseDualEpsSVMRegressor"
  return(BaseDualEpsSVMRegressor)
}

#' Epsilon Insensitive Loss Support Vector Regression
#'
#' \code{eps_svr} is an R implementation of Epsilon Insensitive Loss SVR.
#'
#' @author Zhang Jiaqi.
#' @param X,y dataset and label.
#' @param C plenty term.
#' @param epsilon paramter of epsilon insensitive zone.
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
#' @return return \code{SVMRegressor} object.
#' @import Rcpp
#' @useDynLib SupportVectorLab, .registration=TRUE
#' @export
eps_svr <- function(X, y, C = 1, epsilon = 0.2, kernel = c("linear", "rbf", "poly"),
                    gamma = 1 / ncol(X), degree = 3, coef0 = 0,
                    eps = 1e-5, max.steps = 4000,
                    fit_intercept = TRUE, ...) {
  randx <- 1
  X <- as.matrix(X)
  y <- as.matrix(y)
  kernel <- match.arg(kernel)
  solver <- "dual"
  if (fit_intercept == TRUE) {
    X <- cbind(X, 1)
  }
  kso <- kernel_select_option(X, kernel, solver, randx,
                              gamma, degree, coef0)
  KernelX <- kso$KernelX
  X <- kso$X
  if (solver == "dual") {
    solver.res <- eps_svr_dual_solver(KernelX, y, C, epsilon, eps, max.steps)
  }
  SVMRegressor <- list("X" = X, "y" = y,
                       "C" = C, "kernel" = kernel,
                       "gamma" = gamma, "degree" = degree, "coef0" = coef0,
                       "solver" = solver, "coef" = solver.res$coef,
                       "fit_intercept" = fit_intercept)
  class(SVMRegressor) <- "SVMRegressor"
  return(SVMRegressor)
}

#' Coef Method for Support Vector Regression
#'
#' @author Zhang Jiaqi
#' @param object a fitted object of class inheriting from \code{SVMRegressor}.
#' @param ... unused parameter.
#' @importFrom stats coef
#' @export
coef.SVMRegressor <- function(object, ...) {
  if (object$solver == "dual") {
    return(t(object$X) %*% object$coef)
  } else if (object$solver == "primal") {
    return(object$coef)
  }
}

#' Predict Method for Support Vector Regression
#'
#' @author Zhang Jiaqi
#' @param object a fitted object of class inheriting from \code{SVMRegressor}.
#' @param X new data for predicting.
#' @param ... unused parameter.
#' @importFrom stats predict
#' @export
predict.SVMRegressor <- function(object, X, ...) {
  X <- as.matrix(X)
  if (object$fit_intercept == TRUE) {
    X <- cbind(X, 1)
  }
  if (object$kernel == "linear" & object$solver == "primal") {
    KernelX <- X
  } else {
    KernelX <- kernel_function(X, object$X,
                               kernel.type = object$kernel,
                               gamma = object$gamma,
                               degree = object$degree,
                               coef0 = object$coef0)
  }
  pred <- KernelX %*% object$coef
  return(pred)
}

