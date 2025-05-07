hinge_eps_tsvr_dual_solver <- function(KernelX, y, C1, C2, C3, C4, epsilon1, epsilon2,
                                       eps, max.steps) {
  xn <- nrow(KernelX)
  xp <- ncol(KernelX)
  G <- KernelX
  GramG <- t(G) %*% G
  GTG_C3_inv_GT <- cholsolve(GramG + diag(C3, xp), t(G))
  dualH1 <- G %*% GTG_C3_inv_GT

  if (C3 != C4) {
    GTG_C4_inv_GT <- cholsolve(GramG + diag(C4, xp), t(G))
    dualH2 <- G %*% GTG_C4_inv_GT
  } else {
    GTG_C4_inv_GT <- GTG_C3_inv_GT
    dualH2 <- dualH1
  }
  q1 <- dualH1 %*% y - y - epsilon1
  q2 <- y - epsilon2 - dualH2 %*% y

  lb <- matrix(0, xn, 1)
  ub1 <- matrix(C1, xn, 1)
  ub2 <- matrix(C2, xn, 1)

  x0 <- lb

  alphas <- clip_dcd_optimizer(dualH1, q1, lb, ub1, eps, max.steps, x0)$x
  gammas <- clip_dcd_optimizer(dualH2, q2, lb, ub2, eps, max.steps, x0)$x

  u1 <- GTG_C3_inv_GT %*% (y - alphas)
  u2 <- GTG_C4_inv_GT %*% (y + gammas)

  BaseDualHingeEPSTWSVRRegressor <- list("coef1" = as.matrix(u1),
                                         "coef2" = as.matrix(u2))
}

#' Hinge Epsilon Twin Support Vector Regression
#'
#' \code{hinge_eps_tsvr} is an R implementation of Hinge-EPS-TSVR
#'
#' @author Zhang Jiaqi.
#' @param X,y dataset and label.
#' @param C1,C2 weight of loss term.
#' @param C3,C4 weight of regularization term.
#' @param epsilon1,epsilon2 parameter for epsilon tube.
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
#' @param reduce_set reduce set for reduce SVM, default \code{reduce_set = NULL}.
#' @return return \code{TSVRClassifier} object.
#' @export
hinge_eps_tsvr <- function(X, y, C1 = 1, C2 = C1, C3 = 1, C4 = C3,
                       epsilon1 = 0.1, epsilon2 = epsilon1,
                       kernel = c("linear", "rbf", "poly"),
                       gamma = 1 / ncol(X), degree = 3, coef0 = 0,
                       eps = 1e-7, max.steps = 4000, fit_intercept = TRUE,
                       reduce_set = NULL) {
  X <- as.matrix(X)
  y <- as.matrix(y)
  kernel <- match.arg(kernel)
  if (kernel != "linear") {
    kso <- kernel_select_option(X, kernel, reduce_set, gamma, degree, coef0)
    KernelX <- kso$KernelX
  } else {
    KernelX <- X
  }
  if (fit_intercept == TRUE) {
    KernelX <- cbind(KernelX, 1)
  }
  solver.res <- hinge_eps_tsvr_dual_solver(KernelX, y, C1, C2, C3, C4, epsilon1, epsilon2,
                                           eps, max.steps)
  EPSTSVRRegressor <- list("X" = X, "y" = y,
                         "C1" = C1, "C2" = C2,
                         "epsilon1" = epsilon1, "epsilon2" = epsilon2,
                         "kernel" = kernel,
                         "gamma" = gamma, "degree" = degree, "coef0" = coef0,
                         "coef1" = solver.res$coef1,
                         "coef2" = solver.res$coef2,
                         "fit_intercept" = fit_intercept,
                         "solver.res" = solver.res)
  class(EPSTSVRRegressor) <- "EPSTSVRRegressor"
  return(EPSTSVRRegressor)
}


#' Predict Method for Twin Support Vector Regression
#'
#' @author Zhang Jiaqi
#' @param object a fitted object of class inheriting from \code{SVMClassifier}.
#' @param X new data for predicting.
#' @param ... unused parameter.
#' @importFrom stats predict
#' @export
predict.EPSTSVRRegressor <- function(object, X, ...) {
  X <- as.matrix(X)
  if (object$kernel == "linear") {
    KernelX <- X
  } else {
    KernelX <- kernel_function(X, object$X,
                               kernel.type = object$kernel,
                               gamma = object$gamma,
                               degree = object$degree,
                               coef0 = object$coef0)
  }
  if (object$fit_intercept == TRUE) {
    KernelX <- cbind(KernelX, 1)
  }
  f <- KernelX %*% (object$coef1 + object$coef2) / 2
  return(f)
}


#' Plot Method for Twin Support Vector Regression
#'
#' @author Zhang Jiaqi
#' @param x a fitted object of class inheriting from \code{EPSTSVRRegressor}.
#' @param ... unused parameter.
#' @importFrom graphics abline grid points
#' @export
plot.EPSTSVRRegressor <- function(x, ...) {
  xp <- ncol(x$X)
  xlim_c <- c(min(x$X[,1]), max(x$X[, 1]))
  ylim_c <- c(min(x$y), max(x$y))
  if (xp == 1) {
    plot(x$X, x$y, xlim = xlim_c, ylim = ylim_c,
         xlab = "x", ylab = "y")
    grid(lwd = 2,col = "grey")
    abline(x$coef1[2], x$coef1[1], col = "blue")
    abline(x$coef1[2] - x$epsilon1, x$coef1[1], col = "blue", lty = 2)
    abline(x$coef2[2], x$coef2[1], col = "blue")
    abline(x$coef2[2] + x$epsilon2, x$coef2[1], col = "blue", lty = 2)
    abline((x$coef1[2] + x$coef2[2]) / 2, (x$coef1[1] + x$coef2[1]) / 2, col = "black")
  }
}
