hinge_wtsvr_dual_solver <- function(KernelX, y, D1, D2, C1, C2, epsilon1, epsilon2,
                                    eps, max.steps) {
  xn <- nrow(KernelX)
  xp <- ncol(KernelX)
  G <- KernelX
  GTD1G_inv_GT <- cholsolve(t(G) %*% (D1*G) + diag(1e-7, xp), t(G))
  dualH <- G %*% GTD1G_inv_GT

  f <- as.matrix(y - epsilon1)
  h <- as.matrix(y + epsilon2)

  q1 <- dualH %*% (D1*f) - f
  q2 <- h - dualH %*% (D1*h)

  lb <- matrix(0, xn, 1)
  ub1 <- matrix(C1*D2, xn, 1)
  ub2 <- matrix(C2*D2, xn, 1)

  x0 <- lb

  dual_coef1 <- clip_dcd_optimizer(dualH, q1, lb, ub1, eps, max.steps, x0)$x
  dual_coef2 <- clip_dcd_optimizer(dualH, q2, lb, ub2, eps, max.steps, x0)$x

  coef1      <- GTD1G_inv_GT %*% (D1*f - dual_coef1)
  coef2      <- GTD1G_inv_GT %*% (D1*h + dual_coef2)

  BaseDualHingeTSVRRegressor <- list("coef1" = as.matrix(coef1),
                                     "coef2" = as.matrix(coef2),
                                     "lag1" = dual_coef1,
                                     "lag2" = dual_coef2)
}

#' Hinge Weighted Twin Support Vector Regression
#'
#' \code{hinge_wtsvr} is an R implementation of Hinge Weighted Twin Support Vector Regression
#'
#' @author Zhang Jiaqi.
#' @param X,y dataset and label.
#' @param C1,C2 plenty term.
#' @param epsilon1,epsilon2 parameter for epsilon tube.
#' @param D1,D2 weight vectors for least squares loss and hinge loss, respectively, applied per sample.
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
#' @param weight_f1,weight_f2 optional functions to compute sample weights for the
#'   least squares loss and hinge loss, respectively. Default is NULL, meaning the
#'   corresponding \code{D1} or \code{D2} vector is used directly. When a function is
#'   supplied and return a numeric vector of length \code{nrow(X)}. If provided, the relevant \code{D1} or
#'   \code{D2} argument is ignored.
#' @param weight_option1,weight_option2 optional named lists of additional arguments
#'   to be passed to \code{weight_f1} and \code{weight_f2}, respectively. Default is
#'   NULL (no extra arguments). These are passed via \code{do.call}, for example
#'   \code{weight_option1 = list(k = 2)} if \code{weight_f1} expects a \code{k}
#'   parameter.
#' @return return \code{TSVMRegressor} object.
#' @export
hinge_wtsvr <- function(X, y, C1 = 1, C2 = C1,
                        epsilon1 = 0.1, epsilon2 = epsilon1,
                        D1 = rep(1, nrow(X)), D2 = rep(1, nrow(X)),
                        kernel = c("linear", "rbf", "poly"),
                        gamma = 1 / ncol(X), degree = 3, coef0 = 0,
                        eps = 1e-7, max.steps = 4000, fit_intercept = TRUE,
                        reduce_set = NULL,
                        weight_f1 = NULL, weight_f2 = NULL,
                        weight_option1 = NULL, weight_option2 = NULL) {

  X <- as.matrix(X)
  y <- as.matrix(y)

  kernel  <- match.arg(kernel)
  KernelR <- NULL

  if (kernel != "linear") {
    kso <- kernel_select_option_(X, kernel, reduce_set, gamma, degree, coef0)
    KernelX <- kso$KernelX
    KernelR <- kso$KernelR
  } else {
    KernelX <- X
  }
  kxp <- ncol(KernelX)
  if (fit_intercept == TRUE) {
    KernelX <- cbind(KernelX, 1)
  }

  if (!is.null(weight_f1)) {
    D1 <- do.call(weight_f1, append(list("X" = X, "y" = y), weight_option1))
  }
  if (!is.null(weight_f2)) {
    D2 <- do.call(weight_f2, append(list("X" = X, "y" = y), weight_option2))
  }
  solver.res <- hinge_wtsvr_dual_solver(KernelX, y, D1, D2, C1, C2, epsilon1, epsilon2,
                                        eps, max.steps)
  model_specs   <- list("X" = X, "y" = y,
                        "C1" = C1, "C2" = C2,
                        "epsilon1" = epsilon1, "epsilon2" = epsilon2,
                        "fit_intercept" = fit_intercept)
  model_coef    <- list("coef1" = solver.res$coef1,
                        "coef2" = solver.res$coef2)
  kernel_config <- list("kernel" = kernel,
                        "gamma"  = gamma,
                        "degree" = degree,
                        "coef0" = coef0,
                        "reduce_set" = reduce_set,
                        "KernelR" = KernelR,
                        "KernelX" = KernelX[, 1:kxp, drop = FALSE])

  TSVMRegressor <- structure(list("model_specs" = model_specs,
                                  "model_coef" = model_coef,
                                  "kernel_config" = kernel_config),
                             "class" = "TSVMRegressor")
  return(TSVMRegressor)
}
