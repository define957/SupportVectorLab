hinge_w_eps_tsvr_dual_solver <- function(KernelX, y, D1, D2, C1, C2, C3, C4, epsilon1, epsilon2,
                                       eps, max.steps) {
  xn <- nrow(KernelX)
  xp <- ncol(KernelX)

  G               <- KernelX
  GramGD1         <- t(G) %*% (D1*G)
  GTD1G_C3_inv_GT <- cholsolve(GramGD1 + diag(C3, xp), t(G))
  dualH1          <- G %*% GTD1G_C3_inv_GT

  if (C3 != C4) {
    GTD1G_C4_inv_GT <- cholsolve(GramGD1 + diag(C4, xp), t(G))
    dualH2          <- G %*% GTD1G_C4_inv_GT
  } else {
    GTD1G_C4_inv_GT <- GTD1G_C3_inv_GT
    dualH2          <- dualH1
  }

  q1  <- dualH1 %*% (D1*y) - y - epsilon1
  q2  <- y - epsilon2 - dualH2 %*% (D1*y)
  lb  <- matrix(0, xn, 1)
  ub1 <- matrix(C1*D2, xn, 1)
  ub2 <- matrix(C2*D2, xn, 1)

  u0 <- lb

  dual_coef1 <- clip_dcd_optimizer(dualH1, q1, lb, ub1, eps, max.steps, u0)$x
  dual_coef2 <- clip_dcd_optimizer(dualH2, q2, lb, ub2, eps, max.steps, u0)$x

  coef1      <- GTD1G_C3_inv_GT %*% (D1*y - dual_coef1)
  coef2      <- GTD1G_C4_inv_GT %*% (D1*y + dual_coef2)

  BaseDualHingeEPSTSVRRegressor <- list("coef1" = as.matrix(coef1),
                                        "coef2" = as.matrix(coef2))
}

#' Weighted Hinge Epsilon Twin Support Vector Regression
#'
#' \code{hinge_w_eps_tsvr} is an R implementation of Hinge-EPS-TSVR
#'
#' @author Zhang Jiaqi.
#' @param X,y dataset and label.
#' @param C1,C2 weight of loss term.
#' @param C3,C4 weight of regularization term.
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
#' @return return \code{EPSTSVMRegressor} object.
#' @export
hinge_w_eps_tsvr <- function(X, y, C1 = 1, C2 = C1, C3 = 1e-7, C4 = C3,
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
  solver.res <- hinge_w_eps_tsvr_dual_solver(KernelX, y, D1, D2, C1, C2, C3, C4, epsilon1, epsilon2,
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

  EPSTSVMRegressor <- structure(list("model_specs" = model_specs,
                                     "model_coef" = model_coef,
                                     "kernel_config" = kernel_config),
                                "class" = "EPSTSVMRegressor")
  return(EPSTSVMRegressor)
}
