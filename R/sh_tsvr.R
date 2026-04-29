sh_tsvr_dual_solver <- function(KernelX, y, C1, C2, epsilon1, epsilon2,
                                eps, max.steps) {
  xn <- nrow(KernelX)
  xp <- ncol(KernelX)

  G            <- KernelX
  GTG_inv_G    <- cholsolve(t(G) %*% G + diag(1e-7, xp), t(G))
  dualH        <- G %*% GTG_inv_G
  dualH1       <- dualH
  dualH2       <- dualH
  diag(dualH1) <- diag(dualH1) + 1/C1
  diag(dualH2) <- diag(dualH2) + 1/C2

  f <- as.matrix(y - epsilon1)
  h <- as.matrix(y + epsilon2)

  q1 <- dualH %*% f - f
  q2 <- h - dualH %*% h

  lb <- matrix(0, xn, 1)
  ub <- matrix(Inf, xn, 1)

  x0 <- lb

  dual_coef1 <- clip_dcd_optimizer(dualH1, q1, lb, ub, eps, max.steps, x0)$x
  dual_coef2 <- clip_dcd_optimizer(dualH2, q2, lb, ub, eps, max.steps, x0)$x

  coef1      <- GTG_inv_G %*% (f - dual_coef1)
  coef2      <- GTG_inv_G %*% (h + dual_coef2)

  BaseDualTSVRRegressor <- list("coef1" = as.matrix(coef1),
                                "coef2" = as.matrix(coef2),
                                "lag1" = dual_coef1,
                                "lag2" = dual_coef2)
}

#' Squared Hinge Twin Support Vector Regression
#'
#' \code{sh_tsvr} is an R implementation of SH-TSVR
#'
#' @author Zhang Jiaqi.
#' @param X,y dataset and label.
#' @param C1,C2 plenty term.
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
#' @return return \code{TSVMRegressor} object.
#' @export
sh_tsvr <- function(X, y, C1 = 1, C2 = C1,
                    epsilon1 = 0.1, epsilon2 = epsilon1,
                    kernel = c("linear", "rbf", "poly"),
                    gamma = 1 / ncol(X), degree = 3, coef0 = 0,
                    eps = 1e-7, max.steps = 4000, fit_intercept = TRUE,
                    reduce_set = NULL) {

  X <- as.matrix(X)
  y <- as.matrix(y)

  kernel  <- match.arg(kernel)
  KernelR <- NULL

  kernel <- match.arg(kernel)
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
  solver.res <- sh_tsvr_dual_solver(KernelX, y, C1, C2, epsilon1, epsilon2,
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
