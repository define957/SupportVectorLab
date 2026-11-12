rhinge_eps_tsvr_dual_solver <- function(KernelX, y, C1, C2, C3, C4, epsilon1, epsilon2,
                                        eta1, eta2,
                                        eps, eps.hq,
                                        max.steps, hq.steps, method_ls) {

  ls_loss_function <- function(z) {
    return(sum(z^2))
  }
  rescaled_hinge_loss_function <- function(z, eta, beta) {
    return(sum(-beta * expm1(-eta * pmax(z, 0))))
  }

  update_weight <- function(z, eta, beta) {
    weight <- matrix(0, length(z), 1)
    weight <- beta * eta * exp(-eta*pmax(z, 0))
    return(weight)
  }

  compute_objective1 <- function(coef, residuals_f, C1, C3, epsilon1, eta1, beta1) {
    obj <- 0.5*(C3 * crossprod(coef) + ls_loss_function(residuals_f)) +
      C1 * rescaled_hinge_loss_function(-residuals_f - epsilon1, eta1, beta1)
    return(obj)
  }

  compute_objective2 <- function(coef, residuals_f, C2, C4, epsilon2, eta2, beta2) {
    obj <- 0.5*(C4 * crossprod(coef) + ls_loss_function(residuals_f)) +
      C2 * rescaled_hinge_loss_function(residuals_f - epsilon2, eta2, beta2)
    return(obj)
  }

  beta1 <- 1 / (1 - exp(-eta1))
  beta2 <- 1 / (1 - exp(-eta2))
  xn    <- nrow(KernelX)
  xp    <- ncol(KernelX)

  G             <- KernelX
  if (method_ls == "qr") {
    GTG_C3_inv_GT <- stack_qrsolve(G, t(G), C3)
  } else if (method_ls == "normal-chol") {
    GramG         <- crossprod(G)
    GTG_C3_inv_GT <- cholsolve(GramG + diag(C3, xp), t(G))
  }
  dualH1 <- G %*% GTG_C3_inv_GT

  if (C3 != C4) {
    if (method_ls == "qr") {
      GTG_C4_inv_GT <- stack_qrsolve(G, t(G), C4)
    } else if (method_ls == "normal-chol") {
      GTG_C4_inv_GT <- cholsolve(GramG + diag(C4, xp), t(G))
    }
    dualH2        <- G %*% GTG_C4_inv_GT
  } else {
    GTG_C4_inv_GT <- GTG_C3_inv_GT
    dualH2        <- dualH1
  }

  q1    <- dualH1 %*% y - y - epsilon1
  q2    <- y - epsilon2 - dualH2 %*% y

  lb    <- matrix(0, xn, 1)

  ub1   <- matrix(C1, xn, 1)
  ub2   <- matrix(C2, xn, 1)

  a01   <- lb
  a02   <- lb

  coef1 <- matrix(0, xp)
  coef2 <- matrix(0, xp)

  obj1 <- NULL
  for (i in 1:hq.steps) {
    a1            <- clip_dcd_optimizer(dualH1, q1, lb, ub1, eps, max.steps, a01)$x
    coef1_new     <- GTG_C3_inv_GT %*% (y - a1)
    residuals_f   <- y - KernelX %*% coef1_new
    obj1_new      <- compute_objective1(coef1_new, residuals_f, C1, C3, epsilon1, eta1, beta1)
    if (!is.null(obj1) && function_tolerance(obj1_new, obj1, eps.hq)) {
      coef1 <- coef1_new
      break
    }
    ub1            <- C1*update_weight(-residuals_f - epsilon1, eta1, beta1)
    coef1          <- coef1_new
    obj1           <- obj1_new
    a01            <- a1
    a01[a01 > ub1] <- ub1[a01 > ub1]
    a01[a01 < lb]  <- lb[a01 < lb]
  }

  obj2 <- NULL
  for (i in 1:hq.steps) {
    a2            <- clip_dcd_optimizer(dualH2, q2, lb, ub2, eps, max.steps, a02)$x
    coef2_new     <- GTG_C4_inv_GT %*% (y + a2)
    residuals_f   <- y - KernelX %*% coef2_new
    obj2_new      <- compute_objective2(coef2_new, residuals_f, C2, C4, epsilon2, eta2, beta2)
    if (!is.null(obj2) && function_tolerance(obj2_new, obj2, eps.hq)) {
      coef2 <- coef2_new
      break
    }
    ub2            <- C2*update_weight(residuals_f - epsilon2, eta2, beta2)
    coef2          <- coef2_new
    obj2           <- obj2_new
    a02            <- a2
    a02[a02 > ub2] <- ub2[a02 > ub2]
    a02[a02 < lb]  <- lb[a02 < lb]
  }

  BaseDualRHingeEPSTSVRRegressor <- list("coef1" = as.matrix(coef1),
                                         "coef2" = as.matrix(coef2))
}

#' Rescaled Hinge Epsilon Twin Support Vector Regression
#'
#' \code{rhinge_eps_tsvr} is an R implementation of Res-TSVR, a robust
#' twin support vector regression that replaces the Hinge loss of
#' \code{\link{hinge_eps_tsvr}} with the rescaled Hinge loss
#' \eqn{P_{\text{rhinge}}(\xi) = \beta(1-\exp(-\hat\eta \xi))}{P_rhinge(xi) = beta*(1-exp(-eta_hat*xi))},
#' where \eqn{\beta = 1/(1-\exp(-\hat\eta))}{beta = 1/(1-exp(-eta_hat))} and
#' \eqn{\xi = \max(0, \cdot)}{xi = max(0, .)} is the epsilon-insensitive hinge
#' slack. The loss is bounded, so outliers are down-weighted instead of
#' dominating the solution.
#'
#' The non-convex problem is solved by the half-quadratic technique: each
#' outer iteration solves the weighted dual quadratic programs with
#' \code{\link{clip_dcd_optimizer}} and updates the sample weights
#' \eqn{\beta\hat\eta\exp(-\hat\eta \xi_i)}{beta*eta_hat*exp(-eta_hat*xi_i)}.
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
#' @param eta1,eta2 rescaled Hinge loss constants \eqn{\hat\eta}{eta-hat}. As
#' \eqn{\hat\eta \to 0}{eta-hat -> 0} the loss approaches the conventional
#' Hinge loss.
#' @param eps the precision of the optimization algorithm.
#' @param eps.hq the tolerance of the stopping criterion: the half-quadratic
#'             iterations stop when the relative change of the objective
#'             function value is below \code{eps.hq}.
#' @param max.steps the number of iterations to solve the optimization problem.
#' @param hq.steps the maximum number of half-quadratic iterations. The solver
#'             first solves the unit-weight dual problems as an initialization,
#'             then performs at most \code{hq.steps} reweighted iterations.
#' @param fit_intercept if set \code{fit_intercept = TRUE},
#'                      the function will evaluates intercept.
#' @param reduce_set reduce set for reduce SVM, default \code{reduce_set = NULL}.
#' @param method_ls solver for the internal linear systems: \code{"qr"} (stacked
#'                  QR decomposition, never forms Gram matrices) or
#'                  \code{"normal-chol"} (Cholesky on the normal equations).
#' @return return \code{EPSTSVMRegressor} object.
#' @references Manisha Singla, Debdas Ghosh, K.K. Shukla, Witold Pedrycz (2020).
#'             Robust twin support vector regression based on rescaled Hinge loss.
#'             \emph{Pattern Recognition}, 105, 107395.
#' @export
rhinge_eps_tsvr <- function(X, y, C1 = 1, C2 = C1, C3 = 1e-5, C4 = C3,
                            epsilon1 = 0.1, epsilon2 = epsilon1,
                            kernel = c("linear", "rbf", "poly"),
                            gamma = 1 / ncol(X), degree = 3, coef0 = 0,
                            eta1 = 1, eta2 = eta1,
                            eps = .Machine$double.eps, eps.hq = 1e-4,
                            max.steps = 1e+8, hq.steps = 500,
                            fit_intercept = TRUE,
                            reduce_set = NULL,
                            method_ls = c("qr", "normal-chol")) {

  method_ls <- match.arg(method_ls)
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
  solver.res <- rhinge_eps_tsvr_dual_solver(KernelX, y, C1, C2, C3, C4, epsilon1, epsilon2,
                                            eta1, eta2,
                                            eps, eps.hq,
                                            max.steps, hq.steps, method_ls)

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
