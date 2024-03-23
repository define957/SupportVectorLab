als_svr_dual_solver <- function(KernelX, y, C = 1, p = 0.5,
                                eps = 1e-5, max.steps = 4000) {
  n <- nrow(KernelX)
  H <- KernelX
  H <- cbind(H, -H)
  H <- rbind(H, -H)
  I1 <- diag(1/(p), n, n)
  I12 <- matrix(0, n, n)
  I4 <- diag(1/(1 - p), n, n)
  In <- rbind(cbind(I1, I12), cbind(I12, I4))
  H <- H + In/C
  q <- rbind(y, -y)
  lb <- matrix(0, 2*n)
  ub <- matrix(Inf, 2*n)
  u0 <- lb
  u <- clip_dcd_optimizer(H, q, lb, ub, eps, max.steps, u0)$x
  coef <- u[1:n] - u[(n + 1):(2*n)]
  BaseDualALSSVRRegressor <- list(coef = as.matrix(coef))
  class(BaseDualALSSVRRegressor) <- "BaseDualALSSVRRegressor"
  return(BaseDualALSSVRRegressor)
}


als_svr_primal_solver <- function(KernelX, y, C = 1, p = 0.5,
                                  max.steps = 80, batch_size = nrow(KernelX) / 10,
                                  optimizer = pegasos, ...) {
  gALS <- function(KernelX, y, w, pars,...) {
    C <- pars$C
    p <- pars$p
    xn <- pars$xn
    xmn <- nrow(KernelX)
    xmp <- ncol(KernelX)
    sg <- matrix(0, xmp, 1)
    u <- matrix(1, xmn)
    e <- y -  KernelX %*% w
    idx <- which(e >= 0)
    u[idx] <- p
    u[-idx] <- (1 - p)
    g <- w - (C*xn/xmn) * t(KernelX) %*% (u*e)
    print(g)
    return(g)
  }
  xn <- nrow(KernelX)
  xp <- ncol(KernelX)
  w0 <- matrix(0, nrow = xp, ncol = 1)
  pars <- list("C" = C, "xn" = xn, "p" = p)
  wt <- optimizer(KernelX, y, w0, batch_size, max.steps, gALS, pars, ...)
  BasePrimalALSSVRRegressor <- list(coef = as.matrix(wt[1:xp]))
  class(BasePrimalALSSVRRegressor) <- "BasePrimalALSSVRRegressor"
  return(BasePrimalALSSVRRegressor)
}



#' Asymmetric Least Squares Support Vector Regression
#'
#' \code{als_svr} is an R implementation of ALS-SVR
#'
#' @author Zhang Jiaqi.
#' @param X,y dataset and label.
#' @param C plenty term.
#' @param kernel kernel function. The definitions of various kernel functions are as follows:
#' \describe{
#'     \item{linear:}{\eqn{u'v}{u'*v}}
#'     \item{poly:}{\eqn{(\gamma u'v + coef0)^{degree}}{(gamma*u'*v + coef0)^degree}}
#'     \item{rbf:}{\eqn{e^{(-\gamma |u-v|^2)}}{exp(-gamma*|u-v|^2)}}
#' }
#' @param p parameter for asymmetric least squares loss.
#' @param gamma parameter for \code{'rbf'} and \code{'poly'} kernel. Default \code{gamma = 1/ncol(X)}.
#' @param degree parameter for polynomial kernel, default: \code{degree = 3}.
#' @param coef0 parameter for polynomial kernel,  default: \code{coef0 = 0}.
#' @param eps the precision of the optimization algorithm.
#' @param max.steps the number of iterations to solve the optimization problem.
#' @param batch_size mini-batch size for primal solver.
#' @param solver \code{"dual"} and \code{"primal"} are available.
#' @param fit_intercept if set \code{fit_intercept = TRUE},
#'                      the function will evaluates intercept.
#' @param optimizer default primal optimizer pegasos.
#' @param randx parameter for reduce SVM, default \code{randx = 0.1}.
#' @param ... unused parameters.
#' @return return \code{SVMRegressor} object.
#' @export
als_svr <- function(X, y, C = 1, kernel = c("linear", "rbf", "poly"),
                    gamma = 1 / ncol(X), degree = 3, coef0 = 0,
                    p = 0.5,
                    eps = 1e-5, max.steps = 4000, batch_size = nrow(X) / 10,
                    solver = c("dual", "primal"),
                    fit_intercept = TRUE, optimizer = pegasos, randx = 0.1, ...) {
  if (p >= 1 || p <= 0) {
    stop("Parameter p should satisfy 0 < p < 1")
  }
  X <- as.matrix(X)
  y <- as.matrix(y)
  kernel <- match.arg(kernel)
  solver <- match.arg(solver)
  if (fit_intercept == TRUE) {
    X <- cbind(X, 1)
  }
  kso <- kernel_select_option(X, kernel, solver, randx,
                              gamma, degree, coef0)
  KernelX <- kso$KernelX
  X <- kso$X
  if (solver == "primal") {
    solver.res <- als_svr_primal_solver(KernelX, y, C, p,
                                        max.steps, batch_size,
                                        optimizer, ...)
  } else if (solver == "dual") {
    solver.res <- als_svr_dual_solver(KernelX, y, C,
                                      p, eps, max.steps)
  }
  SVMRegressor <- list("X" = X, "y" = y,
                        "C" = C, "kernel" = kernel,
                        "gamma" = gamma, "degree" = degree, "coef0" = coef0,
                        "solver" = solver, "coef" = solver.res$coef,
                        "fit_intercept" = fit_intercept)
  class(SVMRegressor) <- "SVMRegressor"
  return(SVMRegressor)
}
