ls_svr_dual_solver <- function(KernelX, y, C = 1) {
  coef <- cholsolve(KernelX + diag(1/C, nrow = ncol(KernelX)), y)
  BaseDualLSSVMRegressor <- list(coef = as.matrix(coef))
  class(BaseDualLSSVMRegressor) <- "BaseDualLSSVMRegressor"
  return(BaseDualLSSVMRegressor)
}

ls_svr_primal_solver <- function(KernelX, y, C, max.steps, batch_size,
                                 optimizer, ...) {
  gLeastSquares <- function(KernelX, y, w, pars,...) {
    # gradient of Least Squares loss function
    C <- pars$C
    xn <- pars$xn
    xmn <- nrow(KernelX)
    xmp <- ncol(KernelX)
    g <- matrix(xmp, xp, 1)
    g <- w - (C*xn/xmn) * t(KernelX) %*% (y - KernelX %*% w)
    return(g)
  }
  xn <- nrow(KernelX)
  xp <- ncol(KernelX)
  w0 <- matrix(0, xp, 1)
  pars <- list("C" = C, "xn" = xn)
  wt <- optimizer(KernelX, y, w0, batch_size, max.steps, gLeastSquares, pars, ...)
  BasePrimalLeastSquaresSVMRegressor <- list(coef = as.matrix(wt[1:xp]))
  class(BasePrimalLeastSquaresSVMRegressor) <- "BasePrimalLeastSquaresSVMRegressor"
  return(BasePrimalLeastSquaresSVMRegressor)
}

#' Least Squares Support Vector Regression (Ridge/Kernel Ridge)
#'
#' \code{ls_svr} is an R implementation of Least Squares SVR.
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
ls_svr <- function(X, y, C = 1, kernel = c("linear", "rbf", "poly"),
                   gamma = 1 / ncol(X), degree = 3, coef0 = 0,
                   eps = 1e-5, max.steps = 80, batch_size = nrow(X) / 10,
                   solver = c("dual", "primal"),
                   fit_intercept = TRUE, optimizer = pegasos, randx = 0.1, ...) {
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
    solver.res <- ls_svr_primal_solver(KernelX, y, C,
                                       max.steps, batch_size,
                                       optimizer, ...)
  } else if (solver == "dual") {
    solver.res <- ls_svr_dual_solver(KernelX, y, C)
  }
  SVMRegressor <- list("X" = X, "y" = y,
                       "C" = C, "kernel" = kernel,
                       "gamma" = gamma, "degree" = degree, "coef0" = coef0,
                       "solver" = solver, "coef" = solver.res$coef,
                       "fit_intercept" = fit_intercept)
  class(SVMRegressor) <- "SVMRegressor"
  return(SVMRegressor)
}
