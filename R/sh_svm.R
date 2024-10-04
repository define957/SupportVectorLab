sh_svm_dual_solver <- function(KernelX, y, C = 1,
                                  eps = 1e-5, max.steps = 80) {
  n <- nrow(KernelX)
  H <- calculate_svm_H(KernelX, y)
  H <- H + diag(1/C, n)
  e <- matrix(1, n)
  lb <- matrix(0, n)
  ub <- matrix(Inf, n)
  u0 <- lb
  alphas <- clip_dcd_optimizer(H, e, lb, ub, eps, max.steps, u0)$x
  coef <- y*alphas
  BaseDualSquaredHingeSVMClassifier <- list(coef = as.matrix(coef))
  class(BaseDualSquaredHingeSVMClassifier) <- "BaseDualSquaredHingeSVMClassifier"
  return(BaseDualSquaredHingeSVMClassifier)
}


sh_svm_primal_solver <- function(KernelX, y, C = 1,
                                 max.steps = 80, batch_size = nrow(KernelX) / 10,
                                 optimizer = pegasos, ...) {
  gSquaredHinge <- function(KernelX, y, w, pars, ...) {
    C <- pars$C
    xn <- pars$xn
    xmn <- nrow(KernelX)
    xmp <- ncol(KernelX)
    sg <- matrix(0, nrow = xmp, ncol = 1)
    u <- 1 - y * (KernelX %*% w)
    u[u <  0] <- 0
    u[u >= 0] <- u
    g <- w - (C*xn/xmn) * t(KernelX) %*% (u*y)
    return(g)
  }
  xn <- nrow(KernelX)
  xp <- ncol(KernelX)
  w0 <- matrix(0, nrow = xp, ncol = 1)
  pars <- list("C" = C, "xn" = xn)
  wt <- optimizer(KernelX, y, w0, batch_size, max.steps, gSquaredHinge, pars, ...)
  BasePrimalSquaredHingeSVMClassifier <- list(coef = as.matrix(wt[1:xp]))
  class(BasePrimalSquaredHingeSVMClassifier) <- "BasePrimalSquaredHingeSVMClassifier"
  return(BasePrimalSquaredHingeSVMClassifier)
}


#' Squared Hinge Loss Support Vector Machine
#'
#' \code{sh_svm} is an R implementation of Squared Hinge-SVM
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
#' @return return \code{HingeSVMClassifier} object.
#' @export
sh_svm <- function(X, y, C = 1, kernel = c("linear", "rbf", "poly"),
                   gamma = 1 / ncol(X), degree = 3, coef0 = 0,
                   eps = 1e-5, max.steps = 80, batch_size = nrow(X) / 10,
                   solver = c("dual", "primal"),
                   fit_intercept = TRUE, optimizer = pegasos, randx = 0.1, ...) {
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
  solver <- match.arg(solver)
  if (fit_intercept == TRUE) {
    X <- cbind(X, 1)
  }
  kso <- kernel_select_option(X, kernel, solver, randx,
                              gamma, degree, coef0)
  KernelX <- kso$KernelX
  X <- kso$X
  if (solver == "primal") {
    solver.res <- sh_svm_primal_solver(KernelX, y, C,
                                       max.steps, batch_size,
                                       optimizer, ...)
  } else if (solver == "dual") {
    solver.res <- sh_svm_dual_solver(KernelX, y, C, eps, max.steps)
  }
  SVMClassifier <- list("X" = X, "y" = y, "class_set" = class_set,
                        "C" = C, "kernel" = kernel,
                        "gamma" = gamma, "degree" = degree, "coef0" = coef0,
                        "solver" = solver, "coef" = solver.res$coef,
                        "fit_intercept" = fit_intercept)
  class(SVMClassifier) <- "SVMClassifier"
  return(SVMClassifier)
}
