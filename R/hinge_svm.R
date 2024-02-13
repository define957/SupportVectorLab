hinge_svm_dual_solver <- function(KernelX, y, C = 1,
                                  eps = 1e-5, max.steps = 80) {
  n <- nrow(KernelX)
  H <- calculate_svm_H(KernelX, y)
  e <- matrix(1, nrow = n)
  lb <- matrix(0, nrow = n)
  ub <- matrix(C, nrow = n)
  u0 <- lb
  alphas <- clip_dcd_optimizer(H, e, lb, ub, eps, max.steps, u0)$x #
  coef <- y*alphas
  BaseDualHingeSVMClassifier <- list(coef = as.matrix(coef))
  class(BaseDualHingeSVMClassifier) <- "BaseDualHingeSVMClassifier"
  return(BaseDualHingeSVMClassifier)
}


hinge_svm_primal_solver <- function(KernelX, y, C = 1,
                                    max.steps = 80, batch_size = nrow(KernelX) / 10,
                                    optimizer = pegasos, ...) {
  sgHinge <- function(KernelX, y, w, pars, ...) { # sub-gradient of hinge loss function
    C <- pars$C
    xn <- pars$xn
    xmn <- nrow(KernelX)
    xmp <- ncol(KernelX)
    sg <- matrix(0, nrow = xmp, ncol = 1)
    u <- 1 - y * (KernelX %*% w)
    u[u <= 0] <- 0
    u[u > 0] <- 1
    sg <- w - (C*xn/xmn) * t(KernelX) %*% (u*y)
    return(sg)
  }
  xn <- nrow(KernelX)
  xp <- ncol(KernelX)
  w0 <- matrix(0, nrow = xp, ncol = 1)
  pars <- list("C" = C, "xn" = xn)
  wt <- optimizer(KernelX, y, w0, batch_size, max.steps, sgHinge, pars, ...)
  BasePrimalHingeSVMClassifier <- list(coef = as.matrix(wt[1:xp]))
  class(BasePrimalHingeSVMClassifier) <- "BasePrimalHingeSVMClassifier"
  return(BasePrimalHingeSVMClassifier)
}


#' Hinge Support Vector Machine
#'
#' \code{hinge_svm} is an R implementation of Hinge-SVM
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
hinge_svm <- function(X, y, C = 1, kernel = c("linear", "rbf", "poly"),
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
    solver.res <- hinge_svm_primal_solver(KernelX, y, C,
                                          max.steps, batch_size,
                                          optimizer, ...)
  } else if (solver == "dual") {
    solver.res <- hinge_svm_dual_solver(KernelX, y, C, eps,
                                        max.steps)
  }
  SVMClassifier <- list("X" = X, "y" = y, "class_set" = class_set,
                        "C" = C, "kernel" = kernel,
                        "gamma" = gamma, "degree" = degree, "coef0" = coef0,
                        "solver" = solver, "coef" = solver.res$coef,
                        "fit_intercept" = fit_intercept)
  class(SVMClassifier) <- "SVMClassifier"
  return(SVMClassifier)
}


#' Plot Method for Support Vector Machine
#'
#' @author Zhang Jiaqi
#' @param x a fitted object of class inheriting from \code{SVMClassifier}.
#' @param ... unused parameter.
#' @importFrom stats coef
#' @importFrom graphics abline grid points
#' @export
plot.SVMClassifier <- function(x, ...) {
  coefs <- coef(x)
  idx <- which(x$y == 1)
  xlim_c <- c(min(x$X[,1]), max(x$X[, 1]))
  ylim_c <- c(min(x$X[,2]), max(x$X[, 2]))
  if (ncol(x$X) == 3) {
    plot(x$X[idx, 1], x$X[idx, 2], col = "red", xlim = xlim_c, ylim = ylim_c,
         xlab = "", ylab = "")
    grid(10, 10, lwd = 2,col = "grey")
    points(x$X[-idx, 1], x$X[-idx, 2], col = "blue")
    if (x$kernel == "linear") {
      abline(a = -coefs[3]/coefs[2], b = -coefs[1]/coefs[2],
             lty = 2)
      abline(a = -(coefs[3] + 1)/coefs[2], b = -coefs[1]/coefs[2],
             lty = 2)
      abline(a = -(coefs[3] - 1)/coefs[2], b = -coefs[1]/coefs[2],
             lty = 2)
    }
  }
}

#' Coef Method for Support Vector Machine
#'
#' @author Zhang Jiaqi
#' @param object a fitted object of class inheriting from \code{SVMClassifier}.
#' @param ... unused parameter.
#' @importFrom stats coef
#' @export
coef.SVMClassifier <- function(object, ...) {
  if (object$solver == "dual") {
    return(t(object$X) %*% object$coef)
  } else if (object$solver == "primal") {
    return(object$coef)
  }
}

#' Predict Method for Support Vector Machine
#'
#' @author Zhang Jiaqi
#' @param object a fitted object of class inheriting from \code{SVMClassifier}.
#' @param X new data for predicting.
#' @param values if set \code{values = TRUE}, this function will return predict
#'               values: f = wx+b.
#' @param ... unused parameter.
#' @importFrom stats predict
#' @export
predict.SVMClassifier <- function(object, X, values = FALSE, ...) {
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
  fx <- KernelX %*% object$coef
  if (values == FALSE) {
    decf <- sign(fx)
    idx_pos <- which(decf >= 0)
    idx_neg <- which(decf < 0)
    decf[idx_pos] <- object$class_set[1]
    decf[idx_neg] <- object$class_set[2]
  } else {
    decf <- fx
  }
  return(decf)
}

calculate_svm_H <- function(KernelX, y) {
  H <- (y %*% t(y))*KernelX
  return(H)
}
