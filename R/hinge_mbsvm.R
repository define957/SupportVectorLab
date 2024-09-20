hinge_mbsvm_dual_solver <- function(KernelX, y, C, class_set, class_num,
                                    eps, max.steps) {
  coefk <- matrix(0, ncol(KernelX), class_num)
  xn <- nrow(KernelX)
  xp <- ncol(KernelX)
  for (i in 1:class_num) {
    class_idx <- which(y == class_set[i])
    Hk <- KernelX[-class_idx, ]
    Gk <- KernelX[class_idx, ]
    nGk <- length(class_idx)
    dim(Hk) <- c(xn - nGk, xp)
    dim(Gk) <- c(nGk, xp)
    HTH_inv_Gt <- chol2inv(chol(t(Hk) %*% Hk + diag(1e-7, ncol(Hk)))) %*% t(Gk)
    G_HTH_inv_Gt <- Gk %*% HTH_inv_Gt
    lb <- matrix(0, nGk)
    ub <- matrix(C[i], nGk)
    ek <- matrix(1, nGk)
    alphak <- clip_dcd_optimizer(G_HTH_inv_Gt, ek, lb, ub, eps, max.steps, lb)$x
    coefk[, i] <- HTH_inv_Gt %*% alphak
  }
  BaseDualHingeMBSVMClassifier <- list("coef" = coefk)
  class(BaseDualHingeMBSVMClassifier) <- "BaseDualHingeMBSVMClassifier"
  return(BaseDualHingeMBSVMClassifier)
}

#' Hinge Multiple Birth  Support Vector Machine
#'
#' \code{hinge_mbsvm} is an R implementation of Hinge-MBSVM
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
#' @param solver \code{"dual"} is available.
#' @param fit_intercept if set \code{fit_intercept = TRUE},
#'                      the function will evaluates intercept.
#' @param randx parameter for reduce SVM, default \code{randx = 0.1}.
#' @param ... unused parameters.
#' @return return \code{MBSVMClassifier} object.
#' @export
hinge_mbsvm <- function(X, y, C = 1,
                        kernel = c("linear", "rbf", "poly"),
                        gamma = 1 / ncol(X), degree = 3, coef0 = 0,
                        eps = 1e-5, max.steps = 4000,
                        solver = c("dual"), fit_intercept = TRUE,
                        randx = 1, ...) {
  X <- as.matrix(X)
  y <- as.matrix(y)
  class_set <- sort(unique(y))
  class_num <- length(class_set)
  if (length(C) == 1) {
    C <- matrix(C, class_num)
  }
  kernel <- match.arg(kernel)
  solver <- match.arg(solver)
  kso <- kernel_select_option(X, kernel, "primal", randx,
                              gamma, degree, coef0)
  KernelX <- kso$KernelX
  Kw <- kso$KernelX[kso$sample_idx, ]
  X <- kso$X
  if (fit_intercept == TRUE) {
    KernelX <- cbind(KernelX, 1)
  }
  if (solver == "dual") {
    solver.res <- hinge_mbsvm_dual_solver(KernelX, y, C,
                                          class_set, class_num, eps, max.steps)
  }
  MBSVMClassifier <- list("X" = X, "y" = y, "class_set" = class_set,
                          "class_num" = class_num,
                          "C" = C, "kernel" = kernel,
                          "gamma" = gamma, "degree" = degree, "coef0" = coef0,
                          "solver" = solver, "coef" = solver.res$coef,
                          "fit_intercept" = fit_intercept,
                          "Kw" = Kw)
  class(MBSVMClassifier) <- "MBSVMClassifier"
  return(MBSVMClassifier)
}

#' Predict Method for Multiple Birth Support Vector Machine
#'
#' @author Zhang Jiaqi
#' @param object a fitted object of class inheriting from \code{SVMClassifier}.
#' @param X new data for predicting.
#' @param values if set \code{values = TRUE}, this function will return predict
#'               values: f = <w,x>+b.
#' @param ... unused parameter.
#' @importFrom stats predict
#' @export
predict.MBSVMClassifier <- function(object, X, values = FALSE, ...) {
  X <- as.matrix(X)
  if (object$kernel == "linear") {
    KernelX <- X
  } else {
    KernelX <- kernel_function(X, object$X,
                               kernel.type = object$kernel,
                               gamma = object$gamma,
                               degree = object$degree,
                               coef0 = object$coef0,
                               symmetric = FALSE)
  }
  xp <- ncol(KernelX)
  if (object$fit_intercept == TRUE) {
    KernelX <- cbind(KernelX, 1)
  }
  fx <- KernelX %*% object$coef
  if (object$kernel == "linear") {
    coef_norm <- apply(object$coef[1:xp, ], 2, norm, type = "2")
  } else {
    num_class <- ncol(object$coef)
    coef_norm <- matrix(0, 1, num_class)
    vK <- t(object$coef[1:xp, ]) %*% object$Kw
    for (i in 1:num_class) {
      coef_norm[i] <-  vK[i, ] %*% object$coef[1:xp, i]
    }
    coef_norm <- as.numeric(coef_norm)
  }
  coef_norm[coef_norm == 0] <- 1e-7
  decf_values <- abs(fx) / coef_norm
  if (values == FALSE) {
    class_idx <- apply(decf_values, 1, which.max)
    decf <- object$class_set[class_idx]
    return(decf)
  } else if (values == TRUE) {
    return(decf_values)
  }
}

#' Plot Method for Support Vector Machine
#'
#' @author Zhang Jiaqi
#' @param x a fitted object of class inheriting from \code{SVMClassifier}.
#' @param ... unused parameter.
#' @importFrom stats coef
#' @importFrom graphics abline grid points
#' @export
plot.MBSVMClassifier <- function(x, ...) {
  if (ncol(x$X) == 2) {
    xlim_c <- c(min(x$X[,1]), max(x$X[, 1]))
    ylim_c <- c(min(x$X[,2]), max(x$X[, 2]))
    class_idx <- which(x$y == x$class_set[1])
    slopes <- - x$coef[1, ]/x$coef[2, ]
    intercepts <- - x$coef[3, ]/x$coef[2, ]
    plot(x$X[class_idx, 1], x$X[class_idx, 2], col = 1,
         xlim = xlim_c, ylim = ylim_c,
         xlab = "", ylab = "")
    grid(10, 10, lwd = 2,col = "grey")
    for (i in 2:x$class_num) {
      class_idx <- which(x$y == x$class_set[i])
      points(x$X[class_idx, 1], x$X[class_idx, 2], col = i)
    }
    for (i in 1:x$class_num) {
      abline(intercepts[i], slopes[i], col = i)
    }
  }
}
