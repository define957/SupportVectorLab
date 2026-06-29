hinge_svm_dual_solver <- function(KernelX, y, C,
                                  dual_optimizer, dual_optimizer_option){
  n <- nrow(KernelX)
  H <- calculate_svm_H(KernelX, y)
  e <- matrix(1, nrow = n)
  lb <- matrix(0, nrow = n)
  ub <- matrix(C, nrow = n)
  u0 <- lb
  dual_optimizer_option <- append(list("H" = H,
                                       "q" = e,
                                       "lb" = lb,
                                       "ub" = ub,
                                       "u" = u0),
                                  dual_optimizer_option)
  solver.info <- do.call("dual_optimizer", dual_optimizer_option)
  alphas <- solver.info$x
  coef <- y*alphas
  BaseDualHingeSVMClassifier <- list(coef = as.matrix(coef),
                                     "solver.info" = solver.info)
  class(BaseDualHingeSVMClassifier) <- "BaseDualHingeSVMClassifier"
  return(BaseDualHingeSVMClassifier)
}


hinge_svm_primal_solver <- function(KernelX, X, y, C,
                                    max.steps, batch_size,
                                    optimizer, kernel, reduce_flag,
                                    reduce_set,
                                    ...) {
  sgHinge <- function(batch_KernelX, y, w, pars, ...) { # sub-gradient of hinge loss function
    C <- pars$C
    xn <- pars$xn
    KernelX <- pars$KernelX
    xmn <- nrow(batch_KernelX)
    xmp <- ncol(batch_KernelX)
    sg <- matrix(0, nrow = xmp, ncol = 1)
    u <- 1 - y * (batch_KernelX %*% w)
    u[u <= 0] <- 0
    u[u > 0] <- 1
    if (pars$kernel == "linear" || pars$reduce_flag) {
      sg <- w / xn - C * t(batch_KernelX) %*% (u*y) / xmn
    } else if (pars$kernel != "linear") {
      sg <- KernelX %*% w / xn - C * t(batch_KernelX) %*% (u*y) / xmn
    }
    return(sg)
  }
  xn <- nrow(X)
  if (kernel == "linear") { xp <- ncol(X) } else { xp <- ncol(KernelX)}
  w0 <- matrix(0, nrow = xp, ncol = 1)
  pars <- list("C" = C, "xn" = xn, "kernel" = kernel, "KernelX" = KernelX,
               "reduce_flag" = reduce_flag)
  if (kernel == "linear") {
    wt <- optimizer(X, y, w0, batch_size, max.steps, sgHinge, pars, ...)
  } else {
    wt <- optimizer(KernelX, y, w0, batch_size, max.steps, sgHinge, pars, ...)
  }
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
#' @param reduce_set reduce set for reduce SVM, default \code{reduce_set = NULL}.
#' @param dual_optimizer default optimizer is \code{clip_dcd_optimizer}.
#' @param dual_optimizer_option optimizer options.
#' @param ... unused parameters.
#' @return return \code{SVMClassifier} object.
#' @export
hinge_svm <- function(X, y, C = 1, kernel = c("linear", "rbf", "poly", "precomputed"),
                      gamma = 1 / ncol(X), degree = 3, coef0 = 0,
                      eps = 1e-5, max.steps = 4000, batch_size = nrow(X) / 10,
                      solver = c("dual", "primal"),
                      fit_intercept = TRUE, optimizer = pegasos,
                      reduce_set = NULL, dual_optimizer = clip_dcd_optimizer,
                      dual_optimizer_option = NULL, ...) {

  X <- as.matrix(X)
  y <- as.matrix(y)

  kernel <- match.arg(kernel)
  solver <- match.arg(solver)

  labels <- encode_binary_labels(y)
  y      <- labels$y
  class_set <- labels$class_set

  reduce_res   <- resolve_reduce_set(reduce_set, solver)
  reduce_set   <- reduce_res$reduce_set
  reduce_flag  <- reduce_res$reduce_flag

  if (kernel == "precomputed" && isTRUE(fit_intercept)) {
    # message("fit_intercept is ignored with precomputed kernel.")
    fit_intercept <- FALSE
  }

  X <- handle_intercept(X, fit_intercept)

  KernelX <- resolve_kernel_matrix(X, kernel, reduce_set,
                                   gamma, degree, coef0)

  batch_size <- resolve_batch_size(batch_size, nrow(X), solver)

  if (solver == "primal") {
    solver.res <- hinge_svm_primal_solver(KernelX, X, y, C,
                                          max.steps, batch_size,
                                          optimizer, kernel, reduce_flag,
                                          reduce_set,
                                          ...)
  } else if (solver == "dual") {
    if (is.null(dual_optimizer_option)) {
      dual_optimizer_option <- list("max.steps" = max.steps, "eps" = eps)
    }
    solver.res <- hinge_svm_dual_solver(KernelX, y, C,
                                        dual_optimizer, dual_optimizer_option)
  }

  if (kernel == "precomputed") {
    X <- NULL
  }

  SVMClassifier <- list("X" = X, "y" = y,
                        "reduce_flag" = reduce_flag,
                        "reduce_set" = reduce_set,
                        "class_set" = class_set,
                        "C" = C, "kernel" = kernel,
                        "gamma" = gamma, "degree" = degree, "coef0" = coef0,
                        "solver" = solver, "coef" = solver.res$coef,
                        "fit_intercept" = fit_intercept,
                        "solver.res" = solver.res)
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
    grid(lwd = 2,col = "grey")
    points(x$X[-idx, 1], x$X[-idx, 2], col = "blue")
    if (x$kernel == "linear") {
      abline(a = -coefs[3]/coefs[2], b = -coefs[1]/coefs[2],
             lty = 1)
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
    if (object$solver == "primal") {
      return(object$coef)
    }
    # dual solver
    if (object$kernel == "linear") {
      return(crossprod(object$X, object$coef))
    }
    return(object$coef)
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
  # print(coef(object))
  X <- as.matrix(X)
  if (object$fit_intercept == TRUE) {
    X <- cbind(X, 1)
  }
  if (object$kernel == "precomputed") {
    KernelX <- X
  } else if (object$kernel == "linear" && object$solver == "primal") {
    KernelX <- X
  } else {
    if (is.null(object$reduce_flag) == FALSE && object$reduce_flag == TRUE) {
      mapX <- object$reduce_set
    } else {
      mapX <- object$X
    }
    KernelX <- kernel_function(X, mapX,
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
  H <- KernelX
  idx_pos <- which(y == 1L)
  H[idx_pos, -idx_pos] <- -H[idx_pos, -idx_pos]
  H[-idx_pos, idx_pos] <- -H[-idx_pos, idx_pos]
  return(H)
}
