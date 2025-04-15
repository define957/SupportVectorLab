hinge_nhsvm_svm_dual_solver <- function(KernelX, idx, y, C1, C2,
                                        eps, max.steps) {
  X_pos <- KernelX[-idx, , drop = FALSE]
  X_neg <- KernelX[idx, , drop = FALSE]
  xp <- ncol(KernelX)
  xn <- nrow(KernelX)
  I <- diag(1, xp)
  DK <- diag(as.numeric(y)) %*% KernelX
  PosInv_DK <- cholsolve(I + C1*t(X_pos) %*% X_pos, t(DK))
  NegInv_DK <- cholsolve(I + C1*t(X_neg) %*% X_neg, t(DK))
  dualH <- DK %*% (PosInv_DK + NegInv_DK)
  e <- matrix(1, xn)
  lb <- matrix(0, xn)
  ub <- matrix(C2, xn)
  alpha <- clip_dcd_optimizer(dualH, e, lb, ub, eps, max.steps, lb)$x
  w_pos <-  PosInv_DK %*% alpha
  w_neg <- -NegInv_DK %*% alpha
  BaseDualHingeNHSVMClassifier <- list("coef1" = as.matrix(w_pos),
                                       "coef2" = as.matrix(w_neg))
  class(BaseDualHingeNHSVMClassifier) <- "BaseDualHingeNHSVMClassifier"
  return(BaseDualHingeNHSVMClassifier)
}

#' Hinge Nonparallel Hyperplane Support Vector Machine
#'
#' \code{hinge_nhsvm} is an R implementation of Hinge-NHSVM
#'
#' @author Zhang Jiaqi.
#' @param X,y dataset and label.
#' @param C1,C2 plenty term.
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
#' @return return \code{NHSVMClassifier} object.
#' @export
hinge_nhsvm <- function(X, y, C1 = 1, C2 = 1,
                        kernel = c("linear", "rbf", "poly"),
                        gamma = 1 / ncol(X), degree = 3, coef0 = 0,
                        eps = 1e-5, max.steps = 4000,
                        fit_intercept = TRUE, reduce_set = NULL) {
  X <- as.matrix(X)
  y <- as.matrix(y)

  class_set <- sort(unique(y))
  idx <- which(y == class_set[1])
  y[idx] <- -1
  y[-idx] <- 1
  y <- as.matrix(as.numeric(y))
  if (length(class_set) > 2) {
    stop("The number of class should less 2!")
  }

  kernel <- match.arg(kernel)
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

  solver.res <- hinge_nhsvm_svm_dual_solver(KernelX, idx, y, C1, C2,
                                            eps, max.steps)

  model_specs = list("X" = X, "y" = y,
                     "C1" = C1, "C2" = C2,
                     "fit_intercept" = fit_intercept,
                     "class_set" = class_set)
  model_coef = list("coef1" = solver.res$coef1,
                    "coef2" = solver.res$coef2)
  kernel_config = list("kernel" = kernel,
                       "gamma" = gamma,
                       "degree" = degree,
                       "coef0" = coef0,
                       "reduce_set" = reduce_set,
                       "KernelR" = KernelR,
                       "KernelX" = KernelX[, 1:kxp, drop = FALSE])
  w_norm <- get_coef_norm(kernel_config, model_coef)
  kernel_config$w_norm <- w_norm
  NHSVMClassifier <- structure(
    list("model_specs" = model_specs,
         "model_coef" = model_coef,
         "kernel_config" = kernel_config),
    "class" = "NHSVMClassifier")
  return(NHSVMClassifier)
}

#' Coef Method for Nonparallel Hyperplane Support Vector Machine
#'
#' @author Zhang Jiaqi
#' @param object a fitted object of class inheriting from \code{NHSVMClassifier}.
#' @param ... unused parameter.
#' @importFrom stats coef
#' @export
coef.NHSVMClassifier <- function(object, ...) {
  model_coef <- object$model_coef
  return(model_coef)
}

#' Predict Method for Nonparallel Hyperplane Support Vector Machine
#'
#' @author Zhang Jiaqi
#' @param object a fitted object of class inheriting from \code{NHSVMClassifier}.
#' @param X new data for predicting.
#' @param values if set \code{values = TRUE}, this function will return predict
#'               values: f = abs(wx+b).
#' @param ... unused parameter.
#' @importFrom stats predict
#' @export
predict.NHSVMClassifier <- function(object, X, values = FALSE, ...) {
  X <- as.matrix(X)
  model_specs <- object$model_specs
  kernel_config <- object$kernel_config
  model_coef <- object$model_coef
  reduce_flag <- ifelse(is.null(kernel_config$reduce_set), 0, 1)
  if (kernel_config$kernel == "linear") {
    KernelX <- X
  } else {
    if (reduce_flag) { Xmap <- kernel_config$reduce_set } else {Xmap <- model_specs$X}
    KernelX <- kernel_function(X, Xmap,
                               kernel.type = kernel_config$kernel,
                               gamma = kernel_config$gamma,
                               degree = kernel_config$degree,
                               coef0 = kernel_config$coef0)
  }
  if (model_specs$fit_intercept == TRUE) {
    KernelX <- cbind(KernelX, 1)
  }
  w_norm <- kernel_config$w_norm
  fx1 <- abs(KernelX %*% model_coef$coef1) / w_norm$w1_norm
  fx2 <- abs(KernelX %*% model_coef$coef2) / w_norm$w2_norm
  if (values == FALSE) {
    decf <- apply(cbind(fx1, fx2), 1, which.min)
    idx_pos <- which(decf == 1)
    idx_neg <- which(decf == 2)
    decf[idx_pos] <- model_specs$class_set[2]
    decf[idx_neg] <- model_specs$class_set[1]
  } else {
    dec_values1 <- fx1
    dec_values2 <- fx2
    return(cbind(dec_values1, dec_values2))
  }
  return(decf)
}

get_coef_norm <- function(kernel_config, model_coef) {
  reduce_flag <- ifelse(is.null(kernel_config$reduce_set), 0, 1)
  if (kernel_config$kernel != "linear") {
    if (reduce_flag) { K <- kernel_config$KernelR } else { K <- kernel_config$KernelX }
    xp <- ncol(K)
    w1_norm <- sqrt(t(model_coef$coef1[1:xp,]) %*% K %*% model_coef$coef1[1:xp,])
    w2_norm <- sqrt(t(model_coef$coef2[1:xp,]) %*% K %*% model_coef$coef2[1:xp,])
  } else {
    w1_norm <- norm(model_coef$coef1, type = "2")
    w2_norm <- norm(model_coef$coef2, type = "2")
  }
  return(list("w1_norm" = as.numeric(w1_norm),
              "w2_norm" = as.numeric(w2_norm)))
}
