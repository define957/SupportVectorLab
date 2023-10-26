hinge_nhsvm_svm_dual_solver <- function(KernelX, idx, y, C1, C2,
                                        eps, max.steps) {
  X_pos <- KernelX[idx, ]
  X_neg <- KernelX[-idx,]
  barX <- rbind(X_pos, -X_neg)
  np <- ncol(KernelX)
  nx <- nrow(KernelX)
  I <- diag(1, np)
  PosInv <- chol2inv(chol(I + C1*t(X_pos) %*% X_pos))
  NegInv <- chol2inv(chol(I + C1*t(X_neg) %*% X_neg))
  PosInv_barXT <- PosInv %*% t(barX)
  NegInv_barXT <- NegInv %*% t(barX)
  H <- barX %*% (PosInv_barXT + NegInv_barXT)
  e <- matrix(1, nx)
  lb <- matrix(0, nx)
  ub <- matrix(C2, nx)
  alpha <- clip_dcd_optimizer(H, e, lb, ub, eps, max.steps, lb)$x
  w_pos <-   PosInv_barXT %*% alpha
  w_neg <- - NegInv_barXT %*% alpha
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
#' @param batch_size mini-batch size for primal solver.
#' @param solver \code{"dual"} and \code{"primal"} are available.
#' @param fit_intercept if set \code{fit_intercept = TRUE},
#'                      the function will evaluates intercept.
#' @param randx parameter for reduce SVM, default \code{randx = 0.1}.
#' @param ... unused parameters.
#' @return return \code{NHSVMClassifier} object.
#' @export
hinge_nhsvm <- function(X, y, C1 = 1, C2 = 1,
                        kernel = c("linear", "rbf", "poly"),
                        gamma = 1 / ncol(X), degree = 3, coef0 = 0,
                        eps = 1e-5, max.steps = 80, batch_size = nrow(X) / 10,
                        solver = c("dual", "primal"),
                        fit_intercept = TRUE, randx = 1, ...) {
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
  kso <- kernel_select_option(X, kernel, solver = "primal", randx,
                              gamma, degree, coef0)
  KernelX <- kso$KernelX
  X <- kso$X
  if (fit_intercept == TRUE) {
    KernelX <- cbind(KernelX, 1)
  }
  if (solver == "dual") {
    solver.res <- hinge_nhsvm_svm_dual_solver(KernelX, idx, y, C1, C2,
                                              eps, max.steps)
  }
  NHSVMClassifier <- list("X" = X, "y" = y, "KernelX" = KernelX,
                          "class_set" = class_set,
                          "C1" = C1, "C2" = C2, "kernel" = kernel,
                          "gamma" = gamma, "degree" = degree, "coef0" = coef0,
                          "solver" = solver, "coef1" = solver.res$coef1,
                          "coef2" = solver.res$coef2,
                          "fit_intercept" = fit_intercept)
  class(NHSVMClassifier) <- "NHSVMClassifier"
  return(NHSVMClassifier)
}

#' Predict Method for Twin Support Vector Machine
#'
#' @author Zhang Jiaqi
#' @param object a fitted object of class inheriting from \code{SVMClassifier}.
#' @param X new data for predicting.
#' @param values if set \code{values = TRUE}, this function will return predict
#'               values: f = abs(wx+b).
#' @param ... unused parameter.
#' @importFrom stats predict
#' @export
predict.NHSVMClassifier <- function(object, X, values = FALSE, ...) {
  X <- as.matrix(X)
  if (object$fit_intercept == TRUE) {
    X <- cbind(X, 1)
  }
  if (object$kernel == "linear") {
    KernelX <- X
  } else {
    KernelX <- kernel_function(X, object$X,
                               kernel.type = object$kernel,
                               gamma = object$gamma,
                               degree = object$degree,
                               coef0 = object$coef0)
  }
  if (object$fit_intercept == TRUE) {
    KernelX <- cbind(KernelX, 1)
  }
  np <- ncol(KernelX) - 1
  if (object$kernel == "rbf") {
    w <- t(object$coef1[1:np,]) %*% object$KernelX %*% object$coef1[1:np,]
  }
  fx1 <- abs(KernelX %*% object$coef1)/norm(object$coef1, type = "2")
  fx2 <- abs(KernelX %*% object$coef2)/norm(object$coef1, type = "2")
  if (values == FALSE) {
    decf <- apply(cbind(fx1, fx2), 1, which.min)
    idx_pos <- which(decf == 1)
    idx_neg <- which(decf == 2)
    decf[idx_pos] <- object$class_set[1]
    decf[idx_neg] <- object$class_set[2]
  } else {
    dec_values1 <- fx1
    dec_values2 <- fx2
    return(cbind(dec_values1, dec_values2))
  }
  return(decf)
}
