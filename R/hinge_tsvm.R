hinge_tsvm_dual_solver <- function(KernelX, idx, C1, C2, eps, max.steps) {
  H <- KernelX[-idx, ]
  G <- KernelX[idx, ]
  xn <- nrow(KernelX)
  xp <- ncol(KernelX)
  Hn <- nrow(H)
  Gn <- xn - Hn

  invHTH_GT <- cholsolve(t(H) %*% H + diag(1e-7, xp), t(G))
  dualH <- G %*% invHTH_GT
  dualq1 <- matrix(1, Gn)
  duallb1 <- matrix(0, Gn)
  dualub1 <- matrix(C1, Gn)
  u0 <- matrix(0, Gn)
  a <- clip_dcd_optimizer(dualH, dualq1, duallb1, dualub1,
                          eps, max.steps, u0)$x
  u <- -invHTH_GT %*% a

  invGTG_HT <- cholsolve(t(G) %*% G + diag(1e-7, xp), t(H))
  dualH <- H %*% invGTG_HT
  dualq2 <- matrix(1, Hn)
  duallb2 <- matrix(0, Hn)
  dualub2 <- matrix(C2, Hn)
  u0 <- matrix(0, Hn)
  g <- clip_dcd_optimizer(dualH, dualq2, duallb2, dualub2,
                          eps, max.steps, u0)$x
  v <- invGTG_HT %*% g
  BaseDualHingeTSVMClassifier <- list("coef1" = as.matrix(u),
                                      "coef2" = as.matrix(v))
  return(BaseDualHingeTSVMClassifier)
}

#' Hinge Twin Support Vector Machine
#'
#' \code{hinge_tsvm} is an R implementation of Hinge-TSVM
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
#' @return return \code{TSVMClassifier} object.
#' @export
hinge_tsvm <- function(X, y, C1 = 1, C2 = 1,
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
  solver.res <- hinge_tsvm_dual_solver(KernelX, idx, C1, C2, eps, max.steps)

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
  TSVMClassifier <- structure(
    list("model_specs" = model_specs,
         "model_coef" = model_coef,
         "kernel_config" = kernel_config),
    "class" = "TSVMClassifier")
  return(TSVMClassifier)
}

#' Predict Method for Twin Support Vector Machine
#'
#' @author Zhang Jiaqi
#' @param object a fitted object of class inheriting from \code{TSVMClassifier}.
#' @param X new data for predicting.
#' @param values if set \code{values = TRUE}, this function will return predict
#'               values: f = abs(wx+b).
#' @param ... unused parameter.
#' @importFrom stats predict
#' @export
predict.TSVMClassifier <- function(object, X, values = FALSE, ...) {
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

#' Plot Method for Support Vector Machine
#'
#' @author Zhang Jiaqi
#' @param x a fitted object of class inheriting from \code{SVMClassifier}.
#' @param ... unused parameter.
#' @importFrom stats coef
#' @importFrom graphics abline grid points
#' @export
plot.TSVMClassifier <- function(x, ...) {
  model_specs <- x$model_specs
  model_coef <- x$model_coef
  kernel_config <- x$kernel_config
  y <- model_specs$y
  coef1 <- model_coef$coef1
  coef2 <- model_coef$coef2
  class_set <- model_specs$class_set
  idx <- which(y == class_set[1])
  y[idx] <- -1
  y[-idx] <- 1
  xlim_c <- c(min(model_specs$X[,1]), max(model_specs$X[, 1]))
  ylim_c <- c(min(model_specs$X[,2]), max(model_specs$X[, 2]))
  if (length(coef1) == 3 && length(coef2) == 3) {
    plot(model_specs$X[idx, 1], model_specs$X[idx, 2], col = "red",
         xlim = xlim_c, ylim = ylim_c,
         xlab = "", ylab = "")
    grid(10, 10, lwd = 2,col = "grey")
    points(model_specs$X[-idx, 1], model_specs$X[-idx, 2], col = "blue")
    if (kernel_config$kernel == "linear") {
      abline(a = -coef1[3]/coef1[2], b = -coef1[1]/coef1[2],
             lty = 1, col = "red")
      abline(a = -(coef1[3] + 1)/coef1[2], b = -coef1[1]/coef1[2],
             lty = 2, col = "red")
      abline(a = -coef2[3]/coef2[2], b = -coef2[1]/coef2[2],
             lty = 1, col = "blue")
      abline(a = -(coef2[3] - 1)/coef2[2], b = -coef2[1]/coef2[2],
             lty = 2, col = "blue")
    }
  }
}
