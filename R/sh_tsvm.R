sh_tsvm_dual_solver <- function(KernelX, idx, C1, C2, eps, max.steps) {
  H <- KernelX[-idx, ]
  G <- KernelX[idx, ]
  xn <- nrow(KernelX)
  xp <- ncol(KernelX)
  Hn <- nrow(H)
  Gn <- xn - Hn

  invHTH_GT <- cholsolve(t(H) %*% H + diag(1e-7, xp), t(G))
  dualH <- G %*% invHTH_GT + diag(1/C1, Gn)
  dualq1 <- matrix(1, Gn)
  duallb1 <- matrix(0, Gn)
  dualub1 <- matrix(Inf, Gn)
  u0 <- matrix(0, Gn)
  a <- clip_dcd_optimizer(dualH, dualq1, duallb1, dualub1,
                          eps, max.steps, u0)$x
  u <- -invHTH_GT %*% a

  invGTG_HT <- cholsolve(t(G) %*% G + diag(1e-7, xp), t(H))
  dualH <- H %*% invGTG_HT + diag(1/C2, Hn)
  dualq2 <- matrix(1, Hn)
  duallb2 <- matrix(0, Hn)
  dualub2 <- matrix(Inf, Hn)
  u0 <- matrix(0, Hn)
  g <- clip_dcd_optimizer(dualH, dualq2, duallb2, dualub2,
                          eps, max.steps, u0)$x
  v <- invGTG_HT %*% g
  BaseDualSquaredHingeTSVMClassifier <- list("coef1" = as.matrix(u),
                                             "coef2" = as.matrix(v))
  return(BaseDualSquaredHingeTSVMClassifier)
}

#' Squared Hinge Twin Support Vector Machine
#'
#' \code{sh_tsvm} is an R implementation of SH-TSVM
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
sh_tsvm <- function(X, y, C1 = 1, C2 = 1,
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
  solver.res <- sh_tsvm_dual_solver(KernelX, idx, C1, C2, eps, max.steps)

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
