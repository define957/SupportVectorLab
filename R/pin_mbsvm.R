pin_mbsvm_dual_solver <- function(KernelX, y, C, tau, class_set, class_num,
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
    lb <- matrix(-tau*C[i], nGk)
    ub <- matrix(C[i], nGk)
    ek <- matrix(1, nGk)
    if (tau == 0) {
      u0 <- lb
    } else {
      u0 <- (lb + ub)/2
    }
    alphak <- clip_dcd_optimizer(G_HTH_inv_Gt, ek, lb, ub, eps, max.steps, u0)$x
    coefk[, i] <- HTH_inv_Gt %*% alphak
  }
  BaseDualHingeMBSVMClassifier <- list("coef" = coefk)
  class(BaseDualHingeMBSVMClassifier) <- "BaseDualHingeMBSVMClassifier"
  return(BaseDualHingeMBSVMClassifier)
}

#' Pinball Loss Multiple Birth Support Vector Machine
#'
#' \code{pin_mbsvm} is an R implementation of Pin-MBSVM
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
#' @param tau parameter for pinball loss.
#' @param eps the precision of the optimization algorithm.
#' @param max.steps the number of iterations to solve the optimization problem.
#' @param solver \code{"dual"} is available.
#' @param fit_intercept if set \code{fit_intercept = TRUE},
#'                      the function will evaluates intercept.
#' @param randx parameter for reduce SVM, default \code{randx = 1}.
#' @param ... unused parameters.
#' @return return \code{MBSVMClassifier} object.
#' @export
pin_mbsvm <- function(X, y, C = 1,
                      kernel = c("linear", "rbf", "poly"),
                      gamma = 1 / ncol(X), degree = 3, coef0 = 0,
                      tau = 1,
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
    solver.res <- pin_mbsvm_dual_solver(KernelX, y, C, tau,
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
