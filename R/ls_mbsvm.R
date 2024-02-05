ls_mbsvm_dual_solver <- function(KernelX, y, C, class_set, class_num) {
  coefk <- matrix(0, ncol(KernelX), class_num)
  xn <- nrow(KernelX)
  xp <- ncol(KernelX)
  Ie <- diag(1e-7, xp)
  for (i in 1:class_num) {
    class_idx <- which(y == class_set[i])
    Hk <- KernelX[-class_idx, ]
    Gk <- KernelX[class_idx, ]
    nGk <- length(class_idx)
    dim(Hk) <- c(xn - nGk, xp)
    dim(Gk) <- c(nGk, xp)
    GramHK <- t(Hk) %*% Hk
    GramGK <- t(Gk) %*% Gk
    ek <- matrix(1, nGk)
    coefk[, i] <- C[i] * solve(GramHK + Ie + GramGK, t(Gk) %*% ek)
  }
  BaseDualLeastSquaresMBSVMClassifier <- list("coef" = coefk)
  class(BaseDualLeastSquaresMBSVMClassifier) <- "BaseDualLeastSquaresMBSVMClassifier"
  return(BaseDualLeastSquaresMBSVMClassifier)
}

#' Least Squares Multiple Birth Support Vector Machine
#'
#' \code{ls_mbsvm} is an R implementation of LS-MBSVM
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
#' @param solver \code{"dual"} is available.
#' @param fit_intercept if set \code{fit_intercept = TRUE},
#'                      the function will evaluates intercept.
#' @param randx parameter for reduce SVM, default \code{randx = 0.1}.
#' @param ... unused parameters.
#' @return return \code{MBSVMClassifier} object.
#' @export
ls_mbsvm <- function(X, y, C = 1,
                     kernel = c("linear", "rbf", "poly"),
                     gamma = 1 / ncol(X), degree = 3, coef0 = 0,
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
    solver.res <- ls_mbsvm_dual_solver(KernelX, y, C, class_set, class_num)
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
