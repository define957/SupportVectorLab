ls_tsvm_dual_solver <- function(KernelX, idx, C1, C2) {
  EMat <- KernelX[idx, ]
  FMat <- KernelX[-idx, ]
  xn <- nrow(KernelX)
  xp <- ncol(KernelX)
  Fn <- nrow(FMat)
  Mn <- xn - Fn
  GramF <- t(FMat) %*% FMat
  GramE <- t(EMat) %*% EMat
  e1 <- matrix(1, Fn)
  e2 <- matrix(1, Mn)
  u <- -cholsolve(GramF + GramE/C1 + diag(1e-7, xp), t(FMat) %*% e1)
  v <-  cholsolve(GramE + GramF/C2 + diag(1e-7, xp), t(EMat) %*% e2)
  BaseDualHingeTSVMClassifier <- list("coef1" = as.matrix(u),
                                      "coef2" = as.matrix(v))
  return(BaseDualHingeTSVMClassifier)
}

#' Least Squares Twin Support Vector Machine
#'
#' \code{ls_tsvm} is an R implementation of Hinge-TSVM
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
#' @param solver \code{"dual"} is available.
#' @param fit_intercept if set \code{fit_intercept = TRUE},
#'                      the function will evaluates intercept.
#' @param randx parameter for reduce SVM, default \code{randx = 0.1}.
#' @param ... unused parameters.
#' @return return \code{TSVMClassifier} object.
#' @export
ls_tsvm <- function(X, y, C1 = 1, C2 = C1,
                    kernel = c("linear", "rbf", "poly"),
                    gamma = 1 / ncol(X), degree = 3, coef0 = 0,
                    solver = c("dual"), fit_intercept = TRUE,
                    randx = 0.1, ...) {
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
    solver.res <- ls_tsvm_dual_solver(KernelX, idx, C1, C2)
  }
  TSVMClassifier <- list("X" = X, "y" = y, "class_set" = class_set,
                         "C1" = C1, "C2" = C2, "kernel" = kernel,
                         "gamma" = gamma, "degree" = degree, "coef0" = coef0,
                         "solver" = solver, "coef1" = solver.res$coef1,
                         "coef2" = solver.res$coef2,
                         "fit_intercept" = fit_intercept,
                         "Kw" = Kw)
  class(TSVMClassifier) <- "TSVMClassifier"
  return(TSVMClassifier)
}
