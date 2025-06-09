pin_svm_dual_solver <- function(KernelX, y, C, tau,
                                dual_optimizer, dual_optimizer_option) {
  n <- nrow(KernelX)
  H <- calculate_svm_H(KernelX, y)
  e <- matrix(1, nrow = n)
  lb <- matrix(-tau*C, nrow = n)
  ub <- matrix(C, nrow = n)

  if (tau == 0) {
    u0 <- matrix(0, n)
  } else {
    u0 <- (lb + ub)/2
  }
  dual_optimizer_option <- append(list("H" = H,
                                       "q" = e,
                                       "lb" = lb,
                                       "ub" = ub,
                                       "u" = u0),
                                  dual_optimizer_option)
  u <- do.call("dual_optimizer", dual_optimizer_option)$x
  coef <- y*u
  BaseDualPinSVMClassifier <- list(coef = as.matrix(coef))
  class(BaseDualPinSVMClassifier) <- "BaseDualPinSVMClassifier"
  return(BaseDualPinSVMClassifier)
}


pin_svm_primal_solver <- function(KernelX, X, y, C, tau,
                                  max.steps, batch_size,
                                  optimizer, kernel, reduce_flag,
                                  reduce_set, ...) {
  sgpinball <- function(batch_KernelX, y, w, pars, ...) { # sub-gradient of Pinball loss function
    C <- pars$C
    xn <- pars$xn
    KernelX <- pars$KernelX
    xmn <- nrow(batch_KernelX)
    xmp <- ncol(batch_KernelX)
    sg <- matrix(0, nrow = xmp, ncol = 1)
    u <- 1 - y * (batch_KernelX %*% w)
    u[u <= 0] <- -tau
    u[u > 0] <- 1
    if (pars$kernel == "linear" || pars$reduce_flag) {
      sg <- w / xn - C * t(batch_KernelX) %*% (u*y) / xmn
    } else if (pars$kernel != "linear") {
      sg <- KernelX %*% w / xn - C * t(batch_KernelX) %*% (u*y) / xmn
    }
  }
  xn <- nrow(X)
  if (kernel == "linear") { xp <- ncol(X) } else { xp <- ncol(KernelX)}
  w0 <- matrix(0, nrow = xp, ncol = 1)
  pars <- list("C" = C, "tau" = tau, "xn" = xn, "kernel" = kernel,
               "KernelX" = KernelX, "reduce_flag" = reduce_flag)
  if (kernel == "linear") {
    wt <- optimizer(X, y, w0, batch_size, max.steps, sgpinball, pars, ...)
  } else {
    wt <- optimizer(KernelX, y, w0, batch_size, max.steps, sgpinball, pars, ...)
  }
  BasePrimalPinSVMClassifier <- list(coef = as.matrix(wt[1:xp]))
  class(BasePrimalPinSVMClassifier) <- "BasePrimalPinSVMClassifier"
  return(BasePrimalPinSVMClassifier)
}


#' Pin Support Vector Machine
#'
#' \code{pin_svm} is an R implementation of Pin-SVM
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
#' @param tau parameter for pinball loss
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
pin_svm <- function(X, y, C = 1, kernel = c("linear", "rbf", "poly"),
                    tau = 0.5,
                    gamma = 1 / ncol(X), degree = 3, coef0 = 0,
                    eps = 1e-5, max.steps = 4000, batch_size = nrow(X) / 10,
                    solver = c("dual", "primal"),
                    fit_intercept = TRUE, optimizer = pegasos,
                    reduce_set = NULL, dual_optimizer = clip_dcd_optimizer,
                    dual_optimizer_option = NULL, ...) {
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
  reduce_flag <- is.null(reduce_set) == FALSE
  if (solver == "dual" && reduce_flag == TRUE) {
    reduce_flag <- FALSE
    reduce_set <- NULL
    cat("The dual solver does not support the reduce set; it has been set to NULL.\n")
  }
  kso <- kernel_select_option_(X, kernel, reduce_set, gamma, degree, coef0)
  KernelX <- kso$KernelX
  if (solver == "primal") {
    solver.res <- pin_svm_primal_solver(KernelX, X, y, C, tau,
                                        max.steps, batch_size,
                                        optimizer, kernel, reduce_flag,
                                        reduce_set,
                                        ...)
  } else if (solver == "dual") {
    if (is.null(dual_optimizer_option)) {
      dual_optimizer_option <- list("max.steps" = max.steps, "eps" = eps)
    }
    solver.res <- pin_svm_dual_solver(KernelX, y, C, tau,
                                      dual_optimizer, dual_optimizer_option)
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
