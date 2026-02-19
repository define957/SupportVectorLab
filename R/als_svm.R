als_svm_dual_solver <- function(KernelX, y, C, p,
                                dual_optimizer, dual_optimizer_option) {
  n <- nrow(KernelX)
  H <- calculate_svm_H(KernelX, y)
  H <- cbind(H, -H)
  H <- rbind(H, -H)
  diagelem <- rep(c(1/(C*(p)), 1/(C*(1 - p))), rep(n, 2))
  diag(H) <- diag(H) + diagelem
  q <- matrix(-1, nrow = 2*n)
  q[1:n] <- 1
  lb <- matrix(0, nrow = 2*n)
  ub <- matrix(Inf, nrow = 2*n)
  u0 <- lb
  dual_optimizer_option <- append(list("H" = H,
                                       "q" = q,
                                       "lb" = lb,
                                       "ub" = ub,
                                       "u" = u0),
                                  dual_optimizer_option)
  solver.info <- do.call("dual_optimizer", dual_optimizer_option)
  u <- solver.info$x
  coef <- y*(u[1:n] - u[(n + 1):(2*n)])
  BaseDualALSSVMClassifier <- list(coef = as.matrix(coef),
                                   "solver.info" = solver.info)
  class(BaseDualALSSVMClassifier) <- "BaseDualALSSVMClassifier"
  return(BaseDualALSSVMClassifier)
}


als_svm_primal_solver <- function(KernelX, y, C = 1, p = 0.5,
                                  max.steps = 80, batch_size = nrow(KernelX) / 10,
                                  optimizer = pegasos, ...) {
  gLeastSquares <- function(KernelX, y, w, pars,...) { # gradient of Least Squares loss function
    C <- pars$C
    xn <- pars$xn
    xmn <- nrow(KernelX)
    xmp <- ncol(KernelX)
    sg <- matrix(0, xmp, 1)
    u <- 1 - y * (KernelX %*% w)
    idx <- which(u >= 0)
    u[idx] <- p
    u[-idx] <- (1 - p)
    sg <- w - (C*xn/xmn) * t(KernelX) %*% (u*y)
    return(sg)
  }
  xn <- nrow(KernelX)
  xp <- ncol(KernelX)
  w0 <- matrix(0, nrow = xp, ncol = 1)
  pars <- list("C" = C, "xn" = xn)
  wt <- optimizer(KernelX, y, w0, batch_size, max.steps, gLeastSquares, pars, ...)
  BasePrimalALSSVMClassifier <- list(coef = as.matrix(wt[1:xp]))
  class(BasePrimalALSSVMClassifier) <- "BasePrimalALSSVMClassifier"
  return(BasePrimalALSSVMClassifier)
}


#' Asymmetric Least Squares Support Vector Machine
#'
#' \code{als_svm} is an R implementation of ALS-SVM
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
#' @param p parameter for asymmetric least squares loss.
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
als_svm <- function(X, y, C = 1, kernel = c("linear", "rbf", "poly"),
                    gamma = 1 / ncol(X), degree = 3, coef0 = 0,
                    p = 0.5,
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
    solver.res <- als_svm_primal_solver(KernelX, y, C, p,
                                        max.steps, batch_size,
                                        optimizer, ...)
  } else if (solver == "dual") {
    if (is.null(dual_optimizer_option)) {
      dual_optimizer_option <- list(max.steps = max.steps,
                                    eps = eps)
    }
    solver.res <- als_svm_dual_solver(KernelX, y, C, p,
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
