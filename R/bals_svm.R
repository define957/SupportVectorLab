bals_svm_dual_solver <- function(KernelX, y, C, lambda, p,
                                 eps.hq, hq.steps,
                                 dual_optimizer, dual_optimizer_option) {
  v_update <- function(f, lambda, p) {
    lam_alsf <- rep(0, length(f))
    idx <- which(f >= 0)
    lam_alsf[idx]  <- (lambda*p)*f[idx]^2
    lam_alsf[-idx] <- (lambda*((1 - p)))*f[-idx]^2
    v <- -1 / (1 + lam_alsf)^2
    return(v)
  }
  n      <- nrow(KernelX)
  H0     <- SupportVectorLab:::calculate_svm_H(KernelX, y)
  H      <- cbind(H0, -H0)
  H      <- rbind(H, -H)
  H_R <- H
  diag_H <- diag(H)
  q      <- matrix(-1, 2*n)
  q[1:n] <- 1
  lb     <- matrix(0, 2*n)
  ub     <- matrix(Inf, 2*n)
  u0     <- lb
  v  <- matrix(-1, n)
  Omega <- as.numeric((C*lambda)*(-v))
  diagelem <- c(1/(Omega*(p)), 1/(Omega*(1 - p)))
  diag(H_R) <- diag_H + diagelem

  dual_optimizer_option <- append(list("H" = H_R,
                                       "q" = q,
                                       "lb" = lb,
                                       "ub" = ub,
                                       "u" = u0),
                                  dual_optimizer_option)

  for (i in 1:hq.steps) {
    u <- do.call("clip_dcd_optimizer", dual_optimizer_option)$x
    if (norm(u - u0, type = "2") < eps.hq) {
      break
    } else {
      u0 <- u
    }
    f <- 1 - H0 %*% (u0[1:n] - u0[(n + 1):(2 * n)])
    v <- v_update(f, lambda, p)
    Omega <- as.numeric((C*lambda)*(-v))
    diagelem <- c(1/(Omega*(p)), 1/(Omega*(1 - p)))
    diag(dual_optimizer_option$H) <- diag_H + diagelem
    dual_optimizer_option$u <- u
  }
  coef <- y*(u[1:n] - u[(n + 1):(2*n)])
  BaseDualBALSSVMClassifier <- list(coef = as.matrix(coef))
  class(BaseDualBALSSVMClassifier) <- "BaseDualBALSSVMClassifier"
  return(BaseDualBALSSVMClassifier)
}

#' Bounded Asymmetric Least Squares Loss Support Vector Machine
#'
#' \code{bals_svm} is an R implementation of BALS-SVM
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
#' @param p parameter for bounded asymmetric least squares loss.
#' @param lambda parameter for bq loss (loss increase speed).
#' @param gamma parameter for \code{'rbf'} and \code{'poly'} kernel. Default \code{gamma = 1/ncol(X)}.
#' @param degree parameter for polynomial kernel, default: \code{degree = 3}.
#' @param coef0 parameter for polynomial kernel,  default: \code{coef0 = 0}.
#' @param eps the precision of the optimization algorithm.
#' @param eps.hq the precision of the optimization algorithm.
#' @param max.steps the number of iterations to solve the optimization problem.
#' @param hq.steps the number of iterations of HQ
#' @param batch_size mini-batch size for primal solver.
#' @param fit_intercept if set \code{fit_intercept = TRUE},
#'                      the function will evaluates intercept.
#' @param ... unused parameters.
#' @return return \code{SVMClassifier} object.
#' @export
bals_svm <- function(X, y, C = 1, kernel = c("linear", "rbf", "poly"),
                     gamma = 1 / ncol(X), degree = 3, coef0 = 0,
                     p = 0.5, lambda = 1,
                     eps = 1e-5, eps.hq = 1e-2, max.steps = 4000, hq.steps = 10,
                     batch_size = nrow(X) / 10,
                     fit_intercept = TRUE,
                     dual_optimizer = clip_dcd_optimizer,
                     dual_optimizer_option = NULL) {

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
  solver <- "dual"
  reduce_set <- NULL
  if (fit_intercept == TRUE) {
    X <- cbind(X, 1)
  }
  reduce_flag <- is.null(reduce_set) == FALSE
  kso <- kernel_select_option_(X, kernel, reduce_set, gamma, degree, coef0)
  KernelX <- kso$KernelX
  if (solver == "dual") {
    if (is.null(dual_optimizer_option)) {
      dual_optimizer_option <- list("max.steps" = max.steps, "eps" = eps)
    }
    solver.res <- bals_svm_dual_solver(KernelX, y, C,
                                       lambda, p,
                                       eps.hq, hq.steps,
                                       dual_optimizer, dual_optimizer_option)
  }
  SVMClassifier <- list("X" = X, "y" = y,
                        "class_set" = class_set,
                        "C" = C, "kernel" = kernel,
                        "gamma" = gamma, "degree" = degree, "coef0" = coef0,
                        "solver" = solver, "coef" = solver.res$coef,
                        "fit_intercept" = fit_intercept,
                        "solver.res" = solver.res)
  class(SVMClassifier) <- "SVMClassifier"
  return(SVMClassifier)
}
