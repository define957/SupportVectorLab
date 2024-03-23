rq_svm_dual_solver <- function(KernelX, y, C = 1, update_deltak,
                               tau = 0.5, lambda = 1,
                               eps = 1e-5, eps.cccp = 1e-2, max.steps = 80, cccp.steps = 10) {
  eta <- 1/(1 - exp(-1/lambda))
  n <- nrow(KernelX)
  H <- calculate_svm_H(KernelX, y)
  u0 <- matrix(0, nrow = n, ncol = 1)
  e <- matrix(1, nrow = n, ncol = 1)
  if (tau == 0) {
    u0 <- matrix(0, nrow = n, ncol = 1)
  } else {
    u0 <- matrix((1 - tau)*lambda*C/2, nrow = n, ncol = 1)
  }
  for (i in 1:cccp.steps) {
    f <- 1 - H %*% u0
    delta_k <- update_deltak(f, u0, tau, lambda)
    lb <- -C*delta_k - C*eta*tau/lambda
    ub <- -C*delta_k + C*eta/lambda
    u0[u0 > ub] <- ub[u0 > ub]
    u0[u0 < lb] <- lb[u0 < lb]
    u <- clip_dcd_optimizer(H, e, lb, ub, eps, max.steps, u0)$x
    if (norm(u - u0, type = "2") < eps.cccp) {
      break
    } else {
      u0 <- u
    }
  }
  coef <- y*u
  BaseDualRQSVMClassifier <- list(coef = as.matrix(coef))
  class(BaseDualRQSVMClassifier) <- "BaseDualRQSVMClassifier"
  return(BaseDualRQSVMClassifier)
}


rq_svm_primal_solver <- function(KernelX, y, C = 1, update_deltak,
                                 tau = 0, lambda = 1, eps.cccp = 1e-2,
                                 max.steps = 80, cccp.steps = 10, batch_size = nrow(KernelX) / 10,
                                 optimizer = pegasos, ...) {
  sgRq <- function(KernelX, y, w, pars, At, ...) { # sub-gradient of RQ loss function
    lambda <- pars$lambda
    tau <- pars$tau
    C <- pars$C
    xn <- pars$xn
    deltak <- pars$deltak
    eta <- 1/(1 - exp(-1/lambda))
    xmn <- nrow(KernelX)
    xmp <- ncol(KernelX)
    sg <- matrix(0, xmp, 1)
    f <- 1 - y*(KernelX %*% w)
    u <- matrix(0, xmn)
    u[f < 0] <- -tau
    u[f >= 0] <- 1
    sg <- w - eta*(C*xn/xmn) * t(KernelX) %*% (u*y)/lambda +
      (C*xn/xmn)*t(KernelX) %*% (y*deltak[At, ])
    return(sg)
  }
  xn <- nrow(KernelX)
  xp <- ncol(KernelX)
  w0 <- matrix(0, xp, 1)
  pars <- list("C" = C, "lambda" = lambda, "tau" = tau, "xn" = xn)
  for (i in 1:cccp.steps) {
    f <- 1 - y*(KernelX %*% w0)
    deltak <- update_deltak(f, w0, tau, lambda)
    pars$deltak <- deltak
    wt <- optimizer(KernelX, y, w0, batch_size, max.steps, sgRq, pars, ...)
    if (norm(wt - w0, type = "2") < eps.cccp) {
      break
    } else {
      w0 <- wt
    }
  }
  BasePrimalRqSVMClassifier <- list(coef = as.matrix(wt[1:xp]))
  class(BasePrimalRqSVMClassifier) <- "BasePrimalRqSVMClassifier"
  return(BasePrimalRqSVMClassifier)
}


#' Rescaled Quantile Loss Support Vector Machine
#'
#' \code{rq_svm} is an R implementation of RQ-SVM
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
#' @param tau parameter for RQ loss (quantile parameter).
#' @param lambda parameter for RQ loss (loss increase speed).
#' @param gamma parameter for \code{'rbf'} and \code{'poly'} kernel. Default \code{gamma = 1/ncol(X)}.
#' @param degree parameter for polynomial kernel, default: \code{degree = 3}.
#' @param coef0 parameter for polynomial kernel,  default: \code{coef0 = 0}.
#' @param eps the precision of the optimization algorithm.
#' @param eps.cccp the precision of the optimization algorithm.
#' @param max.steps the number of iterations to solve the optimization problem.
#' @param cccp.steps the number of iterations of CCCP.
#' @param batch_size mini-batch size for primal solver.
#' @param solver \code{"dual"} and \code{"primal"} are available.
#' @param fit_intercept if set \code{fit_intercept = TRUE},
#'                      the function will evaluates intercept.
#' @param optimizer default primal optimizer pegasos.
#' @param randx parameter for reduce SVM, default \code{randx = 0.1}.
#' @param ... unused parameters.
#' @return return \code{SVMClassifier} object.
#' @export
rq_svm <- function(X, y, C = 1, kernel = c("linear", "rbf", "poly"),
                   gamma = 1 / ncol(X), degree = 3, coef0 = 0,
                   tau = 0.5, lambda = 1,
                   eps = 1e-5, eps.cccp = 1e-2, max.steps = 80, cccp.steps = 10,
                   batch_size = nrow(X) / 10,
                   solver = c("dual", "primal"),
                   fit_intercept = TRUE, optimizer = pegasos, randx = 0.1, ...) {
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
  kso <- kernel_select_option(X, kernel, solver, randx,
                              gamma, degree, coef0)
  KernelX <- kso$KernelX
  X <- kso$X
  update_deltak <- function(f, u, tau, lambda) {
    eta <- 1/(1 - exp(-1/lambda))
    idx <- which(f >= 0)
    delta_k <- matrix(0, nrow = nrow(KernelX))
    delta_k[idx] <- eta*(1 - exp(-f[idx]/lambda))/lambda
    delta_k[-idx] <- -tau*eta*(1 - exp(tau*f[-idx]/lambda))/lambda
    return(delta_k)
  }
  if (solver == "primal") {
    solver.res <- rq_svm_primal_solver(KernelX, y, C, update_deltak,
                                       tau, lambda, eps.cccp, max.steps, cccp.steps, batch_size,
                                       optimizer, ...)
  } else if (solver == "dual") {
    solver.res <- rq_svm_dual_solver(KernelX, y, C, update_deltak,
                                     tau, lambda, eps, eps.cccp,
                                     max.steps, cccp.steps)
  }
  SVMClassifier <- list("X" = X, "y" = y, "class_set" = class_set,
                        "C" = C, "kernel" = kernel,
                        "gamma" = gamma, "degree" = degree, "coef0" = coef0,
                        "solver" = solver, "coef" = solver.res$coef,
                        "fit_intercept" = fit_intercept)
  class(SVMClassifier) <- "SVMClassifier"
  return(SVMClassifier)
}
