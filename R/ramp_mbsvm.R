ramp_mbsvm_dual_solver <- function(KernelX, y, C, s, class_set, class_num,
                                   eps, eps.cccp, max.steps, cccp.steps) {
  update_deltak <- function(f, C, s) {
    deltak          <- matrix(0, length(f))
    deltak[f >= s]  <- C
    return(deltak)
  }
  xn <- nrow(KernelX)
  xp <- ncol(KernelX)
  coefk <- matrix(0, xp, class_num)

  for (i in 1:class_num) {
    class_idx <- which(y == class_set[i])
    Hk <- KernelX[-class_idx, ]
    Gk <- KernelX[class_idx, ]
    nGk <- length(class_idx)
    dim(Hk) <- c(xn - nGk, xp)
    dim(Gk) <- c(nGk, xp)
    HTH_inv_Gt <- chol2inv(chol(t(Hk) %*% Hk + diag(1e-7, ncol(Hk)))) %*% t(Gk)
    G_HTH_inv_Gt <- Gk %*% HTH_inv_Gt
    ek <- matrix(1, nGk)
    u0 <- matrix(0, nGk)
    deltak <- matrix(0, nGk)
    lb <- matrix(0, nGk)
    ub <- matrix(C[i], nGk)
    thetak <- matrix(0, xp)
    for (j in 1:cccp.steps) {
      u0[u0 > ub] <- ub[u0 > ub]
      u0[u0 < lb] <- lb[u0 < lb]
      q <- G_HTH_inv_Gt %*% deltak + ek
      u <- clip_dcd_optimizer(G_HTH_inv_Gt, q, lb, ub, eps, max.steps, u0)$x
      thetak <- HTH_inv_Gt %*% (u - deltak)
      if (norm(u - u0, type = "2") < eps.cccp) {
        break
      } else {
        u0 <- u
      }
      f <- 1 - Gk %*% thetak
      deltak <- update_deltak(f, C[i], s)
    }
    coefk[, i] <- thetak
  }
  BaseDualRampMBSVMClassifier <- list("coef" = coefk)
  class(BaseDualRampMBSVMClassifier) <- "BaseDualRampMBSVMClassifier"
  return(BaseDualRampMBSVMClassifier)
}

#' Ramp Multiple Birth Support Vector Machine
#'
#' \code{ramp_mbsvm} is an R implementation of Ramp-MBSVM
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
#' @param s parameter for ramp loss.
#' @param gamma parameter for \code{'rbf'} and \code{'poly'} kernel. Default \code{gamma = 1/ncol(X)}.
#' @param degree parameter for polynomial kernel, default: \code{degree = 3}.
#' @param coef0 parameter for polynomial kernel,  default: \code{coef0 = 0}.
#' @param eps the precision of the optimization algorithm.
#' @param eps.cccp the precision of the optimization algorithm.
#' @param max.steps the number of iterations to solve the optimization problem.
#' @param cccp.steps the number of iterations of CCCP.
#' @param solver \code{"dual"} is available.
#' @param fit_intercept if set \code{fit_intercept = TRUE},
#'                      the function will evaluates intercept.
#' @param randx parameter for reduce SVM, default \code{randx = 0.1}.
#' @param ... unused parameters.
#' @return return \code{MBSVMClassifier} object.
#' @export
ramp_mbsvm <- function(X, y, C = 1,
                       kernel = c("linear", "rbf", "poly"),
                       s = 1,
                       gamma = 1 / ncol(X), degree = 3, coef0 = 0,
                       eps = 1e-5, eps.cccp = 1e-2,
                       max.steps = 4000, cccp.steps = 10,
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
    solver.res <- ramp_mbsvm_dual_solver(KernelX, y, C, s,
                                         class_set, class_num, eps, eps.cccp,
                                         max.steps, cccp.steps)
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
