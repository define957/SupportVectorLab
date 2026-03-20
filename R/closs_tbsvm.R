closs_tbsvm_dual_solver <- function(KernelX, idx, C1, C2, C3, C4,
                                    sigma1, sigma2,
                                    eps.hq, hq.steps) {

  update_omega <- function(f, sigma_sq) {
    return(as.vector(exp(-f^2/(2*sigma_sq))))
  }

  sigma1sq <- sigma1^2
  sigma2sq <- sigma2^2

  C1_ <- C1/(2*sigma1sq)
  C2_ <- C2/(2*sigma2sq)

  H  <- KernelX[-idx, , drop = FALSE]
  G  <- KernelX[idx, , drop = FALSE]
  Hn <- nrow(H)
  Gn <- nrow(G)
  xp    <- ncol(KernelX)

  Omegavec1 <- rep(1, Gn)
  Omegavec2 <- rep(1, Hn)

  e1 <- matrix(1, Hn)
  e2 <- matrix(1, Gn)

  coef1 <- matrix(0, xp)
  coef2 <- matrix(0, xp)

  HTH <- t(H) %*% H
  GTG <- t(G) %*% G

  for (i in 1:hq.steps) {
    G_omega_lam_T <-  t(G*(Omegavec1*C1_))
    A             <-  HTH + G_omega_lam_T %*% G
    diag(A)       <-  diag(A) + C3
    coef1_new     <- -cholsolve((A), G_omega_lam_T %*% e2)
    if (norm(coef1_new - coef1, type = "2") < eps.hq) {
      coef1 <- coef1_new
      break
    }
    f1        <-  G %*% coef1_new
    Omegavec1 <-  update_omega(1 + f1, sigma1sq)
    coef1     <- coef1_new
  }

  for (i in 1:hq.steps) {
    H_omega_lam_T <- t(H*(Omegavec2*C2_))
    B             <-  GTG + H_omega_lam_T %*% H
    diag(B)       <-  diag(B) + C4
    coef2_new     <- cholsolve((B), H_omega_lam_T %*% e1)
    if (norm(coef2_new - coef2, type = "2") < eps.hq) {
      coef2 <- coef2_new
      break
    }
    f2        <- H %*% coef2_new
    Omegavec2 <- update_omega(1 - f2, sigma2sq)
    coef2     <- coef2_new
  }

  BaseDualClossTBSVMClassifier <- list("coef1" = as.matrix(coef1),
                                       "coef2" = as.matrix(coef2))
  return(BaseDualClossTBSVMClassifier)
}

#' C-loss Least Squares Twin Bounded Support Vector Machine
#'
#' \code{closs_tbsvm} is an R implementation of CL2,p-LS-TSVM
#'
#' @author Zhang Jiaqi.
#' @param X,y dataset and label.
#' @param C1,C2,C3,C4 plenty term.
#' @param kernel kernel function. The definitions of various kernel functions are as follows:
#' \describe{
#'     \item{linear:}{\eqn{u'v}{u'*v}}
#'     \item{poly:}{\eqn{(\gamma u'v + coef0)^{degree}}{(gamma*u'*v + coef0)^degree}}
#'     \item{rbf:}{\eqn{e^{(-\gamma |u-v|^2)}}{exp(-gamma*|u-v|^2)}}
#' }
#' @param gamma parameter for \code{'rbf'} and \code{'poly'} kernel. Default \code{gamma = 1/ncol(X)}.
#' @param degree parameter for polynomial kernel, default: \code{degree = 3}.
#' @param coef0 parameter for polynomial kernel,  default: \code{coef0 = 0}.
#' @param sigma1,sigma2 parameters for C-loss function.
#' @param eps.hq the precision of the optimization algorithm.
#' @param hq.steps the number of iterations to solve the optimization problem.
#' @param fit_intercept if set \code{fit_intercept = TRUE},
#'                      the function will evaluates intercept.
#' @param reduce_set reduce set for reduce SVM, default \code{reduce_set = NULL}.
#' @return return \code{TSVMClassifier} object.
#' @export
closs_tbsvm <- function(X, y, C1 = 1, C2 = C1, C3 = 1e-7, C4 = C3,
                        kernel = c("linear", "rbf", "poly"),
                        gamma = 1 / ncol(X), degree = 3, coef0 = 0,
                        sigma1 = 1, sigma2 = sigma1,
                        eps.hq = 1e-2, hq.steps = 10,
                        fit_intercept = TRUE, reduce_set = NULL) {

  X         <- as.matrix(X)
  y         <- as.matrix(y)
  class_set <- sort(unique(y))
  idx       <- which(y == class_set[1])
  y[idx]    <- -1
  y[-idx]   <- 1
  y         <- as.matrix(as.numeric(y))

  if (length(class_set) != 2) {
    stop("This model only supports binary classification; y must have exactly two classes")
  }

  kernel  <- match.arg(kernel)
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
  solver.res <- closs_tbsvm_dual_solver(KernelX, idx, C1, C2, C3, C4,
                                        sigma1, sigma2,
                                        eps.hq, hq.steps)

  model_specs   <- list("X" = X, "y" = y,
                        "C1" = C1, "C2" = C2,
                        "fit_intercept" = fit_intercept,
                        "class_set" = class_set)
  model_coef    <- list("coef1" = solver.res$coef1,
                        "coef2" = solver.res$coef2)
  kernel_config <- list("kernel" = kernel,
                        "gamma"  = gamma,
                        "degree" = degree,
                        "coef0" = coef0,
                        "reduce_set" = reduce_set,
                        "KernelR" = KernelR,
                        "KernelX" = KernelX[, 1:kxp, drop = FALSE])
  w_norm        <- get_coef_norm(kernel_config, model_coef)
  kernel_config$w_norm <- w_norm

  TSVMClassifier <- structure(list("model_specs" = model_specs,
                                   "model_coef" = model_coef,
                                   "kernel_config" = kernel_config),
                              "class" = "TSVMClassifier")
  return(TSVMClassifier)
}
