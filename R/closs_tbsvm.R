closs_tbsvm_dual_solver <- function(KernelX, idx, C1, C2, C3, C4,
                                    sigma1, sigma2,
                                    eps.hq, hq.steps, method_ls) {
  c_loss_function <- function(z, sigmasq, beta) {
    return(sum(-beta*expm1(-z^2 / (2*sigmasq))))
  }
  ls_loss_function <- function(z) {
    return(sum(z^2))
  }

  update_omega <- function(f, sigma_sq) {
    return(as.vector(exp(-f^2/(2*sigma_sq))))
  }

  compute_objective1 <- function(coef, f, idx, C_reg, C_between, sigmasq, beta) {
    obj <- 0.5*(C_reg * crossprod(coef) + ls_loss_function(f[idx]) +
      C_between * c_loss_function(1 + f[-idx], sigmasq, beta))
    return(obj)
  }

  compute_objective2 <- function(coef, f, idx, C_reg, C_between, sigmasq, beta) {
    obj <- 0.5*(C_reg * crossprod(coef) + ls_loss_function(f[-idx]) +
      C_between * c_loss_function(1 - f[idx], sigmasq, beta))
    return(obj)
  }

  sigma1sq <- sigma1^2
  sigma2sq <- sigma2^2

  beta1 <- -1/expm1(-1/(2*sigma1sq))
  beta2 <- -1/expm1(-1/(2*sigma2sq))

  C1_ <- C1*beta1/(2*sigma1sq)
  C2_ <- C2*beta2/(2*sigma2sq)

  H  <- KernelX[idx, , drop = FALSE]
  G  <- KernelX[-idx, , drop = FALSE]
  Hn <- nrow(H)
  Gn <- nrow(G)
  xp    <- ncol(KernelX)

  e1 <- matrix(1, Hn)
  e2 <- matrix(1, Gn)

  coef1 <- matrix(0, xp)
  coef2 <- matrix(0, xp)

  if (method_ls == "normal-chol") {
    HTH <- crossprod(H)
    GTG <- crossprod(G)
  }

  Omegavec1 <- rep(1, Gn)
  obj1      <- NULL
  for (i in 0:hq.steps) {
    G_omega_lam   <- G*(Omegavec1*C1_)
    if (method_ls == "qr") {
      coef1_new     <- -stack_qrsolve(rbind(H, sqrt(Omegavec1*C1_) * G),
                                      crossprod(G_omega_lam, e2), C3)
    } else if (method_ls == "normal-chol") {
      A             <- HTH + crossprod(G_omega_lam, G)
      diag(A)       <- diag(A) + C3
      coef1_new     <- -cholsolve(A, crossprod(G_omega_lam, e2))
    } else if (method_ls == "svd") {
      coef1_new     <- -svd_solve(rbind(H, sqrt(Omegavec1*C1_) * G),
                                      crossprod(G_omega_lam, e2), C3)
    }
    f         <- KernelX %*% coef1_new
    Omegavec1 <- update_omega(1 + f[-idx], sigma1sq)
    obj1_new  <- compute_objective1(coef1_new, f, idx, C3, C1, sigma1sq, beta1)
    if (!is.null(obj1) && function_tolerance(obj1_new, obj1, eps.hq)) {
      coef1 <- coef1_new
      break
    }
    coef1     <- coef1_new
    obj1      <- obj1_new
  }

  Omegavec2 <- rep(1, Hn)
  obj2      <- NULL
  for (i in 0:hq.steps) {
    H_omega_lam   <- H*(Omegavec2*C2_)
    if (method_ls == "qr") {
      coef2_new     <- stack_qrsolve(rbind(G, sqrt(Omegavec2*C2_) * H),
                                     crossprod(H_omega_lam, e1), C4)
    } else if (method_ls == "normal-chol") {
      B             <- GTG + crossprod(H_omega_lam, H)
      diag(B)       <- diag(B) + C4
      coef2_new     <- cholsolve(B, crossprod(H_omega_lam, e1))
    } else if (method_ls == "svd") {
      coef2_new     <- svd_solve(rbind(G, sqrt(Omegavec2*C2_) * H),
                                     crossprod(H_omega_lam, e1), C4)
    }
    f         <- KernelX %*% coef2_new
    Omegavec2 <- update_omega(1 - f[idx], sigma2sq)
    obj2_new  <- compute_objective2(coef2_new, f, idx, C4, C2, sigma2sq, beta2)
    if (!is.null(obj2) && function_tolerance(obj2_new, obj2, eps.hq)) {
      coef2 <- coef2_new
      break
    }
    coef2     <- coef2_new
    obj2      <- obj2_new
  }

  BaseDualClossTBSVMClassifier <- list("coef1" = as.matrix(coef1),
                                       "coef2" = as.matrix(coef2))
  return(BaseDualClossTBSVMClassifier)
}

#' C-loss Least Squares Twin Bounded Support Vector Machine
#'
#' \code{closs_tbsvm} is an R implementation of CTSVM
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
#' @param eps.hq the tolerance of the stopping criterion: the half-quadratic
#'             iterations stop when the relative change of the objective
#'             function value is below \code{eps.hq}.
#' @param hq.steps the maximum number of half-quadratic iterations. The solver
#'             first solves the unit-weight least squares system as an
#'             initialization, then performs at most \code{hq.steps}
#'             reweighted iterations.
#' @param fit_intercept if set \code{fit_intercept = TRUE},
#'                      the function will evaluates intercept.
#' @param reduce_set reduce set for reduce SVM, default \code{reduce_set = NULL}.
#' @param method_ls solver for the internal linear systems: \code{"svd"}
#'                  (singular value decomposition of the stacked system with
#'                  ridge damping, never forms Gram matrices, the default),
#'                  \code{"qr"} (stacked QR decomposition, never forms Gram
#'                  matrices) or \code{"normal-chol"} (Cholesky on the normal
#'                  equations).
#' @return return \code{TSVMClassifier} object.
#' @references Xiaohan Zheng, Li Zhang, Leilei Yan (2021).
#'             CTSVM: A robust twin support vector machine with correntropy-induced
#'             loss function for binary classification problems.
#'             \emph{Information Sciences}, 559, 22-45.
#'             \doi{10.1016/j.ins.2021.01.006}
#' @export
closs_tbsvm <- function(X, y, C1 = 1, C2 = C1, C3 = 1e-7, C4 = C3,
                        kernel = c("linear", "rbf", "poly"),
                        gamma = 1 / ncol(X), degree = 3, coef0 = 0,
                        sigma1 = 1, sigma2 = sigma1,
                        eps.hq = 1e-2, hq.steps = 10,
                        fit_intercept = TRUE, reduce_set = NULL,
                        method_ls = c("svd", "qr", "normal-chol")) {
  method_ls <- match.arg(method_ls)
  X         <- as.matrix(X)
  y         <- as.matrix(y)
  class_set <- sort(unique(y))
  idx       <- which(y == class_set[1])
  y[idx]    <- 1
  y[-idx]   <- -1
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
                                        eps.hq, hq.steps,
                                        method_ls)

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
