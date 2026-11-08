cL2p_ls_tsvm_dual_solver <- function(KernelX, idx, C1, C2,
                                     p, epsilon1, epsilon2,
                                     eps.irls, irls.steps, irls.reg,
                                     method_ls) {

  cl2p_loss_function <- function(z, p, epsilon) {
    return(sum(pmin(abs(z)^p, epsilon)))
  }
  ls_loss_function <- function(z) {
    return(sum(z^2))
  }

  update_weight <- function(f, p, epsilon, irls.reg) {
    abs_f <- abs(f)
    abs_f <- pmax(abs_f, irls.reg)
    fL2p <- abs_f^p
    weight_elem <- (p/2)*(abs_f)^(p - 2)
    weight_elem[fL2p >= epsilon] <- 0
    return(as.vector(weight_elem))
  }

  compute_objective1 <- function(coef, f, idx, C_between, p, epsilon) {
    obj <- 0.5*((1e-7) * crossprod(coef) + cl2p_loss_function(f[idx], p, epsilon) +
           C_between * ls_loss_function(1 + f[-idx]))
    return(obj)
  }

  compute_objective2 <- function(coef, f, idx, C_between, p, epsilon) {
    obj <- 0.5*((1e-7) * crossprod(coef) + cl2p_loss_function(f[-idx], p, epsilon) +
           C_between * ls_loss_function(1 - f[idx]))
    return(obj)
  }

  H     <- KernelX[idx, , drop = FALSE]
  G     <- KernelX[-idx, , drop = FALSE]
  Hn    <- nrow(H)
  Gn    <- nrow(G)
  xp    <- ncol(KernelX)

  e1    <- matrix(1, Hn)
  e2    <- matrix(1, Gn)

  coef1 <- matrix(0, xp)
  coef2 <- matrix(0, xp)

  if (method_ls == "normal-chol") {
    HTH   <- crossprod(H)
    GTG   <- crossprod(G)
  }
  HT_e1 <- crossprod(H, e1)
  GT_e2 <- crossprod(G, e2)

  Fvec      <- rep(1, Hn)
  obj1      <- NULL
  for (i in 0:irls.steps) {
    if (method_ls == "qr") {
      coef1_new <- -stack_qrsolve(rbind(sqrt(Fvec/C1) * H, G), GT_e2, 1e-7)
    } else if (method_ls == "normal-chol") {
      A         <- crossprod(H*(Fvec/C1), H) + GTG
      diag(A)   <- diag(A) + 1e-7
      coef1_new <- -cholsolve(A, GT_e2)
    } else if (method_ls == "svd") {
      coef1_new <- -svd_solve(rbind(sqrt(Fvec/C1) * H, G), GT_e2, 1e-7)
    }
    f         <- KernelX %*% coef1_new
    Fvec      <- update_weight(f[idx], p, epsilon1, irls.reg)
    obj1_new  <- compute_objective1(coef1_new, f, idx, C1, p, epsilon1)
    if (!is.null(obj1) && function_tolerance(obj1_new, obj1, eps.irls)) {
      coef1 <- coef1_new
      break
    }
    coef1     <- coef1_new
    obj1      <- obj1_new
  }

  Dvec      <- rep(1, Gn)
  obj2      <- NULL
  for (i in 0:irls.steps) {
    if (method_ls == "qr") {
      coef2_new <- stack_qrsolve(rbind(sqrt(Dvec / C2) * G, H), HT_e1, 1e-7)
    } else if (method_ls == "normal-chol") {
      B         <- crossprod(G*(Dvec/C2), G) + HTH
      diag(B)   <- diag(B) + 1e-7
      coef2_new <- cholsolve(B, HT_e1)
    } else if (method_ls == "svd") {
      coef2_new <- svd_solve(rbind(sqrt(Dvec / C2) * G, H), HT_e1, 1e-7)
    }
    f         <- KernelX %*% coef2_new
    Dvec      <- update_weight(f[-idx], p, epsilon2, irls.reg)
    obj2_new  <- compute_objective2(coef2_new, f, idx, C2, p, epsilon2)
    if (!is.null(obj2) && function_tolerance(obj2_new, obj2, eps.irls)) {
      coef2 <- coef2_new
      break
    }
    coef2     <- coef2_new
    obj2      <- obj2_new
  }

  BaseDualCL2PLSTSVMClassifier <- list("coef1" = as.matrix(coef1),
                                       "coef2" = as.matrix(coef2))
  return(BaseDualCL2PLSTSVMClassifier)
}

#' Capped L2,p Least Squares Twin Support Vector Machine
#'
#' \code{cL2p_ls_tsvm} is an R implementation of CL2,p-LS-TSVM
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
#' @param p parameter for L2,p norm.
#' @param epsilon1,epsilon2 truncation parameters for the capped L2,p-norm.
#' @param eps.irls the tolerance of the stopping criterion: the IRLS
#'             iterations stop when the relative change of the objective
#'             function value is below \code{eps.irls}.
#' @param irls.steps the maximum number of IRLS iterations. The solver first
#'                  solves the unit-weight least squares system as an
#'                  initialization, then performs at most \code{irls.steps}
#'                  reweighted iterations.
#' @param fit_intercept if set \code{fit_intercept = TRUE},
#'                      the function will evaluates intercept.
#' @param irls.reg \code{irls.reg} is used to avoid division-by-zero issues when updating weights.
#' @param reduce_set reduce set for reduce SVM, default \code{reduce_set = NULL}.
#' @param method_ls solver for the internal linear systems: \code{"svd"}
#'                  (singular value decomposition of the stacked system with
#'                  ridge damping, never forms Gram matrices, the default),
#'                  \code{"qr"} (stacked QR decomposition, never forms Gram
#'                  matrices) or \code{"normal-chol"} (Cholesky on the normal
#'                  equations).
#' @return return \code{TSVMClassifier} object.
#' @references Chao Yuan, Liming Yang (2021).
#'             Capped L2,p-norm metric based robust least squares twin support
#'             vector machine for pattern classification.
#'             \emph{Neural Networks}, 142, 457-478.
#'             \doi{10.1016/j.neunet.2021.06.028}
#' @export
cL2p_ls_tsvm <- function(X, y, C1 = 1, C2 = C1,
                         kernel = c("linear", "rbf", "poly"),
                         gamma = 1 / ncol(X), degree = 3, coef0 = 0,
                         p = 1, epsilon1 = 1, epsilon2 = epsilon1,
                         eps.irls = 1e-2, irls.steps = 10, irls.reg = 1e-4,
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
  solver.res <- cL2p_ls_tsvm_dual_solver(KernelX, idx, C1, C2,
                                         p, epsilon1, epsilon2,
                                         eps.irls, irls.steps, irls.reg,
                                         method_ls)
  model_specs = list("X" = X, "y" = y,
                     "C1" = C1, "C2" = C2,
                     "fit_intercept" = fit_intercept,
                     "class_set" = class_set)
  model_coef = list("coef1" = solver.res$coef1,
                    "coef2" = solver.res$coef2)
  kernel_config = list("kernel" = kernel,
                       "gamma" = gamma,
                       "degree" = degree,
                       "coef0" = coef0,
                       "reduce_set" = reduce_set,
                       "KernelR" = KernelR,
                       "KernelX" = KernelX[, 1:kxp, drop = FALSE])
  w_norm <- get_coef_norm(kernel_config, model_coef)
  kernel_config$w_norm <- w_norm
  TSVMClassifier <- structure(list("model_specs" = model_specs,
                                   "model_coef" = model_coef,
                                   "kernel_config" = kernel_config),
                              "class" = "TSVMClassifier")
  return(TSVMClassifier)
}
