cL2p_ls_tsvm_dual_solver <- function(KernelX, idx, C1, C2,
                                     p, epsilon1, epsilon2,
                                     eps.irls, irls.steps) {

  update_weight <- function(X, w, p, epsilon) {
    f <- X %*% w
    abs_f <- abs(f)
    fL2p <- abs_f^p
    weight_elem <- (p/2)*(abs_f + 1e-5)^(p - 2)
    weight_elem[fL2p >= epsilon] <- 0
    return(as.vector(weight_elem))
  }

  H     <- KernelX[-idx, , drop = FALSE]
  G     <- KernelX[idx, , drop = FALSE]
  Hn    <- nrow(H)
  Gn    <- nrow(G)
  xp    <- ncol(KernelX)

  Fvec  <- rep(1, Hn)
  Dvec  <- rep(1, Gn)

  e1    <- matrix(1, Hn)
  e2    <- matrix(1, Gn)

  coef1 <- matrix(0, xp)
  coef2 <- matrix(0, xp)

  HTH   <- t(H) %*% H
  GTG   <- t(G) %*% G
  HT_e1 <- t(H) %*% e1
  GT_e2 <- t(G) %*% e2

  for (i in 1:irls.steps) {
    A         <- t(H*(Fvec/C1)) %*% H + GTG
    diag(A)   <- diag(A) + 1e-7
    coef1_new <- -cholsolve(A, GT_e2)
    if (norm(coef1_new - coef1, type = "2") < eps.irls) {
      coef1 <- coef1_new
      break
    }
    Fvec  <- update_weight(H, coef1_new, p, epsilon1)
    coef1 <- coef1_new
  }

  for (i in 1:irls.steps) {
    B         <- t(G*(Dvec/C2)) %*% G + HTH
    diag(B)   <- diag(B) + 1e-7
    coef2_new <- cholsolve(B, HT_e1)
    if (norm(coef2_new - coef2, type = "2") < eps.irls) {
      coef2 <- coef2_new
      break
    }
    Dvec  <- update_weight(G, coef2_new, p, epsilon2)
    coef2 <- coef2_new
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
#' @param epsilon1,epsilon2 truncation parameter.
#' @param eps.irls the precision of the optimization algorithm.
#' @param irls.steps the number of iterations to solve the optimization problem.
#' @param fit_intercept if set \code{fit_intercept = TRUE},
#'                      the function will evaluates intercept.
#' @param reduce_set reduce set for reduce SVM, default \code{reduce_set = NULL}.
#' @return return \code{TSVMClassifier} object.
#' @export
cL2p_ls_tsvm <- function(X, y, C1 = 1, C2 = C1,
                         kernel = c("linear", "rbf", "poly"),
                         gamma = 1 / ncol(X), degree = 3, coef0 = 0,
                         p = 1, epsilon1 = 1, epsilon2 = epsilon1,
                         eps.irls = 1e-2, irls.steps = 10,
                         fit_intercept = TRUE, reduce_set = NULL) {

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
                                         eps.irls, irls.steps)
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
