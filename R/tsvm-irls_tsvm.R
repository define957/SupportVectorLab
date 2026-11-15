update_weight <- function(weight_func, r, weight_params) {
  return(do.call(weight_func, c(list(r = r), as.list(weight_params))))
}

tsvm_irls_default_loss <- function(r_within, r_between, C_between) {
  return(0.5*sum(r_within^2) + 0.5*C_between * sum(r_between^2))
}

irls_tsvm_dual_solver <- function(
    KernelX,
    idx,
    C1,
    C2,
    C3,
    C4,
    weight_func1,
    weight_params1,
    weight_func2,
    weight_params2,
    obj_func,
    obj_params,
    eps.irls,
    irls.steps,
    method_ls
  ) {

  compute_objective <- function(coef, r_within, r_between, C_reg, C_between) {
    loss <- do.call(obj_func,
                    c(list(r_within = r_within, r_between = r_between,
                           C_between = C_between), obj_params))
    obj <- 0.5*C_reg * crossprod(coef) + loss
    return(obj)
  }

  H     <- KernelX[idx, , drop = FALSE]
  G     <- KernelX[-idx, , drop = FALSE]
  Hn    <- nrow(H)
  Gn    <- nrow(G)

  e1    <- matrix(1, Hn)
  e2    <- matrix(1, Gn)

  w1    <- rep(1, Hn)
  w2    <- rep(1, Gn)
  obj1  <- NULL
  for (i in 0:irls.steps) {
    sw1 <- sqrt(w1)
    sw2 <- sqrt(w2)
    GW  <- sw2 * G
    if (method_ls == "qr") {
      coef1_new <- -stack_qrsolve(rbind(sw1/sqrt(C1) * H, GW),
                                  crossprod(GW, sw2 * e2), C3)
    } else if (method_ls == "normal-chol") {
      A         <- crossprod(H*(w1/C1), H) + crossprod(G*(w2), G)
      diag(A)   <- diag(A) + C3
      coef1_new <- -cholsolve(A, crossprod(G*(w2), e2))
    }
    f         <- KernelX %*% coef1_new
    w1        <- update_weight(weight_func1, f[idx], weight_params1)
    w2        <- update_weight(weight_func2, 1 + f[-idx], weight_params2)
    obj1_new  <- compute_objective(coef1_new, f[idx], 1 + f[-idx], C3, C1)
    if (!is.null(obj1) && function_tolerance(obj1_new, obj1, eps.irls)) {
      coef1 <- coef1_new
      break
    }
    coef1     <- coef1_new
    obj1      <- obj1_new
  }

  w1    <- rep(1, Gn)
  w2    <- rep(1, Hn)
  obj2  <- NULL
  for (i in 0:irls.steps) {
    sw1 <- sqrt(w1)
    sw2 <- sqrt(w2)
    HW  <- sw2 * H
    if (method_ls == "qr") {
      coef2_new <- stack_qrsolve(rbind(sw1/sqrt(C2) * G, HW),
                                 crossprod(HW, sw2 * e1), C4)
    } else if (method_ls == "normal-chol") {
      B         <- crossprod(G*(w1/C2), G) + crossprod(H*(w2), H)
      diag(B)   <- diag(B) + C4
      coef2_new <- cholsolve(B, crossprod(H*(w2), e1))
    }
    f         <- KernelX %*% coef2_new
    w1        <- update_weight(weight_func1, f[-idx], weight_params1)
    w2        <- update_weight(weight_func2, 1 - f[idx], weight_params2)
    obj2_new  <- compute_objective(coef2_new, f[-idx], 1 - f[idx], C4, C2)
    if (!is.null(obj2) && function_tolerance(obj2_new, obj2, eps.irls)) {
      coef2 <- coef2_new
      break
    }
    coef2     <- coef2_new
    obj2      <- obj2_new
  }

  BaseDualIRLSTSVMClassifier <- list("coef1" = as.matrix(coef1),
                                     "coef2" = as.matrix(coef2))
  return(BaseDualIRLSTSVMClassifier)
}

#' Iteratively Reweighted Least Squares Twin Support Vector Machine
#'
#' \code{irls_tsvm} is an R implementation of IRWLS-LS-TSVM. Both the
#' within-class (own-class) and the between-class (other-class) least
#' squares losses are reweighted at every iteration by user-supplied
#' weight functions.
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
#' @param weight_func1 weight function for the within-class (own-class) loss.
#'                     It must accept the residual vector as its first
#'                     argument \code{r} and return a non-negative weight
#'                     vector of the same length, see the internal \code{capped_lp_weight}.
#' @param weight_params1 additional parameters passed to \code{weight_func1}
#'                       together with the residual \code{r} via \code{do.call}:
#'                       \code{do.call(weight_func1, c(list(r = r), weight_params1))}.
#'                       Numerical stability of the weight function (e.g.
#'                       division-by-zero floors) is the user's responsibility.
#' @param weight_func2 weight function for the between-class (other-class)
#'                     loss, with the same interface as \code{weight_func1}.
#' @param weight_params2 additional parameters passed to \code{weight_func2},
#'                       see \code{weight_params1}.
#' @param obj_func objective (loss) function used to monitor the IRLS
#'                 iterations. It is called as
#'                 \code{do.call(obj_func, c(list(r_within = r_within, r_between = r_between, C_between = C_between), obj_params))}
#'                 and must return the scalar data-fitting loss (the
#'                 regularization term and the outer factor \eqn{1/2} are
#'                 added by the solver). The default \code{squared_loss} is
#'                 the plain least squares objective, matching the default
#'                 unit weights; when custom weight functions implement the
#'                 IRLS weights of some loss \eqn{\rho}{rho}, i.e.
#'                 \eqn{w(r) = \rho'(r)/(2r)}{w(r) = rho'(r)/(2r)}, pass the
#'                 corresponding loss here for a theoretically sound stopping
#'                 criterion.
#' @param obj_params additional parameters passed to \code{obj_func},
#'                   see the internal \code{squared_loss} for an example.
#' @param eps.irls the tolerance of the stopping criterion: the IRLS
#'             iterations stop when the relative change of the objective
#'             function value is below \code{eps.irls}.
#' @param irls.steps the maximum number of IRLS iterations. The solver first
#'                  solves the unit-weight least squares system as an
#'                  initialization, then performs at most \code{irls.steps}
#'                  reweighted iterations.
#' @param fit_intercept if set \code{fit_intercept = TRUE},
#'                      the function will evaluates intercept.
#' @param reduce_set reduce set for reduce SVM, default \code{reduce_set = NULL}.
#' @param method_ls solver for the internal linear systems: \code{"qr"} (stacked
#'                  QR decomposition, never forms Gram matrices) or
#'                  \code{"normal-chol"} (Cholesky on the normal equations).
#' @return return \code{TSVMClassifier} object.
#' @references Yuan-Hai Shao, Chun-Hua Zhang, Zhen Wang, Nai-Yang Deng (2011).
#'             Improvements on twin support vector machines.
#'             \emph{IEEE Transactions on Neural Networks}, 22(6), 962-968.
#' @export
irls_tsvm <- function(
    X,
    y,
    C1             = 1,
    C2             = C1,
    C3             = 1e-5,
    C4             = C3,
    kernel         = c("linear", "rbf", "poly"),
    gamma          = 1 / ncol(X),
    degree         = 3,
    coef0          = 0,
    weight_func1   = function(r) { rep(1, length(r)) },
    weight_params1 = NULL,
    weight_func2   = function(r) { rep(1, length(r)) },
    weight_params2 = NULL,
    obj_func       = tsvm_irls_default_loss,
    obj_params     = list(),
    eps.irls       = 1e-4,
    irls.steps     = 100,
    fit_intercept  = TRUE,
    reduce_set     = NULL,
    method_ls      = c("qr", "normal-chol")) {

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
  solver.res <- irls_tsvm_dual_solver(
    KernelX,
    idx,
    C1,
    C2,
    C3,
    C4,
    weight_func1,
    weight_params1,
    weight_func2,
    weight_params2,
    obj_func,
    obj_params,
    eps.irls,
    irls.steps,
    method_ls
  )

  model_specs = list("X" = X, "y" = y,
                     "C1" = C1, "C2" = C2, "C3" = C3, "C4" = C4,
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
  TSVMClassifier <- structure(
    list("model_specs" = model_specs,
         "model_coef" = model_coef,
         "kernel_config" = kernel_config),
    "class" = "TSVMClassifier")
  return(TSVMClassifier)
}
