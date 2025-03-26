#' Bounded Asymmetric Least Squares Loss Support Vector Regression
#'
#' \code{bals_svr} is an R implementation of BALS-SVR
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
#' @param fit_intercept if set \code{fit_intercept = TRUE},
#'                      the function will evaluates intercept.
#' @return return \code{SVMRegressor} object.
#' @export
bals_svr <- function(X, y, C = 1, kernel = c("linear", "rbf", "poly"),
                     gamma = 1 / ncol(X), degree = 3, coef0 = 0,
                     p = 0.5, lambda = 1,
                     eps = 1e-5, eps.hq = 1e-2, max.steps = 4000, hq.steps = 10,
                     fit_intercept = TRUE) {
  bals_svr_dual_solver <- function(KernelX, y, C, lambda, p, eps, eps.hq,
                                   max.steps, hq.steps) {
    v_update <- function(f, lambda, p) {
      lam_alsf <- rep(0, length(f))
      idx <- which(f >= 0)
      lam_alsf[idx]  <- (lambda*p)*f[idx]^2
      lam_alsf[-idx] <- (lambda*((1 - p)))*f[-idx]^2
      v <- -1 / (1 + lam_alsf)^2
      return(v)
    }
    n      <- nrow(KernelX)
    H0     <- KernelX
    H      <- cbind(H0, -H0)
    H      <- rbind(H, -H)
    q      <- rbind(y, -y)
    lb     <- matrix(0, 2*n)
    ub     <- matrix(Inf, 2*n)
    u0     <- lb
    v  <- matrix(-1, n)
    Omega_mat <- matrix(0, n, hq.steps + 1)
    Omega <- as.numeric((C*lambda)*(-v))
    Omega_mat[, 1] <- Omega
    I1 <- diag(1/(Omega*(p)), n)
    I12 <- diag(0, n)
    I4 <- diag(1/(Omega*(1 - p)), n)
    R <- rbind(cbind(I1, I12), cbind(I12, I4))
    H_R <- H + R
    for (i in 1:hq.steps) {
      u <- SupportVectorLab::clip_dcd_optimizer(H_R, q, lb, ub, eps, max.steps, u0)$x
      if (norm(u - u0, type = "2") < eps.hq) {
        break
      } else {
        u0 <- u
      }
      f <- y - H0 %*% (u0[1:n] - u0[(n + 1):(2 * n)])
      v <- v_update(f, lambda, p)
      Omega <- as.numeric((C*lambda)*(-v))
      Omega_mat[, (i + 1)] <- Omega
      I1 <- diag(1/(Omega*(p)), n)
      I4 <- diag(1/(Omega*(1 - p)), n)
      R <- rbind(cbind(I1, I12), cbind(I12, I4))
      H_R <- H + R
    }
    coef <- u[1:n] - u[(n + 1):(2*n)]
    BaseDualBALSSVRRegressor <- list("coef" = as.matrix(coef),
                                     "Omega_mat" = Omega_mat)
    class(BaseDualBALSSVRRegressor) <- "BaseDualBALSSVRRegressor"
    return(BaseDualBALSSVRRegressor)
  }
  X <- as.matrix(X)
  y <- as.matrix(y)
  kernel <- match.arg(kernel)
  solver <- "dual"
  if (fit_intercept == TRUE) {
    X <- cbind(X, 1)
  }
  kso <- kernel_select_option(X, kernel, solver, 1, gamma, degree, coef0)
  KernelX <- kso$KernelX
  X <- kso$X
  if (solver == "dual") {
    solver.res <- bals_svr_dual_solver(KernelX, y, C,
                                       lambda, p, eps, eps.hq,
                                       max.steps, hq.steps)
  }
  SVMRegressor <- list("X" = X, "y" = y,
                       "C" = C, "kernel" = kernel,
                       "gamma" = gamma, "degree" = degree, "coef0" = coef0,
                       "solver" = solver, "coef" = solver.res$coef,
                       "fit_intercept" = fit_intercept,
                       "Omega_mat" = solver.res$Omega_mat)
  class(SVMRegressor) <- "SVMRegressor"
  return(SVMRegressor)
}
