#'  Primal Estimated sub-Gradient SOlver for SVM (Pegasos)
#'
#' @author Zhang Jiaqi.
#' @param X,y dataset and label.
#' @param w initial point.
#' @param m mini-batch size for pegasos solver.
#' @param max.steps the number of iterations to solve the optimization problem.
#' @param fx sub-gradient of objective function.
#' @param pars parameters list for the sub-gradient.
#' @param projection projection option.
#' @param ... additional settings for the sub-gradient.
#' @return return optimal solution.
#' @references ${1:Pegasos: Primal Estimated sub-GrAdient SOlver for SVM}
#' @export
pegasos <- function(X, y, w, m, max.steps, fx, pars,
                    projection = TRUE, ...) {
  C <- pars$C
  sample_seed <- list(...)$sample_seed
  if (is.null(sample_seed) == FALSE) {
    set.seed(sample_seed)
  }
  nx <- nrow(X)
  px <- ncol(X)
  for (t in 1:max.steps) {
    At <- sample(nx, m)
    xm <- X[At, ]
    xm <- X[At, ]
    dim(xm) <- c(m, px)
    ym <- as.matrix(y[At])
    # update parameter
    dF <- fx(xm, ym, w, pars, At = At)
    w <- w - (1/t)*dF
    if (projection == TRUE) {
      w <- min(1, sqrt(C)/norm(w, type = "2"))*w
    }
  }
  return(w)
}

#' Nesterov's Accelerated Gradient Descent
#'
#' @author Zhang Jiaqi.
#' @param X,y dataset and label.
#' @param w initial point.
#' @param m mini-batch size for pegasos solver.
#' @param max.steps the number of iterations to solve the optimization problem.
#' @param fx sub-gradient of objective function.
#' @param pars parameters list for the sub-gradient.
#' @param v initial velocity.
#' @param lr initial learning rate.
#' @param gam momentum parameter.
#' @param decay_option decay option for learning rate.
#' @param ... additional settings for the sub-gradient.
#' @return return optimal solution.
#' @export
nesterov <- function(X, y, w, m, max.steps, fx, pars,
                     v = matrix(0, nrow(w)), lr = 1, gam = 0.5,
                     decay_option = exp_decay, ...) {
  sample_seed <- list(...)$sample_seed
  if (is.null(sample_seed) == FALSE) {
    set.seed(sample_seed)
  }
  nx <- nrow(X)
  px <- ncol(X)
  for (t in 1:max.steps) {
    At <- sample(nx, m)
    xm <- X[At, ]
    dim(xm) <- c(m, px)
    ym <- as.matrix(y[At])
    v <- gam*v - lr*fx(xm, ym, w + gam*v, pars, At = At)
    w <- w + v
    if (is.null(decay_option) == FALSE) {
      lr <- decay_option(lr, steps = t, ...)
    }
  }
  return(w)
}

#' Rmsprop Optimizer
#'
#' @author Zhang Jiaqi.
#' @param X,y dataset and label.
#' @param w initial point.
#' @param m mini-batch size for pegasos solver.
#' @param max.steps the number of iterations to solve the optimization problem.
#' @param fx sub-gradient of objective function.
#' @param pars parameters list for the sub-gradient.
#' @param lr initial stepsize.
#' @param rho momentum parameter.
#' @param delta avoid division by 0.
#' @param ... additional settings for the sub-gradient.
#' @return return optimal solution.
#' @export
rmsprop <- function(X, y, w, m, max.steps, fx, pars,
                    lr = 0.001, rho = 0.9, delta = 1e-5, ...) {
  sample_seed <- list(...)$sample_seed
  if (is.null(sample_seed) == FALSE) {
    set.seed(sample_seed)
  }
  xn <- nrow(X)
  xp <- ncol(X)
  r <- matrix(0.1, xp, 1)
  g <- r
  for (t in 1:max.steps) {
    At <- sample(xn, m)
    xm <- X[At, ]
    dim(xm) <- c(m, xp)
    ym <- as.matrix(y[At])
    dF <- fx(xm, ym, w, pars, At = At)
    r <- rho*r + (1 - rho)*dF*dF
    w <- w - lr*dF/(sqrt(r + delta))
  }
  return(w)
}

#'  Conjugate Gradient Method for Solving Linear Equation Ax = b
#'
#' @author Zhang Jiaqi.
#' @param A,b matrix of linear equation Ax = b (Since A is a PSD matrix).
#' @param x initial point.
#' @param max.steps the number of iterations to solve the optimization problem.
#' @param eps the precision of the optimization algorithm.
#' @param ... additional settings for the sub-gradient.
#' @return return optimal solution.
#' @export
conjugate_gradient_method <- function(A, b, x, max.steps, eps = 1e-5, ...) {
  rk <- b - A%*%x
  pk <- rk
  for (t in 1:max.steps) {
    rk2 <- (t(rk) %*% rk)
    alphak <-  as.numeric(rk2 / (t(pk) %*% A %*% pk))
    x <- x + alphak*pk
    rk_1 <- rk - alphak*(A%*%pk)
    if (norm(rk_1, type = "2") < eps) {
      cat("converge after", t, "steps", "\n")
      break
    } else {
      betak <- as.numeric((t(rk_1) %*% rk_1) / rk2)
      pk <- rk_1 + betak*pk
    }
    rk <- rk_1
  }
  return(x)
}


#' Adaptive Moment Estimation Optimizer
#'
#' @author Zhang Jiaqi.
#' @param X,y dataset and label.
#' @param w initial point.
#' @param m mini-batch size for pegasos solver.
#' @param max.steps the number of iterations to solve the optimization problem.
#' @param fx sub-gradient of objective function.
#' @param pars parameters list for the sub-gradient.
#' @param lr initial learning rate.
#' @param beta1 first order moment parameter.
#' @param beta2 second order moment parameter.
#' @param delta avoid division by 0.
#' @param ... additional settings for the sub-gradient.
#' @return return optimal solution.
#' @export
adam <- function(X, y, w, m, max.steps, fx, pars,
                 lr = 0.001, beta1 = 0.9, beta2 = 0.999, delta = 1e-5, ...) {
  sample_seed <- list(...)$sample_seed
  if (is.null(sample_seed) == FALSE) {
    set.seed(sample_seed)
  }
  xn <- nrow(X)
  xp <- ncol(X)
  mt <- matrix(0, xp, 1)
  vt <- mt
  for (t in 1:max.steps) {
    At <- sample(xn, m)
    xm <- X[At, ]
    dim(xm) <- c(m, xp)
    ym <- as.matrix(y[At])
    dF <- fx(xm, ym, w, pars, At = At)
    mt <- beta1*mt + (1 - beta1)*dF
    vt <- beta2*vt + (1 - beta2)*(dF*dF)
    mt_hat <- mt/(1 - beta1^t)
    vt_hat <- vt/(1 - beta2^t)
    w <- w - lr*(mt_hat/(sqrt(vt_hat) + delta)) * dF
  }
  return(w)
}


exp_decay <- function(lr, steps, s = 0.001, ...) {
  lr <- lr*exp(-s*steps)
  return(lr)
}

