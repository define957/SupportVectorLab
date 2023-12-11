hinge_tsvm_dual_solver <- function(KernelX, idx, C1, C2, eps, max.steps) {
  H <- KernelX[idx, ]
  G <- KernelX[-idx, ]
  xn <- nrow(KernelX)
  xp <- ncol(KernelX)
  Hn <- nrow(H)
  Gn <- xn - Hn
  invHTH <- chol2inv(chol(t(H) %*% H + diag(1e-7, xp)))
  dualH <- G %*% invHTH %*% t(G)
  dualq1 <- matrix(1, Gn)
  duallb1 <- matrix(0, Gn)
  dualub1 <- matrix(C1, Gn)
  u0 <- matrix(0, Gn)
  a <- clip_dcd_optimizer(dualH, dualq1, duallb1, dualub1,
                          eps, max.steps, u0)$x
  u <- - invHTH %*% t(G) %*% a

  invGTG <- chol2inv(chol(t(G) %*% G + diag(1e-7, xp)))
  dualH <- H %*% invGTG %*% t(H)
  dualq2 <- matrix(1, Hn)
  duallb2 <- matrix(0, Hn)
  dualub2 <- matrix(C2, Hn)
  u0 <- matrix(0, Hn)
  g <- clip_dcd_optimizer(dualH, dualq2, duallb2, dualub2,
                          eps, max.steps, u0)$x
  v <- invGTG %*% t(H) %*% g
  BaseDualHingeTSVMClassifier <- list("coef1" = as.matrix(u),
                                      "coef2" = as.matrix(v))
}

#' Hinge Twin Support Vector Machine
#'
#' \code{hinge_tsvm} is an R implementation of Hinge-TSVM
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
#' @param eps the precision of the optimization algorithm.
#' @param max.steps the number of iterations to solve the optimization problem.
#' @param solver \code{"dual"} is available.
#' @param fit_intercept if set \code{fit_intercept = TRUE},
#'                      the function will evaluates intercept.
#' @param randx parameter for reduce SVM, default \code{randx = 0.1}.
#' @param ... unused parameters.
#' @return return \code{HingeSVMClassifier} object.
#' @export
hinge_tsvm <- function(X, y, C1 = 1, C2 = 2,
                       kernel = c("linear", "rbf", "poly"),
                       gamma = 1 / ncol(X), degree = 3, coef0 = 0,
                       eps = 1e-14, max.steps = 4000,
                       solver = c("dual"), fit_intercept = TRUE,
                       randx = 0.1, ...) {
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
    solver.res <- hinge_tsvm_dual_solver(KernelX, idx, C1, C2, eps, max.steps)
  }
  TSVMClassifier <- list("X" = X, "y" = y, "class_set" = class_set,
                         "C1" = C1, "C2" = C2, "kernel" = kernel,
                         "gamma" = gamma, "degree" = degree, "coef0" = coef0,
                         "solver" = solver, "coef1" = solver.res$coef1,
                         "coef2" = solver.res$coef2,
                         "fit_intercept" = fit_intercept,
                         "Kw" = Kw)
  class(TSVMClassifier) <- "TSVMClassifier"
  return(TSVMClassifier)
}

#' Predict Method for Twin Support Vector Machine
#'
#' @author Zhang Jiaqi
#' @param object a fitted object of class inheriting from \code{SVMClassifier}.
#' @param X new data for predicting.
#' @param values if set \code{values = TRUE}, this function will return predict
#'               values: f = abs(wx+b).
#' @param ... unused parameter.
#' @importFrom stats predict
#' @export
predict.TSVMClassifier <- function(object, X, values = FALSE, ...) {
  X <- as.matrix(X)
  if (object$kernel == "linear") {
    KernelX <- X
  } else {
    KernelX <- kernel_function(X, object$X,
                               kernel.type = object$kernel,
                               gamma = object$gamma,
                               degree = object$degree,
                               coef0 = object$coef0)
  }
  xp <- ncol(KernelX)
  if (object$fit_intercept == TRUE) {
    KernelXe <- cbind(KernelX, 1)
  }
  f1 <- abs(KernelXe %*% object$coef1)
  f2 <- abs(KernelXe %*% object$coef2)
  if (object$kernel == "linear") {
    norm1 <- norm(object$coef1[1:xp], type = "2")
    norm2 <- norm(object$coef2[1:xp], type = "2")
  } else {
    norm1 <- t(object$coef1[1:xp]) %*% object$Kw %*% object$coef1[1:xp]
    norm2 <- t(object$coef2[1:xp]) %*% object$Kw %*% object$coef2[1:xp]
    norm1 <- as.numeric(norm1)
    norm2 <- as.numeric(norm2)
  }
  if (norm1 == 0 || norm2 == 0) {
    norm1 <- norm1 + 1e-7
    norm2 <- norm2 + 1e-7
  }
  fx1 <- f1 / norm1
  fx2 <- f2 / norm2
  if (values == FALSE) {
    decf <- apply(cbind(fx1, fx2), 1, which.min)
    idx_pos <- which(decf == 1)
    idx_neg <- which(decf == 2)
    decf[idx_pos] <- object$class_set[1]
    decf[idx_neg] <- object$class_set[2]
  }
  if (values == FALSE) {
    decf <- apply(cbind(fx1, fx2), 1, which.min)
    idx_pos <- which(decf == 1)
    idx_neg <- which(decf == 2)
    decf[idx_pos] <- object$class_set[1]
    decf[idx_neg] <- object$class_set[2]
  } else {
    dec_values1 <- fx1
    dec_values2 <- fx2
    return(cbind(dec_values1, dec_values2))
  }
  return(decf)
}


#' Plot Method for Support Vector Machine
#'
#' @author Zhang Jiaqi
#' @param x a fitted object of class inheriting from \code{SVMClassifier}.
#' @param ... unused parameter.
#' @importFrom stats coef
#' @importFrom graphics abline grid points
#' @export
plot.TSVMClassifier <- function(x, ...) {
  y <- x$y
  coef1 <- x$coef1
  coef2 <- x$coef2
  class_set <- sort(unique(y))
  idx <- which(y == class_set[1])
  y[idx] <- -1
  y[-idx] <- 1
  xlim_c <- c(min(x$X[,1]), max(x$X[, 1]))
  ylim_c <- c(min(x$X[,2]), max(x$X[, 2]))
  if (length(coef1) ==3 && length(coef2) == 3) {
    plot(x$X[idx, 1], x$X[idx, 2], col = "red", xlim = xlim_c, ylim = ylim_c,
         xlab = "", ylab = "")
    grid(10, 10, lwd = 2,col = "grey")
    points(x$X[-idx, 1], x$X[-idx, 2], col = "blue")
    if (x$kernel == "linear") {
      abline(a = -coef1[3]/coef1[2], b = -coef1[1]/coef1[2],
             lty = 2, col = "red")
      abline(a = -(coef1[3]+1)/coef1[2], b = -coef1[1]/coef1[2],
             lty = 2, col = "red")
      abline(a = -coef2[3]/coef2[2], b = -coef2[1]/coef2[2],
             lty = 2, col = "blue")
      abline(a = -(coef2[3]-1)/coef2[2], b = -coef2[1]/coef2[2],
             lty = 2, col = "blue")
    }
  }
}
