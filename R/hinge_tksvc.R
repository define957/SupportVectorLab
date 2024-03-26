hinge_tksvc_dual_solver <- function(KernelX, y, C, epsilon, class_set, class_num,
                                    eps, max.steps) {
  xp <- ncol(KernelX)
  num_classifier <- class_num*(class_num - 1)/2
  coef_pos <- matrix(0, xp, num_classifier)
  coef_neg <- matrix(0, xp, num_classifier)

  classifier_idx <- 1
  for (i in 1:(class_num - 1)) {
    for (j in (i + 1):class_num) {
      idxPos <- which(y == class_set[i]) # positive samples
      idxNeg <- which(y == class_set[j]) # negtive samples
      idxRest <- which(y != class_set[i] & y != class_set[j]) # rest samples

      # Hyperplane 1
      H <- KernelX[idxPos, ]
      G <- KernelX[idxNeg, ]
      M <- KernelX[idxRest, ]
      N <- rbind(G, M)

      Hn <- nrow(H)
      Gn <- nrow(G)
      Mn <- nrow(M)

      invHTH <- chol2inv(chol(t(H) %*% H + diag(1e-7, xp)))
      dualH <- N %*% invHTH %*% t(N)
      dualq <- rbind(matrix(1, Gn), matrix(1 - epsilon, Mn))
      duallb <- matrix(0, Gn + Mn)
      dualub <- rbind(matrix(C[1], Gn), matrix(C[2], Mn))
      u0 <- matrix(0, Gn + Mn)
      x <- clip_dcd_optimizer(dualH, dualq, duallb, dualub,
                              eps, max.steps, u0)$x
      coef1 <- -invHTH %*% t(N) %*% x

      # Hyperplane 2
      P <- rbind(H, M)
      invGTG <- chol2inv(chol(t(G) %*% G + diag(1e-7, xp)))
      dualH <- P %*% invGTG %*% t(P)
      dualq <- rbind(matrix(1, Hn), matrix(1 - epsilon, Mn))
      dualub <- rbind(matrix(C[3], Hn), matrix(C[4], Mn))
      duallb <- matrix(0, Hn + Mn)
      u0 <- matrix(0, Hn + Mn)
      x <- clip_dcd_optimizer(dualH, dualq, duallb, dualub,
                              eps, max.steps, u0)$x
      coef2 <- invGTG %*% t(P) %*% x

      coef_pos[, classifier_idx] <- coef1
      coef_neg[, classifier_idx] <- coef2
      classifier_idx <- classifier_idx + 1
    }
  }
  BaseDualHingeTKSVCClassifier <- list("coef_pos" = coef_pos,
                                       "coef_neg" = coef_neg)
  return(BaseDualHingeTKSVCClassifier)
}

#' Hinge Twin Multi-Class Support Vector Machine
#'
#' \code{hinge_tksvc} is an R implementation of Hinge-TKSVC
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
#' @param gamma parameter for \code{'rbf'} and \code{'poly'} kernel. Default \code{gamma = 1/ncol(X)}.
#' @param degree parameter for polynomial kernel, default: \code{degree = 3}.
#' @param coef0 parameter for polynomial kernel,  default: \code{coef0 = 0}.
#' @param epsilon rest class parameter.
#' @param eps the precision of the optimization algorithm.
#' @param max.steps the number of iterations to solve the optimization problem.
#' @param solver \code{"dual"} is available.
#' @param fit_intercept if set \code{fit_intercept = TRUE},
#'                      the function will evaluates intercept.
#' @param randx parameter for reduce SVM, default \code{randx = 0.1}.
#' @return return \code{HingeSVMClassifier} object.
#' @export
hinge_tksvc <- function(X, y, C = 1,
                        kernel = c("linear", "rbf", "poly"),
                        gamma = 1 / ncol(X), degree = 3, coef0 = 0,
                        epsilon = 0.1,
                        eps = 1e-5, max.steps = 4000,
                        solver = c("dual"), fit_intercept = TRUE,
                        randx = 1) {
  C <- as.vector(C)
  X <- as.matrix(X)
  y <- as.matrix(y)
  class_set <- sort(unique(y))
  class_num <- length(class_set)
  if (length(C) == 1) {
    C <- matrix(C, 4)
  } else if (length(C) == 2) {
    C <- c(C, C)
  } else if (length(C) != 4) {
    stop("length(C) should equal to 1, 2 or 4!")
  }
  kernel <- match.arg(kernel)
  solver <- match.arg(solver)
  kso <- kernel_select_option(X, kernel, "primal", randx,
                              gamma, degree, coef0)
  KernelX <- kso$KernelX
  X <- kso$X
  if (fit_intercept == TRUE) {
    KernelX <- cbind(KernelX, 1)
  }
  if (solver == "dual") {
    solver.res <- hinge_tksvc_dual_solver(KernelX, y, C, epsilon,
                                          class_set, class_num, eps, max.steps)
  }
  TKSVCClassifier <- list("X" = X, "y" = y, "class_set" = class_set,
                          "class_num" = class_num,
                          "C" = C, "kernel" = kernel,
                          "gamma" = gamma, "degree" = degree, "coef0" = coef0,
                          "epsilon" = epsilon,
                          "solver" = solver, "coef_pos" = solver.res$coef_pos,
                          "coef_neg" = solver.res$coef_neg,
                          "fit_intercept" = fit_intercept)
  class(TKSVCClassifier) <- "TKSVCClassifier"
  return(TKSVCClassifier)
}

#' Predict Method for Twin Support Vector Machine
#'
#' @author Zhang Jiaqi
#' @param object a fitted object of class inheriting from \code{SVMClassifier}.
#' @param X new data for predicting.
#' @param ... unused parameter.
#' @importFrom stats predict
#' @export
predict.TKSVCClassifier <- function(object, X, ...) {
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
  if (object$fit_intercept == TRUE) {
    KernelX <- cbind(KernelX, 1)
  }
  xn <- nrow(KernelX)
  fx_pos <- KernelX %*% object$coef_pos
  fx_neg <- KernelX %*% object$coef_neg
  vote_mat <- matrix(0, xn, object$class_num)
  classifier_idx <- 1

  for (i in 1:(object$class_num - 1)) {
    for (j in (i + 1):object$class_num) {
      idx_pos <- which(fx_pos[, classifier_idx] > -1 + object$epsilon)
      idx_neg <- which(fx_neg[, classifier_idx] <  1 - object$epsilon)
      vote_mat[idx_pos, i] <- vote_mat[idx_pos, i] + 1
      vote_mat[idx_neg, j] <- vote_mat[idx_neg, i] + 1
      idx <- unique(c(idx_pos, idx_neg))
      vote_mat[-idx, i] <- vote_mat[-idx, i] - 1
      vote_mat[-idx, j] <- vote_mat[-idx, j] - 1
      classifier_idx <- classifier_idx + 1

    }
  }
  decf_idx <- apply(vote_mat, 1, which.max)
  decf <- matrix(0, xn)
  for (i in 1:object$class_num) {
    decf[decf_idx == i] <- object$class_set[i]
  }
  return(decf)
}
