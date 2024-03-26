pin_tksvc_dual_solver <- function(KernelX, y, C1, C2, C3, C4, epsilon,
                                  tau1, tau2, tau3, tau4,
                                  class_set, class_num,
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
      P <- rbind(H, M)

      Hn <- nrow(H)
      Gn <- nrow(G)
      Mn <- nrow(M)

      invHTH <- chol2inv(chol(t(H) %*% H + diag(1e-7, xp)))
      dualH <- N %*% invHTH %*% t(N)
      dualq <- rbind(matrix(1, Gn), matrix(1 - epsilon, Mn))
      duallb <- rbind(matrix(-tau1*C1, Gn), matrix(-tau2*C2, Mn))
      dualub <- rbind(matrix(C1, Gn), matrix(C2, Mn))
      if (tau1 == 0) {
        u0_1 <- matrix(0, Gn)
      } else {
        u0_1 <- matrix((1 - tau1)*C1 / 2, Gn)
      }
      if (tau2 == 0) {
        u0_2 <- matrix(0, Mn)
      } else {
        u0_2 <- matrix((1 - tau2*C2) / 2, Mn)
      }
      u0 <- rbind(u0_1, u0_2)
      x <- clip_dcd_optimizer(dualH, dualq, duallb, dualub, eps, max.steps, u0)$x
      coef1 <- -invHTH %*% t(N) %*% x

      # Hyperplane 2
      P <- rbind(H, M)
      invGTG <- chol2inv(chol(t(G) %*% G + diag(1e-7, xp)))
      dualH <- P %*% invGTG %*% t(P)
      dualq <- rbind(matrix(1, Hn), matrix(1 - epsilon, Mn))
      dualub <- rbind(matrix(C3, Hn), matrix(C4, Mn))
      duallb <- rbind(matrix(-tau3*C3, Hn), matrix(-tau4*C4, Mn))
      if (tau3 == 0) {
        u0_1 <- matrix(0, Hn)
      } else {
        u0_1 <- matrix((1 - tau3)*C3 / 2, Hn)
      }
      if (tau4 == 0) {
        u0_2 <- matrix(0, Mn)
      } else {
        u0_2 <- matrix((1 - tau4)*C4 / 2, Mn)
      }
      u0 <- rbind(u0_1, u0_2)
      x <- clip_dcd_optimizer(dualH, dualq, duallb, dualub, eps, max.steps, u0)$x
      coef2 <- invGTG %*% t(P) %*% x

      coef_pos[, classifier_idx] <- coef1
      coef_neg[, classifier_idx] <- coef2
      classifier_idx <- classifier_idx + 1
    }
  }
  BaseDualPinTKSVCClassifier <- list("coef_pos" = coef_pos,
                                     "coef_neg" = coef_neg)
  return(BaseDualPinTKSVCClassifier)
}

#' Hinge Twin Multi-Class Support Vector Machine
#'
#' \code{hinge_tksvc} is an R implementation of Hinge-TKSVC
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
#' @param epsilon rest class parameter.
#' @param tau1,tau2,tau3,tau4 parameter for pinball loss function.
#' @param eps the precision of the optimization algorithm.
#' @param max.steps the number of iterations to solve the optimization problem.
#' @param solver \code{"dual"} is available.
#' @param fit_intercept if set \code{fit_intercept = TRUE},
#'                      the function will evaluates intercept.
#' @param randx parameter for reduce SVM, default \code{randx = 0.1}.
#' @return return \code{HingeSVMClassifier} object.
#' @export
pin_tksvc <- function(X, y, C1 = 1, C2 = 1, C3 = C1, C4 = C2,
                      kernel = c("linear", "rbf", "poly"),
                      gamma = 1 / ncol(X), degree = 3, coef0 = 0,
                      epsilon = 0.1, tau1 = 1, tau2 = 1, tau3 = tau1, tau4 = tau2,
                      eps = 1e-5, max.steps = 5000,
                      solver = c("dual"), fit_intercept = TRUE,
                      randx = 1) {
  X <- as.matrix(X)
  y <- as.matrix(y)
  class_set <- sort(unique(y))
  class_num <- length(class_set)
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
    solver.res <- pin_tksvc_dual_solver(KernelX, y, C1, C2, C3, C4, epsilon,
                                        tau1, tau2, tau3, tau4,
                                        class_set, class_num, eps, max.steps)
  }
  TKSVCClassifier <- list("X" = X, "y" = y, "class_set" = class_set,
                          "class_num" = class_num,
                          "C1" = C1, "C2" = C2, "C3" = C3, "C4" = C4,
                          "kernel" = kernel,
                          "gamma" = gamma, "degree" = degree, "coef0" = coef0,
                          "epsilon" = epsilon,
                          "solver" = solver, "coef_pos" = solver.res$coef_pos,
                          "coef_neg" = solver.res$coef_neg,
                          "fit_intercept" = fit_intercept)
  class(TKSVCClassifier) <- "TKSVCClassifier"
  return(TKSVCClassifier)
}
