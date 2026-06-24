encode_binary_labels <- function(y) {
  class_set <- sort(unique(y))
  if (length(class_set) != 2) {
    stop("Binary classification requires exactly 2 classes.")
  }
  y <- ifelse(y == class_set[1], 1L, -1L)
  y <- as.numeric(y)
  list(y = as.matrix(y), class_set = class_set)
}

handle_intercept <- function(X, fit_intercept) {
  if (isTRUE(fit_intercept)) { cbind(X, 1) } else { X }
}

resolve_reduce_set <- function(reduce_set, solver) {
  reduce_flag <- !is.null(reduce_set)
  if (reduce_flag && solver == "dual") {
    warning("Dual solver does not support reduce_set; it has been set to NULL.")
    reduce_set <- NULL
    reduce_flag <- FALSE
  }
  list(reduce_set = reduce_set, reduce_flag = reduce_flag)
}

resolve_kernel_matrix <- function(X, kernel, reduce_set,
                                  gamma, degree, coef0) {
  if (kernel == "precomputed") {
    X <- as.matrix(X)
    if (nrow(X) != ncol(X)) {
      stop("Precomputed training kernel matrix must be square.")
    }
    return(X)
  }
  kso <- kernel_select_option_(X, kernel, reduce_set, gamma, degree, coef0)
  kso$KernelX
}

resolve_batch_size <- function(batch_size, n, solver) {
  if (solver != "primal" || is.null(batch_size)) return(NULL)
  batch_size <- as.integer(batch_size)
  max(1L, min(batch_size, n))
}
