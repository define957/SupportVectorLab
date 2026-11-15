#' Least Absolute Deviations (LAD) Loss Function
#'
#' \code{least_absolute_deviations_loss} computes the L1 loss (least
#' absolute deviations), i.e. the sum of absolute residuals
#' \code{sum(abs(r))}.
#'
#' @param r numeric vector of residuals.
#'
#' @return numeric scalar, the total LAD loss.
#' @export
least_absolute_deviations_loss <- function(r) {
  return(sum(abs(r)))
}

#' Least Squares (LS) Loss Function
#'
#' \code{least_squares_loss} computes the squared (L2) loss, i.e. the sum
#' of squared residuals \code{sum(r^2)}.
#'
#' @param r numeric vector of residuals.
#'
#' @return numeric scalar, the total squared loss.
#' @export
least_squares_loss <- function(r) {
  return(sum(r^2))
}
