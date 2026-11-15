#' Least Absolute Deviations (LAD) Weight Function
#'
#' \code{lad_weight} computes the per-sample weight used in the IRLS
#' reweighting step for the L1 loss, i.e. \code{1 / |r|}, where \code{r}
#' is the residual vector.
#'
#' @param r numeric vector of residuals.
#' @param reg numeric scalar; residuals with \code{|r| < reg} are assigned
#'   weight \code{1 / reg} to avoid division-by-zero. Default is \code{1e-4}.
#'
#' @return numeric vector of the same length as \code{r} containing the
#'   weights.
#' @export
lad_weight <- function(r, reg = 1e-4) {
  return(1/pmax(abs(r), reg))
}
