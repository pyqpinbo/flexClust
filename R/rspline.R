#' Generate B-Spline Basis and Derivatives
#'
#' @description
#' Constructs B-spline basis and their first-order derivatives for given input points,
#' using specified order, internal knots, and boundary constraints.
#'
#' @param t A numeric vector of observed data points where the B-spline basis is evaluated.
#' @param t0 A numeric vector of target points where the B-spline basis and derivatives are evaluated.
#' @param m An integer specifying the order of the B-spline basis (degree + 1). For example, `m = 4` represents cubic splines.
#' @param k An integer specifying the number of internal knots to use for the B-spline.
#'
#' @details
#' - The function computes a set of B-spline basis functions and their derivatives evaluated at the points `t0`.
#' - The knot sequence is constructed using quantiles of `t` and boundary knots slightly beyond the range of `t`.
#' - The resulting matrices are reordered to match the original order of `t0`.
#'
#' @return A list containing:
#' \item{Bt}{A matrix representing the B-spline basis evaluated at `t0`.}
#' \item{B1t}{A matrix representing the first-order derivatives of the B-spline basis evaluated at `t0`.}
#'
#' @keywords spline basis, B-spline, spline derivatives
#' @import splines
#' @export
#'
#' @examples
#' t <- seq(0, 1, length.out = 100)
#' t0 <- seq(0, 1, length.out = 10)
#' result <- rspline(t, t0, m = 4, k = 5)
#' Bt <- result$Bt
#' B1t <- result$B1t

rspline <- function(t, t0, m, k) {
  # Generate B-spline basis matrix and its derivative matrix
  # Arguments:
  # t  - sample points (vector)
  # t0 - given points for evaluation (vector or scalar)
  # m  - order of the B-splines (integer)
  # k  - number of internal knots (integer)
  
  # Ensure no replicated points in t
  if (any(duplicated(t))) {
    stop("The input sample points 't' must not contain duplicates.")
  }
  
  # Convert t0 to a vector if it's a scalar
  if (is.numeric(t0) && length(t0) == 1) {
    t0 <- c(t0)
  }
  
  # Size of the evaluation points
  r1 <- length(t0)
  
  # Sort t0 and save original order
  order_idx <- order(t0)
  t0_sorted <- t0[order_idx]
  
  # Generate internal knots
  internal_knots <- quantile(t, probs = seq(1, k) / (k + 1))
  
  # Define boundary knots with a small adjustment to avoid edge effects
  b <- max(t) + 1e-10
  a <- min(t) - 1e-10
  knots <- c(rep(a, m), internal_knots, rep(b, m))
  
  # Generate B-spline basis using splines::splineDesign
  library(splines)
  
  # Generate B-spline collocation matrix
  ttt <- rep(t0_sorted, each = 2)
  collocation_matrix <- splineDesign(knots, ttt, ord = m, derivs = rep(0:1, length.out = length(ttt)))
  
  # Split the collocation matrix into values and derivatives
  A1 <- collocation_matrix[seq(1, nrow(collocation_matrix), by = 2), , drop = FALSE]
  A2 <- collocation_matrix[seq(2, nrow(collocation_matrix), by = 2), , drop = FALSE]
  
  # Reorder rows to match the original order of t0
  Bt <- matrix(0, nrow = r1, ncol = ncol(A1))
  B1t <- matrix(0, nrow = r1, ncol = ncol(A1))
  
  Bt[order_idx, ] <- A1
  B1t[order_idx, ] <- A2
  
  return(list(Bt = Bt, B1t = B1t))
}







