#' Select Optimal Knots for Varying Coefficient Additive Model
#'
#' @description Selects optimal knots for the initial and varying-coefficient
#' estimation steps in a spline-based regression model.
#'
#' @param Sam A data frame of predictors.
#' @param y A vector of responses.
#' @param x01 A vector of prediction points for the first covariate.
#' @param x02 A vector of prediction points for the second covariate.
#' @param indt Index of a reference time point for normalizing estimates.
#' @param ord Integer. Order of the B-spline basis.
#' @param delta Regularization parameter for ridge regression.
#'
#' @return A vector of selected knot values for the model.
#' @import Matrix
#' @examples
#'
#' @export

lsknot <- function(Sam, y, x01, x02, indt, ord, delta) {
  L <- nrow(Sam)
  
  # Potential knots
  k1seq <- seq(ceiling(0.5 * L^(1/6)), ceiling(2 * L^(1/6)))
  k2seq <- seq(ceiling(0.5 * L^(1/5)), ceiling(2 * L^(1/5)))
  
  # Step I: Choose knots for initial estimation
  lsbic1 <- matrix(0, nrow = length(k1seq), ncol = length(k1seq))
  for (i in seq_along(k1seq)) {
    for (j in seq_along(k1seq)) {
      ini_res <- lsestIni(Sam, y, x01, x02, indt, k1seq[i], k1seq[j], ord, delta)
      lsbic1[i, j] <- log(ini_res$rss2) + log(L) * ini_res$paranum / (2 * L)
    }
  }
  
  opt1 <- which(lsbic1 == min(lsbic1), arr.ind = TRUE)
  lshC <- k1seq[opt1[1]]
  lshA <- k1seq[opt1[2]]
  
  lsinialp <- lsestIni(Sam, y, x01, x02, indt, lshC, lshA, ord, delta)$inialp
  
  # Step II: Choose knots for varying-coefficient estimation
  lsbic2 <- matrix(0, nrow = length(k2seq), ncol = length(k2seq))
  for (i in seq_along(k2seq)) {
    for (j in seq_along(k2seq)) {
      step2_res <- lsest2(lsinialp, Sam, y, k2seq[i], k2seq[j], ord, delta)
      lsbic2[i, j] <- log(step2_res$rss2) + log(L) * step2_res$paranum / (2 * L)
    }
  }
  
  opt2 <- which(lsbic2 == min(lsbic2), arr.ind = TRUE)
  lskC <- k2seq[opt2[1]]
  lskA <- k2seq[opt2[2]]
  
  return(c(lshC, lshA, lskC, lskA))
}


#' Initial Estimation of Varying Coefficient Model
#'
#' @description Computes initial estimates for varying-coefficient functions
#' using a spline-based regression approach.
#'
#' @param Sam A data frame of predictors.
#' @param y A vector of responses.
#' @param x01 Prediction points for the first covariate.
#' @param x02 Prediction points for the second covariate.
#' @param indt Index of a reference time point for normalization.
#' @param hC Number of knots for the time component.
#' @param hA Number of knots for covariate components.
#' @param ord Integer. Order of the B-spline basis.
#' @param delta Regularization parameter for ridge regression.
#'
#' @return A list containing initial coefficient estimates and model metrics:
#' \item{inialp}{Matrix of initial estimates.}
#' \item{rss2}{Residual sum of squares.}
#' \item{paranum}{Number of parameters in the model.}
#'
#' @examples
#' Sam <- data.frame(matrix(rnorm(300), ncol = 3))
#' y <- rnorm(100)
#' res <- lsestIni(Sam, y, x01 = rnorm(100), x02 = rnorm(100), indt = 1, hC = 5, hA = 3, ord = 3, delta = 0.01)
#' str(res)
#'
#' @export

lsestIni <- function(Sam, y, x01, x02, indt, hC, hA, ord, delta) {
  L <- nrow(Sam)
  
  # B-spline matrix for rescaled time
  bt1 <- rspline(Sam[, 1], Sam[, 1], ord, hC)$Bt
  
  # B-spline matrix for x1 (centered)
  bx1 <- rspline(Sam[, 2], Sam[, 2], ord, hA)$Bt
  c1 <- colMeans(bx1)
  bx1 <- sweep(bx1, 2, c1, "-")[, -1]
  
  # B-spline matrix for x2 (centered)
  bx2 <- rspline(Sam[, 3], Sam[, 3], ord, hA)$Bt
  c2 <- colMeans(bx2)
  bx2 <- sweep(bx2, 2, c2, "-")[, -1]
  
  
  # Tensor product of B-spline
  tensb1 <- kronecker(bx1, matrix(1, nrow = 1, ncol = ncol(bt1))) * matrix(rep(bt1, ncol(bx1)), nrow = nrow(bt1), byrow = FALSE)
  tensb2 <- kronecker(bx2, matrix(1, nrow = 1, ncol = ncol(bt1))) * matrix(rep(bt1, ncol(bx1)), nrow = nrow(bt1), byrow = FALSE)
  xMat <- cbind(bt1, tensb1, tensb2)
  
  # Ridge regression
  coeff1 <- solve(t(xMat) %*% xMat + delta * diag(ncol(xMat))) %*% t(xMat) %*% y
  l1 <- ncol(tensb1)
  l0 <- ncol(bt1)
  
  # Estimate component functions
  library(Matrix)
  blk_diag <- bdiag(bt1, tensb1, tensb2)
  mfun <- (blk_diag %*% coeff1)@x
  est1 <- mfun[1:L]
  est2 <- mfun[(L + 1):(2 * L)]
  est3 <- mfun[(2 * L + 1):(3 * L)]
  
  # RSS and parameter count
  rss2 <- mean((y - est1 - est2 - est3)^2)
  paranum <- ncol(xMat)
  
  # Step II: Initial estimation of varying-coefficient functions
  bx01 <- rspline(Sam[, 2], x01, ord, hA)$Bt
  bx01 <- sweep(bx01, 2, c1, "-")[, -1]
  tem1 <- kronecker(t(bx01), bt1) %*% coeff1[(l0+1):(l0+l1)]
  xsi1 <- tem1 / tem1[indt]
  xsi1 <- xsi1 / sqrt(mean(xsi1^2))
  
  
  bx02 <- rspline(Sam[, 3], x02, ord, hA)$Bt
  bx02 <- sweep(bx02, 2, c2, "-")[, -1]
  tem2 <- kronecker(t(bx02), bt1) %*% coeff1[(l0 + l1 + 1):(l0+2*l1)]
  xsi2 <- tem2 / tem2[indt]
  xsi2 <- xsi2 / sqrt(mean(xsi2^2))
  
  inialp <- cbind(est1, xsi1, xsi2)
  
  return(list(inialp = inialp, rss2 = rss2, paranum = paranum))
}


#' Varying Coefficient Estimation
#'
#' @description Performs the second-step estimation of varying-coefficient
#' functions using spline regression.
#'
#' @param inialp Initial estimates from the first step.
#' @param Sam A data frame of predictors.
#' @param y A vector of responses.
#' @param kC Number of knots for the time component.
#' @param kA Number of knots for covariate components.
#' @param ord Integer. Order of the B-spline basis.
#' @param delta Regularization parameter for ridge regression.
#'
#' @return A list containing estimated functions and model metrics:
#' \item{estfun}{Matrix of estimated functions.}
#' \item{rss2}{Residual sum of squares.}
#' \item{paranum}{Number of parameters in the model.}
#'
#' @examples
#' Sam <- data.frame(matrix(rnorm(300), ncol = 3))
#' y <- rnorm(100)
#' inialp <- matrix(rnorm(300), ncol = 3)
#' res <- lsest2(inialp, Sam, y, kC = 5, kA = 3, ord = 3, delta = 0.01)
#' str(res)
#'
#' @export


lsest2 <- function(inialp, Sam, y, kC, kA, ord, delta) {
  L <- nrow(Sam)
  
  # B-spline matrix for rescaled time
  bt <- rspline(Sam[, 1], Sam[, 1], ord, kC)$Bt
  
  # B-spline matrix for x1 (centered)
  bx1 <- rspline(Sam[, 2], Sam[, 2], ord, kA)$Bt
  bx1 <- sweep(bx1, 2, colMeans(bx1), "-")[, -1]
  
  # B-spline matrix for x2 (centered)
  bx2 <- rspline(Sam[, 3], Sam[, 3], ord, kA)$Bt
  bx2 <- sweep(bx2, 2, colMeans(bx2), "-")[, -1]
  
  # Step II: Additive component functions
  newY <- y - inialp[, 1]
  B21 <- bx1 * inialp[, 2]
  B22 <- bx2 * inialp[, 3]
  B2 <- cbind(B21, B22)
  
  coeff2 <- solve(t(B2) %*% B2 + delta * diag(ncol(B2))) %*% t(B2) %*% newY
  estbeta1 <- bx1 %*% coeff2[1:ncol(B21)]
  estbeta1 <- estbeta1 - mean(estbeta1)
  estbeta2 <- bx2 %*% coeff2[(ncol(B21) + 1):ncol(B2)]
  estbeta2 <- estbeta2 - mean(estbeta2)
  
  # Step III: Varying-coefficient functions
  B31 <- sweep(bt, 1, estbeta1, `*`)
  B32 <- sweep(bt, 1, estbeta2, `*`)
  B3 <- cbind(bt, B31, B32)
  
  coeff3 <- solve(t(B3) %*% B3 + delta * diag(ncol(B3))) %*% t(B3) %*% y
  alpMat <- bdiag(bt, bt, bt) %*% coeff3
  estalp0 <- alpMat@x[1:L]
  estalp1 <- alpMat@x[(L + 1):(2 * L)]
  estalp2 <- alpMat@x[(2 * L + 1):(3 * L)]
  estalp1 <- estalp1 / sqrt(mean(estalp1^2))
  estalp2 <- estalp2 / sqrt(mean(estalp2^2))
  
  # Residuals and RSS
  res <- y - estalp0 - estalp1 * estbeta1 - estalp2 * estbeta2
  rss2 <- mean(res^2)
  paranum <- ncol(B2) + ncol(B3)
  
  estfun <- cbind(estalp0, estalp1, estalp2, estbeta1, estbeta2)
  
  return(list(estfun = estfun, rss2 = rss2, paranum = paranum))
}






