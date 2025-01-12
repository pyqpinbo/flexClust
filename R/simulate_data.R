#' Simulate Data for Dynamic Varying Coefficient Additive Model (VCAM)
#'
#' @description Simulates data for a dynamic varying coefficient additive model.
#'
#' @param n Integer. Number of units.
#' @param T Integer. Number of time points.
#' @param x_mat Matrix. Covariate values for simulation.
#'
#' @return A list containing simulated data components:
#' \item{ymat}{Matrix of response variable values.}
#' \item{ylag}{Matrix of lagged response variable values.}
#' \item{umat}{Matrix of U values.}
#' \item{xmat}{Matrix of covariate values.}
#' \item{alpha0}{Matrix of intercept terms.}
#' \item{alp1}{Matrix of group-specific varying coefficients (part 1).}
#' \item{alp2}{Matrix of group-specific varying coefficients (part 2).}
#' \item{beta1}{Matrix of additive function beta1 values.}
#' \item{beta2}{Matrix of additive function beta2 values.}
#'
#' @examples
#' x_mat <- matrix(rnorm(100), nrow = 10, ncol = 10)
#' result <- simulate_data(10, 5, x_mat)
#' str(result)
#'
#' @export

simulate_data <- function(n, T, x_mat) {
  # Generate the U matrix
  u <- seq(0, 1, length.out = T + 1)[-1]
  U <- replicate(n, u) 
  
  # Generate error terms
  eta <- rnorm(n, mean = 0, sd = 0.1)
  eta_mat <- matrix(eta, nrow = T, ncol = n, byrow = TRUE) 
  eta_mat <- rbind(eta_mat, eta_mat)
  
  
  eps <- rnorm(n * T, mean = 0, sd = 0.1)
  eps_mat <- matrix(eps, nrow = T, ncol = n)
  eps_mat <- rbind(eps_mat, eps_mat)
  
  # Define varying-coefficient functions
  a0 <- function(u) 1.5 * cos(2 * pi * u)
  a2 <- function(u) 2 * sin(1.5 * pi * u) - 1.2 * (u - 0.5) * (1 - u) + 1
  a11 <- function(u) 1.3 * u * sin(2 * pi * u) + 1
  a12 <- function(u) 1.3 * u * cos(2 * pi * u) - 1
  
  # Compute group-specific varying-coefficient functions
  alpha0 <- rbind(a0(U), a0(U))
  alpha11 <- a11(U[, 1:(n / 2)])
  alpha12 <- a12(U[, (n / 2 + 1):n])
  alpha2 <- a2(U)
  
  # Normalize alpha values
  norm11 <- matrix(colMeans(alpha11), nrow = T, ncol = n / 2, byrow = TRUE)
  norm12 <- matrix(colMeans(alpha12), nrow = T, ncol = n / 2, byrow = TRUE)
  norm2 <- matrix(colMeans(alpha2), nrow = T, ncol = n, byrow = TRUE)
  
  alp11 <- alpha11 / norm11
  alp12 <- alpha12 / norm12
  alp1 <- cbind(alp11, alp12)
  alp2 <- alpha2 / norm2
  alp1 <- rbind(alp1, alp1)
  alp2 <- rbind(alp2, alp2)
  
  # Define additive functions
  b1 <- function(x) -0.8 * (1 - x^2) / (1 + x^2)
  b21 <- function(x) 2 * cos(pi * x / 2) + 1.8 * sin(pi * x / 3)
  b22 <- function(x) 1.5 * sin(pi * x / 4) - 1.2 * cos(pi * x / 3)
  
  # Compute beta values
  beta21 <- b21(x_mat[, 1:(n / 2)])
  beta22 <- b22(x_mat[, (n / 2 + 1):n])
  
  # Center beta values
  beta21 <- beta21 - matrix(colMeans(beta21), nrow = T, ncol = n / 2, byrow = TRUE)
  beta22 <- beta22 - matrix(colMeans(beta22), nrow = T, ncol = n / 2, byrow = TRUE)
  beta2 <- cbind(beta21, beta22)
  beta2 <- rbind(beta2, beta2)
  
  # Generate response data
  Y <- matrix(0, nrow = 2 * T, ncol = n)
  Y[1, ] <- 0.1 * rnorm(n)
  
  for (t in 2:(2 * T)) {
    lag <- Y[t - 1, ]
    Y[t, ] <- alpha0[t, ] + alp1[t, ] * b1(lag) + alp2[t, ] * beta2[t, ] + eta_mat[t, ] + eps_mat[t, ]
  }
  
  # Extract relevant data
  Y <- Y[-(1:T), ]
  ymat <- Y[2:T, ]
  ylag <- Y[1:(T - 1), ]
  umat <- U[2:T, ]
  xmat <- x_mat[2:T, ]
  alpha0 <- alpha0[2:T, ]
  alp1 <- alp1[2:T, ]
  alp2 <- alp2[2:T, ]
  beta2 <- beta2[2:T, ]
  beta1 <- -0.8 * (1 - ylag^2) / (1 + ylag^2)
  #beta1 <- beta1 - matrix(colMeans(beta1), nrow = nrow(beta1), ncol = ncol(beta1), byrow = TRUE)
  
  # Return results as a list
  return(list(
    ymat = ymat,
    ylag = ylag,
    umat = umat,
    xmat = xmat,
    alpha0 = alpha0,
    alp1 = alp1,
    alp2 = alp2,
    beta1 = beta1,
    beta2 = beta2
  ))
}

#' Generate Time-Varying Autoregressive covariate
#'
#' @description Generates a time-varying autoregressive series and discards the first half of the generated series.
#'
#' @param n Integer. Number of values to generate after discarding.
#' @param rho Numeric. Autoregressive coefficient.
#' @param sigma Numeric. Standard deviation of the noise term.
#'
#' @return A numeric vector containing the generated time series.
#'
#' @examples
#' ts <- autoreg(100, 0.8, 0.1)
#' plot(ts, type = "l")
#'
#' @export

autoreg <- function(n, rho, sigma) {
  # Generate tvNAR(1) 2*n data, and discard the first n
  # Initial value ts0 = 0
  ts <- numeric(2 * n)
  tim <- seq(0, 1, length.out = 2 * n + 1)
  t <- tim[-1] # Remove the first value
  
  # Set the initial value
  ts[1] <- sigma * rnorm(1, 0, 1)
  
  # Generate the autoregressive series
  for (i in 2:(2 * n)) {
    ts[i] <- rho * ts[i - 1] + sigma * rnorm(1, 0, 1)
  }
  
  # Discard the first n values
  ts <- ts[-seq_len(n)]
  
  return(ts)
}

