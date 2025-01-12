# flexClust

A flexible and parsimonious modeling strategy for clustered data analysis.

# Installation

To install the `flexClust` package from GitHub, follow these steps:

```R
# Install the devtools package if not already installed
install.packages("devtools")

# Load the devtools library
library(devtools)

# Install the flexClust package from GitHub
install_github("pyqpinbo/flexClust")
```

# References

Tao Huang, Youquan Pei, Jinhong You, Wenyang Zhang.
A flexible and parsimonious modeling strategy for clustered data analysis.

# Usage

```R
# Load the flexClust library
library(flexClust)

# Define the number of clusters and time steps
n <- 20
T <- 100

# Generate autoregressive data
x <- autoreg(n * T, 0.6, 0.5)
x_mat <- matrix(x, nrow = T, ncol = n)

# Simulate data using the generated matrix
data <- simulate_data(n, T, x_mat)

# Extract components from the simulated data
ymat <- data$ymat
ylag <- data$ylag
umat <- data$umat
xmat <- data$xmat

# Initialize the estimation function array
estfun <- array(0, dim = c(T - 1, 5, n))

# Loop through each individual for estimation
for (i in 1:n) {
  # Prepare data for the individual
  xMat <- cbind(ylag[, i], xmat[, i])
  u <- umat[, i]
  y <- ymat[, i]
  
  # Determine baseline points
  ord_idx <- order(u)
  t0 <- u[floor(T / 2)]
  x01 <- sort(xMat[, 1])[floor(T / 2)]
  x02 <- sort(xMat[, 2])[floor(T / 2)]
  indt <- floor(T / 2)
  Sam <- cbind(u, xMat)
  
  # Choose optimal knots for varying-coefficient additive model
  knot1 <- lsknot(Sam, y, x01, x02, indt, ord = 3, delta = 0.01)
  
  # Perform three-step spline estimation
  estfun[, , i] <- vcam(Sam, y, c(x01, x02), indt, ord = 3, knot1, delta = 0.01)
}

# Extract coefficients from the estimation function
alpha0 <- estfun[, 1, ]
alpha1 <- estfun[, 2, ]
alpha2 <- estfun[, 3, ]
# Set up grid parameters
ngrid <- 100
x0 <- seq(-1, 1, length.out = ngrid)
beta1 <- matrix(0, nrow = ngrid, ncol = n)
beta2 <- matrix(0, nrow = ngrid, ncol = n)
# Estimate beta1 and beta2 for each individual
for (i in 1:n) {
  # Interpolate for beta1
  y <- sort(ylag[, i])
  id1 <- order(ylag[, i])
  res4 <- estfun[, 4, i]
  f1 <- res4[id1]
  beta1[, i] <- spline(y, f1, xout = x0)$y
  # Interpolate for beta2
  x <- sort(xmat[, i])
  id2 <- order(xmat[, i])
  res5 <- estfun[, 5, i]
  f2 <- res5[id2]
  beta2[, i] <- spline(x, f2, xout = x0)$y
}
```
# Structure identification for varying coefficient functions
```R
# Combine alpha1 and alpha2 into a single matrix
alpha <- cbind(alpha1, alpha2)

# Convert each column of alpha into a time series list
ts_list <- lapply(seq_len(ncol(alpha)), function(i) alpha[, i])

# Compute the distance matrix
distance_matrix <- dist.parallel(tsList = ts_list, distFunc = tsdiss.euclidean, cores = 1)

# Normalize the distance matrix
distance_matrix <- dist.normalize(distance_matrix)

# Plot the distance matrix as a heatmap
heatmap(distance_matrix, main = "Distance Matrix")

# Build the network using a specified percentile value
net <- network_build(distance_matrix = distance_matrix, percentile_value = 0.4)

# Perform community detection
res <- network_clust(net)

# Plot the detected communities
plot(res$community, net, layout = res$layout, main = "Community Detection")

# Extract community membership
membership <- res$community$membership

# Define true community structure
true <- list(1:(n/2), (n/2 + 1):n, (n + 1):(2 * n))

# Group indices based on detected community membership
groups <- unique(membership)
grouped_list <- lapply(groups, function(g) which(membership == g))

# Compute the normalized mutual information (NMI) score
nmi <- NMI(true, grouped_list, 2 * n)

# Print the NMI score
cat("Normalized Mutual Information (NMI):", nmi, "\n")

# Compute the mean estimates for each identified group
alpha1_est <- rowMeans(alpha[, which(membership == 1), drop = FALSE])
alpha2_est <- rowMeans(alpha[, which(membership == 2), drop = FALSE])
alpha3_est <- rowMeans(alpha[, which(membership == 3), drop = FALSE])

# True values
alpha1_true <- data$alp1[, 1]
alpha2_true <- data$alp1[, 20]
alpha3_true <- data$alp2[, 1]

# Plot the true values with solid lines
plot(alpha1_true, type = "l", col = "blue", lwd = 2, 
     xlab = "Time", ylab = "Alpha Value", 
     main = "Estimates and True Values of Varying Coefficient Functions", 
     ylim = range(c(alpha1_est, alpha2_est, alpha3_est, alpha1_true, alpha2_true, alpha3_true)))

lines(alpha2_true, col = "red", lwd = 2)
lines(alpha3_true, col = "green", lwd = 2)

# Add the estimates with dotted lines
lines(alpha1_est, col = "blue", lwd = 2, lty = 3)
lines(alpha2_est, col = "red", lwd = 2, lty = 3)
lines(alpha3_est, col = "green", lwd = 2, lty = 3)

# Add a legend to distinguish estimates and true values
legend("topright", 
       legend = c("True: Group 1", "True: Group 2", "True: Group 3",
                  "Estimate: Group 1", "Estimate: Group 2", "Estimate: Group 3"), 
       col = c("blue", "red", "green", "blue", "red", "green"), 
       lwd = 2, lty = c(1, 1, 1, 3, 3, 3), bty = "n")
```
# Structure identification for transformation functions
```R
# Combine beta1 and beta2 into one matrix
beta <- cbind(beta1, beta2)

# Convert each column of beta into a time series list
ts_list2 <- lapply(seq_len(ncol(beta)), function(i) beta[, i])

# Compute the distance matrix
distance_matrix2 <- dist.parallel(tsList = ts_list2, distFunc = tsdiss.euclidean, cores = 1)

# Normalize the distance matrix
distance_matrix2 <- dist.normalize(distance_matrix2)

# Plot the distance matrix as a heatmap
heatmap(distance_matrix2, main = "Distance Matrix")
# Build the network using a specified percentile value
net2 <- network_build(distance_matrix = distance_matrix2, percentile_value = 0.3)

# Perform community detection
communities2 <- network_clust(net2)

# Plot the detected communities
plot(communities2$community, net2, layout = communities2$layout, 
     main = "Community Detection", vertex.color = membership(communities2$community) + 1)

# Extract community membership
membership2 <- communities2$community$membership
# Define true community structure
true <- list(1:n, (n + 1):(3/2 * n), (3/2 * n + 1):(2 * n))

# Group indices based on detected community membership
groups <- unique(membership2)
grouped_list <- lapply(groups, function(g) which(membership2 == g))

# Compute the normalized mutual information (NMI) score
nmi <- NMI(true, grouped_list, 2 * n)

# Print the NMI score
cat("Normalized Mutual Information (NMI):", nmi, "\n")
# Compute the mean estimates for each identified group
beta1_est <- rowMeans(beta[, which(membership2 == 1), drop = FALSE])
beta2_est <- rowMeans(beta[, which(membership2 == 2), drop = FALSE])
beta3_est <- rowMeans(beta[, which(membership2 == 3), drop = FALSE])

# Define the true values
beta1_true <- -0.8 * (1 - y0^2) / (1 + y0^2)
beta2_true <- 2 * cos(pi * x0 / 2) + 1.8 * sin(pi * x0 / 3)
beta2_true <- beta2_true - mean(beta2_true)
beta3_true <- 1.5 * sin(pi * x0 / 4) - 1.2 * cos(pi * x0 / 3)
beta3_true <- beta3_true - mean(beta3_true)

# Plot the true values with solid lines
plot(x0, beta1_true, type = "l", col = "blue", lwd = 2, 
     xlab = "Covariate", ylab = "Transformation Function", 
     main = "Estimates and True Values of Transformation Functions", 
     ylim = range(c(beta1_est, beta2_est, beta3_est, beta1_true, beta2_true, beta3_true)))

lines(x0, beta2_true, col = "red", lwd = 2)
lines(x0, beta3_true, col = "green", lwd = 2)

# Add the estimates with dotted lines
lines(x0, beta1_est, col = "blue", lwd = 2, lty = 3)
lines(x0, beta2_est, col = "red", lwd = 2, lty = 3)
lines(x0, beta3_est, col = "green", lwd = 2, lty = 3)

# Add a legend to distinguish estimates and true values
legend("bottomright", 
       legend = c("True: Group 1", "True: Group 2", "True: Group 3",
                  "Estimate: Group 1", "Estimate: Group 2", "Estimate: Group 3"), 
       col = c("blue", "red", "green", "blue", "red", "green"), 
       lwd = 2, lty = c(1, 1, 1, 3, 3, 3), bty = "n")

