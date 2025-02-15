#' Estimation of Varying Coefficient Additive Model for real data(VCAM)
#'
#' @description
#' This function estimates a varying coefficient additive model using a three-step
#' spline-based regression approach. It accounts for both varying coefficients
#' and additive components for clustered data.
#'
#' @param Sam A data frame containing the covariates.
#' @param y A vector containing the response variable.
#' @param x0 A vector of grid points for prediction.
#' @param indt An index representing a reference time point.
#' @param ord An integer specifying the order of the B-spline.
#' @param knot A vector specifying the optimal number of knots for different steps.
#'             (e.g., `c(hC, hA, kC, kA)`).
#' @param delta A small positive value used for ridge regularization.
#'
#' @return A matrix containing the estimated functions for varying coefficients
#' and additive components.
#'
#' @keywords varying coefficient additive model, spline regression
#' @export
#'

dyVCAM <- function(Sam, y, x0, indt, ord, knot, delta) {
    
    # Extract knot values from Mknot
    hC <- knot[1]
    hA <- knot[2]
    kC <- knot[3]
    kA <- knot[4]
    
    L <- nrow(Sam)
    
    # Step I estimation
    # B-spline matrix for rescaled time
    bt1 <- rspline(Sam[,1], Sam[,1], ord, hC)$Bt
    
    # B-spline matrices for covariates (centered)
    bx1 <- rspline(Sam[,2], Sam[,2], ord, hA)$Bt
    bx1 <- scale(bx1, center = colMeans(bx1), scale = FALSE)[, -1]
    
    bx2 <- rspline(Sam[,3], Sam[,3], ord, hA)$Bt
    bx2 <- scale(bx2, center = colMeans(bx2), scale = FALSE)[, -1]
    
    bx3 <- rspline(Sam[,4], Sam[,4], ord, hA)$Bt
    bx3 <- scale(bx3, center = colMeans(bx3), scale = FALSE)[, -1]
    
    # Tensor product of B-spline
    tensb1 <- kronecker(bx1, matrix(1, nrow = 1, ncol = ncol(bt1))) * matrix(rep(bt1, ncol(bx1)), nrow = nrow(bt1), byrow = FALSE)
    tensb2 <- kronecker(bx2, matrix(1, nrow = 1, ncol = ncol(bt1))) * matrix(rep(bt1, ncol(bx1)), nrow = nrow(bt1), byrow = FALSE)
    tensb3 <- kronecker(bx3, matrix(1, nrow = 1, ncol = ncol(bt1))) * matrix(rep(bt1, ncol(bx1)), nrow = nrow(bt1), byrow = FALSE)
    
    # Full design matrix for Step I
    xMat <- cbind(bt1, tensb1, tensb2, tensb3)
    
    # Ridge regression for Step I estimation
    coeff1 <- solve(t(xMat) %*% xMat + delta * diag(ncol(xMat))) %*% t(xMat) %*% y
    l1 <- ncol(tensb1)
    l0 <- ncol(bt1)
    
    # Estimate component functions
    alp0 <- bt1 %*% coeff1[1:l0]
    
    # Initial estimation for varying-coefficient functions
    bx01 <- rspline(Sam[,2], x0[1], ord, hA)$Bt
    bx1 <- rspline(Sam[,2], Sam[,2], ord, hA)$Bt
    bx01 <- scale(bx01, center = colMeans(bx1), scale = FALSE)[, -1]
    tem1 <- kronecker(t(bx01), bt1) %*% coeff1[(l0+1):(l0+l1)]
    xsi1 <- tem1 / tem1[indt]
    xsi1 <- xsi1 / sqrt(mean(xsi1^2))
    
    bx02 <- rspline(Sam[,3], x0[2], ord, hA)$Bt
    bx2 <- rspline(Sam[,3], Sam[,3], ord, hA)$Bt
    bx02 <- scale(bx02, center = colMeans(bx2), scale = FALSE)[, -1]
    tem2 <- kronecker(t(bx02), bt1) %*% coeff1[(l0+l1+1):(l0+2*l1)]
    xsi2 <- tem2 / tem2[indt]
    xsi2 <- xsi2 / sqrt(mean(xsi2^2))
    
    bx03 <- rspline(Sam[,4], x0[3], ord, hA)$Bt
    bx3 <- rspline(Sam[,4], Sam[,4], ord, hA)$Bt
    bx03 <- scale(bx03, center = colMeans(bx3), scale = FALSE)[, -1]
    tem3 <- kronecker(t(bx03), bt1) %*% coeff1[(l0+2*l1+1):length(coeff1)]
    xsi3 <- tem3 / tem3[indt]
    xsi3 <- xsi3 / sqrt(mean(xsi3^2))
    
    # Step II estimation
    # B-spline matrices for centered covariates
    bx21 <- rspline(Sam[,2], Sam[,2], ord, kA)$Bt
    bx21 <- scale(bx21, center = colMeans(bx21), scale = FALSE)[, -1]
    
    bx22 <- rspline(Sam[,3], Sam[,3], ord, kA)$Bt
    bx22 <- scale(bx22, center = colMeans(bx22), scale = FALSE)[, -1]
    
    bx23 <- rspline(Sam[,4], Sam[,4], ord, kA)$Bt
    bx23 <- scale(bx23, center = colMeans(bx23), scale = FALSE)[, -1]
    
    # Adjust response variable
    newY <- y - alp0
    
    # Multiply with varying-coefficient estimates
    B1 <- sweep(bx21, 1, xsi1, "*")
    B2 <- sweep(bx22, 1, xsi2, "*")
    B3 <- sweep(bx23, 1, xsi3, "*")
    B <- cbind(B1, B2, B3)
    
    # Ridge regression for Step II estimation
    coeff2 <- solve(t(B) %*% B + delta * diag(ncol(B))) %*% t(B) %*% newY
    
    # Estimate additive component functions
    estbeta1 <- bx21 %*% coeff2[1:ncol(bx21)] - mean(bx21 %*% coeff2[1:ncol(bx21)])
    estbeta2 <- bx22 %*% coeff2[(ncol(bx21)+1):(2*ncol(bx21))] - mean(bx22 %*% coeff2[(ncol(bx21)+1):(2*ncol(bx21))])
    estbeta3 <- bx23 %*% coeff2[(2*ncol(bx21)+1):ncol(B)] - mean(bx23 %*% coeff2[(2*ncol(bx21)+1):ncol(B)])
    
    # Step III estimation
    bt2 <- rspline(Sam[,1], Sam[,1], ord, kC)$Bt
    
    B21 <- sweep(bt2, 1, estbeta1, "*")
    B22 <- sweep(bt2, 1, estbeta2, "*")
    B23 <- sweep(bt2, 1, estbeta3, "*")
    
    B3 <- cbind(bt2, B21, B22, B23)
    
    # Ridge regression for Step III estimation
    coeff3 <- solve(t(B3) %*% B3 + delta * diag(ncol(B3))) %*% t(B3) %*% y
    
    # Extract final estimates for varying-coefficient functions
    estalp0 <- bt2 %*% coeff3[1:ncol(bt2)]
    estalp1 <- bt2 %*% coeff3[(ncol(bt2)+1):(ncol(bt2)+ncol(B21))]
    estalp2 <- bt2 %*% coeff3[(ncol(bt2)+ncol(B21)+1):(2*ncol(B21)+ncol(bt2))]
    estalp3 <- bt2 %*% coeff3[(2*ncol(B21)+ncol(bt2)+1):ncol(B3)]
    
    # Normalize the estimates
    estalp1 <- estalp1 / sqrt(mean(estalp1^2))
    estalp2 <- estalp2 / sqrt(mean(estalp2^2))
    estalp3 <- estalp3 / sqrt(mean(estalp3^2))
    
    # Return the estimated functions
    estfun <- cbind(estalp0, estalp1, estalp2, estalp3, estbeta1, estbeta2, estbeta3)
    
    return(estfun)
  }