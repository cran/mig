## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.align="center", 
  fig.width = 5, 
  fig.height = 5
)

## -----------------------------------------------------------------------------
library(mig)
set.seed(0)

## -----------------------------------------------------------------------------
# library(statmod)
# penult <- mev::smith.penult("invgauss", method = "pot",
#                   qu = seq(0.9, 0.9995, by = 0.0005),
#                mean = 1, shape = 1 / 1)
# plot(penult$qu,
#      penult$shape,
#      xlab = "quantile level",
#      xlim = c(min(penult$qu), 1),
#      ylim = c(0, 0.25),
#      yaxs = "i",
#      xaxs = "i",
#      pch = 20,
#      ylab = "penultimate shape",
#      type = 'l',
#      bty = "l")

## -----------------------------------------------------------------------------
# Create projection matrix onto the orthogonal complement of beta
d <- 5L # dimension of vector
n <- 1e4L # number of simulations
beta <- rexp(d)
xi <- rexp(d)
Omega <- matrix(0.5, d, d) + diag(d)
# Project onto orthogonal complement of vector
Mbeta <- (diag(d) - tcrossprod(beta)/(sum(beta^2)))
# Matrix is rank-deficient: compute eigendecomposition 
# Shed matrix to remove the eigenvector corresponding to the 0 eigenvalue
Q2 <- t(eigen(Mbeta, symmetric = TRUE)$vectors[,-d])
# Check Q2 satisfies the conditions
all.equal(rep(0, d-1), c(Q2 %*% beta)) # check orthogonality
all.equal(diag(d-1), tcrossprod(Q2)) # check basis is orthonormal
Qmat <- rbind(beta, Q2)
covmat <- solve(Q2 %*% solve(Omega) %*% t(Q2))

# Compute mean and variance for Z1
mu <- sum(beta*xi)
omega <- sum(beta * c(Omega %*% beta))
Z1 <- rmig(n = n, xi = mu, Omega = omega) # uses statmod, with mean = mu and shape mu^2/omega
# Generate Gaussian vectors in two-steps (vectorized operations)
Z2 <- sweep(TruncatedNormal::rtmvnorm(n = n, mu = rep(0, d-1), sigma = covmat), 1, sqrt(Z1), "*")
Z2 <- sweep(Z2, 2, c(Q2 %*% xi), "+") + tcrossprod(Z1 - mu, c(Q2 %*% c(Omega %*% beta)/omega))
# Compute inverse of Q matrix (it is actually invertible)
samp <- t(solve(Qmat) %*% t(cbind(Z1, Z2)))
# Check properties
mle <- mig::fit_mig(samp, beta = beta)
max(abs(mle$xi - xi))
norm(mle$Omega - Omega, type = "f")
max(abs(1 - mle$Omega/Omega))

