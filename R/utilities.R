# Title: utilities.R
# Created by: Caleb Leedy
# Created on: September 07, 2023
# Purpose: This file contains code that is common and used in 
# multiple scripts but is still very general.

# TODO: Test me!
#' Thanks Hengfang Wang for this matrix inversion function!
hw_inv = function(X, eps = 1e-12) {
  eig.X = eigen(X, symmetric = TRUE)
  P = eig.X[[2]] 
  lambda = eig.X[[1]] 
  ind = lambda > eps
  lambda[ind] = 1 / lambda[ind] 
  lambda[!ind] = 0
  ans = P %*% diag(lambda) %*% t(P)
  return(ans)
}

# TODO: Test me!
expit <- function(x) {1 / (1 + exp(-x))}
