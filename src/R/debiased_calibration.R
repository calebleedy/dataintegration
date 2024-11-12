# File: debiased_calibration.R
# Author: Caleb Leedy
# Created on: November 2, 2024
# Purpose:
# This file contains the updated debiased calibration code to be used in the
# scripts for my simulation study.
#
# The functions in this script have been developed and ran in the previous
# exploration documents:
# * ../explore/20240411-tpdcsim.qmd
# * ../explore/20240425-nndcsim.qmd
# * ../explore/20240506-msdcsim.qmd
# * ../explore/20240511-msdcsim2.qmd
# * ../explore/20241025-mcbssim.qmd
# * ../explore/20241101-msdcsim_estalpha.qmd
# * ../explore/proto-dctp.R

# *************
# * Libraries *
# *************

library(nleqslv)

# *************
# * Functions *
# *************

#' This  function does the debiased calibration
#'
#' @param zmat - A n x p matrix
#' @param w1i - A vector of length n
#' @param d2i - A vector of length n
#' @param T1 - A vector of length p
#' @param qi - A scaler or a vector of length n
#' @param entropy - One of "EL", "ET"
solveGE <- function(zmat, w1i, d2i, T1, qi, entropy) {

  f <- function(lam, zmat, w1i, d2i, T1, qi, entropy, returnw = FALSE) {

    if (entropy == "EL") {
      w <- (-d2i) / drop(zmat %*% lam)
    } else if (entropy == "ET") {
      w <- d2i * exp(drop(zmat %*% lam))
    } else {
      stop("We only accept entropy of EL or ET.")
    }

    if (returnw) {
      return(w)
    } else {
      return(T1 - drop((w * w1i * qi) %*% zmat))
    }
  }

  j <- function(lam, zmat, w1i, d2i, T1, qi, entropy) {
    end <- ncol(zmat)
    if (entropy == "EL") {
      res <- (-1) * t(zmat) %*% diag((w1i * qi * d2i) / drop(zmat %*% lam)^2) %*% zmat
    } else if (entropy == "ET") {
      res <- (-1) * t(zmat) %*% 
        diag((w1i * qi * d2i) * exp(drop(zmat %*% lam))) %*% zmat
    }

    return(res)
  }

  # zmat still contains g(di) as the last column
  init <- c(rep(0, ncol(zmat) - 1), 1)
  res <- 
  nleqslv(init, f, jac = j, zmat = zmat, w1i = w1i, d2i = d2i,
          T1 = T1, qi = qi, entropy = entropy,
          method = "Newton", control = list(maxit = 1e5, allowSingular = TRUE))

  resw <- f(res$x, zmat, w1i, d2i, T1, qi, entropy, returnw = TRUE)

  if (!(res$termcd %in% c(1))) {
    return(NA)
  }

  return(resw)

}

#' The solveGE_alpha function estimates alpha and the weights using debiased
#' calibration.
#'
#' This is useful for when we do not have design weights for the entire
#' population.
#'
#' @param alpha A vector of length c.
#' @param zmat A matrix of dimension n x c.
#' @param w1i A vector of length n or a scaler value.
#' @param d2i A vector of length n.
#' @param T1 A vector of length c.
#' @param qi A numeric value. Currently, qi must be 1.
#' @param entropy A string indicating the entropy used. Currently on "EL" is
#' accepted.
#' @param retw A boolean value indicating if the function should return the
#' estimated weights (TRUE) or if it should return the minimum value of the
#' objective function (FALSE).
#'
#' @return If retw is TRUE, then the function will return the estimated weights.
#' Otherwise, the function will return the value of $\sum_{i \in A} G(w_i) - N
#' * alpha$. The second case if used for the optimization procedure.
solveGE_alpha <- function(alpha, zmat, w1i, d2i, T1, qi, entropy, retw = FALSE) {

  # The first element of T1 should be the population size N, or an estimate of
  # it.
  N <- T1[1]
  T1_adj <- T1
  T1_adj[length(T1_adj)] <- alpha * N

  # Estimate the weights given that T1 is correct.
  dc_w <- solveGE(zmat, w1i, d2i, T1_adj, qi, entropy)

  if (retw) {
    return(dc_w)
  }

  # Return the minimum value of sum(G(w)) - K(alpha)N for use in an optimization
  # routine.
  Ka_type <- "alpha"
  if (Ka_type == "alpha") {
    # HACK: G(x) = -log(x) for EL only
    G <- function(x) {-log(x)}
    return(sum(G(dc_w)) - N * alpha)
  }
}

