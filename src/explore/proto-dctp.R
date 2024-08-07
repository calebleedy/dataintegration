# Title: proto-dctp.R
# Created by: Caleb Leedy
# Created on: April 24, 2024
# Purpose: This file contains some of the initial code for an R package to run
# simulations of the debiased calibration algorithm with two-phase sampling.

# *********
# * Goals *
# *********
# 1. [ ] Create consistent functions to generate two-phase samples.
# 2. [ ] Be able to run two-phase sampling debiased calibration (estimation and
# variance) with minimal code repetitions.

# *************
# * Libraries *
# *************

library(rlang)
library(dplyr)
library(nleqslv)

# ******************************
# * Generate Two-Phase Samples *
# ******************************
# The tricky part is that when we use MCMC, we need to generate additional
# samples from a fixed population. But, should the sampling probabilities be
# included in the population data set?
# Answer: Yes, with the update_pop function we can add the probabilities later.

#' This function generates a population data frame with the requisite variables.
#'
#' @param obs - This is the number of observations in the population.
#'
#' @details - This function does not include any sampling probabilities because
#' it does not take into account the sampling type or formula.
gen_pop <- function(obs) {

  x1 <- rnorm(obs, 2, 1)
  x2 <- runif(obs, 0, 4)
  x3 <- rnorm(obs, 0, 1)
  x4 <- runif(obs, 0.1, 0.9)
  z <- rnorm(obs, 0, 1)
  eps <- rnorm(obs)

  y <- 3 * x1 + 2 * x2 + 0.5 * z + eps

  return(tibble(X1 = x1, X2 = x2, X3 = x3, X4 = x4, Z = z, Y = y))
}

#' This function adds sampling probabilites to a population data frame.
#'
#' @details The phase 2 formula will be included for every element. There is no
#' information about if an element is part of the Phase 1 or Phase 2 sample.
update_pop <- function(pop_df, p1_formula, p2_formula) {

  pi1 <- eval(parse_expr(p1_formula))
  pi2 <- eval(parse_expr(p2_formula))

  pop_df %>%
    mutate(pi1 = pi1, pi2 = pi2)
}

#' This function generates Phase 1 and Phase 2 samples from an updated
#' population with selection probabilites.
gen_samps <- function(uppop_df, p1_type, p2_type) {

  if (p1_type == "poisson") {
    del1 <- rbinom(nrow(uppop_df), 1, uppop_df$pi1)
  } else if (p1_type == "srs") {
    ind <- sample(1:nrow(uppop_df), size = uppop_df$pi1[1] * nrow(uppop_df))
    del1 <- as.numeric(1:nrow(uppop_df) %in% ind)
  } else {
    stop("We have only implemented poisson and srs for Phase 1.")
  }

  p1_df <- mutate(uppop_df, del1 = del1) %>% filter(del1 == 1)

  if (p2_type == "poisson") {
    del2 <- rbinom(nrow(p1_df), 1, p1_df$pi2)
  } else if (p2_type == "srs") {
    ind <- sample(1:nrow(p1_df), size = round(p1_df$pi2[1] * nrow(p1_df)))
    del2 <- as.numeric(1:nrow(p1_df) %in% ind)
  } else {
    stop("We have only implemented poisson and srs for Phase 2.")
  }

  p2_df <- mutate(p1_df, del2 = del2) %>% filter(del2 == 1)

  return(list(p1_df, p2_df))
}

# *********************
# * Estimate \bar Y_N *
# *********************

## Alternative Estimators

#' This function implements the $\pi^*$-estimator for two-phase sampling.
pi_star <- function(p1_df, p2_df, pop_obs, sampling = "poisson-poisson") {

  theta <- sum(p2_df$Y / (p2_df$pi1 * p2_df$pi2)) / pop_obs

  if (sampling == "srs-srs") {
    
    vhat <- (1 / nrow(p2_df) - 1 / pop_obs) * var(p2_df$Y)

  } else if (sampling == "poisson-poisson") {

    pis <- p2_df$pi1 * p2_df$pi2
    vhat <- sum((1 / pis^2 - 1 / pis) * p2_df$Y^2) / pop_obs^2

  } else if (sampling == "srs-poisson") {

    # var_y <- var(p2_df$Y)
    # v1 <- (1 / nrow(p1_df) - 1 / pop_obs) * var_y
    # v2 <- sum((1 - p2_df$pi2) / (p2_df$pi1^2 * p2_df$pi2^2) * p2_df$Y^2) / pop_obs^2
    # vhat <- v1 + v2

    tmp <- matrix(NA, nrow = nrow(p2_df), ncol = nrow(p2_df))
    for (i in 1:nrow(p2_df)) {
      for (j in 1:nrow(p2_df)) {
        if (i == j) {
          pi2i <- p2_df$pi2[i] * p2_df$pi1[i]
          tmp[i, j] <- ((pi2i - (pi2i)^2) / pi2i) * p2_df$Y[i]^2 / (pi2i^2)
        } else {
          pi2i <- p2_df$pi2[i] * p2_df$pi1[i]
          pi2j <- p2_df$pi2[j] * p2_df$pi1[j]
          pi2ijp1ij <- 
            p2_df$pi2[i] * p2_df$pi2[j] * p2_df$pi1[i] * 
            (nrow(p1_df) - 1) / (nrow(pop_df) - 1)
          tmp[i, j] <- 
            ((pi2ijp1ij - pi2i * pi2j) / pi2ijp1ij) * 
            p2_df$Y[i] / pi2i * p2_df$Y[j] / pi2j  
        }
      }
    }

    vhat <- sum(as.numeric(tmp)) / nrow(pop_df)^2


  } else {
    vhat <- NA
  }

  return(list(theta = theta, Var = vhat))
}

#' This function implements the two-phase regression estimator
tp_reg <- function(p1_df, p2_df, pop_df, sampling = "poisson-poisson") {

  # Step 1: Get beta
  
  mod <- lm(Y ~ X1 + X2, data = p2_df)

  theta <- 
    (sum(predict(mod, p1_df) / p1_df$pi1) + 
    sum(1 / (p2_df$pi1 * p2_df$pi2) * (p2_df$Y - predict(mod, p2_df)))) / nrow(pop_df)

  p12_df <- 
    left_join(p1_df, p2_df, by = join_by(X1, X2, X3, X4, Y, pi1, pi2, del1)) %>%
    mutate(del2 = ifelse(is.na(del2), 0, 1))
  pred_12 <- predict(mod, p12_df)
  eta <- pred_12 + p12_df$del2 / p12_df$pi2 * (p12_df$Y - pred_12)

  # Step 2: Predict and Get Variance
  if (sampling == "srs-srs") {

    theta <- mean(predict(mod, p1_df))

    vhat <-
      (1 / n1 - 1 / nrow(pop_df)) * (mod$coefficients[2]^2 * var(p1_df$X1) + 
      mod$coefficients[3]^2 * var(p1_df$X2)) +
      (1 / n2 - 1 / n1) * var(p2_df$Y - predict(mod, p2_df))

  } else if (sampling == "poisson-poisson") {

    v1 <- sum((1 - p12_df$pi1) / p12_df$pi1^2 * eta^2)
    v2 <- sum(1 / (p2_df$pi1 * p2_df$pi2) * 
             (1 / p2_df$pi2 - 1) * (p2_df$Y - predict(mod, p2_df))^2)

    vhat <- (v1 + v2) / nrow(pop_df)^2

  } else if (sampling == "srs-poisson") {

    var_eta <- var(eta)
    v1 <- (1 / nrow(p1_df) - 1 / nrow(pop_df)) * var_eta
    v2 <- sum(1 / (p2_df$pi1 * p2_df$pi2) * 
             (1 / p2_df$pi2 - 1) * (p2_df$Y - predict(mod, p2_df))^2) / nrow(pop_df)^2

    vhat <- (v1 + v2)

  } else {

    theta <- NA
    vhat <- NA
  }
   
  return(list(theta = theta, Var = vhat))
}

## Debiased Calibration Estimators

#' This  function does the debiased calibration
#'
#' @param zmat - A n x p matrix
#' @param w1i - A vector of length n
#' @param d2i - A vector of length n
#' @param T1 - A vector of length p
#' @param qi - A scaler or a vector of length n
#' @param entropy - One of "EL", "ET"
solveGE <- function(zmat, w1i, d2i, T1, qi, entropy) {

  f <- function(lam, zmat, w1i, T1, qi, entropy, returnw = FALSE) {

    if (entropy == "EL") {
      w <- -1 / drop(zmat %*% lam)
    } else if (entropy == "ET") {
      w <- exp(drop(zmat %*% lam))
    } else {
      stop("We only accept entropy of EL or ET.")
    }

    if (returnw) {
      return(w)
    } else {
      return(drop((w * w1i * qi) %*% zmat) - T1)
    }
  }

  j <- function(lam, zmat, w1i, T1, qi, entropy) {
    if (entropy == "EL") {
      res <- t(zmat) %*% diag((w1i * qi) / drop(zmat %*% lam)^2) %*% zmat
    } else if (entropy == "ET") {
      res <- t(zmat) %*% diag((w1i * qi) * exp(drop(zmat %*% lam))) %*% zmat
    }

    return(res)
  }

  init <- c(rep(0, ncol(zmat) - 1), 1)
  res <- 
  nleqslv(init, f, jac = j, zmat = zmat, w1i = w1i,
          T1 = T1, qi = qi, entropy = entropy,
          method = "Newton", control = list(maxit = 1e5, allowSingular = TRUE))

  resw <- f(res$x, zmat, w1i, T1, qi, entropy, returnw = TRUE)

  if (!(res$termcd %in% c(1))) {
    return(NA)
  }

  return(resw)

}

#' This function estimates the mean of Y_N from a two-phase sample using
#' debiased calibration.
dc_ybar <- 
  function(p1_df,
           p2_df,
           pop_df,
           qi = 1,
           entropy = "EL",
           sampling = "poisson-poisson",
           estT1 = TRUE,
           linearized = FALSE) {

  d_vec <- 1 / (p2_df$pi2)
  if (entropy == "EL") {
    gdi <- -1 / d_vec
    gdi1 <- -1 / (1 / p1_df$pi2)
    gdip <- -1 / (1 / pop_df$pi2)
    gpinv <- d_vec^2
    gdi12 <- "-1 / d2i"
  } else if (entropy == "ET") {
    gdi <- log(d_vec)
    gdi1 <- log(1 / p1_df$pi2)
    gdip <- log(1 / pop_df$pi2)
    gpinv <- d_vec
    gdi12 <- "log(d2i)"
  } else {
    stop("This type of entropy has not yet been implemented.")
  }

  samp_X <- cbind(rep(1, nrow(p2_df)), as.matrix(select(p2_df, X1, X2)))
  zmat <- cbind(samp_X, gdi)
  w1i <- 1 / p2_df$pi1
  d2i <- 1 / p2_df$pi2

  if (estT1) {
    T1 <- c(nrow(pop_df),
            sum(p1_df$X1 * qi / p1_df$pi1),
            sum(p1_df$X2 * qi / p1_df$pi1),
            sum(gdip * qi)) # We do NOT want to estimate gdi
  } else {
    T1 <- c(nrow(pop_df),
            sum(pop_df$X1 * qi),
            sum(pop_df$X2 * qi),
            sum(gdip * qi))
  }

  dc_w <- solveGE(zmat, w1i, d2i, T1, qi, entropy)

  # Mean Estimation
  theta <- sum(p2_df$Y * dc_w * w1i) / nrow(pop_df)

  # Variance Estimation
  # We ignore qi for now
  p2_df <- mutate(p2_df, gdi = gdi, weight_v = gpinv * w1i)
  mod <- lm(Y ~ X1 + X2 + gdi, weights = weight_v, data = p2_df)
  # z_dc <- model.matrix(~1 + p2_df$X1 + p2_df$X2 + gdi)
  # gam_dc <- 
  #   solve(t(z_dc) %*% diag(w1i * gpinv) %*% z_dc,
  #         t(z_dc) %*% diag(w1i * gpinv) %*% p2_df$Y)
  gam_dc <- matrix(mod$coefficients)
  pred2 <- predict(mod, p2_df)

  p12_df <- 
    left_join(p1_df,
              mutate(p2_df, dc_w = dc_w),
              by = join_by(X1, X2, X3, X4, Y, pi1, pi2, del1)) %>%
    mutate(del2 = ifelse(is.na(del2), 0, 1)) %>%
    mutate(dc_w = ifelse(is.na(dc_w), 0, dc_w)) %>%
    mutate(d2i = 1 / pi2) %>%
    mutate(gdi = eval(parse_expr(gdi12)))

  # z12_dc <- model.matrix(~1 + X1 + X2 + gdi, data = p12_df)
  # pred_12 <- as.numeric(z12_dc %*% gam_dc)
  pred_12 <- predict(mod, p12_df)
  eps <- p12_df$Y - pred_12
  eps2 <- p2_df$Y - pred2

  if (estT1) {

    eta <- pred_12 + p12_df$del2 * p12_df$dc_w * eps

    if (sampling == "poisson-poisson") {

      v1 <- sum((1 - p12_df$pi1) / p12_df$pi1^2 * eta^2)
      v2 <- sum(1 / (p2_df$pi1 * p2_df$pi2^2) * (1 - p2_df$pi2) * (eps2)^2)

      vhat <- (v1 + v2) / nrow(pop_df)^2

    } else if (sampling == "srs-poisson") {

      v1 <- (1 / nrow(p1_df) - 1 / nrow(pop_df)) * var(eta)
      v2 <- 
      sum(1 / (p2_df$pi1 * p2_df$pi2^2) * (1 - p2_df$pi2) * (eps2)^2) / nrow(pop_df)^2

      vhat <- (v1 + v2)
    }

    if (linearized) {
      # eta <- pred_12 + p12_df$del2 / p12_df$pi2 * eps
      theta <- sum((1 / p1_df$pi1) * eta) / nrow(pop_df)
    }
  } else {

    p12_df <- 
      left_join(pop_df,
                mutate(p2_df, dc_w = dc_w),
                by = join_by(X1, X2, X3, X4, Y, pi1, pi2)) %>%
      mutate(del1 = ifelse(is.na(del1), 0, 1)) %>%
      mutate(del2 = ifelse(is.na(del2), 0, 1)) %>%
      mutate(dc_w = ifelse(is.na(dc_w), 0, dc_w)) %>%
      mutate(d2i = 1 / pi2) %>%
      mutate(w1i = 1 / pi1) %>%
      mutate(gdi = eval(parse_expr(gdi12)))

    pred_12 <- predict(mod, p12_df)
    eps <- p12_df$Y - pred_12
    eps2 <- p2_df$Y - pred2

    eta <- pred_12 + p12_df$del2 * p12_df$w1i * p12_df$dc_w * eps

    if (sampling == "poisson-poisson") {

      pis <- p2_df$pi1 * (1 / dc_w)
      vhat <- sum((1 / pis^2 - 1 / pis) * eps2^2) / nrow(pop_df)^2

    } else if (sampling == "srs-poisson") {

      # # This works but why 1 / nrow(p2_df) in v1 instead of 1 / nrow(p1_df)?
      # Because, we essentially do not use the first stage.
       var_eps <- var(eps2)
       v1 <- (1 / nrow(p2_df) - 1 / nrow(pop_df)) * var_eps
       v2 <- 
         sum((1 - p2_df$pi2) / (p2_df$pi1^2 * (1 / dc_w)^2) * eps2^2) / nrow(pop_df)^2
       vhat <- v1 + v2

      # tmp <- matrix(NA, nrow = nrow(p2_df), ncol = nrow(p2_df))
      # for (i in 1:nrow(p2_df)) {
      #   for (j in 1:nrow(p2_df)) {
      #     if (i == j) {
      #       pi2i <- p2_df$pi2[i] * p2_df$pi1[i]
      #       tmp[i, j] <- ((pi2i - (pi2i)^2) / pi2i) * eps2[i]^2 / (pi2i^2)
      #     } else {
      #       pi2i <- p2_df$pi2[i] * p2_df$pi1[i]
      #       pi2j <- p2_df$pi2[j] * p2_df$pi1[j]
      #       pi2ijp1ij <- 
      #         p2_df$pi2[i] * p2_df$pi2[j] * p2_df$pi1[i] * 
      #         (nrow(p1_df) - 1) / (nrow(pop_df) - 1)
      #       tmp[i, j] <- 
      #         ((pi2ijp1ij - pi2i * pi2j) / pi2ijp1ij) * 
      #         eps2[i] / pi2i * eps2[j] / pi2j  
      #     }
      #   }
      # }

      # vhat <- sum(as.numeric(tmp)) / nrow(pop_df)^2

    }

    if (linearized) {
      theta <- sum(eta) / nrow(pop_df)
    }
  }

  return(list(theta = theta, Var = as.numeric(vhat)))
}

