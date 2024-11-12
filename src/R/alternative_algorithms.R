# File: alternative_algorithms.R
# Created by: Caleb Leedy
# Created on: November 4, 2024
# Purpose:
# This file contains all of the alternative algorithms that are used in the
# paper simulations procedure found in ../final_sims.R.

# *********************************
# * Two-Phase Sampling Algorithms *
# *********************************
# Exploration file: ../explore/proto-dctp.R

#' This function implements the $\pi^*$-estimator for two-phase sampling.
tp_pistar <- function(vars) {

  p1_df <- vars[[1]]
  p2_df <- vars[[2]]
  pop_df <- vars[[3]]
  sampling = "srs-poisson"

  pop_obs <- nrow(pop_df)
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

  return(list(theta = theta, vhat = vhat))
}

#' This function implements the two-phase regression estimator
tp_reg <- function(vars) {

  p1_df <- vars[[1]]
  p2_df <- vars[[2]]
  pop_df <- vars[[3]]
  sampling = "srs-poisson"

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
   
  return(list(theta = theta, vhat = vhat))
}

# *******************************************
# * Nonnested Two-Phase Sampling Algorithms *
# *******************************************
# Exploration file: ../explore/20240425-nndcsim.qmd

nn_ht <- function(vars) {

  p1_df <- vars[[1]]
  p2_df <- vars[[2]]
  pop_df <- vars[[3]]
  sampling = "srs-poisson"

  pop_obs <- nrow(pop_df)
  theta <- sum(p2_df$Y / p2_df$pi2) / pop_obs

  if (sampling == "srs-poisson") {

    v_hat <- sum((1 - p2_df$pi2) / p2_df$pi2^2 * p2_df$Y^2) / pop_obs^2

  } else {

    v_hat <- NA

  }

  return(list(theta = theta, vhat = v_hat))
}

nn_reg <- function(vars) {

  p1_df <- vars[[1]]
  p2_df <- vars[[2]]
  pop_df <- vars[[3]]
  sampling = "srs-poisson"

  pop_obs <- nrow(pop_df)
  v_poisson <- sum((1 - p2_df$pi2) / p2_df$pi2^2 * p2_df$Y^2) / pop_obs^2
  v_srs <- (1 / nrow(p2_df) - 1 / pop_obs) * var(p2_df$Y)
  d_eff <- v_poisson / v_srs
  n2_eff <- nrow(p2_df) / d_eff
  n1_eff <- nrow(p1_df)

  x1_hat <- 
    c(sum(p1_df$X1 / p1_df$pi1) / pop_obs,
      sum(p1_df$X2 / p1_df$pi1) / pop_obs)

  x2_hat <- 
    c(sum(p2_df$X1 / p2_df$pi2) / pop_obs,
      sum(p2_df$X2 / p2_df$pi2) / pop_obs)

  w <- n1_eff / (n1_eff + n2_eff)
  xc_hat <- (n1_eff * x1_hat + n2_eff * x2_hat) / (n1_eff + n2_eff)

  
  # Estimate beta
  mod <- lm(Y ~ X1 + X2, data = p2_df)

  # Regression Estimation
  y_ht <- sum(p2_df$Y / p2_df$pi2) / pop_obs
  theta <- y_ht + sum((xc_hat - x2_hat) * mod$coefficients[2:3])

  # Variance Estimation
  coefs <- mod$coefficients
  eps <- (p2_df$Y - coefs[1] * w - coefs[2] * w * p2_df$X1 - coefs[3] * w * p2_df$X2)
  var_x <- 
  (1 / nrow(p1_df) - 1 / pop_obs) * var(matrix(c(p1_df$X1, p1_df$X2), ncol = 2))

  v1 <- w^2 * matrix(coefs[2:3], nrow = 1) %*% var_x %*% matrix(coefs[2:3], ncol = 1)
  v2 <- sum((1 - p2_df$pi2) / p2_df$pi2^2 * eps^2) / pop_obs^2

  cov_xy <- c(sum((1 - p2_df$pi2) / p2_df$pi2^2 * coefs[2] * 
                 (n2_eff) / (n1_eff + n2_eff) * p2_df$X1 * eps) / pop_obs^2,
              sum((1 - p2_df$pi2) / p2_df$pi2^2 * coefs[3] * 
                 (n2_eff) / (n1_eff + n2_eff) * p2_df$X3 * eps) / pop_obs^2)
  v3 <- sum(cov_xy)

  v_hat <- v1 + v2 + v3

  return(list(theta = theta, vhat = as.numeric(v_hat)))
}


# ***********************************
# * Multisource Sampling Algorithms *
# ***********************************
# Exploration file: ../explore/20240506-msdcsim.qmd, ../explore/20240511-msdcsim2.qmd

ms_ht <- function(vars) {

  p0_df <- vars[[1]]
  p1_df <- vars[[2]]
  p2_df <- vars[[3]]
  pop_df <- vars[[4]]

  N_obs <- nrow(pop_df)
  p0_type = "poisson"

  theta <- sum(p0_df$Y / p0_df$pi0) / N_obs

  if (p0_type == "poisson") {

    v_hat <- sum((1 - p0_df$pi0) / p0_df$pi0^2 * p0_df$Y^2) / N_obs^2

  } else {

    v_hat <- NA

  }

  return(list(theta = theta, vhat = v_hat))
}

ms_reg <- function(vars) {

  p0_df <- vars[[1]]
  p1_df <- vars[[2]]
  p2_df <- vars[[3]]
  pop_df <- vars[[4]]

  N_obs <- nrow(pop_df)
  p0_samp <- "poisson"
  p1_samp <- "srs"

  v_poisson <- sum((1 - p0_df$pi0) / p0_df$pi0^2 * p0_df$Y^2) / N_obs^2
  v_srs <- (1 / nrow(p0_df) - 1 / N_obs) * var(p0_df$Y)
  d_eff <- v_poisson / v_srs
  n0_eff <- nrow(p0_df) / d_eff
  n1_eff <- nrow(p1_df)

  x1_hat <- 
    c(sum(p1_df$X1 / p1_df$pi1) / N_obs,
      sum(p1_df$X2 / p1_df$pi1) / N_obs)

  x0_hat <- 
    c(sum(p0_df$X1 / p0_df$pi0) / N_obs,
      sum(p0_df$X2 / p0_df$pi0) / N_obs)

  w <- n1_eff / (n1_eff + n0_eff)
  xc_hat <- (n1_eff * x1_hat + n0_eff * x0_hat) / (n1_eff + n0_eff)
  
  # Estimate beta
  mod <- lm(Y ~ X1 + X2, data = p0_df)

  # Regression Estimation
  y_ht <- sum(p0_df$Y / p0_df$pi0) / N_obs
  theta <- y_ht + sum((xc_hat - x0_hat) * mod$coefficients[2:3])

  # Variance Estimation
  coefs <- mod$coefficients
  eps <- (p0_df$Y - coefs[1] * w - coefs[2] * w * p0_df$X1 - coefs[3] * w * p0_df$X2)
  var_x <- 
  (1 / nrow(p1_df) - 1 / N_obs) * var(matrix(c(p1_df$X1, p1_df$X2), ncol = 2))

  v1 <- w^2 * matrix(coefs[2:3], nrow = 1) %*% var_x %*% matrix(coefs[2:3], ncol = 1)
  v2 <- sum((1 - p0_df$pi0) / p0_df$pi0^2 * eps^2) / N_obs^2

  cov_xy <- c(sum((1 - p0_df$pi0) / p0_df$pi0^2 * coefs[2] * 
                 (n0_eff) / (n1_eff + n0_eff) * p0_df$X1 * eps) / N_obs^2,
              sum((1 - p0_df$pi0) / p0_df$pi0^2 * coefs[3] * 
                 (n0_eff) / (n1_eff + n0_eff) * p0_df$X3 * eps) / N_obs^2)
  v3 <- sum(cov_xy)

  v_hat <- v1 + v2 + v3

  return(list(theta = theta, vhat = as.numeric(v_hat)))
}


# ********************************************************
# * Multisource with Estimated Alpha Sampling Algorithms *
# ********************************************************
# Exploration file: ../explore/20241101-msdcsim_estalpha.qmd

msa_ht <- function(vars) {

  p0_df <- vars[[1]]
  p1_df <- vars[[2]]
  p2_df <- vars[[3]]
  pop_df <- vars[[4]]

  N_obs <- nrow(pop_df)
  p0_type <- "poisson"

  theta <- sum(p0_df$Y / p0_df$pi0) / N_obs

  if (p0_type == "poisson") {

    v_hat <- sum((1 - p0_df$pi0) / p0_df$pi0^2 * p0_df$Y^2) / N_obs^2

  } else {

    v_hat <- NA

  }

  return(list(theta = theta, vhat = v_hat))
}

msa_reg <- function(vars) {

  p0_df <- vars[[1]]
  p1_df <- vars[[2]]
  p2_df <- vars[[3]]
  pop_df <- vars[[4]]

  N_obs <- nrow(pop_df)
  p0_samp <- "poisson"
  p1_samp <- "srs"

  v_poisson <- sum((1 - p0_df$pi0) / p0_df$pi0^2 * p0_df$Y^2) / N_obs^2
  v_srs <- (1 / nrow(p0_df) - 1 / N_obs) * var(p0_df$Y)
  d_eff <- v_poisson / v_srs
  n0_eff <- nrow(p0_df) / d_eff
  n1_eff <- nrow(p1_df)

  x1_hat <- 
    c(sum(p1_df$X1 / p1_df$pi1) / N_obs,
      sum(p1_df$X2 / p1_df$pi1) / N_obs)

  x0_hat <- 
    c(sum(p0_df$X1 / p0_df$pi0) / N_obs,
      sum(p0_df$X2 / p0_df$pi0) / N_obs)

  w <- n1_eff / (n1_eff + n0_eff)
  xc_hat <- (n1_eff * x1_hat + n0_eff * x0_hat) / (n1_eff + n0_eff)
  
  # Estimate beta
  mod <- lm(Y ~ X1 + X2, data = p0_df)

  # Regression Estimation
  y_ht <- sum(p0_df$Y / p0_df$pi0) / N_obs
  theta <- y_ht + sum((xc_hat - x0_hat) * mod$coefficients[2:3])

  # Variance Estimation
  coefs <- mod$coefficients
  eps <- (p0_df$Y - coefs[1] * w - coefs[2] * w * p0_df$X1 - coefs[3] * w * p0_df$X2)
  var_x <- 
  (1 / nrow(p1_df) - 1 / N_obs) * var(matrix(c(p1_df$X1, p1_df$X2), ncol = 2))

  v1 <- w^2 * matrix(coefs[2:3], nrow = 1) %*% var_x %*% matrix(coefs[2:3], ncol = 1)
  v2 <- sum((1 - p0_df$pi0) / p0_df$pi0^2 * eps^2) / N_obs^2

  cov_xy <- c(sum((1 - p0_df$pi0) / p0_df$pi0^2 * coefs[2] * 
                 (n0_eff) / (n1_eff + n0_eff) * p0_df$X1 * eps) / N_obs^2,
              sum((1 - p0_df$pi0) / p0_df$pi0^2 * coefs[3] * 
                 (n0_eff) / (n1_eff + n0_eff) * p0_df$X3 * eps) / N_obs^2)
  v3 <- sum(cov_xy)

  v_hat <- v1 + v2 + v3

  return(list(theta = theta, vhat = as.numeric(v_hat)))
}


