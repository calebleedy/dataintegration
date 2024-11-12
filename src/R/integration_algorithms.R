# File: integration_algorithms.R
# Created by: Caleb Leedy
# Created on: November 4, 2024
# Purpose:
# This file contains all of the data integration algorithms that are used in the
# paper simulations procedure found in ../final_sims.R.

# *********************************
# * Two-Phase Sampling Algorithms *
# *********************************
# Exploration file: ../explore/proto-dctp.R

#' This function estimates the mean of Y_N from a two-phase sample using
#' debiased calibration.
tp_est <- function(vars) {
  p1_df <- vars[[1]]
  p2_df <- vars[[2]]
  pop_df <- vars[[3]]
  estT1 <- vars[[4]]

  qi = 1
  entropy = "EL"
  sampling = "srs-poisson"
  linearized = FALSE

  d_vec <- 1 / (p2_df$pi2)
  if (entropy == "EL") {
    gdi <- (-1) / d_vec
    gdi1 <- (-1) / (1 / p1_df$pi2)
    gdip <- (-1) / (1 / pop_df$pi2)
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
    mutate(dc_w = ifelse(is.na(dc_w), 0, dc_w)) %>% # NOTE: We use pi2 now.
    mutate(d2i = 1 / pi2) %>%
    mutate(gdi = eval(rlang::parse_expr(gdi12)))

  # z12_dc <- model.matrix(~1 + X1 + X2 + gdi, data = p12_df)
  # pred_12 <- as.numeric(z12_dc %*% gam_dc)
  pred_12 <- predict(mod, p12_df)
  eps <- p12_df$Y - pred_12
  eps2 <- p2_df$Y - pred2

  if (estT1) {

    # NOTE: Changed dc_w to pi2
    eta <- pred_12 + p12_df$del2 / (p12_df$pi2) * eps

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
      mutate(gdi = eval(rlang::parse_expr(gdi12)))

    pred_12 <- predict(mod, p12_df)
    eps <- p12_df$Y - pred_12
    eps2 <- p2_df$Y - pred2

    # NOTE: We removed w1i
    # NOTE: We changed dc_w to 1 / pi2
    eta <- pred_12 + p12_df$del2 / p12_df$pi2 * eps

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

    }

    if (linearized) {
      theta <- sum(eta) / nrow(pop_df)
    }
  }

  return(list(theta = theta, vhat = as.numeric(vhat)))
}

# *******************************************
# * Nonnested Two-Phase Sampling Algorithms *
# *******************************************
# Exploration file: ../explore/20240425-nndcsim.qmd

nn_est <- 
  function(vars) {
           
  p1_df <- vars[[1]]
  p2_df <- vars[[2]]
  pop_df <- vars[[3]]
  estT1 <- vars[[4]]

  qi = 1
  entropy = "EL"
  sampling = "srs-poisson"

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
  w1i <- 1
  d2i <- 1 / p2_df$pi2

  if (estT1) {
    v_poisson <- sum((1 - p2_df$pi2) / p2_df$pi2^2 * p2_df$Y^2) / nrow(pop_df)^2
    v_srs <- (1 / nrow(p2_df) - 1 / nrow(pop_df)) * var(p2_df$Y)
    d_eff <- v_poisson / v_srs
    n2_eff <- nrow(p2_df) / d_eff
    n1_eff <- nrow(p1_df)

    x1_hat <- 
      c(sum(p1_df$X1 * qi / p1_df$pi1),
        sum(p1_df$X2 * qi / p1_df$pi1))

    x2_hat <- 
      c(sum(p2_df$X1 * qi / p2_df$pi2),
        sum(p2_df$X2 * qi / p2_df$pi2))

    w <- n1_eff / (n1_eff + n2_eff)
    xc_hat <- (n1_eff * x1_hat + n2_eff * x2_hat) / (n1_eff + n2_eff)

    T1 <- c(nrow(pop_df),
            xc_hat[1],
            xc_hat[2],
            sum(gdip * qi))
  } else {
    w <- 1
    T1 <- c(nrow(pop_df),
            sum(pop_df$X1 * qi),
            sum(pop_df$X2 * qi),
            sum(gdip * qi))
  }

  dc_w <- solveGE(zmat, w1i, d2i, T1, qi, entropy)

  # Mean Estimation
  theta <- sum(p2_df$Y * dc_w * w1i) / nrow(pop_df)

  z_dc <- model.matrix(~1 + p2_df$X1 + p2_df$X2 + gdi)
  gam_dc <- 
    solve(t(z_dc) %*% diag(d2i * gpinv * qi) %*% z_dc,
          t(z_dc) %*% diag(d2i * gpinv) %*% p2_df$Y)

  if (sampling == "srs-poisson") {

    if (estT1) {
      # Variance Estimation
      coefs <- as.numeric(gam_dc)
      u_mat <- 
      matrix(c(rep(1, nrow(p2_df)), p2_df$X1 * w, p2_df$X2 * w, gdi * qi), ncol = 4)
      eps <- (p2_df$Y - drop(u_mat %*% coefs))
      var_x <- 
      (1/nrow(p1_df) - 1/nrow(pop_df)) * var(matrix(c(p1_df$X1, p1_df$X2), ncol = 2))

      v1 <- 
      w^2 * matrix(coefs[2:3], nrow = 1) %*% var_x %*% matrix(coefs[2:3], ncol = 1)
      v2 <- sum((1 - p2_df$pi2) / p2_df$pi2^2 * eps^2) / nrow(pop_df)^2

      cov_xy <- c(sum((1 - p2_df$pi2) / p2_df$pi2^2 * coefs[2] * 
                     (n2_eff) / (n1_eff + n2_eff) * p2_df$X1 * eps) / nrow(pop_df)^2,
                  sum((1 - p2_df$pi2) / p2_df$pi2^2 * coefs[3] * 
                     (n2_eff) / (n1_eff + n2_eff) * p2_df$X3 * eps) / nrow(pop_df)^2)
      v3 <- sum(cov_xy)

      v_hat <- v1 + v2 + v3

    } else {
      coefs <- as.numeric(gam_dc)
      u_mat <- 
      matrix(c(rep(1, nrow(p2_df)), p2_df$X1, p2_df$X2, gdi * qi), ncol = 4)
      eps <- (p2_df$Y - drop(u_mat %*% coefs))
      v2 <- sum((1 - p2_df$pi2) / p2_df$pi2^2 * eps^2) / nrow(pop_df)^2

      v_hat <- v2

    }

  } else {
    v_hat <- NA
  }

  return(list(theta = theta, vhat = as.numeric(v_hat)))
}

# ***********************************
# * Multisource Sampling Algorithms *
# ***********************************
# Exploration file: ../explore/20240506-msdcsim.qmd, ../explore/20240511-msdcsim2.qmd

ms_est <- function(vars) {

  p0_df <- vars[[1]]
  p1_df <- vars[[2]]
  p2_df <- vars[[3]]
  pop_df <- vars[[4]]
  estT1 <- vars[[5]]

  qi = 1
  entropy = "EL"
  reg = FALSE

  d_vec <- 1 / (p0_df$pi0)
  if (entropy == "EL") {
    gdi <- -1 / d_vec
    gdip <- -1 / (1 / pop_df$pi0)
    gpinv <- d_vec^2
  } else if (entropy == "ET") {
    gdi <- log(d_vec)
    gdip <- log(1 / pop_df$pi0)
    gpinv <- d_vec
  } else {
    stop("This type of entropy has not yet been implemented.")
  }

  samp_X <- cbind(rep(1, nrow(p0_df)), as.matrix(select(p0_df, X1, X2, X3)))
  zmat <- cbind(samp_X, gdi)
  w1i <- 1
  d0i <- 1 / p0_df$pi0

  if (estT1) {
    x1_0 <- sum(p0_df$X1 / p0_df$pi0) / N_obs
    x2_0 <- sum(p0_df$X2 / p0_df$pi0) / N_obs
    x3_0 <- sum(p0_df$X3 / p0_df$pi0) / N_obs

    x1_1 <- sum(p1_df$X1 / p1_df$pi1) / N_obs
    x3_1 <- sum(p1_df$X3 / p1_df$pi1) / N_obs

    x1_2 <- sum(p2_df$X1 / p2_df$pi2) / N_obs
    x2_2 <- sum(p2_df$X2 / p2_df$pi2) / N_obs

    xhat <- c(x1_0, x2_0, x3_0, x1_1, x3_1, x1_2, x2_2)
    x_mat <- 
    matrix(c(1, 0, 0,
             0, 1, 0,
             0, 0, 1,
             1, 0, 0,
             0, 0, 1,
             1, 0, 0,
             0, 1, 0), ncol = 3, byrow = TRUE)
    vmat <- matrix(0, nrow = 7, ncol = 7)
    vmat[1, 1] <- sum((1 / p0_df$pi0^2 - 1 / p0_df$pi0) * p0_df$X1^2) / nrow(pop_df)^2
    vmat[1, 2] <- 
      sum((1 / p0_df$pi0^2 - 1 / p0_df$pi0) * p0_df$X1 * p0_df$X2) / nrow(pop_df)^2
    vmat[1, 3] <-
      sum((1 / p0_df$pi0^2 - 1 / p0_df$pi0) * p0_df$X1 * p0_df$X3) / nrow(pop_df)^2
    vmat[2, 1] <- vmat[1, 2]
    vmat[3, 1] <- vmat[1, 3]
    vmat[2, 2] <- sum((1 / p0_df$pi0^2 - 1 / p0_df$pi0) * p0_df$X2^2) / nrow(pop_df)^2
    vmat[2, 3] <-
      sum((1 / p0_df$pi0^2 - 1 / p0_df$pi0) * p0_df$X2 * p0_df$X3) / nrow(pop_df)^2
    vmat[3, 2] <- vmat[2, 3]
    vmat[3, 3] <- sum((1 / p0_df$pi0^2 - 1 / p0_df$pi0) * p0_df$X3^2) / nrow(pop_df)^2

    vmat[4, 4] <- (1 / nrow(p1_df) - 1 / nrow(pop_df)) * var(p1_df$X1)
    vmat[5, 5] <- (1 / nrow(p1_df) - 1 / nrow(pop_df)) * var(p1_df$X3)
    vmat[4, 5] <- (1 / nrow(p1_df) - 1 / nrow(pop_df)) * cov(p1_df$X1, p1_df$X3)
    vmat[5, 4] <- vmat[4, 5]

    vmat[6, 6] <- sum((1 / p2_df$pi2^2 - 1 / p2_df$pi2) * p2_df$X1^2) / nrow(pop_df)^2
    vmat[7, 7] <- sum((1 / p2_df$pi2^2 - 1 / p2_df$pi2) * p2_df$X2^2) / nrow(pop_df)^2
    vmat[6, 7] <- 
      sum((1 / p2_df$pi2^2 - 1 / p2_df$pi2) * p2_df$X1 * p2_df$X2) / nrow(pop_df)^2
    vmat[7, 6] <- vmat[6, 7]

    xgls <- drop(MASS::ginv(t(x_mat) %*% solve(vmat) %*% x_mat) %*% 
                (t(x_mat) %*% solve(vmat) %*% xhat))

    T1 <- c(nrow(pop_df),
            xgls[1] * nrow(pop_df),
            xgls[2] * nrow(pop_df),
            xgls[3] * nrow(pop_df),
            sum(gdip * qi))
  } else {
    T1 <- c(nrow(pop_df),
            sum(pop_df$X1),
            sum(pop_df$X2),
            sum(pop_df$X3),
            sum(gdip * qi))
  }

  dc_w <- solveGE(zmat, w1i, d0i, T1, qi, entropy)

  # Mean Estimation
  theta <- sum(p0_df$Y * dc_w * w1i) / nrow(pop_df)
  z_dc <- model.matrix(~1 + p0_df$X1 + p0_df$X2 + p0_df$X3 + gdi)
  gam_dc <- 
    solve(t(z_dc) %*% diag(d0i * gpinv * qi) %*% z_dc,
          t(z_dc) %*% diag(d0i * gpinv) %*% p0_df$Y)
  coefs <- as.numeric(gam_dc)
  u_mat <- 
  matrix(c(rep(1, nrow(p0_df)), p0_df$X1, p0_df$X2, p0_df$X3, gdi * qi), ncol = 5)
  eps <- (p0_df$Y - drop(u_mat %*% coefs))
  if (estT1) {
    # Variance Estimation
    # This is the variance of bar x gls not x gls total.
    var_x <- MASS::ginv(t(x_mat) %*% MASS::ginv(vmat) %*% x_mat)

    v1 <- matrix(coefs[2:4], nrow = 1) %*% var_x %*% matrix(coefs[2:4], ncol = 1)
    v2 <- sum((1 - p0_df$pi0) / p0_df$pi0^2 * eps^2) / nrow(pop_df)^2
    v_hat <- v1 + v2 
  } else {
    v2 <- sum((1 - p0_df$pi0) / p0_df$pi0^2 * eps^2) / nrow(pop_df)^2
    v_hat <- v2
  }
  if (reg) {
    theta <- 
    (drop(T1 %*% gam_dc) + 
      sum(1 / p0_df$pi0 * (p0_df$Y - drop(z_dc %*% gam_dc)))) / nrow(pop_df) 
  }

  return(list(theta = theta, vhat = as.numeric(v_hat)))
}

# ********************************************************
# * Multisource with Estimated Alpha Sampling Algorithms *
# ********************************************************
# Exploration file: ../explore/20241101-msdcsim_estalpha.qmd

msa_est <- function(vars) {

  p0_df <- vars[[1]]
  p1_df <- vars[[2]]
  p2_df <- vars[[3]]
  pop_df <- vars[[4]]

  qi = 1
  entropy = "EL"
  reg = FALSE

  # Since we have a SRS, it is fine to use this info.
  N_obs <- nrow(pop_df)
  d_vec <- 1 / (p0_df$pi0)
  if (entropy == "EL") {
    gdi <- -1 / d_vec
    gpinv <- d_vec^2
  } else if (entropy == "ET") {
    gdi <- log(d_vec)
    gdip <- log(1 / pop_df$pi0)
    gpinv <- d_vec
  } else {
    stop("This type of entropy has not yet been implemented.")
  }

  samp_X <- model.matrix(~X1 + X2 + X3, data = p0_df)
  zmat <- cbind(samp_X, gdi)
  w1i <- 1

  x1_0 <- sum(p0_df$X1 / p0_df$pi0) / N_obs
  x2_0 <- sum(p0_df$X2 / p0_df$pi0) / N_obs
  x3_0 <- sum(p0_df$X3 / p0_df$pi0) / N_obs

  x1_1 <- sum(p1_df$X1 / p1_df$pi1) / N_obs
  x3_1 <- sum(p1_df$X3 / p1_df$pi1) / N_obs

  x1_2 <- sum(p2_df$X1 / p2_df$pi2) / N_obs
  x2_2 <- sum(p2_df$X2 / p2_df$pi2) / N_obs

  xhat <- c(x1_0, x2_0, x3_0, x1_1, x3_1, x1_2, x2_2)
  x_mat <- 
  matrix(c(1, 0, 0,
           0, 1, 0,
           0, 0, 1,
           1, 0, 0,
           0, 0, 1,
           1, 0, 0,
           0, 1, 0), ncol = 3, byrow = TRUE)
  vmat <- matrix(0, nrow = 7, ncol = 7)
  vmat[1, 1] <- sum((1 / p0_df$pi0^2 - 1 / p0_df$pi0) * p0_df$X1^2) / nrow(pop_df)^2
  vmat[1, 2] <- 
    sum((1 / p0_df$pi0^2 - 1 / p0_df$pi0) * p0_df$X1 * p0_df$X2) / nrow(pop_df)^2
  vmat[1, 3] <-
    sum((1 / p0_df$pi0^2 - 1 / p0_df$pi0) * p0_df$X1 * p0_df$X3) / nrow(pop_df)^2
  vmat[2, 1] <- vmat[1, 2]
  vmat[3, 1] <- vmat[1, 3]
  vmat[2, 2] <- sum((1 / p0_df$pi0^2 - 1 / p0_df$pi0) * p0_df$X2^2) / nrow(pop_df)^2
  vmat[2, 3] <-
    sum((1 / p0_df$pi0^2 - 1 / p0_df$pi0) * p0_df$X2 * p0_df$X3) / nrow(pop_df)^2
  vmat[3, 2] <- vmat[2, 3]
  vmat[3, 3] <- sum((1 / p0_df$pi0^2 - 1 / p0_df$pi0) * p0_df$X3^2) / nrow(pop_df)^2

  vmat[4, 4] <- (1 / nrow(p1_df) - 1 / nrow(pop_df)) * var(p1_df$X1)
  vmat[5, 5] <- (1 / nrow(p1_df) - 1 / nrow(pop_df)) * var(p1_df$X3)
  vmat[4, 5] <- (1 / nrow(p1_df) - 1 / nrow(pop_df)) * cov(p1_df$X1, p1_df$X3)
  vmat[5, 4] <- vmat[4, 5]

  vmat[6, 6] <- sum((1 / p2_df$pi2^2 - 1 / p2_df$pi2) * p2_df$X1^2) / nrow(pop_df)^2
  vmat[7, 7] <- sum((1 / p2_df$pi2^2 - 1 / p2_df$pi2) * p2_df$X2^2) / nrow(pop_df)^2
  vmat[6, 7] <- 
    sum((1 / p2_df$pi2^2 - 1 / p2_df$pi2) * p2_df$X1 * p2_df$X2) / nrow(pop_df)^2
  vmat[7, 6] <- vmat[6, 7]

  xgls <- drop(MASS::ginv(t(x_mat) %*% solve(vmat) %*% x_mat) %*% 
              (t(x_mat) %*% solve(vmat) %*% xhat))

  T1 <- c(N_obs,
          xgls[1] * N_obs,
          xgls[2] * N_obs,
          xgls[3] * N_obs,
          NA)

  # Estimating Alpha
  alpha_init <- (-1) * nrow(p0_df) / N_obs
  a_res <- 
    optim(par = alpha_init, fn = solveGE_alpha, 
        zmat = zmat, w1i = w1i, d2i = d_vec, T1 = T1, qi = qi, entropy = entropy)

  if (a_res$convergence != 0) {
    warning(paste0("a_res has convergence", a_res$convergence))
  }

  dc_w <- 
    solveGE_alpha(a_res$par, zmat, w1i, d_vec, T1, qi, entropy, retw = TRUE)

  # Mean Estimation
  # Do we want N_obs? Or an estimated N_hat? 
  # N_obs is probably fine since we have SRS.
  theta <- sum(p0_df$Y * dc_w * w1i) / N_obs
  z_dc <- samp_X
  gam_dc <- 
    solve(t(z_dc) %*% diag(d_vec * gpinv * qi) %*% z_dc,
          t(z_dc) %*% diag(d_vec * gpinv) %*% p0_df$Y)
  coefs <- as.numeric(gam_dc)
  u_mat <- 
  matrix(c(rep(1, nrow(p0_df)), p0_df$X1, p0_df$X2, p0_df$X3), ncol = 4)
  eps <- (p0_df$Y - drop(u_mat %*% coefs))
  # Variance Estimation
  # This is the variance of bar x gls not x gls total.
  var_x <- MASS::ginv(t(x_mat) %*% MASS::ginv(vmat) %*% x_mat)

  v1 <- matrix(coefs[2:4], nrow = 1) %*% var_x %*% matrix(coefs[2:4], ncol = 1)
  v2 <- sum((1 - p0_df$pi0) / p0_df$pi0^2 * eps^2) / nrow(pop_df)^2
  v_hat <- v1 + v2 
  if (reg) {
    theta <- 
    (drop(T1 %*% gam_dc) + 
      sum(1 / p0_df$pi0 * (p0_df$Y - drop(z_dc %*% gam_dc)))) / nrow(pop_df) 
  }

  return(list(theta = theta, vhat = as.numeric(v_hat)))
}


