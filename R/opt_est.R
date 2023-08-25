# Title: opt_est.R
# Created by: Caleb Leedy
# Created on: August 19, 2023
# Updated on: August 22, 2023
# Purpose: This file contains the functions needed for computation of the
# optimal estimator in the non-monotone case.

# ***********
# * Outline *
# ***********
# 1. Generate data.
# 2. EM Algorithm for Optimal Estimate
# 3. Proposed Estimator

# *************
# * Libraries *
# *************

library(dplyr)
library(stringr)

# ************************
# * Additional Functions *
# ************************

expit <- function(x) {1 / (1 + exp(-x))}

# *****************
# * Generate Data *
# *****************

#' This function generates a multivariate normal distribution with missingness.
#'
#' @param n - An integer for the number of observations
#'
#' @return df - A data frame.
gen_optsim_data <- function(n, theta = 0, cor_xe1 = 0, cor_xe2 = 0, cor_e1e2 = 0) {

  # 1. Generate X, Y1, Y2
  mean_vec <- c(0, 0, 0)
  sigma_mat <-
    matrix(c(1, cor_xe1, cor_xe2, cor_xe1, 1, cor_e1e2, cor_xe2, cor_e1e2, 1),
           nrow = 3)

  data <- MASS::mvrnorm(n, mean_vec, sigma_mat)

  x_vec <- data[, 1]
  y1 <- x_vec + data[, 2]
  y2 <- theta + x_vec + data[, 3]

  # 2. Generate the missing patterns:
  # | Pattern | Prob |
  # |---------|------|
  # | 11      | 0.4  |
  # | 10      | 0.2  |
  # | 01      | 0.2  |
  # | 00      | 0.2  |
  u_vec <- runif(n)
  patt_vec <-
    case_when(
      u_vec < 0.4 ~ "11",
      u_vec < 0.6 ~ "10",
      u_vec < 0.8 ~ "01",
      u_vec < 1.0 ~ "00"
    )
  delta_y1 <-
    case_when(
      patt_vec == "11" | patt_vec == "10" ~ 1,
      patt_vec == "01" | patt_vec == "00" ~ 0
    )
  delta_y2 <-
    case_when(
      patt_vec == "11" | patt_vec == "01" ~ 1,
      patt_vec == "10" | patt_vec == "00" ~ 0
    )
  prob_1 <- 0.6
  prob_2 <- 0.6
  prob_11 <- 0.4

  tibble(X = x_vec, Y1 = y1, Y2 = y2,
         delta_1 = delta_y1, delta_2 = delta_y2,
         prob_11 = prob_11, prob_1 = prob_1, prob_2 = prob_2)

}

# ****************
# * EM Algorithm *
# ****************

#' The get_cond_exp function gets the conditional expectation from a conditional
#' normal distribution. This is the expectation of Y1 given Y2.
#'
#' @param obs_y2 - A n x c2 matrix
#' @param mean_y1 - A n x c1 matrix
#' @param mean_y2 - A n x c2 matrix
#' @param sigma_mat - A 2 x 2 list. Each element could be a matrix or scaler.
#'
#' @return ret - A n x c1 matrix
#'
#' @details
#' Preconditions:
#' - dim(sigma_mat[[1]][[1]]) = c(c1, c1)
#' - dim(sigma_mat[[1]][[2]]) = c(c1, c2)
#' - dim(sigma_mat[[2]][[1]]) = c(c2, c1)
#' - dim(sigma_mat[[2]][[2]]) = c(c2, c2)
#' - sigma_mat[[2]][[2]] must be invertible
#'
#' Postconditions:
#'
get_cond_exp <- function(obs_y2, mean_y1, mean_y2, sigma_mat) {

  if (length(mean_y1) == 1) {
    mean_y1 <- matrix(rep(mean_y1, nrow(obs_y2)), ncol = 1)
  }
  if (length(mean_y2) == 1) {
    mean_y2 <- matrix(rep(mean_y2, nrow(obs_y2)), ncol = ncol(obs_y2))
  }

  # FIXME: The dimensions on this are wrong. It does not work for Y1 = (Y1, Y2), Y2 = X.
  out <- 
  purrr::map(1:nrow(obs_y2), function(i) {
    t(mean_y1[i, ]) + 
    sigma_mat[[1]][[2]] %*% solve(sigma_mat[[2]][[2]]) %*% 
      matrix(obs_y2[i, ] - mean_y2[i, ], ncol = ncol(mean_y1))
  }) |> rbind()

  return(t(out))

}


# The EM algorithm consists of two steps: expectation and maximization.
# The E-Step consists of:
# - Imputing missing values for Y1 and Y2. This means that we need to compute the
#   conditional distribution of Y1 | X, Y2 | X, Y1 | X, Y2, and Y2 | X, Y1.
# - We impute based on these distributions.
#
# The M-Step consists of:
# - Estimating parameters from the imputed data from the distributions used in
#   the E-Step.
#
# We assume that the covariance matrix is known.
# The true covariance matrix for (X, Y1, Y2) is
# [1, 1, 1]
# [1, 2, 1]
# [1, 1, 2]
#
# Recall that:
# Z1 | Z2 ~ N(mu_1 + Sigma_12 Sigma_22^(-1) (z_2 - my_2),
#             Sigma_11 - Sigma_12 Sigma_22^(-1) Sigma_21)

em_optsim <- function(df) {

  sigma_mat <- matrix(c(1, 1, 1, 1, 2, 1, 1, 1, 2), nrow = 3)
  s_xy1 <- matrix(c(1, 1, 1, 2), nrow = 2)
  s_xy2 <- s_xy1
  s_y1y2 <- matrix(c(2, 1, 1, 2), nrow = 2)

  # Steps:
  # 1. Initialization
  # 2. E-Step
  # 3. M-Step

  df_11 <- filter(df, delta_1 == 1, delta_2 == 1)
  df_0 <- filter(df, delta_1 == 0, delta_2 == 0)
  df_1 <- filter(df, delta_1 == 1)
  df_2 <- filter(df, delta_2 == 1)

  # 1. Initialization
  # For the initialization we need the following parameters:
  # mean_y1, mean_y2,
  mean_x <- mean(df$X)
  mean_y1 <- mean(df_11$Y1)
  mean_y2 <- mean(df_11$Y2)

  max_iter <- 100
  tol <- 1e-6
  for (iter in 1:max_iter) {

    old_my1 <- mean_y1
    old_my2 <- mean_y2

    # 2. E-Step
    # Impute missing Y1
    obs <- dplyr::select(df_2, X, Y2) |> as.matrix()
    mis_y1 <- 
    get_cond_exp(obs, mean_y1,
                 matrix(c(mean_x, mean_y2), ncol = 2, nrow = nrow(obs)),
                 list(list(2, matrix(c(1, 1), nrow = 1)),
                      list(matrix(c(1, 1), nrow = 2), matrix(c(1, 1, 1, 2), nrow = 2))))

    # Impute missing Y2
    obs <- dplyr::select(df_1, X, Y1) |> as.matrix()
    mis_y2 <-
    get_cond_exp(obs, mean_y2,
                 matrix(c(mean_x, mean_y1), ncol = 2, nrow = nrow(obs)),
                 list(list(2, matrix(c(1, 1), nrow = 1)),
                      list(matrix(c(1, 1), nrow = 2), matrix(c(1, 1, 1, 2), nrow = 2))))

    # Impute missing Y1 and Y2
    # FIXME:
    mis_y1y2 <-
    get_cond_exp(matrix(df_0$X, ncol = 1),
                 matrix(c(mean_y1, mean_y2), ncol = 2, nrow = nrow(df_0)),
                 mean_x,
                 list(list(matrix(c(2, 1, 1, 2), nrow = 2), matrix(c(1, 1), nrow = 2)),
                      list(matrix(c(1, 1), nrow = 1), 1)))

    # 3. M-Step
    mean_y1 <-
      (sum(df_11$Y1) + sum(df_1$Y1) + sum(mis_y1) + sum(mis_y1y2[1, ])) / nrow(df)
    mean_y2 <-
      (sum(df_11$Y2) + sum(df_2$Y2) + sum(mis_y2) + sum(mis_y1y2[2, ])) / nrow(df)

    if ((mean_y1 - old_my1)^2 + (mean_y2 - old_my2)^2 < tol) {
      # We have g = E[Y2].
      return(mean_y2)
    }
  }

  warning("Max iters reached.")
  return(mean_y2)

}


# **********************
# * Proposed Estimator *
# **********************

# #' This function computes the probability that a given piece of data is
# #' contained in the sample.
# #'
# #' @param df - A data frame.
# #' @param obs - A vector of columns
# #'
# #' @return ret - A vector
# #'
# #' @details:
# #' Preconditions:
# #' - obs is a vector containing an element of the power set of \{"Y1", "Y2"\}.
# #' Postconditions:
# #' - ret should not depend on Y1 or Y2 but it can depend on X.
# get_pi_w <- function(df, obs) {
# 
# }

#' This computes the expectation of g conditional on the variables in the 
#' conditional vector.
#'
#' @param df - A data frame
#' @param g - A string
#' @param cond - A vector of variables to condition on.
#' @param pred_df - A data frame.
#' @param d_vec - A vector of length n.
#'
#' @return ret - A numeric vector
#'
#' @details 
#' Preconditions:
#' - g in names(df)
#' - all(cond in names(df))
#' - all(cond in names(pred_df))
#' - length(d_vec) == n
#' - d_vec contains only 0s and 1s.
#' - sum(d_vec) = nrow(pred_df)'
#' Postconditions:
#' - length(ret) == n
expect_g <- function(df, g, cond, pred_df, d_vec) {

  ret <- numeric(length(d_vec))

  if (g %in% cond) {
    res <- pred_df[[g]]
    ret[d_vec == 1] <- res
    return(ret)
  } else {
    mod <- lm(as.formula(str_c(g, " ~ ", str_flatten(cond, " + "))), data = df)
    res <- predict(mod, pred_df)
    ret[d_vec == 1] <- res
    return(ret)
  }
}

#' The proposed non-monotone estimator. This function gives the result based
#' on an input data frame.
#'
#' @param df - A data frame
#' @param oracle - A boolean value
#'
#' @return ret - A numeric value
#'
#' @notes
#' Preconditions:
#'  - df has the following columns: X, Y1, Y2, delta_y1, delta_y2, prob_vec
#'  - If delta_yi = 1 then we assume that Y_i is observed if delta_yi = 0 then
#'  - we assume that Y_i is missing.
prop_nmono_est <- function(df, oracle = TRUE) {

  # HACK: We only test the case where we have g(X, Y1, Y2) = Y2.
  gfun <- "Y2"

  df_11 <- filter(df, delta_1 == 1, delta_2 == 1)
  df_1 <- filter(df, delta_1 == 1)
  df_2 <- filter(df, delta_2 == 1)

  # Steps:
  # 1. Compute pi weights.
  # 2. Compute expectations.
  # 3. Get estimate.

  # 1. Compute pi weights.
  # pi_1+
  # pi_2+
  # pi_11
  pi1_w <- df$prob_1
  pi2_w <- df$prob_2
  pi11_w <- df$prob_11

  # 2. Compute expectations.
  # E[g | X]
  # E[g | X, Y1]
  # E[g | X, Y2] = Y2
  eg_x <- expect_g(df, g = gfun, cond = c("X"), pred_df = df, d_vec = rep(1, nrow(df)))
  eg_xy1 <-
    expect_g(df, g = gfun, cond = c("X", "Y1"), pred_df = df_1, d_vec = df$delta_1)
  eg_xy2 <-
    expect_g(df, g = gfun, cond = c("X", "Y2"), pred_df = df_2, d_vec = df$delta_2)

  # 3. Get estimate.
  # This is correct because eg_xyi both estimate E[g | X, Yi] using all of the 
  # observed Y_i values. Not just the values in A_{01} for example.
  mean(eg_x) + 
    mean(df$delta_1 / pi1_w * (eg_xy1 - eg_x)) +
    mean(df$delta_2 / pi2_w * (eg_xy2 - eg_x)) +
    mean(df$delta_1 * df$delta_2 / pi11_w * (df[[gfun]] - eg_xy1 - eg_xy2 + eg_x))
}

prop2_nmono_est <- function(df, oracle = TRUE) {

  # HACK: We only test the case where we have g(X, Y1, Y2) = Y2.
  gfun <- "Y2"

  df_11 <- filter(df, delta_1 == 1, delta_2 == 1)
  df_1 <- filter(df, delta_1 == 1)
  df_2 <- filter(df, delta_2 == 1)

  # Steps:
  # 1. Compute pi weights.
  # 2. Compute expectations.
  # 3. Get estimate.

  # 1. Compute pi weights.
  # pi_1+
  # pi_2+
  # pi_11
  pi1_w <- df$prob_1
  pi2_w <- df$prob_2
  pi11_w <- df$prob_11

  # 2. Compute expectations.
  # E[g | X]
  # E[g | X, Y1]
  # E[g | X, Y2] = Y2
  eg_x <- expect_g(df, g = gfun, cond = c("X"), pred_df = df, d_vec = rep(1, nrow(df)))
  eg_xy1 <-
    expect_g(df, g = gfun, cond = c("X", "Y1"), pred_df = df_1, d_vec = df$delta_1)
  eg_xy2 <-
    expect_g(df, g = gfun, cond = c("X", "Y2"), pred_df = df_2, d_vec = df$delta_2)

  # 3. Get estimate.
  # This is correct because eg_xyi both estimate E[g | X, Yi] using all of the 
  # observed Y_i values. Not just the values in A_{01} for example.
  mean(eg_x) + 
    mean(df$delta_1 / pi1_w * (eg_xy1 - eg_x)) +
    mean(df$delta_2 / pi2_w * (eg_xy2 - eg_x)) 
}

# ********************
# * Linear Estimator *
# ********************

# FIXME: Do we need to account for survey weights?
opt_lin_est <- function(df,
                        mean_x = 0, var_x = 1,
                        cov_xy1 = 1, cov_xy2 = 1,
                        var_y1 = 1, var_y2 = 1) {

  if (mean_x != 0) {
    df$X <- df$X - mean_x
  }

  df$W1 <- df$Y1 - (cov_xy1 / var_x) * df$X
  df$W2 <- df$Y2 - (cov_xy1 / var_x) * df$X

  df_11 <- filter(df, delta_1 == 1, delta_2 == 1)
  df_1 <- filter(df, delta_1 == 1)
  df_2 <- filter(df, delta_2 == 1)

  w_mat <- 
    matrix(c(mean(df_11$W1), mean(df_11$W2), mean(df_1$W1), mean(df_2$W2)),
           ncol = 1)

  m_mat <- matrix(c(1, 0, 1, 0, 0, 1, 0, 1), ncol = 2)
  v_mat <- diag(c(var_y1 / nrow(df_11),
                  var_y2 / nrow(df_11),
                  var_y1 / nrow(df_1),
                  var_y2 / nrow(df_2)))

  # WLS
  mu_hat <- 
    solve(t(m_mat) %*% solve(v_mat) %*% m_mat) %*% 
    t(m_mat) %*% solve(v_mat) %*% w_mat

  # G function is Y2
  return(mu_hat[2, ])

}
