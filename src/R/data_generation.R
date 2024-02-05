# Title: data_generation.R
# Created by: Caleb Leedy
# Created on: February 03, 2024
# Purpose: This file contains the code for all of the functions used to generate
# data within simulation studies. I needed to store them here so that they could
# be reused.

#' A data generation function to create a simple random sample of equally sized
#' segments of X_1, X_2, and Y. There variables are correlated but all normal.
#'
#' @param n_obs_seg - An integer greater than zero. This is the number of
#' observations in each segment
#' @param mu - The true mean of Y
#' @param rho - The correlation between (X_1 - X_2) and (Y - X_1).
#'
#' @details We always observe X1 but X2 and Y are only observed depending on the
#' value of \delta_{ij}. If i = 1 then X2 is observed. If j = 1 then Y is
#' observed.
gen_simple_data <- function(n_obs_seg, mu, rho) {

  # This code is modified from the case where we had X, Y_1, and Y_2. It only
  # changes at the end to X_1, X_2, and Y.
  # A_11
  x <- rnorm(n_obs_seg)
  e_mat <- MASS::mvrnorm(n_obs_seg, c(0, 0), matrix(c(1, rho, rho, 1), nrow = 2))
  e1 <- e_mat[, 1]
  e2 <- e_mat[, 2]
  y1 <- x + e1
  y2 <- mu + x + e2

  df_11 <- 
    tibble(X1 = x, X2 = y1, Y = y2,
           delta_00 = 0, delta_10 = 0, delta_01 = 0, delta_11 = 1,
           pi_00 = 0.25, pi_10 = 0.25, pi_01 = 0.25, pi_11 = 0.25)
  # A_10
  x <- rnorm(n_obs_seg)
  e_mat <- MASS::mvrnorm(n_obs_seg, c(0, 0), matrix(c(1, rho, rho, 1), nrow = 2))
  e1 <- e_mat[, 1]
  e2 <- e_mat[, 2]
  y1 <- x + e1
  y2 <- mu + x + e2

  df_10 <- 
    tibble(X1 = x, X2 = y1, Y = y2,
           delta_00 = 0, delta_10 = 1, delta_01 = 0, delta_11 = 0,
           pi_00 = 0.25, pi_10 = 0.25, pi_01 = 0.25, pi_11 = 0.25)
  # A_01
  x <- rnorm(n_obs_seg)
  e_mat <- MASS::mvrnorm(n_obs_seg, c(0, 0), matrix(c(1, rho, rho, 1), nrow = 2))
  e1 <- e_mat[, 1]
  e2 <- e_mat[, 2]
  y1 <- x + e1
  y2 <- mu + x + e2

  df_01 <- 
    tibble(X1 = x, X2 = y1, Y = y2,
           delta_00 = 0, delta_10 = 0, delta_01 = 1, delta_11 = 0,
           pi_00 = 0.25, pi_10 = 0.25, pi_01 = 0.25, pi_11 = 0.25)
  # A_00
  x <- rnorm(n_obs_seg)
  e_mat <- MASS::mvrnorm(n_obs_seg, c(0, 0), matrix(c(1, rho, rho, 1), nrow = 2))
  e1 <- e_mat[, 1]
  e2 <- e_mat[, 2]
  y1 <- x + e1
  y2 <- mu + x + e2

  df_00 <- 
    tibble(X1 = x, X2 = y1, Y = y2,
           delta_00 = 1, delta_10 = 0, delta_01 = 0, delta_11 = 0,
           pi_00 = 0.25, pi_10 = 0.25, pi_01 = 0.25, pi_11 = 0.25)

  df <- bind_rows(df_11, df_10, df_01, df_00)

  # Return
  # X, Y_1, Y_2,
  # delta_00, delta_10, delta_01, delta_11,
  # pi_00, pi_10, pi_01, pi_11
  return(df)
  
}
