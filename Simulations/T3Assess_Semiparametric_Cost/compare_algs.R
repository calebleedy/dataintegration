# Title: compare_algs.R
# Created by: Caleb Leedy
# Created on: January 09, 2024
# Purpose: This file contains the code to compare the parametric and 
# semiparametric estimators.

# *************
# * Libraries *
# *************

library(dplyr)
library(CVXR)

source("R/opt_est.R")

# **************
# * Algorithms *
# **************

#' The min_var_alg is the algorithm that construct constants that minimize 
#' the variance of the independent estimators subject to given contraints.
#'
#' @param df - A data frame generated from R/opt_est.R::gen_optsim_data.
#' @param gfun - A string
#' @param mod - A string. Currently, this only supports one of
#' "parametric", "sumzero", "outcome_robust", "response_robust", or
#' "double_robust".
min_var_alg <- 
  function(df,
           gfun = "Y2",
           mod = "para",
           theta = NA,
           mean_x = 0,
           cov_e1e2 = 0) {

  df <- mutate(df, g_i = eval(rlang::parse_expr(gfun)))
  df_11 <- filter(df, delta_1 == 1, delta_2 == 1)
  df_10 <- filter(df, delta_1 == 1, delta_2 == 0)
  df_01 <- filter(df, delta_1 == 0, delta_2 == 1)
  df_00 <- filter(df, delta_1 == 0, delta_2 == 0)

  if (is.na(theta)) {
    theta <- opt_lin_est(df, gfun = "Y2", mean_x = mean_x, cov_y1y2 = cov_e1e2)
  }

  # We use the population functional forms for gamma_hat and v_gamma
  if (gfun == "Y2") {

    Eg <- theta
    Egx <- theta + df$X 
    Egxy1 <- theta + df$X + cov_e1e2 * (df$Y1 - df$X)
    Egxy2 <- df$Y2
    EEg2 <- theta^2 + 2
    EEgx2 <- theta^2 + 1
    EEgxy12 <- theta^2 + 1 + cov_e1e2^2
    EEgxy22 <- theta^2 + 2
    EEgxy1Egxy2 <- theta^2 + 1 + cov_e1e2^2

  } else if (gfun == "Y1^2 * Y2") {

    Eg <- 2 * theta
    Egx <- df$X^3 + theta * df$X^2 + theta + 2 * df$X * cov_e1e2 + df$X
    Egxy1 <- df$Y1^2 * (theta + df$X + cov_e1e2 * (df$Y1 - df$X))
    Egxy2 <- df$Y2 * (df$X^2 + 2 * cov_e1e2 * df$X * (df$Y2 - df$X - theta) + 
                      1 + cov_e1e2^2 * (df$Y2 - df$X - theta)^2 - cov_e1e2^2)
    EEg2 <- 12 * (theta^2 + 2 * cov_e1e2^2 + 4 * cov_e1e2 + 4)
    EEgx2 <- 6 * theta^2 + 4 * cov_e1e2^2 + 16 * cov_e1e2 + 22
    EEgxy12 <- 12 * (theta^2 + 3 * cov_e1e2^2 + 4 * cov_e1e2 + 3)
    EEgxy22 <- 
      6 * theta^2 + 4 * theta^2 * cov_e1e2^2 + 2 * theta^2 * cov_e1e2^4 + 34 +
      32 * cov_e1e2 + 20 * cov_e1e2^2 + 16 * cov_e1e2^3 + 18 * cov_e1e2^4
    EEgxy1Egxy2 <- 
      theta^2 * (6 + 4 * cov_e1e2^2 + 2 * cov_e1e2^4) + 22 + 24 * cov_e1e2 +
      26 * cov_e1e2^2 + 20 * cov_e1e2^3 + 18 * cov_e1e2^4 + 4 * cov_e1e2^5 +
      6 * cov_e1e2^6

  } else {
    stop("We have only implemented two g-functions.")
  }

  g_11 <- mean(df$delta_00 / df$prob_00 * (Egx))
  g_21 <- mean(df$delta_10 / df$prob_10 * (Egx))
  g_22 <- mean(df$delta_10 / df$prob_10 * (Egxy1))
  g_31 <- mean(df$delta_01 / df$prob_01 * (Egx))
  g_33 <- mean(df$delta_01 / df$prob_01 * (Egxy2))
  g_41 <- mean(df$delta_11 / df$prob_11 * (Egx))
  g_42 <- mean(df$delta_11 / df$prob_11 * (Egxy1))
  g_43 <- mean(df$delta_11 / df$prob_11 * (Egxy2))
  g_44 <- mean(df$delta_11 / df$prob_11 * (df$g_i))

  gam_vec <- c(g_11, g_21, g_22, g_31, g_33, g_41, g_42, g_43, g_44)

  v_gam <- matrix(rep(-Eg^2, 81), nrow = 9)
  v_gam[1, 1] <- 1 / df$prob_00[1] * (EEgx2) - Eg^2

  v_gam[2, 2] <- 1 / df$prob_10[1] * (EEgx2) - Eg^2
  v_gam[2, 3] <- 1 / df$prob_10[1] * (EEgx2) - Eg^2
  v_gam[3, 2] <- v_gam[2, 3]
  v_gam[3, 3] <- 1 / df$prob_10[1] * (EEgxy12) - Eg^2 

  v_gam[4, 4] <- 1 / df$prob_01[1] * (EEgx2) - Eg^2
  v_gam[4, 5] <- 1 / df$prob_01[1] * (EEgx2) - Eg^2
  v_gam[5, 4] <- v_gam[4, 5]
  v_gam[5, 5] <- 1 / df$prob_01[1] * (EEgxy22) - Eg^2 

  v_gam[6, 6] <- 1 / df$prob_11[1] * (EEgx2) - Eg^2
  v_gam[6, 7] <- 1 / df$prob_11[1] * (EEgx2) - Eg^2
  v_gam[7, 6] <- v_gam[6, 7] 
  v_gam[6, 8] <- 1 / df$prob_11[1] * (EEgx2) - Eg^2
  v_gam[8, 6] <- v_gam[6, 8]
  v_gam[6, 9] <- 1 / df$prob_11[1] * (EEgx2) - Eg^2
  v_gam[9, 6] <- v_gam[6, 9]
  v_gam[7, 7] <- 1 / df$prob_11[1] * (EEgxy12) - Eg^2
  v_gam[7, 8] <- 1 / df$prob_11[1] * (EEgxy1Egxy2) - Eg^2
  v_gam[8, 7] <- v_gam[7, 8] 
  v_gam[7, 9] <- 1 / df$prob_11[1] * (EEgxy12) - Eg^2
  v_gam[9, 7] <- v_gam[7, 9] 
  v_gam[8, 8] <- 1 / df$prob_11[1] * (EEgxy22) - Eg^2
  v_gam[8, 9] <- 1 / df$prob_11[1] * (EEgxy22) - Eg^2
  v_gam[9, 8] <- v_gam[8, 9]
  v_gam[9, 9] <- 1 / df$prob_11[1] * (EEg2) - Eg^2

  v_gam <- v_gam / nrow(df)

  # CVXR
  c_hat <- Variable(length(gam_vec))

  # TODO: We can change the inequality constraints to equality constraints.
  if (mod == "para") {
    con1 <- matrix(rep(1, length(gam_vec)), nrow = 1) %*% c_hat <= 1
    con2 <- matrix(rep(1, length(gam_vec)), nrow = 1) %*% c_hat >= 1
    con_list <- list(con1, con2)
  } else if (mod == "sumzero") {
    con1 <- matrix(c(rep(0, length(gam_vec) - 1), 1), nrow = 1) %*% c_hat <= 1
    con2 <- matrix(c(rep(0, length(gam_vec) - 1), 1), nrow = 1) %*% c_hat >= 1
    con3 <- matrix(c(rep(1, length(gam_vec) - 1), 0), nrow = 1) %*% c_hat <= 0
    con4 <- matrix(c(rep(1, length(gam_vec) - 1), 0), nrow = 1) %*% c_hat >= 0
    con_list <- list(con1, con2, con3, con4)
  } else if (mod == "out_rob") {
    con1 <- matrix(c(rep(0, length(gam_vec) - 1), 1), nrow = 1) %*% c_hat <= 1
    con2 <- matrix(c(rep(0, length(gam_vec) - 1), 1), nrow = 1) %*% c_hat >= 1
    con3 <- matrix(c(1, 1, 0, 1, 0, 1, 0, 0, 0), nrow = 1) %*% c_hat <= 0
    con4 <- matrix(c(1, 1, 0, 1, 0, 1, 0, 0, 0), nrow = 1) %*% c_hat >= 0
    con5 <- matrix(c(0, 0, 1, 0, 0, 0, 1, 0, 0), nrow = 1) %*% c_hat <= 0
    con6 <- matrix(c(0, 0, 1, 0, 0, 0, 1, 0, 0), nrow = 1) %*% c_hat >= 0
    con7 <- matrix(c(0, 0, 0, 0, 1, 0, 0, 1, 0), nrow = 1) %*% c_hat <= 0
    con8 <- matrix(c(0, 0, 0, 0, 1, 0, 0, 1, 0), nrow = 1) %*% c_hat >= 0
    con_list <- list(con1, con2, con3, con4, con5, con6, con7, con8)
  } else if (mod == "resp_rob") {
    con1 <- matrix(c(1, 0, 0, 0, 0, 0, 0, 0, 0), nrow = 1) %*% c_hat <= df$prob_00[1]
    con2 <- matrix(c(1, 0, 0, 0, 0, 0, 0, 0, 0), nrow = 1) %*% c_hat >= df$prob_00[1]
    con3 <- matrix(c(0, 1, 1, 0, 0, 0, 0, 0, 0), nrow = 1) %*% c_hat <= df$prob_10[1]
    con4 <- matrix(c(0, 1, 1, 0, 0, 0, 0, 0, 0), nrow = 1) %*% c_hat >= df$prob_10[1]
    con5 <- matrix(c(0, 0, 0, 1, 1, 0, 0, 0, 0), nrow = 1) %*% c_hat <= df$prob_01[1]
    con6 <- matrix(c(0, 0, 0, 1, 1, 0, 0, 0, 0), nrow = 1) %*% c_hat >= df$prob_01[1]
    con7 <- matrix(c(0, 0, 0, 0, 0, 1, 1, 1, 1), nrow = 1) %*% c_hat <= df$prob_11[1]
    con8 <- matrix(c(0, 0, 0, 0, 0, 1, 1, 1, 1), nrow = 1) %*% c_hat >= df$prob_11[1]
    con_list <- list(con1, con2, con3, con4, con5, con6, con7, con8)
  } else if (mod == "double_rob") {
    con1 <- matrix(c(rep(0, length(gam_vec) - 1), 1), nrow = 1) %*% c_hat <= 1
    con2 <- matrix(c(rep(0, length(gam_vec) - 1), 1), nrow = 1) %*% c_hat >= 1
    con3 <- matrix(c(1, 1, 0, 1, 0, 1, 0, 0, 0), nrow = 1) %*% c_hat <= 0
    con4 <- matrix(c(1, 1, 0, 1, 0, 1, 0, 0, 0), nrow = 1) %*% c_hat >= 0
    con5 <- matrix(c(0, 0, 1, 0, 0, 0, 1, 0, 0), nrow = 1) %*% c_hat <= 0
    con6 <- matrix(c(0, 0, 1, 0, 0, 0, 1, 0, 0), nrow = 1) %*% c_hat >= 0
    con7 <- matrix(c(0, 0, 0, 0, 1, 0, 0, 1, 0), nrow = 1) %*% c_hat <= 0
    con8 <- matrix(c(0, 0, 0, 0, 1, 0, 0, 1, 0), nrow = 1) %*% c_hat >= 0

    con9 <- matrix(c(1, 0, 0, 0, 0, 0, 0, 0, 0), nrow = 1) %*% c_hat <= df$prob_00[1]
    con10 <- matrix(c(1, 0, 0, 0, 0, 0, 0, 0, 0), nrow = 1) %*% c_hat >= df$prob_00[1]
    con11 <- matrix(c(0, 1, 1, 0, 0, 0, 0, 0, 0), nrow = 1) %*% c_hat <= df$prob_10[1]
    con12 <- matrix(c(0, 1, 1, 0, 0, 0, 0, 0, 0), nrow = 1) %*% c_hat >= df$prob_10[1]
    con13 <- matrix(c(0, 0, 0, 1, 1, 0, 0, 0, 0), nrow = 1) %*% c_hat <= df$prob_01[1]
    con14 <- matrix(c(0, 0, 0, 1, 1, 0, 0, 0, 0), nrow = 1) %*% c_hat >= df$prob_01[1]
    con15 <- matrix(c(0, 0, 0, 0, 0, 1, 1, 1, 1), nrow = 1) %*% c_hat <= df$prob_11[1]
    con16 <- matrix(c(0, 0, 0, 0, 0, 1, 1, 1, 1), nrow = 1) %*% c_hat >= df$prob_11[1]

    con_list <- list(con1, con2, con3, con4, con5, con6, con7, con8,
                     con9, con10, con11, con12, con13, con14, con15, con16)
  } else {
    stop("Other models not supported.")
  }

  obj <- Minimize(quad_form(c_hat, v_gam))
  problem <- Problem(obj, constraints = con_list)
  res <- solve(problem)

  # Output
  c_est <- res$getValue(c_hat)
  return(list(theta_est = as.numeric(gam_vec %*% c_est), c_hat = as.numeric(c_est)))

} 

