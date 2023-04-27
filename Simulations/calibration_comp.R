# Title: Simulations/calibration_comp.R
# Author: Caleb Leedy
# Date Created: April 25, 2023
# Purpose: This script contains the code to compare the calibration estimators
# with the regression estimators. The regression estimators can be originally
# found in Simulations/nonmonotone.R but they have been imported here.

# ******************
# * Script Outline *
# ******************
# A. Monotone Case
#
# B. Nonmonotone Case

# *************
# * Libraries *
# *************

library(dplyr)
library(CVXR)

library(doParallel)
library(doRNG)
library(parallelly)

# *****************************
# * Data Generation Functions *
# *****************************
source("generate_data.R")

# *************************************
# * Non/monotone Estimation Functions *
# *************************************
source("nmono_est_funs.R")

# ************************
# * Monotone Calibration *
# ************************

mono_est_weights <- function(df) {

  # Steps:
  # 1. Calibrate with w1 weights on Y_1
  # 2. Calibrate with w2 weights on Y_2 using w1 weights on Y_1.

  # 1. Calibrate with w1 weights on Y_1
  n1 <- sum(df$r1)
  w1 <- Variable(n1)
  objective <- Minimize(sum_squares(w1))
  Tx <- mean(df$x)
  xr1_vec <-
    filter(df, r1 == 1) %>%
    pull(x)
  constraints <- 
    list(t(w1) %*% matrix(xr1_vec, ncol = 1) == Tx,
         w1 >= 0,
         t(w1) %*% matrix(rep(1, n1), ncol = 1) == 1)
  p1 <- Problem(objective, constraints)
  res1 <- solve(p1)

  df$w1 <- 0
  df$w1[df$r1 == 1] <- res1$getValue(w1)

  # 2. Calibrate with w2 weights on Y_2 using w1 weights on Y_1.
  n2 <- sum(df$r2) # Only works because of monotone.
  w2 <- Variable(n2)
  objective <- Minimize(sum_squares(w2))
  Exw <- sum(df$w1 * df$x)
  Ey1w <- sum(df$w1 * df$y1)
  xr2_vec <- filter(df, r2 == 1) %>% pull(x)
  y1r2_vec <- filter(df, r2 == 1) %>% pull(y1)
  constraints <-
    list(w2 >= 0,
         t(w2) %*% matrix(rep(1, n2), ncol = 1) == 1,
         t(w2) %*% matrix(xr2_vec, ncol = 1) == Exw,
         t(w2) %*% matrix(y1r2_vec, ncol = 1) == Ey1w)
  p2 <- Problem(objective, constraints)
  res2 <- solve(p2)

  df$w2 <- 0
  df$w2[df$r2 == 1] <- res2$getValue(w2)

  # Objective is E[g] = E[y_2]
  sum(df$w2 * df$y2)

}

# *****************************
# * Monotone Simulation Study *
# *****************************

clust <- parallel::makeCluster(min(parallelly::availableCores() - 2, 100))
registerDoParallel(clust)

B <- 1000
true_mean <- 5

mc_theta <- foreach(b = 1:B, 
                  #.options.RNG = 1,
                  .packages = c("dplyr", "CVXR")) %dorng% {

  # Parameters
  n_sim <- 1000

  df <- mono_mar(n_sim, mean_y2 = true_mean)
  tibble(oracle = mean(df$y2),
         ipworacle = ipw_mono_est(df, est = FALSE),
         ipwest = ipw_mono_est(df, est = TRUE),
         semi = est_mono(df),
         calib = mono_est_weights(df))
         
} %>% bind_rows()

stopCluster(clust)

# ****************
# * Analyze Data *
# ****************

mc_theta %>%
 summarize(
   bias_oracle = mean(oracle) - true_mean,
   bias_ipworacle = mean(ipworacle) - true_mean,
   bias_ipwest = mean(ipwest) - true_mean,
   bias_semi = mean(semi) - true_mean,
   bias_calib = mean(calib) - true_mean,
   sd_oracle = sd(oracle),
   sd_ipworacle = sd(ipworacle),
   sd_ipwest = sd(ipwest),
   sd_semi = sd(semi),
   sd_calib = sd(calib),
   tstat_oracle = (mean(oracle) - true_mean) / sqrt(var(oracle) / B),
   tstat_ipworacle = (mean(ipworacle) - true_mean) / sqrt(var(ipworacle) / B),
   tstat_ipwest = (mean(ipwest) - true_mean) / sqrt(var(ipwest) / B),
   tstat_semi = (mean(semi) - true_mean) / sqrt(var(semi) / B), 
   tstat_calib = (mean(calib) - true_mean) / sqrt(var(calib) / B)) %>%
  tidyr::pivot_longer(cols = everything(),
               names_to = c(".value", "algorithm"),
               names_pattern = "(.*)_(.*)") %>%
  mutate(pval = pt(-abs(tstat), df = B)) %>%
  knitr::kable("latex", booktabs = TRUE,
               digits = 3, caption = paste0("True Value is ", true_mean)) %>%
  cat(., file = paste0("../Tables/calimono_t" , true_mean, ".tex"))
  
# Overall, pretty good. It looks like this is working as expected.

# ***************************
# * Nonmonotone Calibration *
# ***************************

nonmono_est_weights <- function(df) {

  gfun <- "y2"

  # Steps:
  # 1. Compute w1--the weights between U and A1.
  # 2. Compute w2--the weights between U and A2.
  # 3. Compute wc--the weights between A1, A2, U, and the core points.

  # 1. Compute w1--the weights between U and A1.
  #   a. Estimate E[g_i \mid X_i].
  #   b. Compute weights on A1.
  Eg_mod <- lm(as.formula(paste0(gfun, " ~ x")), data = df)
  df$Egx <- predict(Eg_mod)

  Mgx <- mean(df$Egx)
  n1 <- sum(df$r1)
  w1 <- Variable(n1)
  Egx_r1 <- filter(df, r1 == 1) %>% pull(Egx)
  const <- list(w1 >= 0,
                t(w1) %*% matrix(rep(1, n1), ncol = 1) == 1,
                t(w1) %*% matrix(Egx_r1, ncol = 1) == Mgx)

  p1 <- Problem(Minimize(sum_squares(w1)), const)
  res1 <- solve(p1)

  df$w1 <- 0
  df$w1[df$r1 == 1] <- res1$getValue(w1)

  # 2. Compute w2--the weights between U and A2.
  #   a. Compute weights on A2.
  n2 <- sum(df$r2)
  w2 <- Variable(n2)
  Egx_r2 <- filter(df, r2 == 1) %>% pull(Egx)
  const <- list(w2 >= 0,
                t(w2) %*% matrix(rep(1, n2), ncol = 1) == 1,
                t(w2) %*% matrix(Egx_r2, ncol = 1) == Mgx)

  p2 <- Problem(Minimize(sum_squares(w2)), const)
  res2 <- solve(p2)

  df$w2 <- 0
  df$w2[df$r2 == 1] <- res2$getValue(w2)

  # 3. Compute wc--the weights between A1, A2, U, and the core points.
  #   a. Estimate E[g_i \mid X_i, Y_{1i}].
  #   b. Estimate E[g_i \mid X_i, Y_{2i}].
  #   c. Compute core weights on each part simultaneously.
  df_1 <- filter(df, r1 == 1)
  Eg1_mod <- lm(as.formula(paste0(gfun, " ~ x + y1")), data = df_1)
  df$Eg1 <- 0
  df$Eg1[df$r1 == 1] <- predict(Eg1_mod)

  # E[y_2 | y_2 ] = y_2
  # df_2 <- filter(df, r2 == 1)
  # Eg2_mod <- lm(as.formula(paste0(gfun, " ~ x + y2")), data = df_2)

  df_11 <- filter(df, r1 == 1, r2 == 1)
  n_11 <- nrow(df_11)
  wc <- Variable(n_11)

  G1 <- sum(df$w1 * df$Eg1)
  G2 <- sum(df$w2 * df$y2)
  G3 <- mean(df$Egx)
  const <- list(wc >= 0,
                t(wc) %*% matrix(rep(1, n_11), ncol = 1) == 1,
                t(wc) %*% matrix(df_11$Eg1, ncol = 1) == G1,
                t(wc) %*% matrix(df_11$y2, ncol = 1) == G2,
                t(wc) %*% matrix(df_11$Egx, ncol = 1) == G3)

  p3 <- Problem(Minimize(sum_squares(wc)), const)
  res3 <- solve(p3)

  sum(df[[gfun]][df$r1 == 1 & df$r2 == 1] * res3$getValue(wc))

}

# **************************
# * Nonmonotone Simulation *
# **************************

run_cali_sims <- function(n_sim, B, true_mean, cor_y1y2,
                          miss_out, est_weights, seed, r_ind_y = FALSE) {


  clust <- parallel::makeCluster(min(parallelly::availableCores() - 2, 100))
  registerDoParallel(clust)

  mc_theta <- 
    foreach(b = 1:B, 
            .options.RNG = seed,
            .export = c("nonmono_mar", "ipw_cc_oracle", "est_nonmono", "expit",
                        "nonmono_est_weights"),
            .packages = c("dplyr", "CVXR")) %dorng% {

    # Parameters
    df <- nonmono_mar(n_sim, mean_y2 = true_mean, cor_y1y2 = cor_y1y2,
                      r_ind_y = r_ind_y, miss_out = miss_out)

    # Estimation
    # We want to compute:
    # 1. Oracle estimate
    # 2. pi_11 estimate (ipw using fully observed cases)
    # 3. New estimate

    tibble(oracle = mean(df$y2),
           ipworacle = ipw_cc_oracle(df),
           proposed = est_nonmono(df, oracle = !est_weights),
           calib = nonmono_est_weights(df))

  } %>% bind_rows()

  stopCluster(clust)

  # ****************
  # * Analyze Data *
  # ****************

  mc_theta %>%
   summarize(
     bias_oracle = mean(oracle) - true_mean,
     bias_ipworacle = mean(ipworacle) - true_mean,
     bias_proposed = mean(proposed) - true_mean,
     bias_calib = mean(calib) - true_mean,
     sd_oracle = sd(oracle),
     sd_ipworacle = sd(ipworacle),
     sd_proposed = sd(proposed),
     sd_calib = sd(calib),
     tstat_oracle = (mean(oracle) - true_mean) / sqrt(var(oracle) / B),
     tstat_ipworacle = (mean(ipworacle) - true_mean) / sqrt(var(ipworacle) / B),
     tstat_proposed = (mean(proposed) - true_mean) / sqrt(var(proposed) / B),
     tstat_calib = (mean(calib) - true_mean) / sqrt(var(calib) / B)) %>%
    tidyr::pivot_longer(cols = everything(),
                 names_to = c(".value", "algorithm"),
                 names_pattern = "(.*)_(.*)") %>%
    mutate(pval = pt(-abs(tstat), df = B)) 

}

# **********************************
# * Simulations: Sim 1 Nonmonotone *
# **********************************

tr_m <- -5
cy1y2 <- 0
run_cali_sims(n_sim = 1000, B = 1000,
              true_mean = tr_m, cor_y1y2 = cy1y2,
              miss_out = FALSE, est_weights = FALSE, seed = 1) %>%
  knitr::kable("latex", booktabs = TRUE, digits = 3,
  caption = paste0("True Value is ", tr_m, ". Cor(Y1, Y2) = ", cy1y2)) %>%
  cat(., file = paste0("../Tables/calinonmono1" , "m", tr_m, ".tex"))

tr_m <- 0
cy1y2 <- 0
run_cali_sims(n_sim = 1000, B = 1000,
              true_mean = tr_m, cor_y1y2 = cy1y2,
              miss_out = FALSE, est_weights = FALSE, seed = 2) %>%
  knitr::kable("latex", booktabs = TRUE, digits = 3,
  caption = paste0("True Value is ", tr_m, ". Cor(Y1, Y2) = ", cy1y2)) %>%
  cat(., file = paste0("../Tables/calinonmono1" , "m", tr_m, ".tex"))

tr_m <- 5
cy1y2 <- 0
run_cali_sims(n_sim = 1000, B = 1000,
              true_mean = tr_m, cor_y1y2 = cy1y2,
              miss_out = FALSE, est_weights = FALSE, seed = 4) %>%
  knitr::kable("latex", booktabs = TRUE, digits = 3,
  caption = paste0("True Value is ", tr_m, ". Cor(Y1, Y2) = ", cy1y2)) %>%
  cat(., file = paste0("../Tables/calinonmono1" , "m", tr_m, ".tex"))










# ****************
# * Clean Tables *
# ****************
# The tables look better if they have the [h!] attribute, so we correct them.
# To be run from Simulations/
source("clean_tables.R")

