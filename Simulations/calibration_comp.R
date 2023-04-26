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

