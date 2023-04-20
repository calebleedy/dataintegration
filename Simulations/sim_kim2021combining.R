# Title: sim_kim2021combining.R
# Author: Caleb Leedy
# Date Created: April 20, 2023
# Purpose: This document contains the code to reproduce the simulation study of
# Kim et al. (2021) (see ../Reference/kim2021combining.pdf).

# *************
# * Libraries *
# *************

library(dplyr)
library(purrr)

library(doParallel)
library(parallelly)
library(doRNG)

# ***************************
# * Data Generating Process *
# ***************************

gen_data <- function(n_B, model_num) {

  # Finite population
  N <- 100000
  x_vec <- rnorm(N, mean = 2, sd = 1)
  e_vec <- rnorm(N)

  if (model_num == 1) {
    y_vec <- 1 + 2 * x_vec + e_vec
  } else if (model_num == 2) {
    y_vec <- 3 + x_vec + 2 * e_vec
  } else if (model_num == 3) {
    y_vec <- 2.5 + 0.5 * x_vec^2 + e_vec
  } else {
    stop("model_num must be 1, 2, or 3")
  }

  pop_df <- tibble(X = x_vec, Y = y_vec)

  # Sample A
  n_A <- 500
  A_df <- slice_sample(pop_df, n = n_A)

  # Sample B
  strat1 <- filter(pop_df, X <= 2)
  strat2 <- filter(pop_df, X > 2)
  B_df <- bind_rows(slice_sample(strat1, n = 0.7 * n_B),
                    slice_sample(strat2, n = 0.3 * n_B))


  return(list(A = A_df, B = B_df))

}

# *************************
# * Estimation Algorithms *
# *************************

mean_A <- function(data) {
  mean(data$A$Y)
}

mean_B <- function(data) {
  mean(data$B$Y)
}

mass_imp <- function(data) {

  imp_mod <- lm(Y ~ X, data = data$B)
  mean(predict(imp_mod, data$A))

}

chen_ipw <- function(data) {

  A_df <- data$A
  B_df <- data$B

  return(1)
  # FIXME:

}

# **************************
# * Monte Carlo Experiment *
# **************************

run_sim <- function(n_B, model_number, seed) {

  clust <- parallel::makeCluster(min(parallelly::availableCores() - 2, 100))
  registerDoParallel(clust)

  B <- 1000

  res_df <- 
    foreach(iter = 1:B,
            .packages = c("dplyr"),
            .export = c("mean_A", "mean_B", "mass_imp", "chen_ipw", "gen_data"),
            .options.RNG = seed) %dorng% {

      data <- gen_data(n_B, model_num = model_number)
      tibble(A_mean = mean_A(data),
             B_mean = mean_B(data),
             MassImp = mass_imp(data),
             IPW_Chen = chen_ipw(data))

  } %>% bind_rows()

  stopCluster(clust)

  res_df %>%
    summarize(bias.A = mean(A_mean) - 5,
              bias.B = mean(B_mean) - 5,
              bias.massimp = mean(MassImp) - 5,
              bias.ipw = mean(IPW_Chen) - 5,
              var.A = var(A_mean),
              var.B = var(B_mean),
              var.massimp = var(MassImp),
              var.ipw = var(IPW_Chen)) %>%
  mutate(mseA = bias.A^2 + var.A) %>%
  mutate(ReMSE.A = mseA / mseA,
         ReMSE.B = (bias.B^2 + var.B) / mseA,
         ReMSE.massimp = (bias.massimp^2 + var.massimp) / mseA,
         ReMSE.ipw = (bias.ipw^2 + var.ipw) / mseA) %>%
  select(-mseA) %>%
    tidyr::pivot_longer(cols = everything(),
                 names_to = c(".value", "algorithm"),
                 names_pattern = "(.*)\\.(.*)") 


}

run_sim(500, 1, seed = 1)
run_sim(500, 2, seed = 1)
run_sim(500, 3, seed = 1)
run_sim(1000, 1, seed = 1)
run_sim(1000, 2, seed = 1)
run_sim(1000, 3, seed = 1)

