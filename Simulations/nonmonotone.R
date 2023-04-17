# Title: nonmonotone.R
# Author: Caleb Leedy
# Date Created: March 29, 2023
# Purpose: This document contains the simulations for computing estimating
# equations with both monotone and nonmonotone missingness. See
# Notes/Non-monotone.pdf for more information.
#
# Questions:
# 1. What estimating functions should we use? 
# - U(\theta, Y) = (\bar y_1 + \bar y_2) / 2 - \theta
# - U(\theta, Y) = \bar y_1 + \bar y_2 - \theta
# - U(\theta, Y) = \bar y_1 - \bar y_2 - \theta
# Other ideas include moments of a distribution of quantiles.
#
# 2. How should nonmonotone missingness data be generated?
# In the monotone setting we think about data being generated conditionally (or
# at least the missingness is conditionally generated). This is because we want
# to have data MAR instead of MCAR. The challenge in this is that we also have
# to think about how to ensure that 
#
# \[p(R_1, R_2 | X) = p(R_1 | R_2, X) p(R_2 | X) = p(R_2 | R_1, X) p(R_1 | X).\]
# 
# A further challenge is that these conditional probabilites can depend on Y_1
# and Y_2 if they have been observed. However, we cannot predetermine if an
# observation should have wlog R_2 depend on Y_1 because we don't know if Y_1 is
# observed for that element.
#
# Answers:
# 1. Currently, we will use U(\theta, Y) = \bar y_2 - \theta.
#
# 2. While one cannot have a global model this is MAR, using the ideas of Robins
# and Gill (1998), we can construct an MAR model where each probability is
# conditional.

# *************
# * Libraries *
# *************

library(MASS)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)

library(doParallel)
library(doRNG)
library(parallelly)

# ************************
# * Additional Functions *
# ************************

expit <- function(x) {1 / (1 + exp(-x))}

# *****************
# * Generate Data *
# *****************

mono_mar <- function(n, mean_y2 = 0) {

  # Steps:
  # 1. Generate X, Y_1, Y_2.
  # 2. Get probabilities conditional on past.
  # 3. Reveal hidden variables.

  # 1. Generate X, Y_1, Y_2.
  x_vec <- rnorm(n)
  y1 <- rnorm(n)
  y2 <- rnorm(n, mean = mean_y2)

  # 2. Get probabilities conditional on past.
  r1 <- rbinom(n, 1, expit(x_vec))
  r2 <- rbinom(n, 1, expit(y1))

  # 3. Reveal hidden variables.
  tibble(x = x_vec, y1 = y1, y2 = y2, r1 = r1, r2 = r2) %>%
    mutate(r2 = ifelse(r1 == 0, 0, r2))

}

nonmono_mar <- 
  function(n, mean_y2 = 0, cor_xy1 = 0, cor_xy2 = 0, cor_y1y2 = 0,
           r_ind_y = FALSE, miss_out = FALSE, miss_resp = FALSE) {

  # Steps:
  # 1. Generate X, Y_1, Y_2.
  # 2. Get probabilities conditional on past.
  # 3. Reveal hidden variables.

  # 1. Generate X, Y_1, Y_2.
  mean_vec <- c(0, 0, mean_y2)
  Sigma_mat <- 
    matrix(c(1, cor_xy1, cor_xy2, cor_xy1, 1, cor_y1y2, cor_xy2, cor_y1y2, 1),
           nrow = 3)
  data <- MASS::mvrnorm(n, mean_vec, Sigma_mat)
  x_vec <- data[, 1]
  y1 <- x_vec + data[, 2]
  y2 <- x_vec + data[, 3]

  # Misspecified outcome model
  if (miss_out) {
    y1 <- x_vec + x_vec^2 + data[, 2]
    y2 <- -x_vec + x_vec^3 + data[, 3]
  }

  # 2. Get probabilities conditional on past.
  if (r_ind_y) {

    # We need:
    # p0, p1, p2, p12, p21, out

    if (miss_resp) {

      p0 <- abs(x_vec^2 - 1)
      p1 <- abs(x_vec^2 - 2)
      p2 <- abs(x_vec^2)

      tot <- p0 + p1 + p2
      p0 <- p0 / tot
      p1 <- p1 / tot
      p2 <- p2 / tot

      p12 <- expit(x_vec^2)
      p21 <- expit(x_vec^2)


    } else {

      p0 <- 0.2
      p1 <- expit(1 + 2 * x_vec)
      p2 <- expit(-1 + -2 * x_vec)
      tot <- p0 + p1 + p2
      p0 <- p0 / tot
      p1 <- p1 / tot
      p2 <- p2 / tot

      p12 <- expit(x_vec)
      p21 <- expit(x_vec)

    }

    rand1 <- runif(n)
    rand2 <- runif(n)

  } else {

    # We want: 
    # p0 = Pr(R1 = 0, R2 = 0)
    # p1 = Pr(R1 = 1) # NOTE: This is NOT true because R1 = 1 if p21 works.
    # p2 = Pr(R2 = 1) # NOTE: This is NOT true because R2 = 1 if p12 works.
    # p12 = Pr(R2 = 1 | R1 == 1)
    # p21 = Pr(R1 = 1 | R2 == 1)
    # Instead, p0, p1, and p2 are Pr(i in Group(r)) which determines the order
    # which we fill missing values. Thus,
    # Pr(R1 = 1) = p1 + p21 * p2
    # Pr(R2 = 1) = p2 + p12 * p1

    if (miss_resp) {

      p0 <- abs(x_vec^2 - 1)
      p1 <- abs(x_vec^2 - 2)
      p2 <- abs(x_vec^2)

      tot <- p0 + p1 + p2
      p0 <- p0 / tot
      p1 <- p1 / tot
      p2 <- p2 / tot

      p12 <- expit(y1^2)
      p21 <- expit(y2^2)


    } else {

      p0 <- 0.2
      # p1 <- expit(1 + 2 * x_vec)
      p1 <- 0.4
      # p2 <- expit(-1 + -2 * x_vec)
      p2 <- 0.4
      tot <- p0 + p1 + p2
      p0 <- p0 / tot
      p1 <- p1 / tot
      p2 <- p2 / tot

      p12 <- expit(y1)
      p21 <- expit(y2)

    }

    rand1 <- runif(n)
    rand2 <- runif(n)

  }

  # 3. Reveal hidden variables.
  tibble(x = x_vec, y1 = y1, y2 = y2,
         p0 = p0, p1 = p1, p2 = p2,
         p12 = p12, p21 = p21,
         rand1 = rand1, rand2 = rand2) %>%
  mutate(out = case_when(
                rand1 < p0 ~ "00",
                rand1 < p0 + p1 & rand2 < 1 - p12 ~ "10",
                rand1 < p0 + p1 & rand2 > 1 - p12 ~ "11a",
                rand1 < p0 + p1 + p2 & rand2 < 1 - p21 ~ "01",
                rand1 < p0 + p1 + p2 & rand2 > 1 - p21 ~ "11b")) %>%
  mutate(r1 = ifelse(out %in% c("10", "11a", "11b"), 1, 0)) %>%
  mutate(r2 = ifelse(out %in% c("01", "11a", "11b"), 1, 0)) %>%
  dplyr::select(-rand1, -rand2)

}

# ***********************************
# * Estimating Monotone Missingness *
# ***********************************

est_mono <- function(df) {

  # g(x, y1, y2) = y2
  gfun <- "y2"

  # Algorithm Steps:
  # 1. Estimate pi_11, pi_1+
  # 2. Estimate b1 and b2
  # 3. Estimate \hat \theta

  df_1 <- filter(df, r1 == 1)  # Contains data from first stage
  df_11 <- filter(df, r2 == 1) # Contains data from first and second stage

  # 1. Estimate pi_11, pi_1+
  pi_11_mod <- glm(r2 ~ x + y1, data = df_1, family = binomial(link = logit))
  pi_11 <- pi_11_mod$fitted.values # NOTE: This is NOT of length n.

  df$pi_11 <- 0
  df$pi_11[df$r1 == 1] <- pi_11

  pi_1_mod <- glm(r1 ~ x, data = df, family = binomial(link = logit))
  pi_1 <- pi_1_mod$fitted.values

  # 2. Estimate b1 and b2
  b1_mod <- lm(as.formula(paste0(gfun, " ~ x")), data = df_11)
  b1 <- predict(b1_mod, newdata = df)
  b2_mod <- lm(as.formula(paste0(gfun, " ~ x + y1")), data = df_11)
  b2 <- predict(b2_mod, newdata = df_1) # NOTE: This is NOT of length n.
  df$b2 <- 0 
  df$b2[df$r1 == 1] <- b2


  # 3. Estimate \hat \theta
  df$r1r2pi_11 <- df$r1 * df$r2 / df$pi_11
  df$r1r2pi_11 <- ifelse(is.na(df$r1r2pi_11), 0, df$r1r2pi_11)
  theta <- 
    mean(df$r1r2pi_11 * df[[gfun]]) -
    mean((df$r1 / pi_1 - 1) * b1) - 
    mean((df$r1r2pi_11 - df$r1 / pi_1) * df$b2)

}

ipw_mono_est <- function(df, est = FALSE) {

  if (est) {
    pi_1_mod <- glm(r1 ~ x, data = df, family = binomial(link = logit))
    pi_1 <- pi_1_mod$fitted.values

    pi_11_mod <- 
      glm(r2~x+y1, data = filter(df, r1 == 1), family = binomial(link = logit))
    pi_11 <- numeric(nrow(df))
    pi_11[df$r1 == 1] <- pi_11_mod$fitted.values
    pi_11[df$r1 != 1] <- Inf

  } else {

    pi_1 <- expit(df$x)
    pi_11 <- expit(df$y1)

  }

  mean(df$y2 * df$r2 / (pi_1 * pi_11))

}

# ****************************************
# * Run Simulation: Monotone Missingness *
# ****************************************

clust <- parallel::makeCluster(min(parallelly::availableCores() - 2, 100))
registerDoParallel(clust)

B <- 1000
true_mean <- -5

mc_theta <- foreach(b = 1:B, 
                  #.options.RNG = 1,
                  .packages = c("dplyr")) %dorng% {

  # Parameters
  n_sim <- 1000

  df <- mono_mar(n_sim, mean_y2 = true_mean)
  tibble(oracle = mean(df$y2),
         ipworacle = ipw_mono_est(df, est = FALSE),
         ipwest = ipw_mono_est(df, est = TRUE),
         semi = est_mono(df))
         
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
   sd_oracle = sd(oracle),
   sd_ipworacle = sd(ipworacle),
   sd_ipwest = sd(ipwest),
   sd_semi = sd(semi),
   tstat_oracle = (mean(oracle) - true_mean) / sqrt(var(oracle) / B),
   tstat_ipworacle = (mean(ipworacle) - true_mean) / sqrt(var(ipworacle) / B),
   tstat_ipwest = (mean(ipwest) - true_mean) / sqrt(var(ipwest) / B),
   tstat_semi = (mean(semi) - true_mean) / sqrt(var(semi) / B)) %>%
  tidyr::pivot_longer(cols = everything(),
               names_to = c(".value", "algorithm"),
               names_pattern = "(.*)_(.*)") %>%
  mutate(pval = pt(-abs(tstat), df = B)) %>%
  knitr::kable("latex", booktabs = TRUE,
               digits = 3, caption = paste0("True Value is ", true_mean)) %>%
  cat(., file = paste0("../Tables/monomar_t" , true_mean, ".tex"))
  

# **************************************
# * Estimating Nonmonotone Missingness *
# **************************************

est_nonmono <- function(df, oracle = FALSE) {

  # g(x, y1, y2) = y2
  gfun <- "y2"

  # Algorithm Steps:
  # 1. Estimate pi_11, pi_1+, pi_2+
  # 2. Estimate b1 and b2, a2
  # 3. Estimate \hat \theta

  df_10 <- filter(df, r1 == 1)          # Contains data from first stage
  df_11 <- filter(df, r1 == 1, r2 == 1) # Contains data from both stages
  df_01 <- filter(df, r2 == 1)          # Contains data from first stage

  # 1. Estimate pi_11, pi_1+, pi_2+
  if (oracle) {
    pi_1 <- df$p1 + df$p2 * df$p21
    pi_2 <- df$p2 + df$p1 * df$p12
    pi_11 <- df$p12 * df$p1 + df$p21 * df$p2
    df$pi_11 <- pi_11

  } else {
    pi_11_mod <- glm(r2 ~ x + y1, data = df_10, family = binomial(link = logit))
    pi_11 <- pi_11_mod$fitted.values # NOTE: This is NOT of length n.

    df$pi_11 <- 0
    df$pi_11[df$r1 == 1] <- pi_11

    # NOTE: We do not estimate
    # glm(r1 ~ x + y2, data = df_01, family = binomial(link = logit))

    pi_1_mod <- glm(r1 ~ x, data = df, family = binomial(link = logit))
    pi_1 <- pi_1_mod$fitted.values

    pi_2_mod <- glm(r2 ~ x, data = df, family = binomial(link = logit))
    pi_2 <- pi_2_mod$fitted.values
  }

  # 2. Estimate b1 and b2
  b1_mod <- lm(as.formula(paste0(gfun, " ~ x")), data = df_11)
  b1 <- predict(b1_mod, newdata = df)
  b2_mod <- lm(as.formula(paste0(gfun, " ~ x + y1")), data = df_11)
  b2 <- predict(b2_mod, newdata = df_10) # NOTE: This is NOT of length n.
  df$b2 <- 0 
  df$b2[df$r1 == 1] <- b2

  # a2_mod <- lm(as.formula(paste0(gfun, " ~ x + y2")), data = df_11)
  # a2 <- predict(a2_mod, newdata = df_01) # NOTE: This is NOT of length n.
  df$a2 <- df$y2
  df$a2[df$r2 == 0] <- 0

  # 3. Estimate \hat \theta
  df$r1r2pi_11 <- df$r1 * df$r2 / df$pi_11
  df$r1r2pi_11 <- ifelse(is.na(df$r1r2pi_11), 0, df$r1r2pi_11)

  theta <-
    mean(b1) +
    mean(df$r1 / pi_1 * (df$b2 - b1)) +
    mean(df$r2 / pi_2 * (df$a2 - b1)) +
    mean(df$r1r2pi_11 * (df[[gfun]] - df$b2 - df$a2 + b1))

}

# **************************
# * Alternative Algorithms *
# **************************

ipw_cc_oracle <- function(df) {

  gfun <- "y2"

  tmp <- df %>%
    filter(r1 == 1, r2 == 1) %>%
    mutate(pi_11 = p12 * p1 + p21 * p2) 

  sum(tmp[[gfun]] / tmp$pi_11) / nrow(df)

}

ipw_cc_est <- function(df) {

  gfun <- "y2"

  # Steps:
  # 1. Likelihood of being in G(1)
  # 2. Likelihood of being in G(2)
  # 3. Relative likelihood of each. This is a weight.
  dist_10 <-
    filter(df, out == "10") %>%
    summarise(mean = mean(y1),
              sd = sd(y1))
  dist_01 <-
    filter(df, out == "01") %>%
    summarise(mean = mean(y2), sd = sd(y2))

  # TODO: This is wrong, but I don't know how to use the weight to estimate the
  # probability.
  tmp <- df %>%
    mutate(like_10 = dnorm(y1, mean = dist_10$mean, sd = dist_10$sd)) %>%
    mutate(like_01 = dnorm(y2, mean = dist_01$mean, sd = dist_01$sd)) %>%
    mutate(w_10 = like_10 / (like_10 + like_01)) %>%
    mutate(w_10 = case_when(r1 == 1 & r2 == 1 ~ w_10,
                            r1 == 1 & r2 == 0 ~ 1,
                            r1 == 0 & r2 == 1 ~ 0, 
                            r1 == 0 & r2 == 0 ~ NA_real_)) 


}

twophase_reg <- function(df) {

  gfun <- "y2"


}

# *******************************************
# * Run Simulation: Nonmonotone Missingness *
# *******************************************

run_init_sims <- function(n_sim, B, true_mean, cor_y1y2,
                          miss_out, oracle, seed) {


  clust <- parallel::makeCluster(min(parallelly::availableCores() - 2, 100))
  registerDoParallel(clust)

  mc_theta <- 
    foreach(b = 1:B, 
            .options.RNG = seed,
            .export = c("nonmono_mar", "ipw_cc_oracle", "est_nonmono", "expit"),
            .packages = c("dplyr")) %dorng% {

    # Parameters
    df <- nonmono_mar(n_sim, mean_y2 = true_mean, cor_y1y2 = cor_y1y2,
                      r_ind_y = FALSE, miss_out = miss_out)

    # Estimation
    # We want to compute:
    # 1. Oracle estimate
    # 2. pi_11 estimate (ipw using fully observed cases)
    # 3. New estimate

    tibble(oracle = mean(df$y2),
           ipworacle = ipw_cc_oracle(df),
           proposed = est_nonmono(df, oracle = oracle))

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
     sd_oracle = sd(oracle),
     sd_ipworacle = sd(ipworacle),
     sd_proposed = sd(proposed),
     tstat_oracle = (mean(oracle) - true_mean) / sqrt(var(oracle) / B),
     tstat_ipworacle = (mean(ipworacle) - true_mean) / sqrt(var(ipworacle) / B),
     tstat_proposed = (mean(proposed) - true_mean) / sqrt(var(proposed) / B))%>%
    tidyr::pivot_longer(cols = everything(),
                 names_to = c(".value", "algorithm"),
                 names_pattern = "(.*)_(.*)") %>%
    mutate(pval = pt(-abs(tstat), df = B)) 

}

# **********************************
# * Simulations: Sim 1 Nonmonotone *
# **********************************

run_init_sims(n_sim = 2000, B = 2000,
              true_mean = -5, cor_y1y2 = 0,
              miss_out = FALSE, oracle = TRUE, seed = 1) %>%
  knitr::kable("latex", booktabs = TRUE, digits = 3,
  caption = "True Value is -5. Cor(Y1, Y2) = 0") %>%
  cat(., file = paste0("../Tables/nonmonosim1" , "m-5", ".tex"))

run_init_sims(n_sim = 2000, B = 2000,
              true_mean = 0, cor_y1y2 = 0,
              miss_out = FALSE, oracle = TRUE, seed = 1) %>%
  knitr::kable("latex", booktabs = TRUE, digits = 3,
  caption = "True Value is 0. Cor(Y1, Y2) = 0") %>%
  cat(., file = paste0("../Tables/nonmonosim1" , "m0", ".tex"))

run_init_sims(n_sim = 2000, B = 2000,
              true_mean = 5, cor_y1y2 = 0,
              miss_out = FALSE, oracle = TRUE, seed = 1) %>%
  knitr::kable("latex", booktabs = TRUE, digits = 3,
  caption = "True Value is 5. Cor(Y1, Y2) = 0") %>%
  cat(., file = paste0("../Tables/nonmonosim1" , "m5", ".tex"))

# **********************************
# * Simulations: Sim 2 Nonmonotone *
# **********************************

run_init_sims(n_sim = 2000, B = 2000,
              true_mean = 0, cor_y1y2 = 0.1,
              miss_out = FALSE, oracle = TRUE, seed = 2) %>%
  knitr::kable("latex", booktabs = TRUE, digits = 3,
  caption = "True Value is 0. Cor(Y1, Y2) = 0.1") %>%
  cat(., file = paste0("../Tables/nonmonosim2" , "c0.1", ".tex"))

run_init_sims(n_sim = 2000, B = 2000,
              true_mean = 0, cor_y1y2 = 0.5,
              miss_out = FALSE, oracle = TRUE, seed = 2) %>%
  knitr::kable("latex", booktabs = TRUE, digits = 3,
  caption = "True Value is 0. Cor(Y1, Y2) = 0.5") %>%
  cat(., file = paste0("../Tables/nonmonosim2" , "c0.5", ".tex"))

run_init_sims(n_sim = 2000, B = 2000,
              true_mean = 0, cor_y1y2 = 0.9,
              miss_out = FALSE, oracle = TRUE, seed = 2) %>%
  knitr::kable("latex", booktabs = TRUE, digits = 3,
  caption = "True Value is 0. Cor(Y1, Y2) = 0.9") %>%
  cat(., file = paste0("../Tables/nonmonosim2" , "c0.9", ".tex"))

# **********************************
# * Simulations: Sim 3 Nonmonotone *
# **********************************

run_init_sims(n_sim = 2000, B = 2000,
              true_mean = 0, cor_y1y2 = 0,
              miss_out = TRUE, oracle = TRUE, seed = 3) %>%
  knitr::kable("latex", booktabs = TRUE, digits = 3,
  caption = "True Value is 0. Cor(Y1, Y2) = 0") %>%
  cat(., file = paste0("../Tables/nonmonosim3" , "c0", ".tex"))

run_init_sims(n_sim = 2000, B = 2000,
              true_mean = 0, cor_y1y2 = 0.1,
              miss_out = TRUE, oracle = TRUE, seed = 3) %>%
  knitr::kable("latex", booktabs = TRUE, digits = 3,
  caption = "True Value is 0. Cor(Y1, Y2) = 0.1") %>%
  cat(., file = paste0("../Tables/nonmonosim3" , "c0.1", ".tex"))

run_init_sims(n_sim = 2000, B = 2000,
              true_mean = 0, cor_y1y2 = 0.5,
              miss_out = TRUE, oracle = TRUE, seed = 3) %>%
  knitr::kable("latex", booktabs = TRUE, digits = 3,
  caption = "True Value is 0. Cor(Y1, Y2) = 0.5") %>%
  cat(., file = paste0("../Tables/nonmonosim3" , "c0.5", ".tex"))

# **********************************
# * Simulations: Sim 4 Nonmonotone *
# **********************************

run_init_sims(n_sim = 2000, B = 2000,
              true_mean = 0, cor_y1y2 = 0,
              miss_out = FALSE, oracle = FALSE, seed = 4) %>%
  knitr::kable("latex", booktabs = TRUE, digits = 3,
  caption = "True Value is 0. Cor(Y1, Y2) = 0") %>%
  cat(., file = paste0("../Tables/nonmonosim4" , "c0", ".tex"))

run_init_sims(n_sim = 2000, B = 2000,
              true_mean = 0, cor_y1y2 = 0.1,
              miss_out = FALSE, oracle = FALSE, seed = 4) %>%
  knitr::kable("latex", booktabs = TRUE, digits = 3,
  caption = "True Value is 0. Cor(Y1, Y2) = 0.1") %>%
  cat(., file = paste0("../Tables/nonmonosim4" , "c0.1", ".tex"))

run_init_sims(n_sim = 2000, B = 2000,
              true_mean = 0, cor_y1y2 = 0.5,
              miss_out = FALSE, oracle = FALSE, seed = 4) %>%
  knitr::kable("latex", booktabs = TRUE, digits = 3,
  caption = "True Value is 0. Cor(Y1, Y2) = 0.5") %>%
  cat(., file = paste0("../Tables/nonmonosim4" , "c0.5", ".tex"))

# %>%
#   knitr::kable("latex", booktabs = TRUE, digits = 3,
#   caption = paste0("True Value is ", true_mean, ". Cor(Y1, Y2) = ", cor_y1y2)) %>%
#   cat(., file = paste0("../Tables/nonmonocmar_t" , true_mean,
#                        "cy1y2", cor_y1y2, ".tex"))
  


# ****************
# * Clean Tables *
# ****************
# The tables look better if they have the [h!] attribute, so we correct them.
# To be run from Simulations/
source("clean_tables.R")

