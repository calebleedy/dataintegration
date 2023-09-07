# Title: Comparing_Estimator_Types/generate_data.R
# Author: Caleb Leedy
# Date Created: April 25, 2023
# Purpose: This file contains the code to generate data for different simulation
# studies. Instead of writing the functions in Simulations/nonmonotone.R and
# then copying them to other files, I thought that it would be a good idea to
# put them all into one script and then source it.

# ------------------------------------------------------------------------------

# ********************
# * Script Functions *
# ********************
# 1. mono_mar -- This generates monotone missingness.
# 2. nonmono_mar -- This generates nonmonotone missing data.
# 3. expit -- A helper function for generating probabilites. 

# ------------------------------------------------------------------------------

# **********************
# * Required Libraries *
# **********************

library(MASS)
library(dplyr)

# ------------------------------------------------------------------------------

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
    mutate(p1 = expit(x_vec)) %>%
    mutate(p12 = expit(y1)) %>%
    mutate(p21 = 0) %>%
    mutate(p2 = 0) %>%
    mutate(pi_11 = p1 * p12) %>%
    mutate(r2 = ifelse(r1 == 0, 0, r2))

}

nonmono_mar <- 
  function(n, mean_y2 = 0, cor_xy1 = 0, cor_xy2 = 0, cor_y1y2 = 0,
           r_ind_y = FALSE, miss_out = FALSE, miss_resp = FALSE, mcar = FALSE) {

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

      if (mcar) {
        p1 <- 0.4 
        p2 <- 0.4
      } else {
        p1 <- expit(2 * x_vec)
        p2 <- expit(-2 * x_vec)
      }
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
      # p1 <- expit(2 * x_vec)
      p1 <- 0.4
      # p2 <- expit(-2 * x_vec)
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

