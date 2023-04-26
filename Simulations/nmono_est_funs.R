# Title: ../Simulations/nmono_est_funs.R
# Author: Caleb Leedy
# Date Created: April 25, 2023
# Purpose: This script contains the functions used to estimate monotone and
# nonmonotone missingness. Like Simulations/generate_data.R, I wanted to
# separate the function estimation from the simulation study itself in order to
# gain more modular code.

# ------------------------------------------------------------------------------

# ********************
# * Script Functions *
# ********************
# A. Monotone estimators
#   1. est_mono -- Classical monotone estimator
#   2. ipw_mono_est -- IPW estimator with estimated weights
# B. Nonmonotone estimators
#   1. est_nonmono -- The proposed estimator in the nonmonotone setting
#   2. ipw_cc_oracle -- IPW estimator with oracle weights
#   3. ipw_cc_est -- IPW estimator with weights estimated
#   4. twophase_reg -- Unfinished

# ------------------------------------------------------------------------------

# **********************
# * Required Libraries *
# **********************

require(dplyr)

# ------------------------------------------------------------------------------

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

# **************************************
# * Estimating Nonmonotone Missingness *
# **************************************

est_nonmono <- function(df, oracle = FALSE, alpha = NA) {

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
    pi_12 <- pi_11_mod$fitted.values # NOTE: This is NOT of length n.

    df$pi_11 <- 0
    df$pi_11[df$r1 == 1] <- pi_12

    if (!is.na(alpha)) {
      pi_21_m <- glm(r1 ~ x + y2, data = df_01, family = binomial(link = logit))
      pi_21 <- pi_21_m$fitted.values
      tmp12_vec <- rep(0, nrow(df))
      tmp21_vec <- rep(0, nrow(df))
      tmp12_vec[df$r1 == 1] <- pi_12
      tmp21_vec[df$r2 == 1] <- pi_21
      df$pi_11 <- (1 - alpha) * tmp12_vec + (alpha) * tmp21_vec
    }

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

ipw_cc_est <- function(df, oracle = FALSE, alpha = NA) {

  gfun <- "y2"


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
    pi_12 <- pi_11_mod$fitted.values # NOTE: This is NOT of length n.

    df$pi_11 <- 0
    df$pi_11[df$r1 == 1] <- pi_12

    if (!is.na(alpha)) {
      pi_21_m <- glm(r1 ~ x + y2, data = df_01, family = binomial(link = logit))
      pi_21 <- pi_21_m$fitted.values
      tmp12_vec <- rep(0, nrow(df))
      tmp21_vec <- rep(0, nrow(df))
      tmp12_vec[df$r1 == 1] <- pi_12
      tmp21_vec[df$r2 == 1] <- pi_21
      df$pi_11 <- (1 - alpha) * tmp12_vec + (alpha) * tmp21_vec
    }
  }

  df$r1r2pi_11 <- df$r1 * df$r2 / df$pi_11
  df$r1r2pi_11 <- ifelse(is.na(df$r1r2pi_11), 0, df$r1r2pi_11)

  mean(df[[gfun]] * df$r1r2pi_11)

}

# png("diffhist.png")
# filter(df, r1 == 1, r2 == 1) %>%
#   dplyr::select(-p0, -p1, -p2, -p12, -p21) %>%
#   mutate(diff = pi_11 - true_pi11) %>%
#   pull(diff) %>% 
#   hist(main = "Histogram of Estimated - True") 
# dev.off()

twophase_reg <- function(df) {

  gfun <- "y2"


}

