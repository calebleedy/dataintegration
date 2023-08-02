# Title: regression_est.R
# Author: Caleb Leedy
# Date Created: May  5, 2023
# Purpose: This file runs some code that I use for better understanding of the
# regression estimator.

# ***********
# * Outline *
# ***********
# 1. Setup
#   a. Libraries
#   b. Additional Code
# 2. Least Squares Estimator
# 3. Least Squares Estimator without Intercept
# 4. Regression Estimator

# ------------------------------------------------------------------------------

# *********
# * Setup *
# *********

# *************
# * Libraries *
# *************

library(dplyr)
library(CVXR)

# *******************
# * Additional Code *
# *******************

source("generate_data.R")
source("nmono_est_funs.R")

df <- mono_mar(n = 1000)
df_2 <- filter(df, r2 == 1)

# ***************************
# * Least Squares Estimator *
# ***************************
ls_est <- lm(y2 ~ x, data = df_2)

# *********************************************
# * Least Squares Estimator without Intercept *
# *********************************************
ls0_est <- lm(y2 ~ 0 + x, data = df_2)

# ************************
# * Regression Estimator *
# ************************
mean(df_2$y2) + (mean(df$x) - mean(df_2$x)) * ls_est$coefficients[2]

# ***********
# * Scratch *
# ***********

lm(y2 ~ x, data = df_2)
lm(y2 ~ I(x - mean(df$x)), data = df_2)
lm(y2 ~ 0 + I(x - mean(df$x)), data = df_2)
lm(y2 ~ 0 + x, data = df_2)

mean(df_2$y2) + (mean(df$x) - mean(df_2$x)) * ls0_est$coefficients[1]




