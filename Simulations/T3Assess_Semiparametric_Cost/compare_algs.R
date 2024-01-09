# Title: compare_algs.R
# Created by: Caleb Leedy
# Created on: January 09, 2024
# Purpose: This file contains the code to compare the parametric and 
# semiparametric estimators.

# *************
# * Libraries *
# *************

library(dplyr)

# **************
# * Algorithms *
# **************

#' The para_alg is the parametric algorithm that should be optimal when the 
#' model is correctly specified.
#'
#' @param df - A data frame generated from R/opt_est.R::gen_optsim_data.
#' @param gfun - A string
para_alg <- function(df, gfun = "Y2") {

  df <- mutate(df, g_i = eval(rlang::parse_expr(gfun)))

  # We use the population functional forms for gamma_hat and v_gamma

} 
out_rob_alg
resp_rob_alg
doub_rob_alg
