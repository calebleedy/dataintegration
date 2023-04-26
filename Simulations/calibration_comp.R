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
