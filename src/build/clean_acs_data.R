# File: clean_acs_data.R
# Created by: Caleb Leedy
# Created on: October 28, 2024
# Purpose:
# This document cleans the ACS data so that it can be more easily used. The data
# set is quite large and even simply tidyverse operations cause R to close.

# *************
# * Libraries *
# *************

library(dplyr)
library(readr)
library(stringr)

# ********
# * Data *
# ********

acs_a <- read_csv("data/ACS/psam_pusa.csv")
acs_b <- read_csv("data/ACS/psam_pusb.csv")

# *********************
# * Cleaning the Data *
# *********************

clean_acs_data <- function(df) {

  df <- df %>%
    mutate(DEM_SEX = as.character(SEX))

  # Clean Education
  df <- df %>%
    mutate(DEM_EDU = case_when(
      SCHL == "bb" ~ "1",
      SCHL == "01" ~ "1", 
      SCHL == "02" ~ "1", 
      SCHL == "03" ~ "1", 
      SCHL == "04" ~ "1", 
      SCHL == "05" ~ "1", 
      SCHL == "06" ~ "1", 
      SCHL == "07" ~ "1", 
      SCHL == "08" ~ "1", 
      SCHL == "09" ~ "1", 
      SCHL == "10" ~ "1", 
      SCHL == "11" ~ "1", 
      SCHL == "12" ~ "1", 
      SCHL == "13" ~ "1", 
      SCHL == "14" ~ "1", 
      SCHL == "15" ~ "1", 
      SCHL == "16" ~ "2", 
      SCHL == "17" ~ "2", 
      SCHL == "18" ~ "2", 
      SCHL == "19" ~ "2", 
      SCHL == "20" ~ "3", 
      SCHL == "21" ~ "3", 
      SCHL == "22" ~ "3", 
      SCHL == "23" ~ "3", 
      SCHL == "24" ~ "3", 
      TRUE ~ "M" # There are some NAs in the data
    ))

  # Clean marital statu
  df <- df %>%
    mutate(DEM_MARSTA = case_when(
      MAR == 1 ~ "1",
      MAR == 2 ~ "2",
      MAR == 3 ~ "3",
      MAR == 4 ~ "3",
      MAR == 5 ~ "4",
    ))

  # Clean race / ethnicity
  # HISP = 01 for non-Hispanic
  df <- df %>%
    mutate(DEM_RACE = case_when(
      (RAC1P == 1 & HISP == "01") ~ 1, # Non-Hispanic white
      (RAC1P == 2 & HISP == "01") ~ 2, # Non-Hispanic Black
      HISP != "01" ~ 3,                # Hispanic
      TRUE ~ 4))

  # Clean Medicare
  # The Medicare ACS variance HINS3 should actually be used to filter this data
  # set since the MCBS data focuses on people on Medicare--not people who are 65
  # years and older.
  
  # Clean Medicaid
  # We actually do not want to balance on Medicaid because the definitions are
  # unclear about how to match.

  return(select(df, SERIALNO, SPORDER, AGEP, HINS3, starts_with("PWGTP"),
                    DEM_SEX, DEM_EDU, DEM_MARSTA, DEM_RACE))
}

acs_a_cl <- clean_acs_data(acs_a)
acs_b_cl <- clean_acs_data(acs_b)

acs_df <- bind_rows(acs_a_cl, acs_b_cl)

# *************
# * Save Data *
# *************

write_csv(acs_df, "cleaned_acs_data.csv")
