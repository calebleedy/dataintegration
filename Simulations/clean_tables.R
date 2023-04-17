# Title: clean_tables.R
# Author: Caleb Leedy
# Date Created: April 17, 2023
# Purpose: This script edits the tables in Tables/ so that I do not have to edit
# them by hand.

# *************
# * Libraries *
# *************

library(stringr)

# Steps:
# 1. Read .tex files from Tables/
# 2. Check to see if "[h!]" is at the end of line 1. 
# 3. If not add it.

tables_to_check <- list.files("../Tables/", full.names = TRUE)

for (tab in tables_to_check) {

  tab_vec <- readLines(tab)

  if (str_sub(tab_vec[1], start = -4) != "[h!]") {
    tab_vec[1] <- str_c(tab_vec[1], "[h!]")
    cat(tab_vec, file = tab, sep = "\n")
  }

}
