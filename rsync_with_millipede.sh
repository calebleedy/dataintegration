#!/bin/bash

# File: rsync_with_millipede.sh
# Created by: Caleb Leedy
# Created on: April 12, 2021
# Updated on: May 4, 2023
# 
# Purpose:
# This file syncs my local repository with my repository on millipede.

echo "Current Directory is `pwd`"

# rsync -rtvz ./ cleedy@prontodtn.las.iastate.edu:~/Data_Integration/
# rsync -rtvz ./ cleedy@millipede.cssm.iastate.edu:~/Data_Integration/
rsync -rtvz Simulations/ cleedy@millipede.cssm.iastate.edu:~/Data_Integration/Simulations

