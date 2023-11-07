#!/bin/bash

# File: rsync_with_nova.sh
# Created by: Caleb Leedy
# Created on: October 25, 2023
# Updated on: October 25, 2023
# 
# Purpose:
# This file syncs my local repository with my repository on millipede.
rsync -rtvz R/ cleedy@nova.its.iastate.edu:~/Data_Integration/R

