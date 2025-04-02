#######################################
# NOTES
  
# This script is meant to run and save output from the Example As-Data model
    # Based on the '2024 Pacific sardine benchmark assessment'
    # Manually modified data.ss file to:
      # Increase number of fleets from 4 to 5
      # Provide definition for new fleet (`3 1 1 2 0 BEUTI`)
      # Specify type of survey index provided by new fleet (`5 31 0 0`)
      # Add unstandardized BEUTI values (2005-2023) to survey data section for new fleet
        # Data obtained from 'beuti_39.csv'
        # SE values were selected arbitrarily to avoid breaking the model
      # Add new lines to length and age composition sections for new fleet 
    # Manually modified control.ss file to:
      # Add catchability information for new fleet
      # Add new lines to size and age selectivity sections for new fleet
# Ensure 'Obtain SS3 exe files.R' has been previously run
# Model currently throws warning messages for not providing WAA for new survey fleet
  # Unsure as to best method for handling this...

# Developed by Alex Jensen, SWFSC, NMFS, NOAA

#######################################

# Universal analysis set-up ----------------

# Load libraries
require(r4ss)
require(ss3diags)
require(dplyr)
require(tidyverse)
require(ggsidekick)
require(ggridges)
require(patchwork)
require(scales)
require(here)
require(reshape2)
require(plyr)

# Select preferred SS3 version for model fitting
sel_SS <- "v3.30.22" # indicate version of SS to run this analysis with
  # 2023 Benchmark used v3.30.22

# Set directory for SS files
ex_AsData_wd <- here::here("scenarioModels", "Example As-Data")

# Run Example As-Data assessment, save output/plots ---------------------------
  # Run model, save results as .rds file ----

r4ss::run(dir = ex_AsData_wd,
          exe = file.path(getwd(), "SS3 exe files", sel_SS, "ss3.exe"),
          skipfinished = FALSE)

ex_AsData_assess <- SS_output(ex_AsData_wd)
# saveRDS(ex_AsData_assess, file.path(ex_AsData_wd, "ex_AsData_assess.rds"))

  # Read model output, make default plots ----

ex_AsData_assess <- readRDS(file.path(ex_AsData_wd, "ex_AsData_assess.rds"))
SS_plots(ex_AsData_assess)
