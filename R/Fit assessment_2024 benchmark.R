#######################################
# NOTES
  
# This script is meant to run and save output from the 
  # 2024 Pacific sardine benchmark assessment
# Ensure 'Obtain SS3 exe files.R' has been previously run

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
  # 2024 Benchmark used v3.30.22

# Set directory for SS files
base_wd <- here::here("scenarioModels", "Pacific sardine 2024 benchmark")

# Run 2024 Pacific sardine benchmark assessment, save output/plots ---------------------------
  # Run model, save results as .rds file ----

r4ss::run(dir = base_wd,
          exe = file.path(getwd(), "SS3 exe files", sel_SS, "ss3.exe"),
          skipfinished = FALSE)

base_assess <- SS_output(base_wd)
saveRDS(base_assess, file.path(base_wd, "base_assess.rds"))

  # Read model output, make default plots ----

base_assess <- readRDS(file.path(base_wd, "base_assess.rds"))
SS_plots(base_assess)

 