
#######################################
# NOTES

# This script is meant to obtain SS3 executable files
# SS3 executable files are included in .gitignore and are not uploaded to GitHub

# Developed by Alex Jensen, SWFSC, NMFS, NOAA

#######################################

# ID all SS3 versions to install .exe file for
SS_vers <- c("v3.30.22.1", "v3.30.22")

# Create local-only directory for storing different versions of SS3 exe files
  # Only creates directory if not already present
if(!("SS3 exe files" %in% substr(list.dirs(), 3, 100))) {
  dir.create(here::here("SS3 exe files"))
}

# Create directories for and install different versions of SS3 exe files
  # Only creates directories if directories are not already present
for(i in 1:length(SS_vers)){
  if(!(SS_vers[i] %in% substr(list.dirs(), 3, 100))) {
    dir.create(here::here("SS3 exe files", SS_vers[i]))
  }
  r4ss::get_ss3_exe(dir = here::here("SS3 exe files", SS_vers[i]), 
                    version = SS_vers[i])
}
