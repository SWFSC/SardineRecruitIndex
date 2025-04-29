#### Code to extract Assessment Ages and merge with Trawl data
# packages 
library(DBI)
library(openxlsx)
library(tidyverse)

#need to be connected to the VPN or Ethernet at building to access databases

# connect to life history database and read in ages #####
cpslh <- DBI::dbConnect(odbc::odbc(),
                        DRIVER="ODBC Driver 18 for SQL Server",
                        Encrypt = "Optional",
                        DATABASE="CPSLifeHistory",
                        Trusted_Connection= "Yes",
                        SERVER="swc-estrella-s")
# read in assessment ages view
ages <- dbReadTable(cpslh, "Assessment_Ages")

# connect to trawl database and read in specimen data ####
trawl <- DBI::dbConnect(odbc::odbc(),
                        DRIVER="ODBC Driver 18 for SQL Server",
                        Encrypt = "Optional",
                        DATABASE="Trawl",
                        Trusted_Connection= "Yes",
                        SERVER="swc-estrella-s")
# specimen data from NOAA ships
specimen.noaa <- dbReadTable(trawl, "Specimen") %>%
  # reduce the number of columns, can always add any you want back in
  select(cruise,ship,haul,collection,species,specimenNumber,otolithNumber,individual_ID, standardLength_mm, forkLength_mm,weightg, sex)

# specimen data from nearshore survey 
specimen.nearshore <- dbReadTable(trawl,"Nearshore_Specimen") %>%
  select(cruise,ship,date,set,species,specimenNumber,individual_ID,standardLength_mm, forkLength_mm,weightg, sex)

# combining both 
specimen <- bind_rows(specimen.noaa,specimen.nearshore)

# getting scientific name 
sci.name <- dbReadTable(trawl, "SpeciesCodes") %>%
  select(species,scientificName)
# merging age data with scientific name 
ages <- left_join(ages, sci.name)

#### merge age data and trawl data #####
assessmentAges <- left_join(ages,specimen)

# getting haul date/time
haul.time <- dbReadTable(trawl, "haul") %>%
  select(collection, equilibriumTime)
assessmentAges2 <- left_join(assessmentAges, haul.time) %>%
  mutate(core_date = as.Date(substr(equilibriumTime, 1, 10), format = "%Y-%m-%d"),
         nearshore_date = as.Date(date, format = "%m/%d/%Y"),
         final_date = as.Date(ifelse(is.na(core_date),
                             nearshore_date,
                             core_date))) %>%
  select(-c(date, equilibriumTime, core_date, nearshore_date))

# Save age/trawl data
write.csv(assessmentAges2, here::here("Data", "Indiv assess age+weight.csv"),
          row.names = FALSE)
