# Model selection of recruitment DFAs based on LFOIC
# Created: 12/14/2023, Robert Wildermuth
# Copied from https://github.com/futureseas/recrmntDFA 1/16/2025

library(tidyverse)
library(MARSS)
source("R/LFOXV.R")
source("R/OSAResids.R")

# read prepped dataset
datDFA <- read_csv("../SardineRecruitIndex/Data/recrDFAdat.csv")
load(file = "Data/indicatorSetNames_LUSI39spawnHabsprSST.RData")
# Assess LFOIC for last 10 years of data
peels <- 10


# Data for full historical dataset ---------------------------------------

# remove unused indicators
allDat <- datDFA %>% filter(year %in% 1985:2021) %>%
            select(all_of(setNames), -anchBioSmrySeas1) 

datNames <- names(allDat)[-1]

# Create a custom R obs error matrix assuming each data source has it's own common error
Rcustom <- matrix(list(0),length(datNames),length(datNames))
diag(Rcustom) <- c("COP", "COP", 
                   "BEUTI", "BEUTI", 
                   "CUTI", "CUTI",
                   "LUSI",
                   "STI",
                   "RREAS", "RREAS", 
                   "WAA", "WAA",
                   "NEMURO", "NEMURO",
                   "sardRec",
                   # "anchBio",
                   "SST",
                   "Transp", "Transp", "Transp", "Transp",
                   "condK",
                   "LUSI",
                   "sardSDM", 
                   "sardlarvSDM", 
                   "SST")



# LFOIC  ------------------------------------------------------------------

# Set up table of model structures to hold LFOIC vals
xvModSel <- tibble(initYr = 0, 
                   Rstructure = "", 
                   mTrends = 0)[0,]

# loop over initial dates
for(y in c(#1980, 1990
            1985)){
  cat("\n")
  print(y)
  
  # subset from 1980, 1985, or 1990 to 2023
  initDat <- allDat %>% filter(year %in% y:2021)
  
  # transpose for MARSS formatting
  itDat <- initDat %>% select(-year) %>% t()
  
  # loop over number of trends
  for(m in 5:1){
    cat("\n Trends: ", m)
    cat("\n Diagonal and equal R matrix")
    itEqRMSE <- LFOXV(dfaDat = itDat,
                       Rstructure = "diagonal and equal",
                       mTrends = m,
                       peels = peels,
                      horizon = 5,
                      colsRMSE = "sardRec")

    itEqRMSE <- itEqRMSE %>% mutate(initYr = y,
                                    mTrends = m,
                                    Rstructure = "diag & equal")

    cat("\n Diagonal and unequal R matrix")

    itUneqRMSE <- LFOXV(dfaDat = itDat,
                      Rstructure = "diagonal and unequal",
                      mTrends = m,
                      peels = peels,
                      horizon = 5,
                      colsRMSE = "sardRec")

    itUneqRMSE <- itUneqRMSE %>% mutate(initYr = y,
                                    mTrends = m,
                                    Rstructure = "diag & unequal")

    cat("\n Custom R matrix")

    itCustRMSE <- LFOXV(dfaDat = itDat,
                        Rstructure = Rcustom,
                        mTrends = m,
                        peels = peels,
                        horizon = 5,
                        colsRMSE = "sardRec")

    itCustRMSE <- itCustRMSE %>% mutate(initYr = y,
                                    mTrends = m,
                                    Rstructure = "custom R by SDM")

    xvModSel <- xvModSel %>% bind_rows(itUneqRMSE, itCustRMSE, itEqRMSE)
    # xvModSel <- xvModSel %>% bind_rows(itCustRMSE)
  } # end trends loop
} # end year loop 

xvModSel <- xvModSel %>% mutate(nIndices = length(datNames),
                                peels = peels)
  



# Calculate Persistence Prediction RMSE -----------------------------------

datLen <- length(1985:2021)

perstResids <- allDat %>% select(year, sardRec) %>% 
                  filter(year >= 1985, year < 2022) %>% 
                  mutate(zscoreSardRec = zscore(sardRec),
                         perst1Sard = c(NA, zscoreSardRec[1:(datLen-1)]),
                         resid1Sard = zscoreSardRec - perst1Sard,
                         perst2Sard = c(NA,NA, zscoreSardRec[1:(datLen-2)]),
                         resid2Sard = zscoreSardRec - perst2Sard,
                         perst3Sard = c(NA,NA,NA, zscoreSardRec[1:(datLen-3)]),
                         resid3Sard = zscoreSardRec - perst3Sard,
                         perst4Sard = c(NA,NA,NA,NA, zscoreSardRec[1:(datLen-4)]),
                         resid4Sard = zscoreSardRec - perst4Sard,
                         perst5Sard = c(NA,NA,NA,NA,NA, zscoreSardRec[1:(datLen-5)]),
                         resid5Sard = zscoreSardRec - perst5Sard,
                         resid0Mean = zscoreSardRec - 0) %>% # residual from mean S-R value is just the rec dev
                  select(year, resid1Sard, resid2Sard, resid3Sard, resid4Sard, resid5Sard, resid0Mean) %>% 
                  pivot_longer(cols = -year, names_prefix = "resid", names_sep = 1, 
                               names_to = c("predHoriz", "variable"), 
                               values_to = "perstResid") %>% 
                  filter(year %in% 2012:2021) %>%
                  group_by(variable, predHoriz) %>%
                  summarize(sosRes = sum(perstResid^2, na.rm = TRUE),
                            peels = sum(!is.na(perstResid))) %>%
                  mutate(RMSE = sqrt(sosRes/peels),
                         resType = "resid.Perst",
                         predHoriz = as.numeric(predHoriz)) %>%
                  select(-sosRes)



# write_csv(xvModSel, file = "out/historicalModelSelection_NoAnch.csv")

xvModSel <- read_csv("out/historicalModelSelection_noAnch.csv") %>%
              mutate(dataset = "noAnch")
testAnch <- read_csv("out/historicalModelSelection_Anch.csv") %>%
              mutate(dataset = "Anch")

xvModSel <- bind_rows(xvModSel, testAnch)
xvModSel <- bind_rows(xvModSel, perstResids)
# plot out best performing model structures over prediction horizons
xvModSel %>% filter(resType %in% "resid.Inf", variable %in% c("Sard", "sardRec")) %>%
  ggplot(aes(x = predHoriz, y = RMSE)) +
  geom_line(aes(color = as.character(mTrends)), linewidth = 1) + #paste(mTrends, Rstructure, sep = "-"))) +
  geom_point(data = perstResids %>% filter(!variable == "Mean")) +
  geom_hline(yintercept = perstResids %>% filter(variable == "Mean") %>% pull(RMSE)) + 
  # scale_color_viridis_d() +
  facet_grid(cols= vars(Rstructure), rows = vars(dataset), scales = "free_y")
  # facet_wrap(~Rstructure)

xvModSel %>% filter(resType %in% c("resid.Inf", "resid.Perst"), 
                    variable %in% c("Sard", "sardRec"),
                    predHoriz == 1) %>%
  arrange(RMSE)
