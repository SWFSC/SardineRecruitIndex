# Code to calculate LFOIC over 10 peels for candidate GAM recruitment index model
# borrowing from LFOXV.R and OSAResids.R
# Created 6/4/2025, Robert Wildermuth

library(r4ss)
library(tidyverse)
library(MARSS)
library(mgcv)
source("R/LFOXV.R") # LFOIC calculation for DFA


# Modified OSHResids fxn --------------------------------------------------

# Modified function from Desiree Tommasi, Robert Wildermuth
#' @param objGAM GAM object fitted outside this function to data up to a specified peel
#' @param fullDat all available data
#' @param p peel used when fitting objMARSS
#' @param horizon integer for steps ahead in the prediction horizon
#' @param colsRMSE # column names to calculate RMSE for. Default: c("sardRec", "anchRec")
#' @return dataframe of residuals for the specified 
#' datasets, scaled observations fed into the model estimate (y.*),
#' original scaled observations (origDat), prediction (est.*), 
#' and timestep (t). 
#' Note only the residuals for the last 'horizon' time steps where y=NA are 
#' out of sample and true innovation residuals

OSAResids.gam <- function(objGAM, fullDat, p, horizon = 1,
                      colsRMSE = c("dev") 
){
  
  
  origDat <- fullDat[1:min(nrow(fullDat), # either full data set or ...
                            nrow(fullDat)-p+horizon), ] # including 'horizon' steps ahead of peel
  
  # naive innovations
  origDat$estDev <- predict(object = objGAM, newdata = origDat)

  osaResids <- origDat %>% dplyr::select(Yr, all_of(colsRMSE), estDev)
  
  osaResids <- osaResids %>% mutate(resid = .data[[colsRMSE]] - estDev)
  
  return(osaResids)
}


# LFO peels ---------------------------------------------------------------

# Read in and prep data
mngtBench2024 <- SS_output("C:/Users/r.wildermuth/Documents/CEFI/SardineRecruitmentESP/SardineRecruitIndex/scenarioModels/Pacific sardine 2024 benchmark",)
mngt2024recdevs <- mngtBench2024$recruit %>% filter(era == "Main")

# read prepped dataset
datDFA <- read_csv("../SardineRecruitIndex/Data/recrDFAdat.csv")
datGAM <- datDFA
datGAM <- mngt2024recdevs %>% dplyr::select(Yr, dev) %>% 
            filter(Yr <=2021) %>%
            left_join(y = datGAM, by = c("Yr"="year")) %>%
            select(Yr, dev, 
                   # top indicators from DFA analyses
                   ZM_NorCal, springSST, NCOPspring, sardNurseHab, BEUTI_39N,
                   CUTI_39N, summerSST, NCOPsummerlag1, ZM_SoCal, OC_LUSI_39N,
                   OC_STI_39N, HCI_30N355N, # high loadings in other indicator set
                   # full indicator set
                   NCOPspring, NCOPsummerlag1, BEUTI_33N, BEUTI_39N, CUTI_33N,
                   CUTI_39N, OC_LUSI_33N, OC_STI_33N, RREAS_Myctophids, 
                   RREAS_YOYsardine, age1SprSardmeanWAA, meanSSBwt, 
                   ZM_NorCal, ZM_SoCal, sardNurseHab, anchBioSmrySeas1, 
                   summerSST, avgOffTransspring, avgOffTranssummer, 
                   sprRelOffTrans, sumRelOffTrans, meanK,
                   OC_STI_39N, sardSpawnHab, springSST,
                    sardRec)
datGAM <-cbind(datGAM[,1:2],
               apply(X = datGAM[, -(1:2)], MARGIN = 2, FUN = zscore))

# Set up table of model structures to hold LFOIC vals
peelRMSE <- tibble(Yr = 0)[0,]

peels <- 10
horizon <- 5

for(i in 1:peels){
  
  cat("\nFitting peel: ", i)
  
  # peelGAM3 <- gam(dev ~ OC_STI_39N + s(BEUTI_39N, k = 4) + s(HCI_30N355N, k = 4) - 1,
  #                data = datGAM[1:(nrow(datGAM)-i),], method = "REML")
  # peelGAM <- gam(dev ~ s(OC_STI_39N, k = 4) + CUTI_33N + s(BEUTI_39N, k = 4) + s(HCI_30N355N, k = 4) - 1,
  #                data = datGAM[1:(nrow(datGAM)-i),], method = "REML")
  # peelGAM <- gam(dev ~ s(meanSSBwt, k = 4) + OC_STI_39N + s(BEUTI_39N, k = 4) + s(HCI_30N355N, k = 4) - 1,
  #                  data = datGAM[1:(nrow(datGAM)-i),], method = "REML")
  # peelGAM <- gam(dev ~ s(meanSSBwt, k = 4) + s(BEUTI_39N, k = 4) + s(HCI_30N355N, k = 4) - 1,
  #                  data = datGAM[1:(nrow(datGAM)-i),], method = "REML")
  # peelGAM <- gam(dev ~ OC_LUSI_36N + OC_STI_39N + s(BEUTI_39N, k = 4) + s(HCI_30N355N, k = 4) - 1,
  #                  data = datGAM[1:(nrow(datGAM)-i),], method = "REML")
  # peelGAM <- gam(dev ~ s(OC_LUSI_39N, k = 4) + s(BEUTI_39N, k = 4) + s(HCI_30N355N, k = 4) - 1,
  #                data = datGAM[1:(nrow(datGAM)-i),], method = "REML", select = TRUE)
  # peelGAM <- gam(dev ~ s(C.pacificus, k = 4) + s(OC_LUSI_39N, k = 4) + s(BEUTI_39N, k = 4) + s(HCI_30N355N, k = 4) - 1,
  #                data = datGAM[1:(nrow(datGAM)-i),], method = "REML", select = TRUE)
  # peelGAM <- gam(dev ~ s(C.pacificus, k = 4) + s(BEUTI_39N, k = 4) + s(HCI_30N355N, k = 4) - 1,
  #                data = datGAM[1:(nrow(datGAM)-i),], method = "REML", select = TRUE)
  peelGAM4 <- gam(dev ~ s(NCOPsummerlag1, k = 4) + s(ZM_NorCal, k = 4) - 1,
                 data = datGAM[1:(nrow(datGAM)-i),], method = "REML")
  # peelGAM3 <- gam(dev ~ s(NCOPsummerlag1, k = 4) + OC_STI_39N + s(BEUTI_39N, k = 4) + s(HCI_30N355N, k = 4) - 1,
  #                    data = datGAM[1:(nrow(datGAM)-i),], method = "REML")  
  peelGAM3 <- gam(dev ~ s(NCOPsummerlag1, k = 4) + s(HCI_30N355N, k = 4) - 1,
                     data = datGAM[1:(nrow(datGAM)-i),], method = "REML")
  # peelGAM3 <- gam(dev ~ s(NCOPsummerlag1, k = 4) + s(sumRelOffTrans, k = 4) - 1,
  #                    data = datGAM[1:(nrow(datGAM)-i),], method = "REML")
  # peelGAM3 <- gam(dev ~ s(NCOPsummerlag1, k = 4) + s(HCI_30N355N, k = 4) - 1,
  #                    data = datGAM[1:(nrow(datGAM)-i),], method = "REML")
  #
  # check convergence
  if(!peelGAM3$converged){
    print("\n Failed to converge for peel ", i)
    return(-99999) # default RMSE for failed convergence
  }
  if(!peelGAM4$converged){
    print("\n Failed to converge for peel ", i)
    return(-99999) # default RMSE for failed convergence
  }
  
  # get one-step-ahead residuals
  resids3 <- OSAResids.gam(objGAM = peelGAM3, fullDat = datGAM, p = i, horizon = horizon,
                      colsRMSE = "dev")
  resids4 <- OSAResids.gam(objGAM = peelGAM4, fullDat = datGAM, p = i, horizon = horizon,
                      colsRMSE = "dev")
  
  # collect residuals for datum of interest
  itPeelRMSE3 <- resids3 %>% mutate(signMatch = sign(dev) * sign(estDev)) %>% 
                  dplyr::select(Yr, resid, signMatch) %>% 
                  filter(Yr %in% (max(datGAM$Yr)-i+1):max(Yr)) %>% 
                  mutate(peel = i)
  itPeelRMSE3 <- itPeelRMSE3 %>% mutate(predHoriz = 1:nrow(itPeelRMSE3),
                                        compMethod = "GAM3")
  itPeelRMSE4 <- resids4 %>% mutate(signMatch = sign(dev) * sign(estDev)) %>% 
                  dplyr::select(Yr, resid, signMatch) %>% 
                  filter(Yr %in% (max(datGAM$Yr)-i+1):max(Yr)) %>% 
                  mutate(peel = i)
  itPeelRMSE4 <- itPeelRMSE4 %>% mutate(predHoriz = 1:nrow(itPeelRMSE4),
                                        compMethod = "GAM4")
  peelRMSE <- bind_rows(peelRMSE, itPeelRMSE3, itPeelRMSE4)
    
                
} # end peel for-loop

overPeelGAM <- peelRMSE %>% 
                # Calc RMSE over all time points
                group_by(predHoriz, compMethod) %>%
                summarize(sosRes = sum(resid^2, na.rm = TRUE),
                          nObs = sum(!is.na(resid)),
                          avgAbsRes = mean(abs(resid), na.rm = TRUE),
                          meanBias = mean(resid, na.rm = TRUE),
                          nSignMatch = sum(signMatch > 0, na.rm = TRUE)) %>%
                mutate(RMSE = sqrt(sosRes/nObs),
                       freqSignMatch = nSignMatch/nObs,
                       variable = "dev",
                       resType = "resid.Inf") %>% arrange(RMSE)


# Re-calculate DFA LFOIC for same time period -----------------------------
load(file = "Data/indicatorSetNames_STI39spawnHabsprSST.RData")

m <- 1

# Add mngt assessment rec devs and remove unused indicators
allDat <- datDFA
allDat <- mngt2024recdevs %>% dplyr::select(Yr, dev) %>% 
            left_join(y = allDat, by = c("Yr"="year")) %>%
            rename(year = Yr)

allDat <- allDat %>% select(year, dev, 
                            # ZM_NorCal, springSST, NCOPspring, sardNurseHab, BEUTI_39N,
                            #        CUTI_39N, summerSST, NCOPsummerlag1, ZM_SoCal, OC_LUSI_39N,
                            #        OC_STI_39N, HCI_30N355N, # high loadings in other indicator set
                            #        year, sardRec) %>% select(-sardRec)
                            all_of(setNames)) %>% select(-sardRec)

datNames <- names(allDat)[-1]

# subset from 2005 to 2023
initDat <- allDat %>% filter(year %in% 2005:2021)

# transpose for MARSS formatting
itDat <- initDat %>% select(-year) %>% t()

itEqRMSE <- LFOXV(dfaDat = itDat,
                  Rstructure = "diagonal and equal",
                  mTrends = m,
                  peels = peels,
                  horizon = horizon,
                  colsRMSE = "dev")
itEqRMSE <- itEqRMSE %>% mutate(nIndices = length(datNames),
                                peels = peels) %>%
                filter(resType == "resid.Inf", variable %in% c("sardRec", "dev"))

# Calculate bias and sign matching

peelRMSE <- tibble(.rownames = '', 
                   t = 0)[0,]


for(i in 1:peels){
  
  cat("\nFitting peel: ", i)
  
  peelDFA <- MARSS(y = itDat[,1:(ncol(itDat)-i)], 
                   form = "dfa",
                   method = "BFGS",
                   inits = list(x0 = matrix(1, 1, 1)),
                   z.score = TRUE,
                   model = list(R = "diagonal and equal",
                                m = m),
                   silent = TRUE)
  
  # check convergence
  if(peelDFA$convergence != 0){
    stop("\n Failed to converge for peel ", i)
    return(-99999) # default RMSE for failed convergence
  }
  
  # get one-step-ahead residuals
  resids <- OSAResids(objMARSS = peelDFA, fullDat = itDat, p = i, horizon = horizon,
                      colsRMSE = "dev")
  
  # collect residuals for datum of interest
  itPeelRMSE <- resids %>% mutate(signMatch = sign(est.Inf) * sign(origDat)) %>%
                  select(.rownames, t, resid.Naiv, resid.Inf, resid.Cont, 
                         resid.Proj, signMatch) %>% 
                  filter(.rownames == "dev", t %in% (ncol(itDat)-i+1):max(t)) %>% 
                  mutate(peel = i)
  peelRMSE <- itPeelRMSE %>% mutate(predHoriz = (t-(ncol(itDat)-i))) %>%
                bind_rows(peelRMSE)
} # end peel for-loop

overPeelDFA <- peelRMSE %>% 
                # Calc RMSE over all time points
                group_by(.rownames, predHoriz) %>%
                summarize(sosRes = sum(resid.Inf^2, na.rm = TRUE),
                          nObs = sum(!is.na(resid.Inf)),
                          avgAbsRes = mean(abs(resid.Inf), na.rm = TRUE),
                          meanBias = mean(resid.Inf, na.rm = TRUE),
                          nSignMatch = sum(signMatch > 0, na.rm = TRUE)) %>%
                mutate(RMSE = sqrt(sosRes/nObs),
                       freqSignMatch = nSignMatch/nObs,) %>% arrange(RMSE)


# Re-do persistence calculations and plot ---------------------------------

datLen <- length(2005:2021)

perstResids <- allDat %>% select(year, dev) %>% 
  filter(year >= 2005, year < 2022) %>% 
  mutate(zscoreSardRec = zscore(dev),
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
         # residual from mean S-R value is just the rec dev
         resid0Mean = zscoreSardRec - 0,
         # sign matching calculations
         signMatch1 = sign(zscoreSardRec) * sign(perst1Sard),
         signMatch2 = sign(zscoreSardRec) * sign(perst2Sard),
         signMatch3 = sign(zscoreSardRec) * sign(perst3Sard),
         signMatch4 = sign(zscoreSardRec) * sign(perst4Sard),
         signMatch5 = sign(zscoreSardRec) * sign(perst5Sard)) 

perstSignMatch <- perstResids %>%
  select(year, signMatch1, signMatch2, signMatch3, signMatch4, signMatch5) %>% 
  pivot_longer(cols = -year, names_prefix = "signMatch", #names_sep = 1, 
               names_to = c("predHoriz"), 
               values_to = "signMatch") %>% 
  filter(year %in% 2012:2021) %>%
  group_by(predHoriz) %>%
  summarize(peels = sum(!is.na(signMatch)),
            nSignMatch = sum(signMatch > 0, na.rm = TRUE)) %>%
  mutate(freqSignMatch = nSignMatch/peels)

perstResids <- perstResids %>%
  select(year, resid1Sard, resid2Sard, resid3Sard, resid4Sard, resid5Sard, resid0Mean) %>% 
  pivot_longer(cols = -year, names_prefix = "resid", names_sep = 1, 
               names_to = c("predHoriz", "variable"), 
               values_to = "perstResid") %>% 
  filter(year %in% 2012:2021) %>%
  group_by(variable, predHoriz) %>%
  summarize(sosRes = sum(perstResid^2, na.rm = TRUE),
            peels = sum(!is.na(perstResid)),
            meanBias = mean(perstResid, na.rm = TRUE)) %>%
  mutate(RMSE = sqrt(sosRes/peels),
         resType = "resid.Perst",
         predHoriz = as.numeric(predHoriz)) %>%
  select(-sosRes)

xvModSel <- bind_rows(itEqRMSE, 
                      perstResids, overPeelGAM) %>%
              mutate(compMethod = case_when(compMethod == "individual" ~ "DFA",
                                            TRUE ~ compMethod))

# plot out best performing model structures over prediction horizons
xvModSel %>% filter(resType %in% "resid.Inf", 
                    variable %in% c("Sard", "sardRec", "dev")) %>%
  ggplot(aes(x = predHoriz, y = RMSE)) +
  geom_line(aes(color = compMethod), linewidth = 2) + 
  geom_point(data = perstResids %>% filter(!variable == "Mean")) +
  geom_hline(yintercept = perstResids %>% filter(variable == "Mean") %>% pull(RMSE)) +
  theme_classic() + ylim(0,3.5)
