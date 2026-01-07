# use DFA fit with recruitment to forecast rec dev
# Created: 4/24/2025, Robert Wildermuth

library(MARSS)
library(r4ss)
library(tidyverse)

# load DFA model fit
load(file = "out/marssFit_1990to2023_noAnch_1trend_EqlVar.RData")

# For now use naive innovations
yExpNaive <- predict(object = sardDFA, interval = "none", 
                       n.ahead = 1, type = 'ytT')

# extract estimated historical trends
yExpHist <- fitted(sardDFA, type = "ytT", interval = "confidence") %>%
  mutate(t = t+1989)
yExpHist %>% filter(.rownames == "sardRec")
yExpNaive <- yExpNaive$pred %>%
                  mutate(t = t+1989) %>% 
                  filter(.rownames == "sardRec")
# For now, these are the same as the 2023 values because no other indicators are 
# informing the sardRec estimate in predict()


# Fit As Forecast SS model ------------------------------------------------

# For now, leave the forecast file alone?

# edit control file to read in the 2024 estimated rec dev
ssCtl <- SS_readctl(file = "scenarioModels/Pacific sardine 2024 benchmark/control.ss")

# check that advanced options for SR turned on
# !!RW: I think this is the right variable to look at?
ssCtl$recdev_adv # should be 1 

# change last year of main rec_devs
ssCtl$MainRdevYrLast <- 2022

# Change forecast recruitment phase to -1 to force input values for late and forecast rec devs
ssCtl$Fcast_recr_phase <- -1

# Change end year for ramp bias adjustment
ssCtl$first_recent_yr_nobias_adj <- 2022.9

# Change number of rec devs to read in
ssCtl$N_Read_recdevs <- 2
# provide the new rec dev for the forecast year
newRecDevs <- yExpctNaive %>% filter(t %in% 2023:2024) %>% 
                select(t, estimate) %>% 
                rename(Year = t,
                       recdev = estimate)
ssCtl$recdev_input <- newRecDevs
  
# save the new control file
SS_writectl(ctllist = ssCtl, outfile = "scenarioModels/benchmarkDFA_AsForecast/control.ss",
            overwrite = TRUE)
# check
# ctlTest <- SS_readctl(file = "scenarioModels/benchmarkDFA_AsForecast/control.ss")

# Run model ----

# Select preferred SS3 version for model fitting
sel_SS <- "v3.30.23" # indicate version of SS to run this analysis with
r4ss::run(dir = "C:/Users/r.wildermuth/Documents/CEFI/SardineRecruitmentESP/SardineRecruitIndex/scenarioModels/benchmarkDFA_AsForecast",
          exe = file.path("C:/Users/r.wildermuth/Documents/SS3.30/ss3_win.exe"),
          skipfinished = FALSE)

dfaAsFcastFit <- SS_output("C:/Users/r.wildermuth/Documents/CEFI/SardineRecruitmentESP/SardineRecruitIndex/scenarioModels/benchmarkDFA_AsForecast",)
SS_plots(dfaAsFcastFit)

# Compare to recruitment pattern from 2024 assessment
mngtDir <- "../SardineRecruitIndex/scenarioModels/Pacific sardine 2024 benchmark"

mngtAssmt2024 <- SS_output(dir = mngtDir, repfile = "Report.sso", printstats = FALSE)

compSmry <- SSsummarize(list(dfaAsFcastFit, 
                             mngtAssmt2024))
dev.off()
SSplotComparisons(compSmry, legendlabels = c("dfaAsFcast", 
                                             "mngtAssmt2024"))

recdevDiffs <- compSmry$recdevs %>% mutate(recdevDiff = model1 - model2) # w/Envt - benchmark
mean(recdevDiffs$recdevDiff)
hist(recdevDiffs$recdevDiff, breaks = 10)

bioDiffs <- compSmry$SmryBio %>% mutate(bioDiff = model1 - model2) # w/Envt - benchmark
mean(bioDiffs$bioDiff)
hist(bioDiffs$bioDiff, breaks = 10)
# Envt-informed model estimates higher final biomass, even though rec devs are smaller?