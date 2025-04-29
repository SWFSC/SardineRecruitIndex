# fit assessment to DFA latent trend "as data"
# Created: 4/21/2025, Robert Wildermuth

library(MARSS)
library(r4ss)
library(tidyverse)

# load DFA model fit
load(file = "out/marssFit_1990to2023_noAnch_noSardRec_1trend_EqlVar.RData")

# extract estimated historical trends
trendsHist <- tsSmooth(sardDFA, type = "xtT", interval = "confidence") %>%
                mutate(model = "Local",
                       t = t+1989)


# Regression of DFA trend on rec devs -------------------------------------

mngtBench2024 <- SS_output("C:/Users/r.wildermuth/Documents/CEFI/SardineRecruitmentESP/SardineRecruitIndex/scenarioModels/Pacific sardine 2024 benchmark",)
mngt2024recdevs <- mngtBench2024$recruit %>% filter(era == "Main")

regrDat <- mngt2024recdevs %>% select(Yr, dev) %>%
              left_join(y = trendsHist, by = c("Yr" = "t"))

regrFit <- lm(dev ~ .estimate, data = regrDat)
summary(regrFit) # no sig relationship

# plot(regrFit)

regrDat <- regrDat %>% mutate(predPts = predict(regrFit, newdata = regrDat))

regrDat %>% ggplot(aes(x = .estimate, y = dev)) +
  geom_point() + 
  geom_line(aes(y = predPts)) +
  theme_classic()

# see what it looks like with research model time series
resAssess2025 <- SS_output("C:/Users/r.wildermuth/Documents/CEFI/SardineRecruitmentESP/SardineRecruitIndex/scenarioModels/2025_research_assessment_LorenzM",)
res2025recdevs <- resAssess2025$recruit %>% filter(era == "Main")

regrDat <- res2025recdevs %>% select(Yr, dev) %>%
              left_join(y = trendsHist, by = c("Yr" = "t"))

regrFit <- lm(dev ~ .estimate, data = regrDat)
summary(regrFit) # sig relationship over longer historical period

# plot(regrFit)

regrDat <- regrDat %>% mutate(predPts = predict(regrFit, newdata = regrDat),
                              # color points by era
                              colr = case_when(Yr >= 2005 ~ "benchmark",
                                               TRUE ~ "DFA"))
              

regrDat %>% ggplot(aes(x = .estimate, y = dev)) +
  geom_point(aes(color = colr)) + 
  geom_line(aes(y = predPts)) +
  theme_classic()

# Fit As Data SS model ----------------------------------------------------

# read in and modify data.ss file
ssDat <- SS_readdat(file = "scenarioModels/Example As-Data/data.ss")

# change fleet name
ssDat$fleetinfo <- ssDat$fleetinfo %>% 
                      mutate(fleetname = case_when(fleetname == "BEUTI" ~ "dfaT1",
                                                   TRUE ~ fleetname))
ssDat$fleetnames <- sub("BEUTI", "dfaT1", ssDat$fleetnames)

# have to change error type to deal with anomalies 
ssDat$CPUEinfo <- ssDat$CPUEinfo %>% 
                    mutate(errtype = case_when(fleet == 5 ~ -1,
                                                 TRUE ~ errtype))

# add in DFA estimates
dfaDat <- ssDat$CPUE %>% filter(index == 5) %>%
            left_join(y = trendsHist, by = c("year" = "t")) %>%
            select(year, month, index, .estimate, .se) %>%
            rename(obs = .estimate,
                   se_log = .se)

ssDat$CPUE <- ssDat$CPUE %>% filter(index != 5) %>% bind_rows(dfaDat)

# save updated data file
SS_writedat(datlist = ssDat, outfile = "scenarioModels/benchmarkDFA_AsData/data.ss",
            overwrite = TRUE)
# datTest <- SS_readdat(file = "scenarioModels/benchmarkDFA_AsData/data.ss")

# Run model ----

# Select preferred SS3 version for model fitting
sel_SS <- "v3.30.23" # indicate version of SS to run this analysis with
r4ss::run(dir = "C:/Users/r.wildermuth/Documents/CEFI/SardineRecruitmentESP/SardineRecruitIndex/scenarioModels/benchmarkDFA_AsData",
          exe = file.path("C:/Users/r.wildermuth/Documents/SS3.30/ss3_win.exe"),
          skipfinished = FALSE)

dfaAsDataFit <- SS_output("C:/Users/r.wildermuth/Documents/CEFI/SardineRecruitmentESP/SardineRecruitIndex/scenarioModels/benchmarkDFA_AsData",)
SS_plots(dfaAsDataFit)

# Compare to recruitment pattern from 2024 assessment
mngtDir <- "../SardineRecruitIndex/scenarioModels/Pacific sardine 2024 benchmark"

mngtAssmt2024 <- SS_output(dir = mngtDir, repfile = "Report.sso", printstats = FALSE)

compSmry <- SSsummarize(list(dfaAsDataFit, 
                             mngtAssmt2024))
dev.off()
SSplotComparisons(compSmry, legendlabels = c("dfaAsData", 
                                             "mngtAssmt2024"))

# diff in estimated rec devs
recdevDiffs <- compSmry$recdevs %>% mutate(recdevDiff = model1 - model2) # w/Envt - benchmark
mean(recdevDiffs$recdevDiff)
hist(recdevDiffs$recdevDiff, breaks = 10)

bioDiffs <- compSmry$SmryBio %>% mutate(bioDiff = model1 - model2) # w/Envt - benchmark
mean(bioDiffs$bioDiff)
hist(bioDiffs$bioDiff, breaks = 10)
