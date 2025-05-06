# fit assessment to final years of DFA latent trend "as nowcast", blending 
# "as data" and "as forecast" methods
# Created: 5/6/2025, Robert Wildermuth

library(MARSS)
library(r4ss)
library(tidyverse)

# load DFA model fit
load(file = "out/marssFit_1990to2023_noAnch_1trend_EqlVar.RData")

# extract estimated historical rec devs
sardRecHist <- tsSmooth(sardDFA, type = "ytT", interval = "confidence") %>%
                  mutate(t = t+1989) %>%
                  filter(.rownames == "sardRec")

# want to emulate making nowcast prediction
# create 'YStar' by removing last observation of 'sardRec' 
YStar <- sardDFA$marss$data
YStar['sardRec', ncol(YStar)] <- NA
# informed innovations
sardRecPred <- predict(object = sardDFA, interval = "confidence", 
                       newdata = list(t = 1:ncol(YStar),
                                      y = YStar),
                       type = 'ytT', x0 = "use.model")
sardRecPred <- sardRecPred$pred %>%
                mutate(t = t+1989) %>%
                filter(.rownames == "sardRec",
                       t %in% 2022:2023) #!!RW will want to see how many years to provide data


# Fit SS model with DFA nowcast ------------------------------------------------

# read in and modify data.ss file
ssDat <- SS_readdat(file = "scenarioModels/Example As-Data/data.ss")

# change fleet name
ssDat$fleetinfo <- ssDat$fleetinfo %>% 
                      mutate(fleetname = case_when(fleetname == "BEUTI" ~ "dfaT1",
                                                   TRUE ~ fleetname))
ssDat$fleetnames <- sub("BEUTI", "dfaT1", ssDat$fleetnames)

# have to change units to select rec dev as inputs 
ssDat$CPUEinfo <- ssDat$CPUEinfo %>% 
                    mutate(errtype = case_when(fleet == 5 ~ -1,
                                                 TRUE ~ errtype),
                           units = case_when(fleet == 5 ~ 36,
                                                 TRUE ~ units))

# add in DFA estimates
dfaDat <- ssDat$CPUE %>% filter(index == 5) %>%
            left_join(y = sardRecPred, by = c("year" = "t")) %>%
            select(year, month, index, estimate, se) %>%
            rename(obs = estimate,
                   se_log = se) %>%
            filter(!is.na(obs))

ssDat$CPUE <- ssDat$CPUE %>% filter(index != 5) %>% bind_rows(dfaDat)

# save updated data file
SS_writedat(datlist = ssDat, outfile = "scenarioModels/benchmarkDFA_AsNowcast/data.ss",
            overwrite = TRUE)
# datTest <- SS_readdat(file = "scenarioModels/benchmarkDFA_AsData/data.ss")

# need to add catchability (Q) parameters for the new index
ssCtl <- SS_readctl(file = "scenarioModels/Example As-Data/control.ss")

# add an extra SE for Q
ssCtl$Q_options <- ssCtl$Q_options %>% 
                    mutate(extra_se = case_when(fleet == 5 ~ 1,
                                                TRUE ~ extra_se))

# fix Q = 1 and estimate extra_se
ssCtl$Q_parms$INIT <- c(0,1)
newQparms <- as.data.frame(t(c(0.001, 1.3, 0.949309, 0, 99, 0, 2, 0, 0, 0, 0, 0, 0, 0)))# extra_se estimate, borrowed from sablefish assessment
names(newQparms) <- names(ssCtl$Q_parms)
ssCtl$Q_parms <- ssCtl$Q_parms %>% bind_rows(newQparms) 

# save the new control file
SS_writectl(ctllist = ssCtl, outfile = "scenarioModels/benchmarkDFA_AsNowcast/control.ss",
            overwrite = TRUE)
# Run model ----

# Select preferred SS3 version for model fitting
sel_SS <- "v3.30.23" # indicate version of SS to run this analysis with
r4ss::run(dir = "C:/Users/r.wildermuth/Documents/CEFI/SardineRecruitmentESP/SardineRecruitIndex/scenarioModels/benchmarkDFA_AsNowcast",
          exe = file.path("C:/Users/r.wildermuth/Documents/SS3.30/ss3_win.exe"),
          skipfinished = FALSE)

dfaAsNowFit <- SS_output("C:/Users/r.wildermuth/Documents/CEFI/SardineRecruitmentESP/SardineRecruitIndex/scenarioModels/benchmarkDFA_AsNowcast",)
SS_plots(dfaAsNowFit)

# Compare to recruitment pattern from 2024 assessment
mngtDir <- "../SardineRecruitIndex/scenarioModels/Pacific sardine 2024 benchmark"

mngtAssmt2024 <- SS_output(dir = mngtDir, repfile = "Report.sso", printstats = FALSE)

compSmry <- SSsummarize(list(dfaAsNowFit, 
                             mngtAssmt2024))
dev.off()
SSplotComparisons(compSmry, legendlabels = c("dfaAsNowcast", 
                                             "mngtAssmt2024"))

# diff in estimated rec devs
recdevDiffs <- compSmry$recdevs %>% mutate(recdevDiff = model1 - model2) # w/Envt - benchmark
mean(recdevDiffs$recdevDiff)
hist(recdevDiffs$recdevDiff, breaks = 10)

bioDiffs <- compSmry$SmryBio %>% mutate(bioDiff = model1 - model2) # w/Envt - benchmark
mean(bioDiffs$bioDiff)
hist(bioDiffs$bioDiff, breaks = 10)
