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

# plot of rec devs with DFA index overlaid
# need points and error bars in same table
recs <- compSmry$recdevs %>% pivot_longer(cols = c(model1, model2), names_to = "Model", values_to = "estRecDev")
recsLo <- compSmry$recdevsLower %>% pivot_longer(cols = c(model1, model2), names_to = "Model", values_to = "recdevLo")
recsHi <- compSmry$recdevsUpper %>% pivot_longer(cols = c(model1, model2), names_to = "Model", values_to = "recdevHi")
recs <- recs %>% full_join(y = recsLo, by = c("Label", "Yr", "Model")) %>%
          full_join(y = recsHi, by = c("Label", "Yr", "Model")) %>%
          mutate(Model = case_when(Model == "model1" ~ "dfaAsData",
                                   Model == "model2" ~ "mngtAssmt2024"),
                 Label = "RecDev")
# add the DFA trend used as index
recs <- trendsHist %>% select(.rownames, t, .estimate, .conf.low, .conf.up) %>%
  rename(Label = .rownames,
         Yr = t,
         estRecDev = .estimate,
         recdevLo = .conf.low,
         recdevHi = .conf.up) %>%
  mutate(Model = "DFAtrend") %>%
  bind_rows(recs)

ggplot(recs, aes(x=Yr, y=estRecDev, group=Model, color=Model)) +
  geom_hline(yintercept = 0, color = "grey") +
  geom_line(linewidth = 1) +
  geom_point(position=position_dodge(0.25), size = 2)+
  geom_errorbar(aes(ymin=recdevLo, ymax=recdevHi), width=.5,
                position=position_dodge(0.25),
                linewidth = 1) +
  labs(y = "Recruitment deviations or Index", x = "Year", cex = 1.5) + xlim(1998.5, 2024) +
  theme_classic() +
  facet_wrap(~Label, ncol = 1) + 
  theme(text = element_text(size = 20))

# diff in estimated rec devs
recdevDiffs <- compSmry$recdevs %>% mutate(recdevDiff = model1 - model2) # w/Envt - benchmark
mean(recdevDiffs$recdevDiff)
hist(recdevDiffs$recdevDiff, breaks = 10)

bioDiffs <- compSmry$SmryBio %>% mutate(bioDiff = model1 - model2) # w/Envt - benchmark
mean(bioDiffs$bioDiff)
hist(bioDiffs$bioDiff, breaks = 10)
