# Code from Alex Jensen to re-create 2024 benchmark retrospective analysis,
# modified to apply to the As-Data DFA implementation

# Necessary libs
require(r4ss)
require(tidyverse)
library(MARSS)
library(ss3diags)

# Fit As Forecast SS model ------------------------------------------------

# Function to modify control file to input late and forecast rec devs from DFA
# and re-run model

ChangeCtlRecDevs <- function(dfDFAExp, # data.frame of DFA expected rec devs
                             rmYr # retro peel directory suffix
                             ){

  # edit control file to read in the 2024 estimated rec dev
  ssCtl <- SS_readctl(file = file.path("scenarioModels/benchmarkDFA_AsForecast/Retro_runs",
                                       paste0("retro", rmYr), "control.ss"))
  
  # last year of main (calibration) rec devs
  lastMainRecs <- ssCtl$MainRdevYrLast
  
  # change last year of main rec_devs
  ssCtl$MainRdevYrLast <- lastMainRecs + rmYr
  
  # Change forecast recruitment phase to -1 to force input values for late and forecast rec devs
  ssCtl$Fcast_recr_phase <- -1
  
  # Change end year for ramp bias adjustment
  ssCtl$first_recent_yr_nobias_adj <- ssCtl$first_recent_yr_nobias_adj + rmYr
  
  # Change number of rec devs to read in
  ssCtl$N_Read_recdevs <- 2
  # provide the new rec dev for the forecast year
  newRecDevs <- dfDFAExp %>% filter(t %in% ((lastMainRecs+rmYr)+(1:2))) %>% 
    select(t, estimate) %>% 
    rename(Year = t,
           recdev = estimate)
  ssCtl$recdev_input <- newRecDevs
  
  # save the new control file
  SS_writectl(ctllist = ssCtl, outfile = file.path("scenarioModels/benchmarkDFA_AsForecast/Retro_runs",
                                                   paste0("retro", rmYr),"control.ss"),
              overwrite = TRUE)
  
  r4ss::run(dir = file.path("C:/Users/r.wildermuth/Documents/CEFI/SardineRecruitmentESP/SardineRecruitIndex/scenarioModels/benchmarkDFA_AsForecast",
                            "Retro_runs", paste0("retro", rmYr)),
            exe = "C:/Users/r.wildermuth/Documents/SS3.30/ss3_win.exe",
            skipfinished = FALSE)
}

# Fit As Nowcast SS model ------------------------------------------------

# Function to modify control file to shorten Main rec devs period and
# input late  rec devs from DFA as index
# and re-run model

ChangeIndexRecDevs <- function(dfDFAExp, # data.frame of DFA expected rec devs
                             rmYr, # retro peel directory suffix
                             lastYr = 2023 # give last year of main rec devs since SS_readctl() doesn't work
                             ){

  # # edit control file to read in the 2024 estimated rec dev
  # ssCtl <- SS_readctl(file = file.path("scenarioModels/benchmarkDFA_AsNowcast/Retro_runs",
  #                                      paste0("retro", rmYr), "control.ss"))
  # 
  # # change reference years for main rec devs
  # # last year of main (calibration) rec devs
  # lastMainRecs <- ssCtl$MainRdevYrLast
  lastMainRecs <- lastYr
  
  # change the input recruitment index
  # read in and modify data.ss file
  ssDat <- SS_readdat(file = file.path("scenarioModels/benchmarkDFA_AsNowcast/Retro_runs",
                                       paste0("retro", rmYr), "data.ss"))
  
  newRecDevs <- dfDFAExp %>% filter(t %in% ((lastMainRecs+rmYr-1):(lastMainRecs+rmYr))) %>% 
                  select(t, estimate, se) %>% 
                  rename(year = t,
                         obs = estimate,
                         se_log = se) %>%
                  mutate(month = 1,
                         index = 5)
  
  ssDat$CPUE <- ssDat$CPUE %>% filter(index != 5) %>% bind_rows(newRecDevs)
  
  # save updated data file
  SS_writedat(datlist = ssDat, outfile = file.path("scenarioModels/benchmarkDFA_AsNowcast/Retro_runs",
                                                   paste0("retro", rmYr), "data.ss"),
              overwrite = TRUE)
  
  
  r4ss::run(dir = file.path("C:/Users/r.wildermuth/Documents/CEFI/SardineRecruitmentESP/SardineRecruitIndex/scenarioModels/benchmarkDFA_AsNowcast",
                            "Retro_runs", paste0("retro", rmYr)),
            exe = "C:/Users/r.wildermuth/Documents/SS3.30/ss3_win.exe",
            skipfinished = FALSE)
}

# -------------------------------------------------------------------------

cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Settings for retrospective analysis
yrs_rm <- c(0:-5)
max_yr <- 2023

# Run the retrospective analyses
# retro(dir = "scenarioModels/benchmarkDFA_AsData/", # ID base directory
# retro(dir = "scenarioModels/Pacific sardine 2024 benchmark/", # ID base directory
retro(dir = "scenarioModels/benchmarkDFA_AsNowcast/", # ID base directory
      oldsubdir = "", # ID input files
      newsubdir = "Retro_runs", # ID new place to store retro runs within dir
      subdirstart = "retro",
      years = yrs_rm, # years relative to ending year of model
      exe = "C:/Users/r.wildermuth/Documents/SS3.30/ss3_win.exe")

# retrospective forecast for As-Forecast application ---------------------------
retro(dir = "scenarioModels/benchmarkDFA_AsForecast/", # ID base directory
      oldsubdir = "", # ID input files
      newsubdir = "Retro_runs", # ID new place to store retro runs within dir
      subdirstart = "retro",
      years = yrs_rm, # years relative to ending year of model
      exe = "C:/Users/r.wildermuth/Documents/SS3.30/ss3_win.exe")

# run retro() first to get model directories set up, but then need to modify 
# input files as in useDFAAsForecast and re-run estimation before continuing
# load DFA model fit
load(file = "out/marssFit_1990to2023_noAnch_1trend_EqlVar.RData")

# For now use naive innovations
yExpNaive <- predict(object = sardDFA, interval = "confidence", 
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

# update assumed recruitments using DFA and re-run models
for(i in yrs_rm) {
  ChangeCtlRecDevs(dfDFAExp = yExpNaive, rmYr = i)
}

# update estimated recruitments using DFA and re-run models
for(i in yrs_rm) {
  ChangeIndexRecDevs(dfDFAExp = yExpNaive, rmYr = i)
}

# Calculate Mohn's rho values ---------------------------------------------

# Extract, process retrospective results
# retroModels <- SSgetoutput(dirvec = file.path("scenarioModels/benchmarkDFA_AsData/Retro_runs",
# retroModels <- SSgetoutput(dirvec = file.path("scenarioModels/benchmarkDFA_AsForecast/Retro_runs",
retroModels <- SSgetoutput(dirvec = file.path("scenarioModels/benchmarkDFA_AsNowcast/Retro_runs",
# retroModels <- SSgetoutput(dirvec = file.path("scenarioModels/Pacific sardine 2024 benchmark/Retro_runs",
                                              paste0("retro", yrs_rm)))
retroSummary <- SSsummarize(retroModels) # summarize the model results
endyrvec <- (retroSummary[["endyrs"]]+1) + yrs_rm # create a vector of the ending year for retrospective forecasts
strtYrVec <- retroSummary[["endyrs"]] - 4 # Start years for retrospective forecasts
# Generate retrospective plots
# Default

# SSplotComparisons(retroSummary, # make plots comparing the 6 models
#                   endyrvec = endyrvec,
#                   legendlabels = paste("Data", yrs_rm, "years"),
#                   print = TRUE, # send plots to PNG file
#                   plot = TRUE, # dont plot to default graphics device
#                   plotdir = "scenarioModels/benchmarkDFA_AsData/Retro_runs")


# Manual using ggplot and custom data processing
biomass_retro <- retroSummary$quants %>%
  rename_with(~as.character(abs(yrs_rm)), 1:6) %>%
  pivot_longer(1:length(yrs_rm), names_to="Years removed") %>%
  mutate(`Years removed`=as.numeric(`Years removed`)) %>%
  filter(substr(Label,1,9)=="SmryBio_2") %>%
  mutate(rm_yr = (Yr > 2023-`Years removed`)) %>% # Identify whether year should be removed
  filter(rm_yr == FALSE)

p1 <- biomass_retro %>%  
  ggplot(aes(x = Yr, y = value, group = `Years removed`,
             color = as.character(`Years removed`))) +
  geom_point() +
  geom_line() + #theme_sleek() + 
  theme(legend.position = c(.8, .8)) +
  # scale_y_continuous(label = comma) +
  ylab("Summary biomass (age-1+ mt)") + xlab("Model year")  +
  scale_color_manual(values = cbPalette) +
  guides(color = guide_legend(title="Years removed"))
p1
# ggsave(file = "custom output/custom_retro.png", width = 7, height = 7)

# Calculate Mohn's rho for retrospective forecasts
rho_output <- SSmohnsrho(summaryoutput = retroSummary,
                         endyrvec = endyrvec,
                         startyr = strtYrVec,
                         verbose = FALSE)
rho_output$AFSC_Hurtado_SSB # Hurtado-Ferro et al. 2015 rho
rho_output$AFSC_Hurtado_Rec # Hurtado-Ferro et al. 2015 rho

# Calculate Mohn's rho for end assessment year
rho_output <- SSmohnsrho(summaryoutput = retroSummary,
                         endyrvec = endyrvec - 1,
                         startyr = strtYrVec - 1,
                         verbose = FALSE)
rho_output$AFSC_Hurtado_SSB # Hurtado-Ferro et al. 2015 rho
rho_output$AFSC_Hurtado_Rec # Hurtado-Ferro et al. 2015 rho

# hindcast cross-validation scores:

SSplotHCxval(retroSummary = retroSummary, Season = 1, subplots = "cpue", indexselect = 1)


# Compare retro to re-forecast --------------------------------------------

# Retrospective run as done above
retro4 <- "../SardineRecruitIndex/scenarioModels/benchmarkDFA_AsForecast/Retro_runs/retro-4"
retro4 <- SS_output(dir = retro4, repfile = "Report.sso", printstats = FALSE)

# Re-forecast by changing endyr, time varying blocks in control file, and benchmarks in forecast file
retro4end2019 <- "../SardineRecruitIndex/scenarioModels/benchmarkDFA_AsForecast/Retro_runs/retro-4_end2019"
retro4end2019 <- SS_output(dir = retro4end2019, repfile = "Report.sso", printstats = FALSE)

compSmry <- SSsummarize(list(retro4, 
                             retro4end2019))
dev.off()
SSplotComparisons(compSmry, legendlabels = c("retro4", 
                                             "retro4end2019"))
