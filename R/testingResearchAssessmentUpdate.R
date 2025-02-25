# Updating 2020 research assessment
# Created: 2/5/2025, Robert Wildermuth

library(tidyverse)
library(r4ss)

# Original research assessment 
resDir <- "C:/Users/r.wildermuth/Documents/CEFI/SardineRecruitmentESP/sardine_research_assessment"

#------------------------------------------------------------------
# look at the data provided

# resDat <- SS_readdat(paste0(resDir, "/data.ss"))

#Making plots
# pull in output
resAssmt <- SS_output(dir = resDir, repfile = "Report.sso", printstats = FALSE)
resAssmt$timeseries
head(resAssmt$parameters )

# save plot info
# SS_plots(resAssmt)

resAssmt$recruit

# read FS DFA dataset
datDFA <- read_csv("C:/Users/r.wildermuth/Documents/FutureSeas/RecruitmentIndex/recrmntDFA/recrDFAdat.csv")


compRecs <- datDFA %>% select(year, sardRec) %>%
              full_join(y = resAssmt$recruit %>% filter(era == "Main") %>% select(Yr, dev),
                        by = c("year" = "Yr")) %>% 
              filter(year > 1980)

# look at recruitment pattern from 2024 assessment
mngtDir <- "../SardineRecruitIndex/scenarioModels/Pacific sardine 2024 benchmark"

mngtAssmt2024 <- SS_output(dir = mngtDir, repfile = "Report.sso", printstats = FALSE)

# Include recruitment pattern from updated (2025) research assessment
res25Dir <- "../SardineRecruitIndex/scenarioModels/2025_research_assessment"

resAssmt2025 <- SS_output(dir = res25Dir, repfile = "Report.sso", printstats = FALSE)


# Compile and plot --------------------------------------------------------

compRecs <- compRecs %>%
              full_join(y = mngtAssmt2024$recruit %>% filter(era == "Main") %>% select(Yr, dev),
                        by = c("year" = "Yr")) %>%
              rename("resAssmt" = "dev.x",
                     "mngtAssmt2024" = "dev.y") %>% 
              full_join(y = resAssmt2025$recruit %>% filter(era == "Main") %>% select(Yr, dev),
                        by = c("year" = "Yr")) %>%
              rename("resAssmt2025" = "dev")
compRecs %>% filter(year > 2004)

compRecs %>% select(-sardRec) %>%
  pivot_longer(cols = -year, names_to = "SSmodel", values_to = "recdev") %>%
  ggplot(aes(x = year, y = recdev, color = SSmodel)) +
  geom_line()

compBio <- resAssmt$timeseries %>% filter(Era == "TIME") %>%
              select(Yr, Seas, Bio_smry) %>%
              full_join(y = mngtAssmt2024$timeseries %>% filter(Era == "TIME") %>%
                              select(Yr, Seas, Bio_smry),
                        by = c("Yr", "Seas")) %>%
              rename("resAssmt" = "Bio_smry.x",
                     "mngtAssmt2024" = "Bio_smry.y") %>% 
              full_join(y = resAssmt2025$timeseries %>% filter(Era == "TIME") %>% 
                              select(Yr, Seas, Bio_smry),
                        by = c("Yr", "Seas")) %>%
              rename("resAssmt2025" = "Bio_smry")

compBio %>% pivot_longer(cols = -c(Yr, Seas), names_to = "SSmodel", values_to = "bioSmry") %>%
  # ggplot(aes(x = Yr, y = log(bioSmry), color = SSmodel)) +
  ggplot(aes(x = Yr, y = bioSmry, color = SSmodel)) +
  geom_line() +
  facet_wrap(~Seas)
