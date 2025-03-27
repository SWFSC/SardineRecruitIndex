# Updating 2020 research assessment
# Created: 2/5/2025, Robert Wildermuth

library(tidyverse)
library(r4ss)
library(reshape2)

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
  ggplot(aes(x = Yr, y = log(bioSmry), color = SSmodel)) +
  # ggplot(aes(x = Yr, y = bioSmry, color = SSmodel)) +
  geom_line() +
  facet_wrap(~Seas)

res25LorenzDir <- "../SardineRecruitIndex/scenarioModels/2025_research_assessment_LorenzM"

res25Lorenz <- SS_output(dir = res25LorenzDir, repfile = "Report.sso", printstats = FALSE)


# Compare estimated parameters --------------------------------------------

# Include recruitment pattern from updated (2025) research assessment, run with 3.30.22.1
res25BloxDir <- "../SardineRecruitIndex/scenarioModels/2025_research_assessment_newBlox"

resAssmt2025Blox <- SS_output(dir = res25BloxDir, repfile = "Report.sso", printstats = FALSE)

resAssmt$parameters

# full_join(rownames_to_column(resAssmt$parameters), rownames_to_column(resAssmt2025.22$parameters), by = "rowname") %>%
full_join(resAssmt$parameters, resAssmt2025.22$parameters, by = "Label") %>%
  filter(grepl("AgeSel", Label)) %>% select(Label, Value.x, Value.y, Active_Cnt.x, Status.x, Active_Cnt.y, Status.y)

full_join(resAssmt$parameters, resAssmt2025$parameters, by = "Label") %>%
  filter(grepl("AgeSel", Label)) %>% select(Label, Value.x, Value.y, Active_Cnt.x, Status.x, Active_Cnt.y, Status.y)


#extract the age selectivity for all fleets and seasons
#since benchmark assessment only starts 2005, select more recent years
S_at_age.y <- resAssmt2025.22$ageselex %>% filter(Yr %in% 2005:2023 & Sex == 1 & Factor=="Asel2")
S_at_age.y <- mngtAssmt2024$ageselex %>% filter(Yr %in% 2005:2023 & Sex == 1 & Factor=="Asel2")
S_at_age.y <- res25Lorenz$ageselex %>% filter(Yr %in% 2005:2023 & Sex == 1 & Factor=="Asel2")
S_at_age.y <- resAssmt$ageselex %>% filter(Yr %in% 2005:2023 & Sex == 1 & Factor=="Asel2")
Saa_long <- melt(S_at_age.y, id.vars=c("Factor", "Fleet", "Yr","Seas","Sex", "Morph", "Label"),
                 variable.name = "Age", value.name = "Sa")
fleet <- 2
ggplot(Saa_long%>%filter(Fleet==fleet), aes(x = as.numeric(Age), y = Sa, color = as.factor(Yr))) +
  geom_line(linewidth = 1) + facet_wrap(~Seas, nrow=2) + theme_bw() + ggtitle(paste0("Fleet",fleet))+
  ylab("Selectivity at age") + xlab("Age") + labs(color = "Year") +
  scale_color_viridis_d(option = "H")


compSmry <- SSsummarize(list(resAssmt, 
                             resAssmt2025, 
                             mngtAssmt2024, 
                             resAssmt2025Blox, 
                             # resAssmt2025yrlyBlox, 
                             res25Lorenz))
dev.off()
SSplotComparisons(compSmry, legendlabels = c("resAssmt", 
                                             "resAssmt2025", 
                                             "mngtAssmt2024", 
                                             "resAssmt2025Blox", 
                                             # "resAssmt2025yrlyBlox", 
                                             "res25Lorenz"))
