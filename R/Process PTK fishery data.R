
# Load required libraries ----------

require(tidyverse)
require(here)
require(r4ss)

# Import data -------------------------

raw_landings <- read.csv(here::here("Data", "Sardine_landings_1981-2017.csv"))
raw_bio <- read.csv(here::here("Data", "Sardine_bio_1981-2015.csv"))

# Work up age compositions ------------
  # MexCal_S1 -------------------------

# Format, filter data, calculate total monthly landings for each model year
landings_m_y <- raw_landings %>% 
  filter(REGION %in% c("OR", "WA", "CA")) %>% # necessary to match landings to sources of age data
  mutate(MY = ifelse(MONTH %in% c(7:12), YEAR, YEAR-1),
         Season = ifelse(MONTH %in% c(7:12), 1, 2),
         `Model Y_S` = paste(MY,Season,sep="-"),
         Fishery = ifelse(REGION %in% c("OR", "WA", "BC"), "PNW", "MexCal")) %>%
  filter(MY >= 1981, SUB_POP == "NSP", 
         Season == 1, Fishery == "MexCal") %>%
  group_by(MY, MONTH) %>% # Removed Port_Area for time being, based on eqns. in briefing book based on month, age, and year only
  dplyr::summarize(L_m_y = sum(MT)) %>% 
  as.data.frame(); head(landings_m_y)

# Initially process, filter bio data
proc_bio <- raw_bio %>%
  filter(AGE != "NULL") %>%
  mutate(MY = ifelse(MONTH %in% c(7:12), YEAR, YEAR-1),
         Season = ifelse(MONTH %in% c(7:12), 1, 2),
         `Model Y_S` = paste(MY,Season,sep="-"),
         Fishery = ifelse(REGION %in% c("OR", "WA", "BC"), "PNW", "MexCal"),
         Age = as.numeric(ifelse(as.numeric(AGE) >= 8, 8, AGE))) %>%
  filter(MY >= 1981, SUB_POP == "NSP",
         Season == 1, Fishery == "MexCal") 
'/ double-check whether I have any weird ages being lumped into age 8+
proc_bio <- raw_bio %>%
  mutate(MY = ifelse(MONTH %in% c(7:12), YEAR, YEAR-1),
         Season = ifelse(MONTH %in% c(7:12), 1, 2),
         `Model Y_S` = paste(MY,Season,sep="-"),
         Fishery = ifelse(REGION %in% c("OR", "WA", "BC"), "PNW", "MexCal")) %>%
         #Age = as.numeric(ifelse(Age >= 8, 8, Age))) %>%
  filter(MY >= 1981, SUB_POP == "NSP",
         Season == 1, Fishery == "MexCal") 
'

# Summarize number and avg weight of individual fish per month/age/year
bio_m_a_y <- proc_bio %>%
  group_by(MY, Season, MONTH, Age) %>%
  dplyr::summarize(N_m_a_y = length(WT_KG),
            avg_wt_m_a_y = sum(as.numeric(WT_KG))/N_m_a_y)

# Summarize biological sample weight
bio_m_y <- proc_bio %>%
  group_by(MY, Season, MONTH) %>%
  dplyr::summarize(B_m_y = sum(as.numeric(WT_KG)),
            N_m_y = length(WT_KG))

# Create expanded framework for all year, month, age combinations, append data, calculate new m_a_y values
age_comp_grid <- expand.grid(Age = c(0:8),
                             MONTH = c(7:12),
                             MY = c(1981:2014)) %>%
  left_join(., bio_m_a_y) %>% # join with age data summarized by m, a, y
  left_join(., bio_m_y) %>% # join with age data summarized by m, y
  left_join(., landings_m_y) %>% # join with landings data summarized by m, y
  mutate(A_prop_m_a_y = (N_m_a_y * avg_wt_m_a_y) / B_m_y, # calculate age proportions by m, a, y
         L_m_a_y = A_prop_m_a_y * L_m_y, # calculate landings (kg) by m, a, y
         F_m_a_y = L_m_a_y / avg_wt_m_a_y) # calculate fish caught by m, a, y

# Summarize total fish caught for each age and year
Fcaught_a_y <- age_comp_grid %>%
  group_by(MY, Age) %>%
  dplyr::summarize(F_a_y = sum(F_m_a_y, na.rm = TRUE))

# Summarize total fish caught for each year
Fcaught_y <- age_comp_grid %>%
  group_by(MY) %>%
  dplyr::summarize(F_y = sum(F_m_a_y, na.rm = TRUE))

# Create final expanded framework for reference
age_comp_grid <- age_comp_grid %>%
  left_join(., Fcaught_a_y) %>%
  left_join(., Fcaught_y) 

# Calculate final proportion at age for each year, convert to long format  
Age_prop_prep <- Fcaught_a_y %>%
  left_join(., Fcaught_y) %>%
  mutate(P_a_y = F_a_y / F_y) %>%
  select(-F_a_y, -F_y) %>%
  pivot_wider(names_from = Age, values_from = P_a_y)

# Calculate adjusted sample size for each year, appended to final proportion at age
Nadj_year <- age_comp_grid %>%
  group_by(MY) %>%
  dplyr::summarize(N_year = sum(N_m_a_y, na.rm = TRUE)) %>%
  mutate(N_adj_year = N_year / 25)
Age_prop_prep <- Age_prop_prep %>%
  left_join(., Nadj_year)

# Format to match structure in data.ss
Age_prop_final <- Age_prop_prep %>%
  drop_na() %>%
  mutate(month = 4, fleet = 1, sex = 0, part = 0,
         Lbin_lo = -1, Lbin_hi = -1) %>%
  left_join(.,data.frame(MY=c(1981:2014), ageerr = rep(1, length(c(1981:2014))))) %>% # ensure I have the right years and ageerr values here!!
  select(MY, month, fleet, sex, part, ageerr, Lbin_lo, Lbin_hi, N_adj_year,
         "0", "1", "2", "3", "4", "5", "6", "7", "8") %>%
  as.data.frame()

# Compare 2005-2014 values to data.ss for 2024 benchmark
bench_agecomp <- r4ss::SS_readdat(here::here("scenarioModels", "Pacific sardine 2024 benchmark", "data.ss"))$agecomp %>%
  filter(fleet == 1)
years <- c(2005:2014)
for(i in 1:10){
  plot(as.numeric(Age_prop_final[which(Age_prop_final$MY == years[i]),10:18]) ~ c(0:8), xlab = "Age", 
       ylab = "Proportion", type="l", ylim = c(0,1), main = paste0("MexCal_S1 ", years[i]))
  lines(as.numeric(bench_agecomp[which(bench_agecomp$year == years[i]),10:18]) ~ c(0:8), lty=2)
  legend(x="topright", legend=c("Data Processing for Research Assessment","2024 Benchmark"), lty=c(1,2))
}

# *Plot comparison of Nsamp for all available years
par(mfrow = c(1,1))
plot(Age_prop_final$N_adj_year ~ Age_prop_final$MY, xlab = "Year",
     ylab = "N", type = "l")
lines(bench_agecomp$Nsamp ~ bench_agecomp$year, lty = 2)
legend(x="topright", legend = c("Data Processing for Research Assessment","2024 Benchmark"), lty = c(1, 2))

# Save estimates of age composition
write.csv(Age_prop_final, here::here("Data", "MexCal_S1_agecomps_research.csv"))

  # MexCal_S2 -------------------------

# Format, filter data, calculate total monthly landings for each model year
landings_m_y <- raw_landings %>% 
  filter(REGION %in% c("OR", "WA", "CA")) %>% # necessary to match landings to sources of age data
  mutate(MY = ifelse(MONTH %in% c(7:12), YEAR, YEAR-1),
         Season = ifelse(MONTH %in% c(7:12), 1, 2),
         `Model Y_S` = paste(MY,Season,sep="-"),
         Fishery = ifelse(REGION %in% c("OR", "WA", "BC"), "PNW", "MexCal")) %>%
  filter(MY >= 1981, SUB_POP == "NSP", 
         Season == 2, Fishery == "MexCal") %>%
  group_by(MY, MONTH) %>% # Removed Port_Area for time being, based on eqns. in briefing book based on month, age, and year only
  dplyr::summarize(L_m_y = sum(MT)) %>% 
  as.data.frame(); head(landings_m_y)

# Initially process, filter bio data
proc_bio <- raw_bio %>%
  filter(AGE != "NULL") %>%
  mutate(MY = ifelse(MONTH %in% c(7:12), YEAR, YEAR-1),
         Season = ifelse(MONTH %in% c(7:12), 1, 2),
         `Model Y_S` = paste(MY,Season,sep="-"),
         Fishery = ifelse(REGION %in% c("OR", "WA", "BC"), "PNW", "MexCal"),
         Age = as.numeric(ifelse(as.numeric(AGE) >= 8, 8, AGE))) %>%
  filter(MY >= 1981, SUB_POP == "NSP",
         Season == 2, Fishery == "MexCal") 

# Summarize number and avg weight of individual fish per month/age/year
bio_m_a_y <- proc_bio %>%
  group_by(MY, Season, MONTH, Age) %>%
  dplyr::summarize(N_m_a_y = length(WT_KG),
            avg_wt_m_a_y = sum(as.numeric(WT_KG))/N_m_a_y)

# Summarize biological sample weight
bio_m_y <- proc_bio %>%
  group_by(MY, Season, MONTH) %>%
  dplyr::summarize(B_m_y = sum(as.numeric(WT_KG)),
            N_m_y = length(WT_KG))

# Create expanded framework for all year, month, age combinations, append data, calculate new m_a_y values
age_comp_grid <- expand.grid(Age = c(0:8),
                             MONTH = c(1:6),
                             MY = c(1981:2014)) %>%
  left_join(., bio_m_a_y) %>% # join with age data summarized by m, a, y
  left_join(., bio_m_y) %>% # join with age data summarized by m, y
  left_join(., landings_m_y) %>% # join with landings data summarized by m, y
  mutate(A_prop_m_a_y = (N_m_a_y * avg_wt_m_a_y) / B_m_y, # calculate age proportions by m, a, y
         L_m_a_y = A_prop_m_a_y * L_m_y, # calculate landings (kg) by m, a, y
         F_m_a_y = L_m_a_y / avg_wt_m_a_y) # calculate fish caught by m, a, y

# Summarize total fish caught for each age and year
Fcaught_a_y <- age_comp_grid %>%
  group_by(MY, Age) %>%
  dplyr::summarize(F_a_y = sum(F_m_a_y, na.rm = TRUE))

# Summarize total fish caught for each year
Fcaught_y <- age_comp_grid %>%
  group_by(MY) %>%
  dplyr::summarize(F_y = sum(F_m_a_y, na.rm = TRUE))

# Create final expanded framework for reference
age_comp_grid <- age_comp_grid %>%
  left_join(., Fcaught_a_y) %>%
  left_join(., Fcaught_y) 

# Calculate final proportion at age for each year, convert to long format  
Age_prop_prep <- Fcaught_a_y %>%
  left_join(., Fcaught_y) %>%
  mutate(P_a_y = F_a_y / F_y) %>%
  select(-F_a_y, -F_y) %>%
  pivot_wider(names_from = Age, values_from = P_a_y)

# Calculate adjusted sample size for each year, appended to final proportion at age
Nadj_year <- age_comp_grid %>%
  group_by(MY) %>%
  dplyr::summarize(N_year = sum(N_m_a_y, na.rm = TRUE)) %>%
  mutate(N_adj_year = N_year / 25)
Age_prop_prep <- Age_prop_prep %>%
  left_join(., Nadj_year)

# Format to match structure in data.ss
Age_prop_final <- Age_prop_prep %>%
  drop_na() %>%
  mutate(month = 10, fleet = 2, sex = 0, part = 0,
         Lbin_lo = -1, Lbin_hi = -1) %>%
  left_join(.,data.frame(MY=c(1981:2014), ageerr = rep(1, length(1981:2014)))) %>% # ensure I have the right years and ageerr values here!!
  select(MY, month, fleet, sex, part, ageerr, Lbin_lo, Lbin_hi, N_adj_year,
         "0", "1", "2", "3", "4", "5", "6", "7", "8") %>%
  as.data.frame()

# Compare 2005-2014 values to data.ss for 2024 benchmark
bench_agecomp <- r4ss::SS_readdat(here::here("scenarioModels", "Pacific sardine 2024 benchmark", "data.ss"))$agecomp %>%
  filter(fleet == 2)
years <- c(2005:2014)
for(i in 1:10){
  plot(as.numeric(Age_prop_final[which(Age_prop_final$MY == years[i]),10:18]) ~ c(0:8), xlab = "Age", 
       ylab = "Proportion", type="l", ylim = c(0,1), main = paste0("MexCal_S1 ", years[i]))
  lines(as.numeric(bench_agecomp[which(bench_agecomp$year == years[i]),10:18]) ~ c(0:8), lty=2)
  legend(x="topright", legend=c("Data Processing for Research Assessment","2024 Benchmark"), lty=c(1,2))
}

# *Plot comparison of Nsamp for all available years
par(mfrow = c(1,1))
plot(Age_prop_final$N_adj_year ~ Age_prop_final$MY, xlab = "Year",
     ylab = "N", type = "l")
lines(bench_agecomp$Nsamp ~ bench_agecomp$year, lty = 2)
legend(x="topright", legend = c("Data Processing for Research Assessment","2024 Benchmark"), lty = c(1, 2))

# Save estimates of age composition
write.csv(Age_prop_final, here::here("Data", "MexCal_S2_agecomps_research.csv"))

# Generate estimates of landings by MY, season
landings_y <- raw_landings %>% 
  #filter(REGION %in% c("OR", "WA", "CA")) %>% # now need to include BC, MexCal catch
  mutate(MY = ifelse(MONTH %in% c(7:12), YEAR, YEAR-1),
         Season = ifelse(MONTH %in% c(7:12), 1, 2),
         `Model Y_S` = paste(MY,Season,sep="-"),
         Fishery = ifelse(REGION %in% c("OR", "WA", "BC"), "PNW", "MexCal")) %>%
  filter(MY >= 1981, SUB_POP == "NSP", 
         Season == 2, Fishery == "MexCal") %>%
  group_by(MY, MONTH) %>% # Removed Port_Area for time being, based on eqns. in briefing book based on month, age, and year only
  dplyr::summarize(L_m_y = sum(MT)) %>% 
  as.data.frame() %>%
  group_by(MY) %>%
  dplyr::summarize(L_y = sum(L_m_y))


  # PNW -------------------------------

# Format, filter data, calculate total monthly landings for each model year
# Only keep fishery ages from S1 based on Briefing Book specification
landings_m_y <- raw_landings %>% 
  filter(REGION %in% c("OR", "WA", "CA")) %>% # necessary to match landings to sources of age data
  mutate(MY = ifelse(MONTH %in% c(7:12), YEAR, YEAR-1),
         Season = ifelse(MONTH %in% c(7:12), 1, 2),
         `Model Y_S` = paste(MY,Season,sep="-"),
         Fishery = ifelse(REGION %in% c("OR", "WA", "BC"), "PNW", "MexCal")) %>%
  filter(MY >= 1981, SUB_POP == "NSP", 
         Season == 1, Fishery == "PNW") %>%
  group_by(MY, MONTH) %>% # Removed Port_Area for time being, based on eqns. in briefing book based on month, age, and year only
  dplyr::summarize(L_m_y = sum(MT)) %>% 
  as.data.frame(); head(landings_m_y)

# Initially process, filter bio data
proc_bio <- raw_bio %>%
  filter(AGE != "NULL", REGION %in% c("OR", "WA", "CA")) %>%
  mutate(MY = ifelse(MONTH %in% c(7:12), YEAR, YEAR-1),
         Season = ifelse(MONTH %in% c(7:12), 1, 2),
         `Model Y_S` = paste(MY,Season,sep="-"),
         Fishery = ifelse(REGION %in% c("OR", "WA", "BC"), "PNW", "MexCal"),
         Age = as.numeric(ifelse(as.numeric(AGE) >= 8, 8, AGE))) %>%
  filter(MY >= 1981, SUB_POP == "NSP",
         Season == 1, Fishery == "PNW") 

# Summarize number and avg weight of individual fish per month/age/year
bio_m_a_y <- proc_bio %>%
  group_by(MY, Season, MONTH, Age) %>%
  dplyr::summarize(N_m_a_y = length(WT_KG),
            avg_wt_m_a_y = sum(as.numeric(WT_KG))/N_m_a_y)

# Summarize biological sample weight
bio_m_y <- proc_bio %>%
  group_by(MY, Season, MONTH) %>%
  dplyr::summarize(B_m_y = sum(as.numeric(WT_KG)),
            N_m_y = length(WT_KG))

# Create expanded framework for all year, month, age combinations, append data, calculate new m_a_y values
age_comp_grid <- expand.grid(Age = c(0:8),
                             MONTH = c(7:12),
                             MY = c(1981:2014)) %>%
  left_join(., bio_m_a_y) %>% # join with age data summarized by m, a, y
  left_join(., bio_m_y) %>% # join with age data summarized by m, y
  left_join(., landings_m_y) %>% # join with landings data summarized by m, y
  mutate(A_prop_m_a_y = (N_m_a_y * avg_wt_m_a_y) / B_m_y, # calculate age proportions by m, a, y
         L_m_a_y = A_prop_m_a_y * L_m_y, # calculate landings (kg) by m, a, y
         F_m_a_y = L_m_a_y / avg_wt_m_a_y) # calculate fish caught by m, a, y

# Summarize total fish caught for each age and year
Fcaught_a_y <- age_comp_grid %>%
  group_by(MY, Age) %>%
  dplyr::summarize(F_a_y = sum(F_m_a_y, na.rm = TRUE))

# Summarize total fish caught for each year
Fcaught_y <- age_comp_grid %>%
  group_by(MY) %>%
  dplyr::summarize(F_y = sum(F_m_a_y, na.rm = TRUE))

# Create final expanded framework for reference
age_comp_grid <- age_comp_grid %>%
  left_join(., Fcaught_a_y) %>%
  left_join(., Fcaught_y) 

# Calculate final proportion at age for each year, convert to long format  
Age_prop_prep <- Fcaught_a_y %>%
  left_join(., Fcaught_y) %>%
  mutate(P_a_y = F_a_y / F_y) %>%
  select(-F_a_y, -F_y) %>%
  pivot_wider(names_from = Age, values_from = P_a_y)

# Calculate adjusted sample size for each year, appended to final proportion at age
Nadj_year <- age_comp_grid %>%
  group_by(MY) %>%
  dplyr::summarize(N_year = sum(N_m_a_y, na.rm = TRUE)) %>%
  mutate(N_adj_year = N_year / 25)
Age_prop_prep <- Age_prop_prep %>%
  left_join(., Nadj_year)

# Format to match structure in data.ss
Age_prop_final <- Age_prop_prep %>%
  drop_na() %>%
  mutate(month = 10, fleet = 2, sex = 0, part = 0,
         Lbin_lo = -1, Lbin_hi = -1) %>%
  left_join(.,data.frame(MY=c(1981:2014), ageerr = rep(5, length(1981:2014)))) %>% # ensure I have the right years and ageerr values here!!
  select(MY, month, fleet, sex, part, ageerr, Lbin_lo, Lbin_hi, N_adj_year,
         "0", "1", "2", "3", "4", "5", "6", "7", "8") %>%
  as.data.frame(); Age_prop_final

# Compare 2005-2014 values to data.ss for 2024 benchmark
bench_agecomp <- r4ss::SS_readdat(here::here("scenarioModels", "Pacific sardine 2024 benchmark", "data.ss"))$agecomp %>%
  filter(fleet == 3)
years <- c(2005:2014)
for(i in 1:10){
  plot(as.numeric(Age_prop_final[which(Age_prop_final$MY == years[i]),10:18]) ~ c(0:8), xlab = "Age", 
       ylab = "Proportion", type="l", ylim = c(0,1), main = paste0("MexCal_S1 ", years[i]))
  lines(as.numeric(bench_agecomp[which(bench_agecomp$year == years[i]),10:18]) ~ c(0:8), lty=2)
  legend(x="topright", legend=c("Data Processing for Research Assessment","2024 Benchmark"), lty=c(1,2))
}

# *Plot comparison of Nsamp for all available years
par(mfrow = c(1,1))
plot(Age_prop_final$N_adj_year ~ Age_prop_final$MY, xlab = "Year",
     ylab = "N", type = "l")
lines(bench_agecomp$Nsamp ~ bench_agecomp$year, lty = 2)
legend(x="topright", legend = c("Data Processing for Research Assessment","2024 Benchmark"), lty = c(1, 2))

# Save estimates of age composition
write.csv(Age_prop_final, here::here("Data", "PNW_agecomps_research.csv"))
