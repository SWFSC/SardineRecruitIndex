
########## NOTES ##########

# NOTES:
#   - Based on methodology documented in Johnson et al. 2016
#     - https://dx.doi.org/10.1016/j.fishres.2016.06.004 


# PREPARED BY: Alexander Jensen, SWFSC

############################

# Load libraries -----------

require(r4ss)
require(here)
require(tidyverse)
require(plyr)
require(ggsidekick)
require(ggridges)
require(patchwork)
require(scales)

# Calculate lag-1 autocorrelation of recdevs from 2024 benchmark -----------

bench_out <- SS_output(here::here("scenarioModels", "Pacific sardine 2024 benchmark"))
recdevs_main <- bench_out$parameters %>%
  slice(grep("Main_RecrDev_", Label)) 
recdevs_main_acf <- acf(recdevs_main$Value) 
lag1_acf <- recdevs_main_acf$acf[2]; lag1_acf
  # Lag 1 autocorrelation = 0.5744528

# Modify the benchmark Pacific sardine model ----------------

# Manual modification, using the following steps:
  # 1. Copy+paste model files from benchmark folder into 'Autocorrelated recruitment' folder
  # 2. Replace the INIT value for SR_autocorr in control.ss with the lag-1 autocorrelation (0.5744528)
  # 3. Run the modified model, extract the S-R bias adjustment parameters from the default HTML output
    # Using ss3 v3.30.22.1
    # Est. bias-adjustment values:
      # 1991.8   #_last_early_yr_nobias_adj_in_MPD (prev. 1992.5)
      # 2002.8   #_first_yr_fullbias_adj_in_MPD  (prev. 2002.9)
      # 2020.8   #_last_yr_fullbias_adj_in_MPD  (prev. 2020.8)
      # 2023.7   #_first_recent_yr_nobias_adj_in_MPD  (prev. 2023.1)
      # 0.9291  #_max_bias_adj_in_MPD (1.0 to mimic pre-2009 models)  (prev. 0.9249) 
  # 4. Replace the bias adjustment parameters in control.ss with the new recommendations
  # 5. Run the final model version
    # Using v3.30.22.1

# Generate HTML output for bias-adj. values in step #3
out <- SS_output(here::here("scenarioModels", "Autocorrelated recruitment"))
SS_plots(out)

# Summarize fit of benchmark model w/ and w/o autocorrelated recruitment --------------
  # Import model output ------------

autocor <- SS_output(here::here("scenarioModels", "Autocorrelated recruitment"))
bench <- SS_output(here::here("scenarioModels", "Pacific sardine 2024 benchmark"))

  # Confirm model convergence w/ autocorrelated recruitment (FINE) ----------

# Summarize likelihoods and convergence based on gradient
likelihoods <- autocor$likelihoods_used
final_gradient <- autocor$maximum_gradient_component
high_param_gradients <- autocor$parameters_with_highest_gradients

likelihoods; final_gradient; high_param_gradients

# Hessian matrix inverted
# If TRUE,  the Hessian is positive definite and the model has converged
autocor$log_det_hessian>0

# Parameters at bounds
# Not currently excluding deviation parameters
params <- autocor$parameters
params %>% # identify parameter values on boundaries
  filter(Phase > 0, Value == Min | Value == Max)

# Highly correlated parameter pairs
corrs <- autocor$CoVar
  # Correlations >= 0.95
corrs %>% 
  mutate(abs_corr = abs(corr)) %>%
  filter(abs_corr >= 0.95, Par..i == "Par", Par..j == "Par") %>%
  arrange(desc(abs_corr)) %>%
  select(label.i, label.j, corr)
  # Correlations >=0.8, <0.95
corrs %>% 
  mutate(abs_corr = abs(corr)) %>%
  filter(abs_corr >= 0.8, abs_corr < 0.95, Par..i == "Par", Par..j == "Par") %>%
  arrange(desc(abs_corr)) %>%
  select(label.i, label.j, corr)

# High parameter variance
stdevs <- autocor$stdtable; stdevs
high_cv_params <- params %>% # ID high parameter variance by calculating CV (%)
  filter(Phase > 0) %>%
  mutate(CV = Parm_StDev/Value*100) %>%
  filter(CV >= 100) %>% # ID parameters where CV>100% (e.g., Carvalho et al. 2021)
  arrange(desc(CV)) %>%
  select(Label, Value, Parm_StDev, CV)
high_cv_params 

  # Summarize model fit, parameter values, derived quantities ----------

# Benchmark

# Likelihoods
like_names <- c("TOTAL", "Catch", "Equil_catch", "Survey", "Length_comp",
                "Age_comp", "Recruitment", "InitEQ_Regime", "Forecast_Recruitment",
                "Parm_priors", "Parm_softbounds", "Parm_devs", "Crash_Pen")
likelihoods <- bench$likelihoods_used
like_df <- data.frame(Type="Likelihoods",
                      Component=rownames(likelihoods)[which(rownames(likelihoods) %in% like_names)],
                      values=round(likelihoods$values[which(rownames(likelihoods) %in% like_names)], 3)) %>%
  arrange(desc(values))

# Parameter estimates
parm_names <- c("NatM_Lorenzen_averageFem_GP_1", "SR_LN(R0)",
                "SR_BH_steep", "SR_sigmaR",
                "SR_regime_BLK3repl_2007", "LnQ_base_AT(2)",
                "LnQ_base_AT(2)_BLK4repl_2016", "LnQ_base_AT(2)_DEVmult_2008",
                "LnQ_base_AT(2)_DEVmult_2012", "LnQ_base_AT(2)_DEVmult_2013",
                "LnQ_base_AT(2)_DEVmult_2014", "LnQ_base_AT(2)_DEVmult_2015")
parm_ests <- bench$parameters
parm_df <- data.frame(Type="Parameters",
                      Component=rownames(parm_ests)[which(rownames(parm_ests) %in% parm_names)],
                      values=round(parm_ests$Value[which(rownames(parm_ests) %in% parm_names)],3))

# Biomass estimates
biomass_names <- c("SmryBio_2020", "SmryBio_2021", 
                   "SmryBio_2022", "SmryBio_2023",
                   "SmryBio_2024", "SmryBio_2025",
                   "SmryBio_2026")
biomass_ests <- bench$derived_quants
biomass_df <- data.frame(Type = "Summary biomass",
                         Component = substr(biomass_ests$Label[which(biomass_ests$Label %in% biomass_names)], 9 , 12),
                         values=round(biomass_ests$Value[which(biomass_ests$Label %in% biomass_names)]))

# Aggregate, save table
summ_table <- rbind(like_df, parm_df, biomass_df)
summ_table

# + Autocorrelated recruitment

summ_table_bench <- summ_table %>% 
  mutate(`2024 Benchmark` = values) %>%
  select(Type, Component, `2024 Benchmark`)

# Likelihoods - autocormark
like_names <- c("TOTAL", "Catch", "Equil_catch", "Survey", "Length_comp",
                "Age_comp", "Recruitment", "InitEQ_Regime", "Forecast_Recruitment",
                "Parm_priors", "Parm_softbounds", "Parm_devs", "Crash_Pen")
likelihoods <- autocor$likelihoods_used
like_df <- data.frame(Type="Likelihoods",
                      Component=rownames(likelihoods)[which(rownames(likelihoods) %in% like_names)],
                      values=round(likelihoods$values[which(rownames(likelihoods) %in% like_names)], 3)) %>%
  arrange(desc(values))

# Parameter estimates - autocormark
parm_names <- c("NatM_Lorenzen_averageFem_GP_1", "SR_LN(R0)",
                "SR_BH_steep", "SR_sigmaR",
                "SR_regime_BLK3repl_2007", "LnQ_base_AT(2)",
                "LnQ_base_AT(2)_BLK4repl_2016", "LnQ_base_AT(2)_DEVmult_2008",
                "LnQ_base_AT(2)_DEVmult_2012", "LnQ_base_AT(2)_DEVmult_2013",
                "LnQ_base_AT(2)_DEVmult_2014", "LnQ_base_AT(2)_DEVmult_2015")
parm_ests <- autocor$parameters
parm_df <- data.frame(Type="Parameters",
                      Component=rownames(parm_ests)[which(rownames(parm_ests) %in% parm_names)],
                      values=round(parm_ests$Value[which(rownames(parm_ests) %in% parm_names)],3))

# Biomass estimates - autocormark
biomass_names <- c("SmryBio_2020", "SmryBio_2021", 
                   "SmryBio_2022", "SmryBio_2023",
                   "SmryBio_2024", "SmryBio_2025",
                   "SmryBio_2026")
biomass_ests <- autocor$derived_quants
biomass_df <- data.frame(Type = "Summary biomass",
                         Component = substr(biomass_ests$Label[which(biomass_ests$Label %in% biomass_names)], 9 , 12),
                         values=round(biomass_ests$Value[which(biomass_ests$Label %in% biomass_names)]))

# Aggregate, save table - benchmark
summ_table_autocor <- rbind(like_df, parm_df, biomass_df) %>%
  mutate(`Autocor Recr` = values) %>%
  select(Type, Component, `Autocor Recr`)
summ_table_autocor

# Combine benchmark, base tables
summ_table_both <- full_join(summ_table_bench, summ_table_autocor)
summ_table_both

write.csv(summ_table_both, here::here("out", "bench_autocor_table_comparison.csv"), 
          row.names = FALSE)

  # Plot derived quantities (smrybio, recr) ----------------

    # Biomass -----------

smrybio_bench <- bench$derived_quants %>% slice(grep("SmryBio_2", Label)) 
smrybio_bench <- smrybio_bench %>% 
  mutate(lo = ifelse(qnorm(.025, mean = Value, sd = StdDev) >= 0,
                     qnorm(.025, mean = Value, sd = StdDev), 0),
         hi = qnorm(.975, mean = Value, sd = StdDev))
smrybio_bench$year <- as.numeric(ldply(strsplit(smrybio_bench$Label, split = "_"))$V2)
smrybio_bench$type <- "2024 Benchmark"

smrybio_autocor <- autocor$derived_quants %>% slice(grep("SmryBio_2", Label)) 
smrybio_autocor <- smrybio_autocor %>% 
  mutate(lo = ifelse(qnorm(.025, mean = Value, sd = StdDev) >= 0, 
                     qnorm(.025, mean = Value, sd = StdDev), 0),
         hi = qnorm(.975, mean = Value, sd = StdDev))
smrybio_autocor$year <- as.numeric(ldply(strsplit(smrybio_autocor$Label, split = "_"))$V2)
smrybio_autocor$type <- "Autocor Recr"

smrybio_all <- rbind(smrybio_bench, smrybio_autocor)

p1 <- ggplot(smrybio_all, aes(x = year, y = Value, group = type, color = type)) + 
  geom_point() + 
  geom_line(aes(color = type), alpha = 0.5, linewidth = 1) + 
  theme_sleek() +
  scale_y_continuous(label = comma) + 
  ylab("Summary biomass (age-1+; mt)") +
  xlab("Year") +
  theme(legend.position = "inside", legend.position.inside = c(.7, .85), 
        legend.title = element_blank()) + 
  geom_line(aes(y = lo), lty = 2, alpha = 0.5) + 
  geom_line(aes(y = hi), lty = 2, alpha = 0.5)
p1

p2 <- ggplot(smrybio_all %>% filter(year >= 2014), 
             aes(x = year, y = Value, group = type, color = type)) + 
  geom_point() + geom_line(aes(color = type), alpha = 0.5, linewidth = 1) + 
  theme_sleek() +
  scale_y_continuous(label = comma) + 
  ylab("Summary biomass (age-1+; mt)") +
  xlab("Year") +
  geom_line(aes(y = lo), lty = 2, alpha = 0.5) + 
  geom_line(aes(y = hi), lty = 2, alpha = 0.5)
p2

p1 / p2

ggsave(file = here::here("out", "bench_autocor_biomass_diff.png"), width = 7, height = 7)

    # Recruitment -----------

recs_bench <- bench$derived_quants %>% slice(grep("Recr_2", Label)) 
recs_bench <- recs_bench %>% 
  mutate(lo = ifelse(qnorm(.025, mean = Value, sd = StdDev) >= 0,
                     qnorm(.025, mean = Value, sd = StdDev), 0),
         hi = qnorm(.975, mean = Value, sd = StdDev))
recs_bench$year <- as.numeric(ldply(strsplit(recs_bench$Label, split = "_"))$V2)
recs_bench$type <- "2023 Benchmark"

recs_autocor <- autocor$derived_quants %>% slice(grep("Recr_2", Label)) 
recs_autocor <- recs_autocor %>% 
  mutate(lo = ifelse(qnorm(.025, mean = Value, sd = StdDev) >= 0,
                     qnorm(.025, mean = Value, sd = StdDev), 0),
         hi = qnorm(.975, mean = Value, sd = StdDev))
recs_autocor$year <- as.numeric(ldply(strsplit(recs_autocor$Label, split = "_"))$V2)
recs_autocor$type <- "Autocor Recr"

recs_all <- rbind(recs_bench, recs_autocor)

p1 <- ggplot(recs_all, aes(x = year, y = Value, group = type, color = type)) + 
  geom_point() + 
  geom_line(aes(color = type), alpha = 0.5, linewidth = 1) + 
  theme_sleek() +
  scale_y_continuous(label = comma) + 
  ylab("Recruits (age-0; thousands)") +
  xlab("Year") +
  theme(legend.position = "inside", legend.position.inside = c(.7, .85), 
        legend.title = element_blank()) + 
  geom_line(aes(y = lo), lty = 2, alpha = 0.5) + 
  geom_line(aes(y = hi), lty = 2, alpha = 0.5)
p1

p2 <- ggplot(recs_all %>% filter(year >= 2014), 
             aes(x = year, y = Value, group = type, color = type)) + 
  geom_point() + 
  geom_line(aes(color = type), alpha = 0.5, linewidth = 1) + 
  theme_sleek() +
  scale_y_continuous(label = comma) + 
  ylab("Recruits (age-0; thousands)") +
  xlab("Year") +
  geom_line(aes(y = lo), lty = 2, alpha = 0.5) + 
  geom_line(aes(y = hi), lty = 2, alpha = 0.5)
p2

p1 / p2

ggsave(file = here::here("out", "bench_autocor_recruits_diff.png"), width = 7, height = 7)
