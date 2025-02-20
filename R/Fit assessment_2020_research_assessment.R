

# 2024 Benchmark used v3.30.22

# Set directory for SS files
base_wd <- here::here("scenarioModels", "Pacific sardine 2024 benchmark")
research_wd <- here::here("scenarioModels", "2020_research_assessment")

# Run 2024 Pacific sardine benchmark assessment, save output/plots ---------------------------
  # Run model, save results as .rds file (already run, so just importing Report.sso) ----

r4ss::run(dir = research_wd,
          exe = file.path(getwd(), "SS3 exe files", sel_SS, "ss.exe"), # need 'ss3.exe' if not using v3.30.15.09
          skipfinished = FALSE)

research_assess <- SS_output(research_wd)
saveRDS(research_assess, file.path(research_wd, "research_assess.rds"))

  # Read model output, make default plots ----

research_assess <- readRDS(file.path(research_wd, "research_assess.rds"))
SS_plots(research_assess)

# Import model results for research, 2024 benchmark, compare outputs ----------

research_assess <- readRDS(file.path(research_wd, "research_assess.rds"))
base_assess <- readRDS(file.path(base_wd, "base_assess.rds"))

  # Plot recdevs -------

research_recdevs <- research_assess$recruit %>%
  mutate(Type = "Research")
base_recdevs <- base_assess$recruit %>%
  mutate(Type = "Base")

recdevs <- rbind(research_recdevs, base_recdevs)

p1 <- ggplot(recdevs, aes(x = Yr, y = dev, group = Type, color = Type)) + 
  geom_point(aes(shape = Type)) + 
  geom_line(aes(color = Type), alpha = 0.5, linewidth = 1) + 
  theme_sleek() + 
  ylab("Recruitment Deviation") +
  xlab("Year") + 
  scale_color_manual(values = cbPalette) +
  theme(legend.position = "inside", legend.position.inside = c(.75, .75), 
        legend.title = element_blank())
p1

  # Plot Recruitment -------

# Process, combine results for summary biomass

# Research assessment
recr_research <- research_assess$derived_quants %>% slice(grep("Recr_", Label)) 
recr_research <- recr_research %>% 
  mutate(lo = qnorm(.025, mean = Value, sd = StdDev),
         hi = qnorm(.975, mean = Value, sd = StdDev))
recr_research$year <- as.numeric(ldply(strsplit(recr_research$Label, split = "_"))$V2)
recr_research[which(recr_research$lo < 0), 'lo'] <- 0
recr_research$type <- "Research"

# Base assessment
recr_base <- base_assess$derived_quants %>% slice(grep("Recr_2", Label)) 
recr_base <- recr_base %>% 
  mutate(lo = qnorm(.025, mean = Value, sd = StdDev),
         hi = qnorm(.975, mean = Value, sd = StdDev))
recr_base$year <- as.numeric(ldply(strsplit(recr_base$Label, split = "_"))$V2)
recr_base[which(recr_base$lo < 0), 'lo'] <- 0
recr_base$type <- "Base"

recr_all <- rbind(recr_research,
                  recr_base)

# Plot results for smrybiomass
p1 <- ggplot(recr_all, aes(x = year, y = Value, group = type, color = type)) + 
  geom_point(aes(shape = type)) + 
  geom_line(aes(color = type), alpha = 0.5, linewidth = 1) + 
  theme_sleek() +
  #scale_y_continuous(label = comma, limits = c(0, 2.2e5)) + 
  ylab("Recruitment (age-0; 1000s)") +
  xlab("Year") + 
  scale_color_manual(values = cbPalette) +
  theme(legend.position = "inside", legend.position.inside = c(.55, .85), 
        legend.title = element_blank())
p1

