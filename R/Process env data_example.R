
###############################
# NOTES

# The Biologically Effective Upwelling Transport Index (BEUTI), estimated at 39N,
  # was selected based on inclusion in Wildermuth et al. (in prep), its relevance
  # in informing sardine recruitment (negative correlation), and public availability
  # (obtained from https://www.integratedecosystemassessment.noaa.gov/regions/california-current/california-current-iea-indicators
  # on 9/24/24)
# Raw data were processed to calculate annual summaries, using the mean BEUTI over the
  # months of May-July, following Wildermuth et al. (in prep)

# Developed by Alex Jensen, SWFSC, NMFS, NOAA
###############################

require(tidyverse)

beuti_39 <- read.csv(here::here("Data", "cciea_OC_BEUTI_b4a8_8c3c_bd5b.csv"),
                     header = TRUE)[-1,]
beuti_39$date <- as.Date(substr(beuti_39$time, 1, 10),
                         "%Y-%m-%d")
beuti_39$yr <- lubridate::year(beuti_39$date)
beuti_39$mo <- lubridate::month(beuti_39$date)

beuti_39_summ <- beuti_39 %>%
  mutate(beuti = as.numeric(beuti)) %>%
  filter(mo %in% c(5, 6, 7), 
         yr %in% c(2005:2023)) %>%
  group_by(yr) %>%
  dplyr::summarize(beuti_yr = mean(beuti)) %>%
  ungroup()
beuti_39_summ$beuti_stand <- (beuti_39_summ$beuti_yr - mean(beuti_39_summ$beuti_yr)) / sd(beuti_39_summ$beuti_yr)
beuti_39_summ$beuti_cent <- beuti_39_summ$beuti_yr - mean(beuti_39_summ$beuti_yr)

write.csv(beuti_39_summ, here::here("Data", "beuti_39.csv"),
          row.names = FALSE)
