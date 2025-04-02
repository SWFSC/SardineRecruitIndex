# Plot example of predictive performance at 5 yr horizon
# Created: 4/1/2025, Robert Wildermuth
# following OSAResids.R, but not using out-of-sample model predictions

library(tidyverse)
library(MARSS)

# read prepped dataset
datDFA <- read_csv("../SardineRecruitIndex/Data/recrDFAdat.csv")

# subset for sardine DFA from 1990 to 2023 - no rec devs for 2024
sardDat <- datDFA %>% filter(year %in% 1990:2023) %>%
  select(-c(NCOPsummer,
            SCOPsummer,
            anchBioSmrySeas1,
            anchBioSmrySeas2,
            # Potential redundant variables:
            # CUTI_39N, 
            # OC_LUSI_36N,
            # OC_STI_36N,
            summerSST,
            SCOPspring,
            SCOPsummerlag1))

datNames <- names(sardDat)[-1]

# transpose for MARSS formatting
sardDat <- sardDat %>% select(-year) %>% t()

load(file = "out/marssFit_1990to2023_noAnch_1trend_EqlVar.RData")

p <- 5
horizon <- 5
colsRMSE <- "sardRec"

origDat <- sardDat[,1:min(ncol(sardDat), # either full data set or ...
                          ncol(sardDat)-p+horizon)] # including 'horizon' steps ahead of peel
# need to scale data for residual calcs and 'newdata' input
origDat <- apply(origDat, 1, scale, simplify = TRUE) %>% t() 
# !RW: scales slightly different from those used to fit the DFA model 'objMARSS'
#      because of added data from 1 step ahead

# create 'YStar' by removing last observation of 'colsRMSE' variables
YStar <- origDat
YStar[colsRMSE, (ncol(YStar)-horizon+1):ncol(YStar)] <- NA

yHistFit <- tsSmooth(sardDFA, type = "ytT", interval = "prediction")

# informed innovations
yExpctInf <- predict(object = sardDFA, interval = "prediction", 
                     newdata = list(t = 1:ncol(origDat),
                                    y = YStar),
                     type = 'ytT', x0 = "use.model")
yExpctInf <- yExpctInf$pred %>% rename(pt.Inf = y)
# informed trend estimates
xExpctInf <- predict(object = sardDFA, interval = "prediction", 
                     newdata = list(t = 1:ncol(origDat),
                                    y = YStar),
                     type = 'xtT', x0 = "use.model")
xExpctInf <- xExpctInf$pred %>% rename(pt.Inf = .x)

expctInf <- bind_rows(yExpctInf, xExpctInf)

predYr <- max(expctInf$t) - horizon + 1989

expctInf %>%  
  mutate(t = t+1989) %>%
  filter(.rownames %in% c("sardRec", "X1")) %>%
  ggplot(aes(x = t, y = estimate, color = .rownames, fill = .rownames)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = `Lo 95`, ymax = `Hi 95`), alpha = 0.3) +
  geom_point(data = tsSmooth(sardDFA, type = "ytT") %>%
                      filter(.rownames == "sardRec") %>%
                      mutate(t = t+1989), 
             aes(x = t, y = y), color = "black") +
  facet_wrap(~.rownames) +
  labs(x= "Year", y = "State") +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = expctInf %>% filter(.rownames == "sardRec", 
                                              t == predYr - 1989) %>%
                            pull(pt.Inf), color = "green") +
  geom_vline(xintercept = predYr, color = "grey") +
  theme_classic()
