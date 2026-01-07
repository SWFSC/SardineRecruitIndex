# MLE dynamic factor analysis of sardine recruitment using MARSS package 
# Created: 4/4/2023, Robert Wildermuth
# Copied from https://github.com/futureseas/recrmntDFA 1/14/2025

library(tidyverse)
library(MARSS)
library(corrplot)

# read prepped dataset
datDFA <- read_csv("../SardineRecruitIndex/Data/recrDFAdat.csv")
# load(file = "Data/indicatorSetNames_LUSI39spawnHabsprSST.RData")
# load(file = "Data/indicatorSetNames_STI39spawnHabHCI.RData")
# load(file = "Data/indicatorSetNames_LUSI39spawnHabHCI.RData")
load(file = "Data/indicatorSetNames_STI39spawnHabsprSST.RData")
# Function to process loadings from MARSS output --------------------------

ProcessLoadings <- function(outMARSS, ...){
  
  # Add CIs to marssMLE object 
  outMARSS <- MARSSparamCIs(outMARSS, ...)
  
  # Look at factor loadings
  # get the inverse of the rotation matrix 
  Z.est <- coef(outMARSS, type = "matrix")$Z
  H.inv <- 1 
  if (ncol(Z.est) > 1){
    H.inv <- varimax(coef(outMARSS, type = "matrix")$Z)$rotmat
  } 
  
  # rotate factor loadings 
  Z.rot <- Z.est %*% H.inv 
  
  # Use coef() to get the upper and lower CIs 
  Z.low <- coef(outMARSS, type = "Z", what = "par.lowCI") 
  Z.up <- coef(outMARSS, type = "Z", what = "par.upCI") 
  Z.rot.up <- Z.up %*% H.inv 
  Z.rot.low <- Z.low %*% H.inv 
  
  # rotate trends 
  trends.rot <- solve(H.inv) %*% outMARSS$states
  # get ts of trends
  ts.trends <- t(trends.rot)
  
  # new df with coordinates
  loadingsDF <- data.frame(est = as.vector(Z.rot), 
                           conf.up = as.vector(Z.rot.up), 
                           conf.low = as.vector(Z.rot.low),
                           trend = rep(1:outMARSS$call$model$m, each = nrow(Z.rot)), 
                           index = rownames(Z.rot),
                           dummy0 = 0)
  
  loadingsDF$isSig <- sign(loadingsDF$conf.up) == sign(loadingsDF$conf.low)
  
  return(list(loadingsDF = loadingsDF, trendTS = ts.trends))
}

# Sardine only model -----------------------------------------------------------

# subset for sardine DFA from 1985 to 2021 - no rec devs for 2024
# leave out 2022 and 2023 b/c rec devs poorly estimated by SS model
sardDat <- datDFA %>% filter(year %in% 1985:2021) %>%
            select(all_of(setNames)) #!!RW: pull 'setNames' from IndicatorSelection.R for now
            
datNames <- names(sardDat)[-1]

# transpose for MARSS formatting
sardDat <- sardDat %>% select(-year) %>% t()

datZscore <- zscore(sardDat)
corrMat <- cor(t(datZscore), use = "pairwise.complete.obs")
pTest <- cor.mtest(t(datZscore), alternative = "two.sided", method = "pearson")
corrplot(corrMat, p.mat = pTest$p, sig.level = 0.05, insig = "blank",
         order = 'hclust', hclust.method = "ward.D2", #"centroid", #"single", #
         tl.col = 'black', type = "lower",
         cl.ratio = 0.1, tl.srt = 45, tl.cex = 0.6, #mar = c(0.1, 0.1, 0.1, 0.1), 
         addrect = 6, rect.col = "green", diag = FALSE)

# Create a custom R obs error matrix assuming each data source has it's own common error
Rcustom <- matrix(list(0),length(datNames),length(datNames))
diag(Rcustom) <- c("COP", "COP", 
                   "BEUTI", "BEUTI", 
                   "CUTI", "CUTI",
                   "LUSI",
                   "STI",
                   "RREAS", "RREAS", 
                   "WAA", "WAA",
                   "NEMURO", "NEMURO",
                   "sardRec",
                   "anchBio",
                   "SST",
                   "Transp", "Transp", "Transp", "Transp",
                   "condK",
                   "LUSI",
                   "sardSDM", 
                   "sardlarvSDM", 
                   "SST")

# number of trends
m <- 1

sardDFA <- MARSS(y = sardDat, 
                    form = "dfa",
                    method = "BFGS",
                 # control = list(maxit = 10000,
                 #                conv.test.slope.tol = 0.1,
                 #                allow.degen = TRUE),
                 inits = list(x0 = matrix(1, 1, 1)),
                 z.score = TRUE,
                 model = list( R = "diagonal and equal", # observation errors are the same
                               # R = "diagonal and unequal", # observation errors independent
                               # R = "equalvarcov", # observation errors equal and covars equal
                               # R = "unconstrained", # all observation errors independent
                               # R = Rcustom,
                               m = m) # number of latent processes
)


# save(sardDFA, file = "out/marssFit_3trendDiagEq1985_2021Anch_SardRec.RData")
# save(sardDFA, file = "out/marssFit_1trendDiagEq1985_2021Anch_SardRec.RData")

load(file = "out/marssFit_3trendDiagEq1985_2021Anch_SardRec.RData")

# calc RMSE
histResids <- residuals(sardDFA, type = "tT")

histRMSE <- histResids %>% filter(name == "model") %>%
  group_by(.rownames) %>%
  summarize(sosRes = sum(.resids^2, na.rm = TRUE),
            nObs = sum(!is.na(.resids))) %>%
  mutate(RMSE = sqrt(sosRes/nObs))
histRMSE %>% filter(.rownames %in% c("sardRec")) %>%
  summarize(totRMSE = sum(RMSE)) %>% pull(totRMSE)

loadingsHist <- ProcessLoadings(sardDFA)

loadingsDF <- loadingsHist$loadingsDF

loadingsDF %>% filter(index %in% c("sardRec")) %>%
  arrange(index, isSig)


# investigate whether loadings are large/significant
# loadingsDF %>% filter(isSig, abs(est) > 0.05) # only 2 variables with moderate significant loadings on trend 5
loadingsDF %>% filter(isSig, abs(est) > 0.2) %>% arrange(abs(est))

# Random indicator threshold loadings for 4 trend model with equal error variances
# index    trend meanLoading
# <chr>    <dbl>       <dbl>
# 1 randTest     1       0.197
# 2 randTest     3       0.197
# 3 randTest     2       0.204
# 4 randTest     4       0.228
loadingsDF %>% filter(isSig, trend == 1 & abs(est) > 0.197 |
                              trend == 2 & abs(est) > 0.197 | # no indicator passes threshold
                              trend == 3 & abs(est) > 0.204 |
                              trend == 4 & abs(est) > 0.228) %>% arrange(trend, abs(est))

# Random indicator threshold loadings for 4 trend model with equal error variances
# index    trend meanLoading
# <chr>    <dbl>       <dbl>
# 1 randTest     2      0.0839
# 2 randTest     1      0.172 
# 3 randTest     3      0.236
loadingsDF %>% filter(isSig, trend == 1 & abs(est) > 0.172 |
                        trend == 2 & abs(est) > 0.0839 | 
                        trend == 3 & abs(est) > 0.236 ) %>% arrange(trend, abs(est))

# Random indicator threshold loadings for 4 trend model with equal error variances
# index    trend meanLoading
# <chr>    <dbl>       <dbl>
#   1 randTest     1       0.288
loadingsDF %>% filter(isSig, trend == 1 & abs(est) > 0.288) %>% arrange(trend, abs(est))

# look at most influential indicators with significant loadings
# significant sardine loadings
loadingsDF %>% filter(trend %in% c(2), isSig, trend == 2 & abs(est) > 0.0839) %>%
  group_by(index) %>%
  summarize(cummLoading = sum(abs(est))) %>%
  arrange(desc(cummLoading)) %>% print(n=45)
# 
# # all strong sardine loadings
loadingsDF %>% filter(trend %in% c(1,2), trend == 1 & abs(est) > 0.172 |
                        trend == 2 & abs(est) > 0.0839) %>%
  group_by(index) %>%
  summarize(cummLoading = sum(abs(est))) %>%
  arrange(desc(cummLoading)) %>% print(n=45)

alpha <- 0.05 
histResids <- histResids %>% mutate(up = qnorm(1- alpha / 2) * .sigma + .fitted,
                                    lo = qnorm(alpha / 2) * .sigma + .fitted,
                                    model = "Local")



# Plots --------------------------------------------------



# plots of model fit, physical variables
histResids %>% filter(name=="model" &
# projResids %>% filter(name=="model" &
                      .rownames %in% c("HCI_30N355N", "BEUTI_33N", "BEUTI_39N",
                                       "CUTI_33N", "CUTI_39N", "OC_LUSI_33N", 
                                       "OC_LUSI_36N", "OC_LUSI_39N", "OC_STI_33N", 
                                       "OC_STI_36N", "OC_STI_39N", "avgSSWIspring", 
                                       "avgSSWIsummer", "sardSpawnHab", 
                                       "daysAbove5pct", "sardNurseHab", 
                                       "springSST", "summerSST", 
                                       "avgNearTransspring", "avgNearTranssummer", 
                                       "avgOffTransspring", "avgOffTranssummer")) %>% 
  mutate(t = t+1984) %>% 
  ggplot() +
  geom_point(aes(t, value)) +
  geom_ribbon(aes(x = t, ymin = lo, ymax = up), linetype = 2, alpha = 0.2) +
  geom_line(aes(t, .fitted), col="blue") +
  facet_wrap(~.rownames) + 
  xlab("Time Step") + ylab("Anomaly") +
  geom_hline(yintercept = 0, color = "black") +
  theme_classic()

# plots of model fit, biological variables
histResids %>% filter(name=="model" &
# projResids %>% filter(name=="model" &
                      .rownames %in% c("NCOPspring", "NCOPsummerlag1", 
                                       "RREAS_Myctophids", "RREAS_YOYsardine",
                                       "sardLarv", "mesopelLarv", 
                                       "anchYoY", "age1SprSardmeanWAA", "meanSSBwt", 
                                       "C.pacificus", "sardRec", "anchBioSmrySeas1", 
                                       "yoySardSL")) %>% 
  mutate(t = t+1984) %>% 
  ggplot() +
  geom_point(aes(t, value)) +
  geom_ribbon(aes(x = t, ymin = lo, ymax = up), linetype = 2, alpha = 0.2) +
  geom_line(aes(t, .fitted), col="blue") +
  facet_wrap(~.rownames) + 
  xlab("Time Step") + ylab("Anomaly") +
  geom_hline(yintercept = 0, color = "black") +
  theme_classic()

# plots of model estimated latent trends
histResids %>% filter(name=="state") %>%
# projResids %>% filter(name=="state") %>% 
  mutate(t = t+1984) %>% 
  ggplot() +
  geom_point(aes(t, value)) +
  geom_ribbon(aes(x = t, ymin = lo, ymax = up), linetype = 2, alpha = 0.2) +
  geom_line(aes(t, .fitted), col="blue") +
  facet_wrap(~.rownames) + 
  xlab("Time Step") + ylab("Anomaly") +
  geom_hline(yintercept = 0, color = "black") +
  theme_classic()

# plot of model fits for variables of interest
histResids %>% filter(name=="model" &
                       .rownames %in% c("sardRec", "RREAS_YOYsardine")) %>% 
  mutate(t = t+1984) %>% 
  ggplot() +
  geom_point(aes(t, value)) +
  geom_ribbon(aes(x = t, ymin = lo, ymax = up), linetype = 2, alpha = 0.2) +
  geom_line(aes(t, .fitted), col="blue") +
  facet_grid(rows = vars(model), cols = vars(.rownames)) + 
  xlab("Time Step") + ylab("Anomaly") +
  geom_hline(yintercept = 0, color = "black") +
  theme_classic()
 

# plot loadings

# order loadings by magnitude and arrange for plotting
varArrang <- loadingsDF %>% filter(trend == 1) %>% arrange(est) %>% pull(index)
# leave out response variables
varArrang <- varArrang[-which(varArrang %in% c("sardRec", "RREAS_YOYsardine", "sardLarv"))]
loadingsDF <- loadingsDF %>% mutate(index = factor(index, 
                                                   level = c("sardRec", 
                                                             "RREAS_YOYsardine", 
                                                             "sardLarv",
                                                      varArrang)))

myCols <- c("#F8766D", "black","#FFB000", "#619CFF", 
            "#00BA38")
names(myCols) <- levels(c("Foraging", "Interest Var", "Preconditioning",  "Predation",
                          "Temperature"
                          ))

test1 <- loadingsDF %>% mutate(est = case_when(abs(est) < 0.05 ~ 0,
                                                TRUE ~ est),
         colCode = case_when(index %in% c("age1SprSardmeanWAA", "meanSSBwt",
                                          "NCOPspring", "NCOPsummerlag1",                  
                                          "SCOPspring", "SCOPsummerlag1",
                                          "ZM_NorCal", "ZM_SoCal") ~"Preconditioning",#  "#FFB000",
                             index %in% c("HCI_30N355N", "sardSpawnHab", "anchSpawnHab",                
                                          "daysAbove5pct", "daysAbove40pct",
                                          "springSST", "summerSST") ~ "Temperature", #"#00BA38",
                             index %in% c("BEUTI_33N", "BEUTI_39N", "CUTI_33N",
                                          "CUTI_39N", "OC_LUSI_33N","OC_LUSI_36N",
                                          "OC_LUSI_39N", "OC_STI_33N", "OC_STI_36N",
                                          "OC_STI_39N", "RREAS_Myctophids",
                                          "avgSSWIspring", "avgSSWIsummer",
                                          "mesopelLarv", "C.pacificus",
                                          "sardNurseHab", "anchNurseHab",
                                          "avgNearTransspring", "avgNearTranssummer",
                                          "avgOffTransspring", "avgOffTranssummer") ~ "Foraging",#"#F8766D",
                             index %in% c("yoySardSL", "anchBioSmrySeas1") ~ "Predation",#"#619CFF",
                             TRUE ~ "Interest Var" ),
         colCode = as.factor(colCode),
              # hypoth =  case_when(trend == 1 ~ "Trend 1",
              #                     trend == 2 ~ "Trend 2",
              #                     trend == 3 ~ "Trend 3",
              #                     trend == 4 ~ "Trend 4",
              #                     trend == 5 ~ "Trend 5"),
         labl = paste0("Trend ", trend)) #, ": ", hypoth)) 

test1 %>%
  ggplot(aes(y = index, color = colCode)) +
  geom_segment(aes(x = dummy0,
                   yend = index,
                   xend = est,
                   linewidth = 4)) +
  scale_color_manual(values = myCols) +
  # scale_color_manual(values = c("#FFB000", "#00BA38", "#F8766D", "#619CFF", "black"),
  #                    labels = c("Preconditioning", "Temperature",
  #                               "Foraging", "Predation", "Interest Var")) +
  labs(x = "Loadings", y = "Index", color = "Hypothesis") +
  geom_vline(xintercept = 0, color = "grey") +
  geom_hline(yintercept = 2.5, color = "black") +
  theme_classic() +
  facet_wrap(~labl, nrow = 1) +
  geom_text(x = .7, color = "black", 
            label = ifelse(test1$isSig & abs(test1$est) > 0.05, "*", "")) +
  guides(linewidth = "none",
         color = guide_legend(override.aes = list(linewidth = 4)))


trendsAll <- tsSmooth(sardDFA, type = "xtT", interval = "confidence") %>%
                mutate(model = "Local")
# trendsAll <- tsSmooth(projectDFA, type = "xtT", interval = "confidence") %>%
#                 mutate(model = "Project") %>%
#                 bind_rows(trendsAll)

trendsAll %>%  
  mutate(t = t+1984) %>%
  ggplot(aes(x = t, y = .estimate, color = model, fill = model)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = .conf.low, ymax = .conf.up), alpha = 0.3) +
  facet_grid(cols = vars(model), rows = vars(.rownames)) +
  labs(x= "Year", y = "State") +
  geom_hline(yintercept = 0) +
  theme_classic()



# Test against random time series -----------------------------------------

# record estimated loading and significance for random variable over 100 iterations
randLoading <- tibble(est = 0, conf.up = 0, conf.low = 0, trend = 0, index = "", 
                      dummy0 = 0, isSig = 0)

for(ii in 1:100){
  # subset for sardine DFA from 1985 to 2021 - years 2022-2024 not well estimated
  sardDat <- datDFA %>% filter(year %in% 1985:2021) %>%
    select(all_of(setNames))
  
  # add random vector to test strength of loadings relative to random var
  sardDat <- sardDat %>% mutate(randTest = rnorm(n = nrow(sardDat)))
  
  datNames <- names(sardDat)[-1]
  
  # transpose for MARSS formatting
  sardDat <- sardDat %>% select(-year) %>% t()
  
  sardDFA <- MARSS(y = sardDat, 
                   form = "dfa",
                   method = "BFGS",
                   # control = list(maxit = 10000,
                   #                conv.test.slope.tol = 0.1,
                   #                allow.degen = TRUE),
                   inits = list(x0 = matrix(1, 1, 1)),
                   z.score = TRUE,
                   model = list( R = "diagonal and equal", # observation errors are the same
                                 m = m) # number of latent processes
  )
  
  loadingsHist <- ProcessLoadings(sardDFA)
  
  loadingsDF <- loadingsHist$loadingsDF
  
  randLoading <- loadingsDF %>% #filter(index == "randTest") %>% 
                    bind_rows(randLoading)
}
randLoading %>% group_by(index, trend) %>% 
  filter(isSig == TRUE, index == "randTest") %>%
  summarize(meanLoading = mean(abs(est))) %>%
  arrange(meanLoading) %>% print(n=34)
# write_csv(randLoading, file = "out/randLoadings_3trendDiagEq1985_2021Anch.csv")

randLoading %>% filter(index == "randTest", isSig == 1) %>% pull(est) %>% abs() %>% summary()
# most absolute loadings < 0.3 no better than random!