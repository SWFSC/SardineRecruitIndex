# MLE dynamic factor analysis of sardine recruitment using MARSS package 
# Created: 4/4/2023, Robert Wildermuth
# Copied from https://github.com/futureseas/recrmntDFA 1/14/2025

library(tidyverse)
library(MARSS)
library(corrplot)

# read prepped dataset
datDFA <- read_csv("C:/Users/r.wildermuth/Documents/FutureSeas/RecruitmentIndex/recrmntDFA/recrDFAdat.csv")


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

# subset for sardine DFA from 1990 to 2019
sardDat <- datDFA %>% filter(year %in% 1990:2019) %>%
            # remove contemporary adult biomass with recruits, should be S2 biomass -> S1 recs
            select(-c(NOI,
                      ENSO,
                      NPGO,
                      PDOsummer, 
                      PDOspring,
                      NCOPsummer,
                      SCOPsummer,
                      sardBioSmrySeas1,
                      sardBioSmrySeas2,
                      anchBioSmrySeas1,
                      anchBioSmrySeas2,
                      copMeanSize,
                      copPropBioSize,
                      naupMeanSize,
                      naupPropBioSize,
                      ZL_NorCal,
                      ZL_SoCal,
                      age1SprAnchmeanWAA, # not enough data in time 
                      anchLarv,
                      anchYoY,
                      anchRec,
                      anchSpawnHab,
                      anchNurseHab,
                      daysAbove40pct,
                      # Potential redundant variables:
                      # CUTI_39N, 
                      # OC_LUSI_36N,
                      # OC_STI_36N,
                      summerSST,
                      HCI_R3,
                      SCOPspring,
                      SCOPsummerlag1))
                      
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
diag(Rcustom) <- c("HCI", 
                   "COP", "COP", 
                   "BEUTI", "BEUTI", 
                   "CUTI", "CUTI",
                   "LUSI", "LUSI", "LUSI",
                   "STI", "STI", "STI",
                   "RREAS",
                   "SSWI", "SSWI",
                   "CalCOFI", "CalCOFI", 
                   "WAA", "WAA",
                   "PRPOOS", "PRPOOS", #"PRPOOS", "PRPOOS", "PRPOOS", #"PRPOOS", "PRPOOS",
                   "NEMURO", "NEMURO",
                   # "SDM", "SDM", "TIME", "TIME", "SDM", "SDM",
                   "sardSDM", "sardSDM", 
                   "sardlarvSDM", 
                   "sardRec",
                   "SST",
                   "Transp", "Transp", "Transp", "Transp",
                   "Alb", "Hake")

sardDFA <- MARSS(y = sardDat, 
                    form = "dfa",
                    method = "BFGS",
                 # control = list(maxit = 10000,
                 #                conv.test.slope.tol = 0.1,
                 #                allow.degen = TRUE),
                 inits = list(x0 = matrix(1, 1, 1)),
                 z.score = TRUE,
                 model = list(# R = "diagonal and equal", # observation errors are the same
                               # R = "diagonal and unequal", # observation errors independent
                               # R = "equalvarcov", # observation errors equal and covars equal
                               # R = "unconstrained", # all observation errors independent
                               R = Rcustom,
                               m = 3) # number of latent processes
)

# save(sardDFA, file = "marssFit_1990to2019_noBio_3trend_Rcustom.RData")

# load(file = "marssFit_1990to2019_noBio_5trend_DiagEql.RData")

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

loadingsDF %>% filter(index %in% c("sardRec"), isSig) %>%
  arrange(index)


# investigate whether loadings are large/significant
# loadingsDF %>% filter(isSig, abs(est) > 0.05) # only 2 variables with moderate significant loadings on trend 5
loadingsDF %>% filter(isSig, abs(est) > 0.2)


# fits to data from pg 137 in MARSS User Guide
alpha <- 0.05 
# d <- residuals(overallDFA, type = "tT") 
# d$up <- qnorm(1- alpha / 2) * d$.sigma + d$.fitted 
# d$lo <- qnorm(alpha / 2) * d$.sigma + d$.fitted 
# ggplot(data = subset(d, name=="model" & 
#                        .rownames %in% c("sardRec", "anchRec", 
#                                         # "sardLarv", "anchLarv", 
#                                         # "anchBioSmrySeas2", "sardBioSmrySeas2",
#                                         "anchYoY"))) + 
#   geom_point(aes(t, value)) + 
#   geom_ribbon(aes(x = t, ymin = lo, ymax = up), linetype = 2, alpha = 0.2) + 
#   geom_line(aes(t, .fitted), col="blue") + 
#   facet_wrap(~.rownames) + xlab("Time Step") + ylab("Count")

histResids <- histResids %>% mutate(up = qnorm(1- alpha / 2) * .sigma + .fitted,
                                    lo = qnorm(alpha / 2) * .sigma + .fitted,
                                    model = "Local")

# Projection Model -----------------------------------------------------------

# subset from 1980 to 2019
projDat <- datDFA %>% filter(year %in% 1990:2019) %>%
  # select only projectable vars and vars of interest
  select("year", "HCI_R3", "HCI_R4", "BEUTI_33N", "BEUTI_39N", "CUTI_33N", 
         "CUTI_39N", "OC_LUSI_33N",
         "OC_LUSI_36N", "OC_LUSI_39N", "OC_STI_33N", "OC_STI_36N", "OC_STI_39N",
         "sardLarv", "anchLarv", "anchYoY", 
         "ZM_NorCal", "ZM_SoCal", "sardSpawnHab", "anchSpawnHab", "daysAbove5pct",
         "daysAbove40pct", "sardNurseHab", "anchNurseHab", "anchRec", "sardRec",
         "springSST", "summerSST", "avgNearTransspring", "avgNearTranssummer",
         "avgOffTransspring", "avgOffTranssummer") 

datNames <- names(projDat)[-1]

# transpose for MARSS formatting
projDat <- projDat %>% select(-year) %>% t()

# Create a custom R obs error matrix assuming each data source has it's own common error
RcustProj <- matrix(list(0),length(datNames),length(datNames)) 
diag(RcustProj) <- c("HCI", "HCI", 
                     "BEUTI", "BEUTI", 
                     "CUTI", "CUTI",
                     "LUSI", "LUSI", "LUSI", 
                     "STI", "STI", "STI",
                     "CalCOFI", "CalCOFI",
                     "RREAS",
                     "NEMURO", "NEMURO",
                     # "SDM", "SDM", "TIME", "TIME", "SDM", "SDM",
                     "sardSDM", "anchSDM", "sardSDM", "anchSDM", 
                     "sardlarvSDM", "anchlarvSDM",
                     # "SDM1", "SDM2", "SDM1", "SDM2", "SDM1", "SDM2",
                     "anchRec",
                     "sardRec",
                     "SST", "SST",
                     "Transp", "Transp", "Transp", 
                     "Transp")

projectDFA <- MARSS(y = projDat, 
                    form = "dfa",
                    method = "BFGS",
                    # control = list(maxit = 10000,
                    #                conv.test.slope.tol = 0.1,
                    #                allow.degen = TRUE),
                    inits = list(x0 = matrix(1, 1, 1)),
                    z.score = TRUE,
                    model = list(R = "diagonal and equal", # observation errors are the same
                      # R = "diagonal and unequal", # observation errors independent
                      # R = "equalvarcov", # observation errors equal and covars equal
                      # R = "unconstrained", # all observation errors independent
                      # R = RcustProj,
                      m = 5) # number of latent processes
)

# save(projectDFA, file = "marssFit_1990to2019_ProjDFA_5trend_Rcustom.RData")
# save(projectDFA, file = "marssFit_1990to2019_ProjDFA_3trend_DiagEql.RData")
load(file = "marssFit_1990to2019_ProjDFA_5trend_Rcustom.RData")

# calc RMSE
projResids <- residuals(projectDFA, type = "tT")

projRMSE <- projResids %>% filter(name == "model") %>%
              group_by(.rownames) %>%
              summarize(sosRes = sum(.resids^2, na.rm = TRUE),
                        nObs = sum(!is.na(.resids))) %>%
              mutate(RMSE = sqrt(sosRes/nObs))
projRMSE %>% filter(.rownames %in% c(#"anchYoY", 
                                     "anchRec", "sardRec")) %>%
  summarize(totRMSE = sum(RMSE)) %>% pull(totRMSE)

loadingsProj <- ProcessLoadings(projectDFA)

loadingsDF <- loadingsProj$loadingsDF

# investigate whether loadings are large/significant
loadingsDF %>% filter(isSig, abs(est) > 0.05) # all trends have variables with moderate significant loadings 
loadingsDF %>% filter(isSig, abs(est) > 0.2) # all trends have variables with strong significant loadings 


projResids <- projResids %>% mutate(up = qnorm(1- alpha / 2) * .sigma + .fitted,
                                    lo = qnorm(alpha / 2) * .sigma + .fitted,
                                    model = "Projection")

# fits to data from pg 137 in MARSS User Guide
# alpha <- 0.05 
# d <- residuals(projectDFA, type = "tT") 
# d$up <- qnorm(1- alpha / 2) * d$.sigma + d$.fitted 
# d$lo <- qnorm(alpha / 2) * d$.sigma + d$.fitted 

# Plots for manuscript --------------------------------------------------

# plots of model fit, physical variables
# histResids %>% filter(name=="model" &
projResids %>% filter(name=="model" &
                      .rownames %in% c("HCI_R3", "HCI_R4", "BEUTI_33N", "BEUTI_39N",
                                       "CUTI_33N", "CUTI_39N", "OC_LUSI_33N", 
                                       "OC_LUSI_36N", "OC_LUSI_39N", "ENSO", "NPGO",
                                       "PDOspring", "PDOsummer", "OC_STI_33N", 
                                       "OC_STI_36N", "OC_STI_39N", "avgSSWIspring", 
                                       "avgSSWIsummer", "sardSpawnHab", 
                                       "anchSpawnHab", "daysAbove5pct", 
                                       "daysAbove40pct", "sardNurseHab", 
                                       "anchNurseHab", "springSST", "summerSST", 
                                       "avgNearTransspring", "avgNearTranssummer", 
                                       "avgOffTransspring", "avgOffTranssummer")) %>% 
  mutate(t = t+1989) %>% 
  ggplot() +
  geom_point(aes(t, value)) +
  geom_ribbon(aes(x = t, ymin = lo, ymax = up), linetype = 2, alpha = 0.2) +
  geom_line(aes(t, .fitted), col="blue") +
  facet_wrap(~.rownames) + 
  xlab("Time Step") + ylab("Anomaly") +
  geom_hline(yintercept = 0, color = "black") +
  theme_classic()

# plots of model fit, biological variables
# histResids %>% filter(name=="model" &
projResids %>% filter(name=="model" &
                      .rownames %in% c("NCOPspring", "NCOPsummerlag1", "SCOPspring", 
                                       "SCOPsummerlag1", "swfscRockfishSurv_Myctophids",
                                       "sardLarv", "anchLarv", "mesopelLarv", 
                                       "anchYoY", "age1SprSardmeanWAA", "meanSSBwt", 
                                       "copBio", "naupBio", "ZM_NorCal", "ZM_SoCal",
                                       "anchRec", "sardRec", "albacore", "hake")) %>% 
  mutate(t = t+1989) %>% 
  ggplot() +
  geom_point(aes(t, value)) +
  geom_ribbon(aes(x = t, ymin = lo, ymax = up), linetype = 2, alpha = 0.2) +
  geom_line(aes(t, .fitted), col="blue") +
  facet_wrap(~.rownames) + 
  xlab("Time Step") + ylab("Anomaly") +
  geom_hline(yintercept = 0, color = "black") +
  theme_classic()

# plots of model estimated latent trends
# histResids %>% filter(name=="state") %>% 
projResids %>% filter(name=="state") %>% 
  mutate(t = t+1989) %>% 
  ggplot() +
  geom_point(aes(t, value)) +
  geom_ribbon(aes(x = t, ymin = lo, ymax = up), linetype = 2, alpha = 0.2) +
  geom_line(aes(t, .fitted), col="blue") +
  facet_wrap(~.rownames) + 
  xlab("Time Step") + ylab("Anomaly") +
  geom_hline(yintercept = 0, color = "black") +
  theme_classic()

# plot of model fits for variables of interest
comResids <- bind_rows(projResids, histResids)

comResids %>% filter(name=="model" &
                       .rownames %in% c("sardRec", "anchRec",
                                        # "sardLarv", "anchLarv",
                                        # "anchBioSmrySeas2", "sardBioSmrySeas2",
                                        "anchYoY")) %>% 
  mutate(t = t+1989) %>% 
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
varArrang <- varArrang[-which(varArrang %in% c("sardRec", "anchRec", "anchYoY",
                                               "sardLarv", "anchLarv"))]
loadingsDF <- loadingsDF %>% mutate(index = factor(index, 
                                                   level = c("sardRec", "anchRec", 
                                                             "sardLarv", "anchLarv",
                                                             "anchYoY",
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
                             index %in% c("HCI_R3", "HCI_R4", "sardSpawnHab", "anchSpawnHab",                
                                          "daysAbove5pct", "daysAbove40pct",
                                          "springSST", "summerSST") ~ "Temperature", #"#00BA38",
                             index %in% c("BEUTI_33N", "BEUTI_39N", "CUTI_33N",
                                          "CUTI_39N", "OC_LUSI_33N","OC_LUSI_36N",
                                          "OC_LUSI_39N", "OC_STI_33N", "OC_STI_36N",
                                          "OC_STI_39N", "swfscRockfishSurv_Myctophids",
                                          "avgSSWIspring", "avgSSWIsummer",
                                          "mesopelLarv", "copBio", "naupBio",
                                          "sardNurseHab", "anchNurseHab",
                                          "avgNearTransspring", "avgNearTranssummer",
                                          "avgOffTransspring", "avgOffTranssummer") ~ "Foraging",#"#F8766D",
                             index %in% c("albacore", "hake") ~ "Predation",#"#619CFF",
                             TRUE ~ "Interest Var" ),
         colCode = as.factor(colCode),
              hypoth =  case_when(trend == 1 ~ "Upwelling Strength",
                                  trend == 2 ~ "Upwelling Timing",
                                  trend == 3 ~ "Preconditioning",
                                  # trend == 3 ~ "LTL Conditions",
                                  trend == 4 ~ "Advection",
                                  trend == 5 ~ "Spawning Conditions"),
         labl = paste0("Trend ", trend, ": ", hypoth)) 

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
  geom_hline(yintercept = 5.5, color = "black") +
  theme_classic() +
  facet_wrap(~labl, nrow = 1) +
  geom_text(x = .8, color = "black", 
            label = ifelse(test1$isSig & abs(test1$est) > 0.05, "*", "")) +
  guides(linewidth = "none",
         color = guide_legend(override.aes = list(linewidth = 4)))

loadingsDF %>% filter(index %in% c("sardRec", "anchRec", "anchYoY",
                                   "sardLarv", "anchLarv"))

# check correlations between zooplankton indicators
loadingsDF %>% filter(index %in% c("NCOPspring", "SCOPspring", "ZM_NorCal"))
loadingsDF %>% filter(index %in% c("copBio", "naupBio", "ZM_SoCal"))

# check correlations between advection indicators
loadingsDF %>% filter(index %in% c("avgSSWIspring", "avgNearTransspring","avgOffTransspring"))
loadingsDF %>% filter(index %in% c("avgSSWIsummer", "avgNearTranssummer", "avgOffTranssummer"))


trendsAll <- tsSmooth(overallDFA, type = "xtT", interval = "confidence") %>%
                mutate(model = "Local")
trendsAll <- tsSmooth(projectDFA, type = "xtT", interval = "confidence") %>%
                mutate(model = "Project") %>%
                bind_rows(trendsAll)

trendsAll %>%  
  mutate(t = t+1989) %>%
  ggplot(aes(x = t, y = .estimate, color = model, fill = model)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = .conf.low, ymax = .conf.up), alpha = 0.3) +
  facet_grid(cols = vars(model), rows = vars(.rownames)) +
  labs(x= "Year", y = "State") +
  geom_hline(yintercept = 0) +
  theme_classic()

trendsAll <- trendsAll %>% mutate(hypoth = case_when(.rownames == "X1" & model == "Local" ~ "Upwelling Strength",
                                        .rownames == "X2" & model == "Local" ~ "Upwelling Timing",
                                        .rownames == "X3" & model == "Local" ~ "Preconditioning",
                                        .rownames == "X4" & model == "Local" ~ "Advection",
                                        .rownames == "X5" & model == "Local" ~ "Spawning Conditions", # new term?
                                        .rownames == "X1" & model == "Project" ~ "Upwelling Strength",
                                        .rownames == "X2" & model == "Project" ~ "Upwelling Timing",
                                        .rownames == "X3" & model == "Project" ~ "LTL Conditions",
                                        .rownames == "X4" & model == "Project" ~ "Advection",
                                        .rownames == "X5" & model == "Project" ~ "Spawning Conditions",
                                        TRUE ~ NA),
                     hypoth = factor(hypoth, 
                          level = c("Upwelling Strength", "Upwelling Timing",
                                    "Preconditioning", "Advection", "Spawning Conditions", "LTL Conditions",
                                    "Trophic Community")))

# invertTrends <- trendsAll %>% 
#                   mutate(invEst = case_when(model == "Project" & .rownames %in% c("X3") ~ -.estimate,
#                                             TRUE ~ .estimate),
#                          invLow = case_when(model == "Project" & .rownames %in% c("X3") ~ -.conf.low,
#                                             TRUE ~ .conf.low),
#                          invHi = case_when(model == "Project" & .rownames %in% c("X3") ~ -.conf.up,
#                                             TRUE ~ .conf.up))

trendsAll %>%  
  mutate(t = t+1989) %>%
  ggplot(aes(x = t, y = .estimate, color = model, fill = model)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = .conf.low, ymax = .conf.up), alpha = 0.3) +
  facet_wrap(~hypoth) +
  labs(x= "Year", y = "State") +
  geom_hline(yintercept = 0) +
  theme_classic()
