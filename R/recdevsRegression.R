# fit GAM to DFA input datasets with high loadings for sardine recruitment time series
# Created: 5/27/2025, Robert Wildermuth

library(MARSS)
library(r4ss)
library(tidyverse)
library(mgcv)
library(corrplot)
# library(Effects)
library(klaR)
library(caret)
library(fuzzySim)

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


# load DFA model fit
load(file = "out/marssFit_1990to2023_noAnch_noSardRec_1trend_EqlVar.RData")


loadingsHist <- ProcessLoadings(sardDFA)

loadingsDF <- loadingsHist$loadingsDF

# investigate whether loadings are large/significant
# loadingsDF %>% filter(isSig, abs(est) > 0.05) # only 2 variables with moderate significant loadings on trend 5
loadingsDF %>% filter(isSig, abs(est) > 0.2) %>% arrange(abs(est))


# extract estimated historical trends
trendsHist <- tsSmooth(sardDFA, type = "xtT", interval = "confidence") %>%
  mutate(model = "Local",
         t = t+1989)

# Simple regression of DFA trend on rec devs -------------------------------------

mngtBench2024 <- SS_output("C:/Users/r.wildermuth/Documents/CEFI/SardineRecruitmentESP/SardineRecruitIndex/scenarioModels/Pacific sardine 2024 benchmark",)
mngt2024recdevs <- mngtBench2024$recruit %>% filter(era == "Main")

regrDat <- mngt2024recdevs %>% select(Yr, dev) %>%
  left_join(y = trendsHist, by = c("Yr" = "t"))

regrFit <- lm(dev ~ .estimate, data = regrDat)
summary(regrFit) # no sig relationship

# plot(regrFit)

regrDat <- regrDat %>% mutate(predPts = predict(regrFit, newdata = regrDat))

regrDat %>% ggplot(aes(x = .estimate, y = dev)) +
  geom_point() + 
  geom_line(aes(y = predPts)) +
  theme_classic()

# see what it looks like with research model time series
resAssess2025 <- SS_output("C:/Users/r.wildermuth/Documents/CEFI/SardineRecruitmentESP/SardineRecruitIndex/scenarioModels/2025_research_assessment_LorenzM",)
res2025recdevs <- resAssess2025$recruit %>% filter(era == "Main")

regrDat <- res2025recdevs %>% select(Yr, dev) %>%
  left_join(y = trendsHist, by = c("Yr" = "t"))

regrFit <- lm(dev ~ .estimate, data = regrDat)
summary(regrFit) # sig relationship over longer historical period

# plot(regrFit)

regrDat <- regrDat %>% mutate(predPts = predict(regrFit, newdata = regrDat),
                              # color points by era
                              colr = case_when(Yr >= 2005 ~ "benchmark",
                                               TRUE ~ "DFA"))


regrDat %>% ggplot(aes(x = .estimate, y = dev)) +
  geom_point(aes(color = colr)) + 
  geom_line(aes(y = predPts)) +
  theme_classic()


# GAM regression with DFA variables ---------------------------------------

# read prepped dataset
datDFA <- read_csv("../SardineRecruitIndex/Data/recrDFAdat.csv")

# use top 10 strong, significant, non-duplicate loading variables
# datGAM <- datDFA %>% dplyr::select(C.pacificus, CUTI_39N, OC_LUSI_39N,
#                             NCOPsummerlag1, OC_STI_39N, BEUTI_39N,
#                             sardNurseHab, HCI_30N355N, NCOPspring,
#                             springSST, year)
datGAM <- datDFA %>% dplyr::select(-c("sardRec", "anchBioSmrySeas1", "anchBioSmrySeas2",
                               "avgSSWIspring", "avgSSWIsummer"))#, # only 6 obs during assessment period (post-2005)
                               # "age1SprSardmeanWAA", "meanSSBwt", # shift in sampling season leads to gap post-2017
                               # "sardLarv", "mesopelLarv", "yoySardSL")) # also remove vars with long sampling lags
datGAM <- mngt2024recdevs %>% dplyr::select(Yr, dev) %>% left_join(y = datGAM,
                               by = c("Yr"="year"))
# test model selection with the longer rec dev dataset instead
# datGAM <- res2025recdevs %>% dplyr::select(Yr, dev) %>% 
#             left_join(y = datGAM, by = c("Yr"="year"))
datGAM <-cbind(datGAM[,1:2],
               apply(X = datGAM[, -(1:2)], MARGIN = 2, FUN = zscore))
# add lagged rec dev for autocorrelation structure
datGAM <- datGAM %>% mutate(devLag1 = c(1, datGAM$dev[1:(nrow(datGAM)-1)]))

# investigate correlations ---------------------------
corrMat <- cor(datGAM, use = "pairwise.complete.obs")
pTest <- cor.mtest(datGAM, alternative = "two.sided", method = "pearson")
corrplot(corrMat, p.mat = pTest$p, sig.level = 0.05, insig = "blank",
         order = 'hclust', hclust.method = "ward.D2", #"centroid", #"single", #
         tl.col = 'black', type = "lower",
         cl.ratio = 0.1, tl.srt = 45, tl.cex = 0.6, #mar = c(0.1, 0.1, 0.1, 0.1), 
         addrect = 6, rect.col = "green", diag = FALSE)

# clusters variables based on correlation structure - members at distant branches are less correlated
corrClusts <- corclust(x = datGAM)
plot(corrClusts)

# identifies which variables to remove to reduce multicolinearity
findCorrelation(x = corrMat, cutoff = 0.8, names = TRUE)

# similar, but can use different ways to determine removal of vars
corSelect(data = datGAM, sp.cols = "dev", var.cols = names(datGAM)[-(1:2)],
          coeff = FALSE) # based on p-value cutoff (0.05)
corSelect(data = datGAM, sp.cols = "dev", var.cols = names(datGAM)[-(1:2)],
          coeff = TRUE) # based on correlation coefficient magnitude (0.8)

# Code to create candidate model structures with low-correlation covariates
candModCovars <- list()
for(ii in 1:500){
  # get names of covariates in 'datGAM'
  allCovarNames <- names(datGAM)
  allCovarNames <- allCovarNames[-which(allCovarNames %in% c("Yr", "dev", "devLag1"))]
  # take sub-sample of covar names
  propNames <- sample(allCovarNames, size = sample(2:5, 1))
  # subDat <- datGAM %>% dplyr::select(all_of(propNames))
  # # find correlation matrix of subset
  # corrMat <- cor(subDat, use = "pairwise.complete.obs")
  # # find and remove highly correlated covars
  # rmNames <- findCorrelation(x = corrMat, cutoff = 0.8)
  # # record remaining combo of low-correlation covars
  # candModCovars[[ii]] <- sort(propNames[-rmNames])

  # could also base off of p-value threshold
  subDat <- datGAM %>% dplyr::select(dev, all_of(propNames))
  candSel <- corSelect(data = subDat, sp.cols = "dev", var.cols = names(subDat)[-1],
                        coeff = FALSE) # based on p-value cutoff (0.05)
  candModCovars[[ii]] <- sort(candSel$selected.vars)
}


# test1 <- unique(candModCovars)
# test2 <- unique(candModCovars)
# test3 <- unique(candModCovars)

list2df_dt <- function(x) {
  tmp <- lapply(x, as.data.frame, stringsAsFactors = FALSE)
  tmp <- data.table::rbindlist(tmp, idcol = "name")
  colnames(tmp)[2] <-  "item"
  tmp
}
test1long <- list2df_dt(test1)
test1long <- as.data.frame(test1long) %>% mutate(inMod = 1) %>% pivot_wider(values_from = inMod, names_from = item)
test2long <-list2df_dt(test2)
test2long <- as.data.frame(test2long) %>% mutate(inMod = 1) %>% pivot_wider(values_from = inMod, names_from = item)
test3long <-list2df_dt(test3)
test3long <- as.data.frame(test3long) %>% mutate(inMod = 1) %>% pivot_wider(values_from = inMod, names_from = item)

candMods <- bind_rows(test1long, test2long, test3long) %>% dplyr::select(-name)
unique(candMods) %>% dim() 
# search over 100 iterations gives 186 models
# search over 500 iterations gives 885 models

# check for models with just top 10 highest loading vars
candMods %>% dplyr::select(C.pacificus, CUTI_39N, OC_LUSI_39N,
                                          NCOPsummerlag1, OC_STI_39N, BEUTI_39N,
                                          sardNurseHab, HCI_30N355N, NCOPspring,
                                          springSST) %>%
              distinct() %>% View()

# template model order
# gamAll <- gam(dev ~ s(C.pacificus, k = 4) + s(CUTI_39N, k = 4) + s(OC_LUSI_39N, k = 4) + s(NCOPsummerlag1, k = 4) + s(OC_STI_39N, k = 4) + s(BEUTI_39N, k = 4) + s(sardNurseHab, k = 4) + s(HCI_30N355N, k = 4) + s(NCOPspring, k = 4) + s(springSST, k = 4) - 1,
#               data = datGAM, method = "REML", select = TRUE)

# additional models to consider from search above
gam1 <- gam(dev ~ s(C.pacificus, k = 4) + s(CUTI_39N, k = 4) + s(OC_LUSI_39N, k = 4) - 1,
              data = datGAM, method = "REML", select = TRUE)
gam2 <- gam(dev ~ s(C.pacificus, k = 4) + s(OC_STI_39N, k = 4) - 1,
              data = datGAM, method = "REML", select = TRUE)
gam3 <- gam(dev ~ s(OC_LUSI_39N, k = 4) + s(springSST, k = 4) - 1,
              data = datGAM, method = "REML", select = TRUE)
gam4 <- gam(dev ~ s(OC_LUSI_39N, k = 4) + s(sardNurseHab, k = 4) + s(NCOPspring, k = 4) - 1,
              data = datGAM, method = "REML", select = TRUE)
gam5 <- gam(dev ~ s(NCOPsummerlag1, k = 4) + s(springSST, k = 4) - 1,
              data = datGAM, method = "REML", select = TRUE)
gam6 <- gam(dev ~ s(OC_STI_39N, k = 4) + s(springSST, k = 4) - 1,
              data = datGAM, method = "REML", select = TRUE)
gam7 <- gam(dev ~ s(CUTI_39N, k = 4) + s(HCI_30N355N, k = 4) - 1,
              data = datGAM, method = "REML", select = TRUE)
gam8 <- gam(dev ~ s(NCOPsummerlag1, k = 4) + s(sardNurseHab, k = 4) - 1,
              data = datGAM, method = "REML", select = TRUE)
gam9 <- gam(dev ~ s(CUTI_39N, k = 4) + s(sardNurseHab, k = 4) + s(HCI_30N355N, k = 4) + s(NCOPspring, k = 4) - 1,
              data = datGAM, method = "REML", select = TRUE)
gam10 <- gam(dev ~ s(OC_LUSI_39N, k = 4) + s(NCOPspring, k = 4) - 1,
              data = datGAM, method = "REML", select = TRUE)
gam11 <- gam(dev ~ s(sardNurseHab, k = 4) + s(springSST, k = 4) - 1,
              data = datGAM, method = "REML", select = TRUE)
gam12 <- gam(dev ~ s(OC_LUSI_39N, k = 4) + s(NCOPsummerlag1, k = 4) - 1,
              data = datGAM, method = "REML", select = TRUE)
gam13 <- gam(dev ~ s(OC_LUSI_39N, k = 4) + s(BEUTI_39N, k = 4) - 1,
              data = datGAM, method = "REML", select = TRUE)
gam14 <- gam(dev ~ s(C.pacificus, k = 4) + s(NCOPsummerlag1, k = 4) - 1,
              data = datGAM, method = "REML", select = TRUE)
gam15 <- gam(dev ~ s(C.pacificus, k = 4) + s(OC_LUSI_39N, k = 4) - 1,
              data = datGAM, method = "REML", select = TRUE)
gam16 <- gam(dev ~ s(OC_LUSI_39N, k = 4) + s(HCI_30N355N, k = 4) - 1,
              data = datGAM, method = "REML", select = TRUE)
gam17 <- gam(dev ~ s(OC_STI_39N, k = 4) + s(BEUTI_39N, k = 4) - 1,
              data = datGAM, method = "REML", select = TRUE)
gam18 <- gam(dev ~ s(BEUTI_39N, k = 4) + s(sardNurseHab, k = 4) - 1,
              data = datGAM, method = "REML", select = TRUE)
gam19 <- gam(dev ~ s(OC_STI_39N, k = 4) + s(sardNurseHab, k = 4) - 1,
              data = datGAM, method = "REML", select = TRUE)
gam20 <- gam(dev ~ s(BEUTI_39N, k = 4) + s(NCOPspring, k = 4) - 1,
              data = datGAM, method = "REML", select = TRUE)
gam21 <- gam(dev ~ s(sardNurseHab, k = 4) + s(NCOPspring, k = 4) - 1,
              data = datGAM, method = "REML", select = TRUE)
gam22 <- gam(dev ~ s(C.pacificus, k = 4) + s(CUTI_39N, k = 4) - 1,
              data = datGAM, method = "REML", select = TRUE)

# Fit GAMs and select covars ----------------------------------------------

# see if previous rec dev helps estimation
# Compare performance against lagged rec devs GAM
gamPersist <- gam(dev ~ devLag1, data = datGAM)
summary(gamPersist)

# try models including only non-correlated variables
# starting with vars with highest loadings
fitSST <- gam(dev ~ s(C.pacificus, k = 4) + s(OC_STI_39N, k = 4) + s(springSST, k = 4) - 1,
              data = datGAM, method = "REML", select = TRUE)
summary(fitSST) # --> just linear model on STI
fitNCOP1 <- gam(dev ~ s(C.pacificus, k = 4) + s(HCI_30N355N, k = 4) + s(NCOPspring, k = 4) - 1,
              data = datGAM, method = "REML", select = TRUE)
summary(fitNCOP1) # --> only smoothed model on HCI
fitNCOP2 <- gam(dev ~ s(CUTI_39N, k = 4) + s(NCOPspring, k = 4) - 1,
              data = datGAM, method = "REML", select = TRUE)
summary(fitNCOP2) # --> only linear model on NCOPspring
fitNCOP3 <- gam(dev ~ s(HCI_30N355N, k = 4) + s(NCOPspring, k = 4) - 1,
              data = datGAM, method = "REML", select = TRUE)
summary(fitNCOP3) # both ok, but linear on NCOPspring
fitHCI1 <- gam(dev ~ s(C.pacificus, k = 4) + s(OC_LUSI_39N, k = 4) + s(BEUTI_39N, k = 4) + s(HCI_30N355N, k = 4) - 1,
              data = datGAM, method = "REML", select = TRUE)
summary(fitHCI1) # no effect of C.pac or LUSI
fitHCI2 <- gam(dev ~ s(C.pacificus, k = 4) + s(BEUTI_39N, k = 4) + s(HCI_30N355N, k = 4) - 1,
              data = datGAM, method = "REML", select = TRUE)
summary(fitHCI2) # no effect of C.pac
fitHCI4 <- gam(dev ~ s(C.pacificus, k = 4) + s(OC_STI_39N, k = 4) + s(BEUTI_39N, k = 4) + s(HCI_30N355N, k = 4) - 1,
              data = datGAM, method = "REML", select = TRUE)
summary(fitHCI4) # no effect of C.pac, only linear STI
fitHCI5 <- gam(dev ~ C.pacificus + OC_STI_39N + s(BEUTI_39N, k = 4) + s(HCI_30N355N, k = 4) - 1,
              data = datGAM, method = "REML")
summary(fitHCI5) # no effect of C.pac, only linear STI
fitNurse1 <- gam(dev ~ s(C.pacificus, k = 4) + s(BEUTI_39N, k = 4) + s(sardNurseHab, k = 4) - 1,
              data = datGAM, method = "REML", select = TRUE)
summary(fitNurse1) # --> only effect on BEUTI
fitNurse2 <- gam(dev ~ s(CUTI_39N, k = 4) + s(sardNurseHab, k = 4) - 1,
              data = datGAM, method = "REML", select = TRUE)
summary(fitNurse2) # no sig effects
fitBEUTI1 <- gam(dev ~ s(C.pacificus, k = 4) + s(OC_LUSI_39N, k = 4) + s(BEUTI_39N, k = 4) + s(sardNurseHab, k = 4) - 1,
              data = datGAM, method = "REML", select = TRUE)
summary(fitBEUTI1) # --> only effect on BEUTI
fitBEUTI6 <- gam(dev ~ s(C.pacificus, k = 4) + s(NCOPsummerlag1, k = 4) + s(BEUTI_39N, k = 4) - 1,
              data = datGAM, method = "REML", select = TRUE)
summary(fitBEUTI6) # --> only effect on BEUTI
fitSTI2 <- gam(dev ~ s(C.pacificus, k = 4) + s(OC_STI_39N, k = 4) + s(HCI_30N355N, k = 4) - 1,
              data = datGAM, method = "REML", select = TRUE)
summary(fitSTI2) # --> only effect on HCI
fitSTI4 <- gam(dev ~ s(CUTI_39N, k = 4) + s(OC_STI_39N, k = 4) - 1,
              data = datGAM, method = "REML", select = TRUE)
summary(fitSTI4) # --> only linear effect on STI
fitSTI5 <- gam(dev ~ s(OC_STI_39N, k = 4) + s(HCI_30N355N, k = 4) - 1,
              data = datGAM, method = "REML", select = TRUE)
summary(fitSTI5) # effects of both, but STI only linear
fitSTI6 <- gam(dev ~ s(OC_STI_39N, k = 4) + s(BEUTI_39N, k = 4) + s(HCI_30N355N, k = 4) - 1,
              data = datGAM, method = "REML", select = TRUE)
summary(fitSTI6) # hardly an effect of STI, no effect of BEUTI
fitSTI6.1 <- gam(dev ~ OC_STI_39N + s(BEUTI_39N, k = 4) + s(HCI_30N355N, k = 4) - 1,
              data = datGAM, method = "REML", select = TRUE)
summary(fitSTI6.1) # removes effect of BEUTI
fitSTI7 <- gam(dev ~ s(OC_STI_39N, k = 4) + s(HCI_30N355N, k = 4) + s(springSST, k = 4) - 1,
              data = datGAM, method = "REML", select = TRUE)
summary(fitSTI7) # no effect of SST, only linear STI
fitNCOPLag2 <- gam(dev ~ s(CUTI_39N, k = 4) + s(NCOPsummerlag1, k = 4) - 1,
              data = datGAM, method = "REML", select = TRUE)
summary(fitNCOPLag2) # no effect of CUTI
fitNCOPLag3 <- gam(dev ~ s(NCOPsummerlag1, k = 4) + s(BEUTI_39N, k = 4) - 1,
              data = datGAM, method = "REML", select = TRUE)
summary(fitNCOPLag3) # no effect of BEUTI
fitLUSI2 <- gam(dev ~ s(CUTI_39N, k = 4) + s(OC_LUSI_39N, k = 4) - 1,
              data = datGAM, method = "REML", select = TRUE)
summary(fitLUSI2) # no sig effects
fitLUSI3 <- gam(dev ~ s(OC_LUSI_39N, k = 4) + s(BEUTI_39N, k = 4) + s(HCI_30N355N, k = 4) - 1,
              data = datGAM, method = "REML", select = TRUE)
summary(fitLUSI3) # no effect of LUSI or BEUTI
fitCpac3 <- gam(dev ~ s(C.pacificus, k = 4) + s(HCI_30N355N, k = 4) - 1,
              data = datGAM, method = "REML", select = TRUE)
summary(fitCpac3) # only effect on HCI
fitCpac4 <- gam(dev ~ s(C.pacificus, k = 4) + s(springSST, k = 4) - 1,
              data = datGAM, method = "REML", select = TRUE)
summary(fitCpac4) # no sig effects
fitCpac6 <- gam(dev ~ s(C.pacificus, k = 4) + s(NCOPspring, k = 4) - 1,
              data = datGAM, method = "REML", select = TRUE)
summary(fitCpac6) # no sig effects
fitCpac8 <- gam(dev ~ s(C.pacificus, k = 4) + s(sardNurseHab, k = 4) - 1,
              data = datGAM, method = "REML", select = TRUE)
summary(fitCpac8) # no sig effects
fitCpac12 <- gam(dev ~ s(C.pacificus, k = 4) + s(BEUTI_39N, k = 4) - 1,
              data = datGAM, method = "REML", select = TRUE)
summary(fitCpac12) # only effect on BEUTI

# candidate models
fitBEUTI.HCI <- gam(dev ~ s(BEUTI_39N, k = 4) + s(HCI_30N355N, k = 4) - 1,
                data = datGAM, method = "REML")
summary(fitBEUTI.HCI) # 63% dev explained
concurvity(fitBEUTI.HCI)
plot(fitBEUTI.HCI, pages = 1)
gam.check(fitBEUTI.HCI)

fitBEUTI.NCOPlag1 <- gam(dev ~ s(BEUTI_39N, k = 4) + s(NCOPsummerlag1, k = 4) - 1,
                    data = datGAM, method = "REML")
summary(fitBEUTI.NCOPlag1) # much lower dev explained (16%)
fitNCOPlag1 <- gam(dev ~ s(NCOPsummerlag1, k = 4) + OC_STI_39N + s(BEUTI_39N, k = 4) + s(HCI_30N355N, k = 4) - 1,
                    data = datGAM, method = "REML")
concurvity(fitNCOPlag1)
summary(fitNCOPlag1)

fitSTI.BEUTI.HCI1 <- gam(dev ~ s(OC_STI_39N, k = 4) + s(BEUTI_39N, k = 4) + s(HCI_30N355N, k = 4) - 1,
                 data = datGAM, method = "REML")
summary(fitSTI.BEUTI.HCI1) # 67% dev expl
fitSTI.BEUTI.HCI2 <- gam(dev ~ OC_STI_39N + s(BEUTI_39N, k = 4) + s(HCI_30N355N, k = 4) - 1,
                 data = datGAM, method = "REML")
summary(fitSTI.BEUTI.HCI2) # 67% dev expl
concurvity(fitSTI.BEUTI.HCI2)
fitSTI.HCI <- gam(dev ~ OC_STI_39N + s(HCI_30N355N, k = 4) - 1,
                        data = datGAM, method = "REML")
summary(fitSTI.HCI) # 35% dev expl

fitSTI.CUTI.BEUTI.HCI1 <- gam(dev ~ s(OC_STI_39N, k = 4) + s(CUTI_33N, k = 4) + s(BEUTI_39N, k = 4) + s(HCI_30N355N, k = 4) - 1,
                        data = datGAM, method = "REML")
summary(fitSTI.CUTI.BEUTI.HCI1) # 89% dev expl
plot(fitSTI.CUTI.BEUTI.HCI1, pages = 1)
concurvity(fitSTI.CUTI.BEUTI.HCI1)
gam.check(fitSTI.CUTI.BEUTI.HCI1)
fitSTI.CUTI.BEUTI.HCI2 <- gam(dev ~ OC_STI_39N + s(CUTI_33N, k = 4) + s(BEUTI_39N, k = 4) + s(HCI_30N355N, k = 4) - 1,
                        data = datGAM, method = "REML")
summary(fitSTI.CUTI.BEUTI.HCI2) # 81.5% dev expl
fitSTI.CUTI.BEUTI.HCI3 <- gam(dev ~ s(OC_STI_39N, k = 4) + CUTI_33N + s(BEUTI_39N, k = 4) + s(HCI_30N355N, k = 4) - 1,
                        data = datGAM, method = "REML")
summary(fitSTI.CUTI.BEUTI.HCI3) # 86% dev expl
fitSTI.CUTI.BEUTI.HCI4 <- gam(dev ~ OC_STI_39N + CUTI_33N + s(BEUTI_39N, k = 4) + s(HCI_30N355N, k = 4) - 1,
                        data = datGAM, method = "REML")
summary(fitSTI.CUTI.BEUTI.HCI4) # 78% dev expl
fitSTI.CUTI.HCI <- gam(dev ~ s(OC_STI_39N, k = 4) + s(CUTI_33N, k = 4) + s(HCI_30N355N, k = 4) - 1,
                        data = datGAM, method = "REML")
summary(fitSTI.CUTI.HCI) # 37.5% dev expl

# try with more vars (lower loading ests) ----------------------------------

# try adding single vars to best model to check for improvement
fitSCOP <- gam(dev ~ s(SCOPspring, k = 4) + s(OC_STI_39N, k = 4) + s(BEUTI_39N, k = 4) + s(HCI_30N355N, k = 4) - 1,
                        data = datGAM, method = "REML")
summary(fitSCOP) # linear SCOP, 78.7% dev expl
concurvity(fitSCOP)
fitSCOPlag <- gam(dev ~ s(SCOPsummerlag1, k = 4) + s(BEUTI_39N, k = 4) + s(HCI_30N355N, k = 4) - 1,
                        data = datGAM, method = "REML")
summary(fitSCOPlag) # linear SCOP, 63% dev expl
concurvity(fitSCOPlag) # removed STI
fitBEUTI33 <- gam(dev ~ s(BEUTI_33N, k = 4) + s(OC_STI_39N, k = 4) + s(BEUTI_39N, k = 4) + s(HCI_30N355N, k = 4) - 1,
                        data = datGAM, method = "REML")
summary(fitBEUTI33) # linear BEUTI and STI, 52.5% dev expl
concurvity(fitBEUTI33)
fitLUSI33 <- gam(dev ~ s(OC_LUSI_33N, k = 4) + s(OC_STI_39N, k = 4) + s(BEUTI_39N, k = 4) - 1,
                        data = datGAM, method = "REML")
summary(fitLUSI33) # mostly just STI kept, 9% dev expl
concurvity(fitLUSI33) # removed HCI
fitLUSI36 <- gam(dev ~ s(OC_LUSI_36N, k = 4) + s(OC_STI_39N, k = 4) + s(BEUTI_39N, k = 4) + s(HCI_30N355N, k = 4) - 1,
                        data = datGAM, method = "REML")
summary(fitLUSI36) # linear LUSI and STI, 79% dev expl
concurvity(fitLUSI36)
fitSTI33 <- gam(dev ~ s(OC_STI_33N, k = 4) + s(BEUTI_39N, k = 4) - 1,
                        data = datGAM, method = "REML")
summary(fitSTI33) # both linear, 3% dev expl
concurvity(fitSTI33) # remove STI39 and HCI
fitSTI36 <- gam(dev ~ s(OC_STI_36N, k = 4) + s(OC_STI_39N, k = 4) + s(BEUTI_39N, k = 4) + s(HCI_30N355N, k = 4) - 1,
                        data = datGAM, method = "REML")
summary(fitSTI36) # linear STI, 68% dev expl
concurvity(fitSTI36)
fitMyct <- gam(dev ~ s(RREAS_Myctophids, k = 4) + s(BEUTI_39N, k = 4) + s(HCI_30N355N, k = 4) - 1,
                        data = datGAM, method = "REML")
summary(fitMyct) # 65% dev expl
concurvity(fitMyct) # remove STI
fitYOYsard <- gam(dev ~ s(RREAS_YOYsardine, k = 4) + s(OC_STI_39N, k = 4) + s(BEUTI_39N, k = 4) - 1,
                        data = datGAM, method = "REML")
summary(fitYOYsard) # YOY linear, 52.5% dev expl
concurvity(fitYOYsard)# remove HCI
fitSardLarv <- gam(dev ~ s(sardLarv, k = 4) + s(OC_STI_39N, k = 4) + s(HCI_30N355N, k = 4) - 1,
                        data = datGAM, method = "REML")
summary(fitSardLarv) # 73% dev expl
concurvity(fitSardLarv)# remove BEUTI
fitMesoLarv <- gam(dev ~ s(mesopelLarv, k = 4) + s(OC_STI_39N, k = 4) + s(BEUTI_39N, k = 4) - 1,
                        data = datGAM, method = "REML")
summary(fitMesoLarv) # all linear, 12.5% dev expl
concurvity(fitMesoLarv)# remove HCI
fitAge1WAA <- gam(dev ~ s(age1SprSardmeanWAA, k = 4) + s(OC_STI_39N, k = 4) + s(BEUTI_39N, k = 4) - 1,
                        data = datGAM, method = "REML")
summary(fitAge1WAA) # 56% dev expl
concurvity(fitAge1WAA) # remove HCI
fitSSBwt1 <- gam(dev ~ s(meanSSBwt, k = 4) + OC_STI_39N + s(BEUTI_39N, k = 4) + s(HCI_30N355N, k = 4) - 1,
                        data = datGAM, method = "REML")
fitSSBwt2 <- gam(dev ~ s(meanSSBwt, k = 4) + s(BEUTI_39N, k = 4) + s(HCI_30N355N, k = 4) - 1,
                        data = datGAM, method = "REML")
fitSSBwt3 <- gam(dev ~ s(meanSSBwt, k = 4) + s(OC_STI_39N, k = 4) + s(BEUTI_39N, k = 4) - 1,
                 data = datGAM, method = "REML")
summary(fitSSBwt1) # 92% dev expl
concurvity(fitSSBwt1)# remove STI
fitSpawnHab <- gam(dev ~ s(sardSpawnHab, k = 4) + s(OC_STI_39N, k = 4) + s(BEUTI_39N, k = 4) - 1,
                        data = datGAM, method = "REML")
summary(fitSpawnHab) # linear SpawnHab and STI
concurvity(fitSpawnHab)# remove HCI
fitSpawnDur <- gam(dev ~ s(daysAbove5pct, k = 4) + s(OC_STI_39N, k = 4) + s(BEUTI_39N, k = 4) - 1,
                        data = datGAM, method = "REML")
summary(fitSpawnDur) # linear daysAbove5pct and STI, 46.5% dev expl
concurvity(fitSpawnDur)# remove HCI
fitSuSST <- gam(dev ~ s(summerSST, k = 4) + s(OC_STI_39N, k = 4) + s(BEUTI_39N, k = 4) - 1,
                        data = datGAM, method = "REML")
summary(fitSuSST) # linear STI and BEUTI, 20% dev expl
concurvity(fitSuSST)# remove HCI
fitSprNear <- gam(dev ~ s(avgNearTransspring, k = 4) + s(OC_STI_39N, k = 4) + s(BEUTI_39N, k = 4) + s(HCI_30N355N, k = 4) - 1,
                        data = datGAM, method = "REML")
summary(fitSprNear) # linear NearTransp, STI and BEUTI, 54% dev expl
concurvity(fitSprNear)
fitSuNear <- gam(dev ~ s(avgNearTranssummer, k = 4) + s(OC_STI_39N, k = 4) + s(BEUTI_39N, k = 4) + s(HCI_30N355N, k = 4) - 1,
                        data = datGAM, method = "REML")
summary(fitSuNear) # linear NearTransp and STI, 78% dev expl
concurvity(fitSuNear)
fitSprOff <- gam(dev ~ s(avgOffTransspring, k = 4) + s(OC_STI_39N, k = 4) + s(HCI_30N355N, k = 4) - 1,
                        data = datGAM, method = "REML")
summary(fitSprOff) # linear STI, 55% dev expl
concurvity(fitSprOff)# remove BEUTI
fitSuOff <- gam(dev ~ s(avgOffTranssummer, k = 4) + s(OC_STI_39N, k = 4) + s(HCI_30N355N, k = 4) - 1,
                        data = datGAM, method = "REML")
summary(fitSuOff) # linear OffTransp and STI, 36% dev expl
concurvity(fitSuOff)# remove BEUTI
fitSLYOY <- gam(dev ~ s(yoySardSL, k = 4) + s(OC_STI_39N, k = 4) + s(BEUTI_39N, k = 4) + s(HCI_30N355N, k = 4) - 1,
                        data = datGAM, method = "REML")
summary(fitSLYOY) # linear STI and BEUTI, 60.5% dev expl
concurvity(fitSLYOY)


# vars with significant correlation with rec devs
# gamCorrDev <- gam(dev ~ s(sardSpawnHab, k = 4) + s(avgNearTransspring, k = 4) + s(avgOffTransspring, k = 4) + s(sardLarv, k = 4) - 1,
gamCorrDev <- gam(dev ~ sardSpawnHab + avgNearTransspring + avgOffTransspring + sardLarv - 1,
                       data = datGAM, method = "REML")
summary(gamCorrDev) # only linear effects
gamNearTransp <- gam(dev ~ avgNearTransspring - 1,
                       data = datGAM, method = "REML")
summary(gamNearTransp)

gamUncorrs <- gam(dev ~ s(CUTI_33N, k = 4) + s(BEUTI_33N, k = 4) + s(yoySardSL, k = 4) + s(OC_LUSI_33N, k = 4) - 1,
                     data = datGAM, method = "REML", select = TRUE)
summary(gamUncorrs)
gamBEUTI.CUTI <- gam(dev ~ s(CUTI_33N, k = 4) + s(BEUTI_33N, k = 4) - 1,
                     data = datGAM, method = "REML", select = TRUE)
summary(gamBEUTI.CUTI)
gamYoySL <- gam(dev ~ s(yoySardSL, k = 4) - 1,
                  data = datGAM, method = "REML", select = TRUE)
summary(gamYoySL)
# ---------------------------
# see if previous rec dev helps estimation
gamPersist <- gam(dev ~ devLag1, data = datGAM, method = "REML")
summary(gamPersist)
# Compare this to other results

# single variable models
gamSTI <- gam(dev ~ + s(OC_STI_39N, k = 4) - 1,
              data = datGAM, method = "REML")
summary(gamSTI) # 8.5% dev expl
gamBEUTI <- gam(dev ~ s(BEUTI_39N, k = 4) - 1,
              data = datGAM, method = "REML")
summary(gamBEUTI) # 22.7% dev expl
gamHCI <- gam(dev ~ s(HCI_30N355N, k = 4) - 1,
              data = datGAM, method = "REML")
summary(gamHCI) # 27.5% dev expl
gamNCOP <- gam(dev ~ s(NCOPspring, k = 4) - 1,
              data = datGAM, method = "REML")
summary(gamNCOP) # 7.15% dev expl
gamNCOPlag <- gam(dev ~ s(NCOPsummerlag1, k = 4) - 1,
              data = datGAM, method = "REML")
gamCpac <- gam(dev ~ s(C.pacificus, k = 4) - 1,
              data = datGAM, method = "REML")
gamCUTI <- gam(dev ~ s(CUTI_39N, k = 4) - 1,
              data = datGAM, method = "REML")
gamLUSI <- gam(dev ~ s(OC_LUSI_39N, k = 4) - 1,
              data = datGAM, method = "REML")
gamNurseHab <- gam(dev ~ s(sardNurseHab, k = 4) - 1,
              data = datGAM, method = "REML")
gamSST <- gam(dev ~ s(springSST, k = 4) - 1,
              data = datGAM, method = "REML")

# Summarize results -------------------------------------------------------


modelObjs <- grep(pattern = "fit", x = ls(), value = TRUE)
modelObjs <- c(modelObjs, grep(pattern = "gam", x = ls(), value = TRUE))
recAIC <- numeric()
recDevExp <- numeric()
recRsq <- numeric()
for(j in 1:length(modelObjs)){
  recAIC[j] <- do.call(what = AIC, args = list(get(modelObjs[[j]])))
  tmp <- do.call(what = summary, list(get(modelObjs[[j]])))
  recDevExp[j] <- tmp$dev.expl
  recRsq[j] <- tmp$r.sq
}
smryTbl <- data.frame(modelName = modelObjs,
                      AIC = recAIC,
                      devExp = recDevExp, 
                      Rsq = recRsq)

smryTbl <- smryTbl %>% mutate(deltaAIC = AIC - min(AIC)) %>% arrange(deltaAIC)
smryTbl %>% arrange(desc(devExp))

# gamARCpac <- gam(dev ~ s(C.pacificus) + devLag1 - 1,
#                  data = datGAM, select = TRUE)
# summary(gamARCpac)
# gamARCUTI <- gam(dev ~ s(CUTI_39N) + devLag1 - 1,
#               data = datGAM, select = TRUE)
# summary(gamARCUTI) # better GCV
# gamARLUSI <- gam(dev ~ s(OC_LUSI_39N) + devLag1 - 1,
#               data = datGAM, select = TRUE)
# summary(gamARLUSI) # better GCV
# gamARNCOPLag1 <- gam(dev ~ s(NCOPsummerlag1) + devLag1 - 1,
#               data = datGAM, select = TRUE)
# summary(gamARNCOPLag1) # better GCV
# gamARSTI <- gam(dev ~ s(OC_STI_39N) + devLag1 - 1,
#               data = datGAM, select = TRUE)
# summary(gamARSTI) # better GCV
# gamARBEUTI <- gam(dev ~ s(BEUTI_39N) + devLag1 - 1,
#               data = datGAM, select = TRUE)
# summary(gamARBEUTI) # better GCV
# gamARnurseHab <- gam(dev ~ s(sardNurseHab) + devLag1 - 1,
#                      data = datGAM, select = TRUE)
# summary(gamARnurseHab) # better GCV
# gamARHCI <- gam(dev ~ s(HCI_30N355N) + devLag1 - 1,
#               data = datGAM, select = TRUE)
# summary(gamARHCI) # better GCV
# gamARNCOP <- gam(dev ~  s(NCOPspring) + devLag1 - 1,
#               data = datGAM, select = TRUE)
# summary(gamARNCOP) # better GCV
# gamARSST <- gam(dev ~ s(springSST) + devLag1 - 1,
#               data = datGAM, select = TRUE)
# summary(gamARSST) # better GCV
# 
# gamARNCOPHCI <- gam(dev ~  s(HCI_30N355N) + s(NCOPspring) + devLag1 - 1,
#                  data = datGAM)
# summary(gamARNCOPHCI)
