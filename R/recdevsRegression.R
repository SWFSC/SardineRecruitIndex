# fit GAM to DFA input datasets with high loadings for sardine recruitment time series
# Created: 5/27/2025, Robert Wildermuth

library(MARSS)
library(r4ss)
library(tidyverse)
library(mgcv)

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
candMods <- read_csv("../SardineRecruitIndex/out/candidateGAMmodels.csv")

#!!RW: for now omit rows with NAs
candMods <- na.omit(candMods)

datGAM <- datDFA %>% select(year, sardRec, #all_of(names(candMods)),
                            NCOPspring, NCOPsummerlag1, BEUTI_33N, BEUTI_39N, CUTI_33N,
                            CUTI_39N, OC_LUSI_33N, OC_STI_33N, RREAS_Myctophids, 
                            RREAS_YOYsardine, age1SprSardmeanWAA, meanSSBwt, 
                            ZM_NorCal, ZM_SoCal, sardNurseHab, anchBioSmrySeas1, 
                            summerSST, avgOffTransspring, avgOffTranssummer, 
                            sprRelOffTrans, sumRelOffTrans, meanK,
                            OC_STI_39N, sardSpawnHab, springSST) %>% 
            filter(year %in% 1985:2021)

datGAM <- mngt2024recdevs %>% dplyr::select(Yr, dev) %>% filter(Yr <= 2021) %>%
              left_join(y = datGAM, by = c("Yr"="year"))

datGAM <-cbind(datGAM[,1:3],
               apply(X = datGAM[, -(1:3)], MARGIN = 2, FUN = zscore))
# add lagged rec dev for autocorrelation structure
datGAM <- datGAM %>% mutate(devLag1 = c(1, datGAM$dev[1:(nrow(datGAM)-1)]))

# Fit GAMs and select covars ----------------------------------------------

indNames <- names(candMods)

# template model order
# gamAll <- gam(dev ~ s(CUTI_39N, k = 4) + s(NCOPspring, k = 4) + s(OC_LUSI_39N, k = 4) + s(ZM_SoCal, k = 4) + s(summerSST, k = 4) + s(ZM_NorCal, k = 4) + s(NCOPsummerlag1, k = 4) + s(springSST, k = 4) + s(sardNurseHab, k = 4) + s(BEUTI_39N, k = 4) - 1,
#               data = datGAM, method = "REML", select = TRUE)


# see if previous rec dev helps estimation
# Compare performance against lagged rec devs GAM
gamPersist <- gam(dev ~ devLag1 - 1, data = datGAM, method = "REML", select = TRUE)
summary(gamPersist)

# fit each model specified in 'candMods'
# use code borrowed from Megan Feddern: https://github.com/mfeddern/YellowtailNorth_EnvIndex/blob/main/StockAssessmentIndex/OceanographicIndexDevelopment2025.R
models <- list() #models list to fill with loop
smryTbl <- tibble(modelN = 0,
                  AIC = 0,
                  devExp = 0, 
                  Rsq = 0)[0,]
for (i in 1:nrow(candMods)) { #loop over each covariate combination i
  # k represent the number of parameters / knots estimating function at, should be small
  smooth_terms <- paste("s(", indNames[as.logical(candMods[i,])], ", k = 4)", collapse = " + ") #generatind smooth term for each combination
  formula_str <- paste("dev ~ ", smooth_terms, "- 1") #generating the formula string for each smooth term

  #fitting full model 
  gam_model <- gam(as.formula(formula_str),
                   data = datGAM, method = "REML", select = TRUE)
  gamSry <- summary(gam_model)
  
  # Store results
  models[[i]] <- gam_model
  smryTbl <- bind_rows(smryTbl,
                       data.frame(modelN = i,
                                  AIC= AIC(gam_model),
                                  devExp = gamSry$dev.expl,
                                  Rsq = gamSry$r.sq))
}

# candidate models
fitBEUTI <- gam(dev ~ s(NCOPsummerlag1 , k = 4) + s(BEUTI_33N  , k = 4) - 1,
                data = datGAM, method = "REML", select = TRUE)
summary(fitBEUTI) # 
concurvity(fitBEUTI)
plot(fitBEUTI, pages = 1)
gam.check(fitBEUTI)

fitNurse <- gam(dev ~ s(sardNurseHab, k = 4) - 1,
                data = datGAM, method = "REML", select = TRUE)
summary(fitNurse)
fitBEUTI.Nurse <- gam(dev ~ s(BEUTI_39N, k = 4) + s(sardNurseHab, k = 4) - 1,
                      data = datGAM, method = "REML", select = TRUE)
summary(fitBEUTI.Nurse)


# try with more vars (lower loading ests) ----------------------------------

# try adding single vars to best model to check for improvement


# vars with significant correlation with rec devs

# Summarize results -------------------------------------------------------

# modelObjs <- grep(pattern = "fit", x = ls(), value = TRUE)
# modelObjs <- c(modelObjs, grep(pattern = "gam", x = ls(), value = TRUE))
# recAIC <- numeric()
# recDevExp <- numeric()
# recRsq <- numeric()
# for(j in 1:length(modelObjs)){
#   recAIC[j] <- do.call(what = AIC, args = list(get(modelObjs[[j]])))
#   tmp <- do.call(what = summary, list(get(modelObjs[[j]])))
#   recDevExp[j] <- tmp$dev.expl
#   recRsq[j] <- tmp$r.sq
# }
# smryTbl <- data.frame(modelName = modelObjs,
#                       AIC = recAIC,
#                       devExp = recDevExp, 
#                       Rsq = recRsq)

smryTbl <- smryTbl %>% mutate(deltaAIC = AIC - min(AIC)) %>% arrange(deltaAIC)
smryTbl %>% arrange(desc(devExp))
# 6 models with best scores: 10, 17, 23, 27, 28, 44
candMods[c(128, 130, 138, 151, 154, 155, 164, 166, 170),]

summary(models[[151]])
concurvity(models[[151]]) # 128, 130, 154, 155, 164, 166, 170 too high of concurvity, 138 borderline
plot(models[[151]], pages = 1)
gam.check(models[[151]])
