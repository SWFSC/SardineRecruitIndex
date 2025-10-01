# Code to evaluate multicollinearity in model variables and select candidate subsets for analyses
# Created: 7/17/2025, Robert Wildermuth

library(tidyverse)
library(MARSS)
library(corrplot)

# read prepped dataset
datDFA <- read_csv("Data/recrDFAdat.csv")

allDat <- datDFA %>% filter(year %in% 1985:2021) %>%
            select(-c(NCOPsummer,
                      SCOPsummer,
                      GCM)) 

datNames <- names(allDat)[-1]

# transpose for MARSS formatting
allDat <- allDat %>% select(-year) %>% t()

datZscore <- zscore(allDat) 
corrMat <- cor(t(datZscore), use = "pairwise.complete.obs")
pTest <- cor.mtest(t(datZscore), alternative = "two.sided", method = "pearson")
corrplot(corrMat, p.mat = pTest$p, sig.level = 0.05, insig = "blank",
         order = 'hclust', hclust.method = "ward.D2", #"centroid", #"single", #
         tl.col = 'black', type = "lower",
         cl.ratio = 0.1, tl.srt = 45, tl.cex = 0.6, #mar = c(0.1, 0.1, 0.1, 0.1), 
         addrect = 6, rect.col = "green", diag = FALSE)


# clusters variables based on correlation structure - members at distant branches are less correlated
corrClusts <- klaR::corclust(x = t(datZscore))
plot(corrClusts)

# try selection based on variance inflation factor
collinear::vif_select(df = t(datZscore), max_vif = 10) # sample size is too small

# ranks variables by cumulative Pearson correlation (low to high)
varsKeep <- collinear::cor_select(df = as.data.frame(t(datZscore)), max_cor = 0.8)
datNames[which(!datNames %in% varsKeep)]

# identifies which variables to remove to reduce multicollinearity
caret::findCorrelation(x = corrMat, cutoff = 0.8, names = TRUE)

# similar, but can use different ways to determine removal of vars
fuzzySim::corSelect(data = t(datZscore),  var.cols = datNames,
                    coeff = FALSE) # based on p-value cutoff (0.05)
hiCorrs <- fuzzySim::corSelect(data = t(datZscore), var.cols = datNames,
                               coeff = TRUE) # based on correlation coefficient magnitude (0.8)
hiCorrs$high.correlations %>% arrange(var1)

# Check for repetitive indicators
corrMat[rownames(corrMat) %in% c("OC_STI_33N", "OC_STI_36N", "OC_STI_39N"),
        colnames(corrMat) %in% c("OC_STI_33N", "OC_STI_36N", "OC_STI_39N")]
# STI at northern and southern extent most similar to STI_36N

corrMat[rownames(corrMat) %in% c("OC_LUSI_33N", "OC_LUSI_36N", "OC_LUSI_39N"),
        colnames(corrMat) %in% c("OC_LUSI_33N", "OC_LUSI_36N", "OC_LUSI_39N")]
# LUSI at northern and southern extent most similar to LUSI_36N

# upwelling indicators
corrMat[rownames(corrMat) %in% c("BEUTI_33N", "BEUTI_39N", "CUTI_33N", "CUTI_39N"),
        colnames(corrMat) %in% c("BEUTI_33N", "BEUTI_39N", "CUTI_33N", "CUTI_39N")]
# BEUTI and CUTI at 39N highly correlated

# temperature indicators
corrMat[rownames(corrMat) %in% c("HCI_30N355N", "springSST", "summerSST"),
        colnames(corrMat) %in% c("HCI_30N355N", "springSST", "summerSST")]
# springSST highly correlated with all
# HCI also correlated with summerSST

# check correlations between zooplankton indicators
corrMat[rownames(corrMat) %in% c("C.pacificus", "ZM_SoCal", "ZL_SoCal"),
        colnames(corrMat) %in% c("C.pacificus", "ZM_SoCal", "ZL_SoCal")]
# C.pacificus and ZM_SoCal not correlated in SoCal Bight

corrMat[rownames(corrMat) %in% c("NCOPspring", "SCOPspring", "ZM_NorCal", "ZL_NorCal"),
        colnames(corrMat) %in% c("NCOPspring", "SCOPspring", "ZM_NorCal", "ZL_NorCal")]
# NCOP and SCOP strongly negatively correlated - can drop SCOP
# both correlated with ZM_NorCal and ZL_NorCal
corrMat[rownames(corrMat) %in% c("NCOPspring", "SCOPspring", "SCOPsummerlag1", "NCOPsummerlag1"),
        colnames(corrMat) %in% c("NCOPspring", "SCOPspring", "SCOPsummerlag1", "NCOPsummerlag1")]
#SCOPsummerlag1 with NCOPspring, or NCOPsummerlag1 with SCOPspring
corrMat[rownames(corrMat) %in% c("ZM_SoCal", "ZM_NorCal"),
        colnames(corrMat) %in% c("ZM_SoCal", "ZM_NorCal")]

# check correlations between advection indicators
corrMat[rownames(corrMat) %in% c("avgSSWIspring", "avgOffTransspring", "avgNearTransspring"),
        colnames(corrMat) %in% c("avgSSWIspring", "avgOffTransspring", "avgNearTransspring")]

corrMat[rownames(corrMat) %in% c("avgSSWIsummer", "avgOffTranssummer", "avgNearTranssummer"),
        colnames(corrMat) %in% c("avgSSWIsummer", "avgOffTranssummer", "avgNearTranssummer")]
# strongest association with SSWI and NearTrans in summer

corrMat[rownames(corrMat) %in% c("avgOffTransspring", "avgNearTransspring", 
                                 "avgOffTranssummer", "avgNearTranssummer",
                                 "sprRelOffTrans", "sumRelOffTrans"),
        colnames(corrMat) %in% c("avgOffTransspring", "avgNearTransspring", 
                                 "avgOffTranssummer", "avgNearTranssummer",
                                 "sprRelOffTrans", "sumRelOffTrans")]
# transport not super correlated with each other
# spring NearTrans highly correlated with sprRelOffTrans

# Check correlations with condition factors
corrMat[rownames(corrMat) %in% c("age1SprSardmeanWAA", "meanSSBwt"),
        colnames(corrMat) %in% c("age1SprSardmeanWAA", "meanSSBwt")]
# not highly correlated - keep both

#### Final Selection of indicators ####

# remove redundant variables
allDat <- datDFA %>% filter(year %in% 1985:2023) %>%
            select(-c(NCOPsummer,
                      SCOPsummer,
                      GCM)) %>% 
            select(-c(avgNearTransspring, avgNearTranssummer,
                      anchBioSmrySeas2,
                      OC_STI_36N, OC_LUSI_36N,
                      PS_NorCal, PS_SoCal, PL_NorCal, PL_SoCal, 
                      ZS_NorCal, ZS_SoCal, ZL_NorCal, ZL_SoCal,
                      meanResid, sdResid, # Could add these as options to sample from also
                      SCOPspring, SCOPsummerlag1,
                      # also remove the ecological (not timely or projectable) indicators
                      avgSSWIspring, avgSSWIsummer,
                      C.pacificus,
                      sardLarv, mesopelLarv,
                      yoySardSL, posCThSk))

# three sets of correlated variables:
# OC_STI_39N    vs    OC_LUSI_39N
# daysAbove5pct    vs   sardSpawnHab
# springSST    vs    HCI_30N355N, sardNurseHab

# number of possible low-mid correlated sets:
2^3

corrMat[rownames(corrMat) %in% c("OC_STI_39N", "OC_LUSI_39N", "daysAbove5pct", "sardSpawnHab", "springSST", "HCI_30N355N", "sardNurseHab"),
        colnames(corrMat) %in% c("OC_STI_39N", "OC_LUSI_39N", "daysAbove5pct", "sardSpawnHab", "springSST", "HCI_30N355N", "sardNurseHab")]

setNames <- names(allDat)[-which(names(allDat) %in% c("OC_STI_39N", "OC_LUSI_39N",
                                                      "daysAbove5pct", "sardSpawnHab",
                                                      "springSST", "HCI_30N355N", "sardNurseHab"))]
#!!RW: For now just work with one combination
setNames <- c(setNames, sample(c("OC_STI_39N", "OC_LUSI_39N"), 1),
              sample(c("daysAbove5pct", "sardSpawnHab"), 1),
              sample(c("springSST", "HCI_30N355N", "sardNurseHab"), 1))
# for now use set with LUSI, spawning habitat, and spring SST
# save(setNames, file = "Data/indicatorSetNames_LUSI39spawnHabsprSST.RData")

# GAM exploration ---------------------------------------------------------

# Top 3 variables with significant trends related to sardRec in DFA fit
# with 3 trends, equal variance, anchovy biomass included
# index          cummLoading
# <chr>                <dbl>
# 1 NCOPspring           0.231
# 2 BEUTI_39N            0.229
# 3 sardRec              0.220
# 4 NCOPsummerlag1       0.100
# Top 10 for significant/strong loadings with sardRec from 3-trend DFA
# index             cummLoading
# <chr>                   <dbl>
# 1 ZM_NorCal              0.921 
# 2 springSST              0.917 
# 3 NCOPspring             0.807 
# 4 sardNurseHab           0.748 
# 5 BEUTI_39N              0.728 
# 6 CUTI_39N               0.650 
# 7 summerSST              0.636 
# 8 NCOPsummerlag1         0.588 
# 9 ZM_SoCal               0.581 
# 10 OC_LUSI_39N            0.574
# All significant strong loadings from model with 1 trend, equal variance, anchovy biomass included
#           est    conf.up    conf.low trend            index dummy0 isSig
# 1   0.3308371  0.6139538  0.04772029     1            meanK      0  TRUE
# 2   0.3314988  0.5494866  0.11351113     1 RREAS_YOYsardine      0  TRUE
# 3   0.3337437  0.6346624  0.03282508     1        meanSSBwt      0  TRUE
# 4  -0.3344017 -0.1187809 -0.55002247     1         CUTI_33N      0  TRUE
# 5   0.3672421  0.5785110  0.15597326     1          sardRec      0  TRUE
# 6  -0.4064455 -0.1889540 -0.62393700     1      OC_LUSI_39N      0  TRUE
# 7  -0.4263450 -0.1713215 -0.68136854     1   NCOPsummerlag1      0  TRUE
# 8   0.4874528  0.7174332  0.25747237     1        summerSST      0  TRUE
# 9  -0.5185126 -0.2752890 -0.76173627     1         CUTI_39N      0  TRUE
# 10 -0.5496007 -0.2789844 -0.82021688     1       NCOPspring      0  TRUE
# 11 -0.5511358 -0.3017253 -0.80054627     1        BEUTI_39N      0  TRUE
# 12  0.5742051  0.8551889  0.29322140     1     sardNurseHab      0  TRUE
# 13 -0.5814889 -0.3314427 -0.83153514     1         ZM_SoCal      0  TRUE
# 14  0.7185261  0.9939870  0.44306519     1        springSST      0  TRUE
# 15 -0.7197106 -0.4413661 -0.99805501     1        ZM_NorCal      0  TRUE

# take top 10 from DFA
datGAM <- datDFA %>% filter(year %in% 1985:2021) %>%
            select(ZM_NorCal, springSST, NCOPspring, sardNurseHab, BEUTI_39N, 
                   CUTI_39N, summerSST, NCOPsummerlag1, ZM_SoCal, OC_LUSI_39N,
                   year, sardRec)

# Code to create candidate model structures with low-correlation covariates
candModCovars <- list()
for(ii in 1:100){
  # get names of covariates in 'datGAM'
  allCovarNames <- names(datGAM)
  allCovarNames <- allCovarNames[-which(allCovarNames %in% c("year", "sardRec"))]
  # take sub-sample of covar names
  propNames <- sample(allCovarNames, size = sample(2:5, 1))
  subDat <- datGAM %>% dplyr::select(all_of(propNames))
  # find correlation matrix of subset
  corrMat <- cor(subDat, use = "pairwise.complete.obs")
  # find and remove highly correlated covars
  rmNames <- caret::findCorrelation(x = corrMat, cutoff = 0.6)
  # record remaining combo of low-correlation covars
  candModCovars[[ii]] <- sort(propNames[-rmNames])
  
  # # could also base off of p-value threshold
  # subDat <- datGAM %>% dplyr::select(sardRec, all_of(propNames))
  # # candSel <- fuzzySim::corSelect(data = subDat, sp.cols = "sardRec", var.cols = names(subDat)[-1],
  # #                      coeff = TRUE, cor.thresh = 0.6) # based on coefficient threshold
  # candSel <- fuzzySim::corSelect(data = subDat, sp.cols = "sardRec", var.cols = names(subDat)[-1],
  #                      coeff = FALSE) # based on p-value cutoff (0.05)
  # candModCovars[[ii]] <- sort(candSel$selected.vars)
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
singles <- diag(10) %>% as_tibble()
names(singles) <- names(candMods)
candMods[is.na(candMods)] <- 0
candMods <- bind_rows(candMods, singles)
candMods <- unique(candMods) 
dim(candMods) 
candMods <- candMods %>% arrange(CUTI_39N, NCOPspring, OC_LUSI_39N, ZM_SoCal, 
                                 summerSST, ZM_NorCal, NCOPsummerlag1, springSST, 
                                 sardNurseHab, BEUTI_39N) 
candMods %>% print(n=101)

write_csv(candMods, file = "out/candidateGAMmodels.csv")
