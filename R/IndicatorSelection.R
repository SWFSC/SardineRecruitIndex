# Code to evaluate multicollinearity in model variables and select candidate subsets for analyses
# Created: 7/17/2025, Robert Wildermuth

library(tidyverse)
library(MARSS)
library(corrplot)

# read prepped dataset
datDFA <- read_csv("C:/Users/r.wildermuth/Documents/FutureSeas/RecruitmentIndex/recrmntDFA/recrDFAdat.csv")

allDat <- datDFA %>% filter(year %in% 1990:2019) %>%
  # remove contemporary adult biomass with recruits, should be S2 biomass -> S1 recs
  select(-c(NOI,
            ENSO,
            NPGO,
            #All_Copepods, # ~same as calanoid copepods
            #euphausiids, # too large for larval mouth gape
            anchBioSmrySeas1,
            sardBioSmrySeas1,
            anchBioSmrySeas2, # leave out biomass since not fit well
            sardBioSmrySeas2,
            # NCOPspring,
            # SCOPspring, # summer copepod index had higher loadings
            PDOsummer, # lower loading than spring, may want to try a lag
            PDOspring,
            NCOPsummer,
            SCOPsummer,
            # BEUTI_33N, # oceanography at 39N had highest loadings
            # OC_LUSI_33N,
            # OC_LUSI_36N,
            # OC_STI_33N,
            # OC_STI_36N,
            copMeanSize,
            copPropBioSize,
            naupMeanSize,
            naupPropBioSize,
            ZL_NorCal,
            ZL_SoCal,
            age1SprAnchmeanWAA)) # not enough data in time windowselect(-c(#sprCalCOFILarvalSardine,

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
fuzzySim::corSelect(data = t(datZscore), var.cols = datNames,
                    coeff = TRUE) # based on correlation coefficient magnitude (0.8)

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
corrMat[rownames(corrMat) %in% c("HCI_R3", "HCI_R4", "springSST", "summerSST"),
        colnames(corrMat) %in% c("HCI_R3", "HCI_R4", "springSST", "summerSST")]
# springSST highly correlated with all
# HCI in R3 and R4 highly correlated
# HCI_R3 highly correlated with summerSST

# check correlations between zooplankton indicators
corrMat[rownames(corrMat) %in% c("copBio", "naupBio", "ZM_SoCal"),
        colnames(corrMat) %in% c("copBio", "naupBio", "ZM_SoCal")]
# copBio and ZM_SoCal somewhat correlated in SoCal Bight

corrMat[rownames(corrMat) %in% c("NCOPspring", "SCOPspring", "ZM_NorCal"),
        colnames(corrMat) %in% c("NCOPspring", "SCOPspring", "ZM_NorCal")]
# NCOP and SCOP strongly negatively correlated - can drop SCOP
# both highly correlated with ZM_NorCal
corrMat[rownames(corrMat) %in% c("NCOPspring", "SCOPspring", "SCOPsummerlag1", "NCOPsummerlag1"),
        colnames(corrMat) %in% c("NCOPspring", "SCOPspring", "SCOPsummerlag1", "NCOPsummerlag1")]
#SCOPsummerlag1 less correlated with others - keep
corrMat[rownames(corrMat) %in% c("ZM_SoCal", "ZM_NorCal"),
        colnames(corrMat) %in% c("ZM_SoCal", "ZM_NorCal")]

# check correlations between advection indicators
corrMat[rownames(corrMat) %in% c("avgSSWIspring", "avgOffTransspring", "avgNearTransspring"),
        colnames(corrMat) %in% c("avgSSWIspring", "avgOffTransspring", "avgNearTransspring")]

corrMat[rownames(corrMat) %in% c("avgSSWIsummer", "avgOffTranssummer", "avgNearTranssummer"),
        colnames(corrMat) %in% c("avgSSWIsummer", "avgOffTranssummer", "avgNearTranssummer")]
# strongest association with SSWI and NearTrans in summer

# Check correlations with condition factors
corrMat[rownames(corrMat) %in% c("age1SprSardmeanWAA", "meanSSBwt"),
        colnames(corrMat) %in% c("age1SprSardmeanWAA", "meanSSBwt")]
# not highly correlated - keep both

#### Final Selection of indicators ####
# keep observation data indicators, remove associated model-derived indicators
localModel <- c("HCI_R4", "NCOPspring", "NCOPsummerlag1", "SCOPsummerlag1", 
                "BEUTI_33N", "BEUTI_39N", "CUTI_33N", "OC_LUSI_33N", "OC_LUSI_39N",
                "OC_STI_33N", "OC_STI_39N", "swfscRockfishSurv_Myctophids",
                "avgSSWIspring", "avgSSWIsummer", "sardLarv", "anchLarv", 
                "mesopelLarv", "anchYoY", "age1SprSardmeanWAA", "meanSSBwt",
                "copBio", "naupBio", "sardSpawnHab", "anchSpawnHab", 
                "daysAbove5pct", "daysAbove40pct", "sardNurseHab", "anchNurseHab",
                "anchRec", "sardRec", "summerSST", "albacore", "hake")

# remove observation indicators, keep associated model-derived indicators
projectModel <- c("HCI_R4", "BEUTI_33N", "BEUTI_39N", "CUTI_33N", "OC_LUSI_33N",
                  "OC_LUSI_39N", "OC_STI_33N", "OC_STI_39N", "ZM_NorCal",
                  "ZM_SoCal", "sardSpawnHab", "anchSpawnHab", "daysAbove5pct",
                  "daysAbove40pct", "sardNurseHab", "anchNurseHab", "anchRec", 
                  "sardRec", "summerSST", "avgNearTransspring", "avgNearTranssummer", 
                  "avgOffTransspring", "avgOffTranssummer")


# GAM exploration ---------------------------------------------------------

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
  # rmNames <- caret::findCorrelation(x = corrMat, cutoff = 0.8)
  # # record remaining combo of low-correlation covars
  # candModCovars[[ii]] <- sort(propNames[-rmNames])
  
  # could also base off of p-value threshold
  subDat <- datGAM %>% dplyr::select(dev, all_of(propNames))
  candSel <- fuzzySim::corSelect(data = subDat, sp.cols = "dev", var.cols = names(subDat)[-1],
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