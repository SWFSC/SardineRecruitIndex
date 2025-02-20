# Leave-Future-Out Cross  validation for MARSS DFA model selection
# Created: 12/13/2023, Robert Wildermuth


# Leave-Future-Out Cross-Validation ---------------------------------------

LFOXV <- function(dfaDat, # data matrix formatted for MARSS input (variables in rows, time in columns)
                  Rstructure, # either MARSS model 'R' matrix character arguments or custom matrix format
                  mTrends, # number of trends as MARSS model 'm' argument
                  peels, # number of years to peel from end of time series
                  colsRMSE = c("sardRec", "anchRec") # column names to calculate RMSE for. Default: c("sardRec", "anchRec")
){
  
  peelRMSE <- tibble(.rownames = '', 
                     t = 0)[0,]
  
  
  for(i in 1:peels){
    
    cat("\nFitting peel: ", i)
    
    peelDFA <- MARSS(y = dfaDat[,1:(ncol(dfaDat)-i)], 
                     form = "dfa",
                     method = "BFGS",
                     inits = list(x0 = matrix(1, 1, 1)),
                     z.score = TRUE,
                     model = list(R = Rstructure,
                                  m = mTrends),
                     silent = TRUE)
    
    # check convergence
    if(peelDFA$convergence != 0){
      stop("\n Failed to converge for peel ", i)
      return(-99999) # default RMSE for failed convergence
    }
    
    # get one-step-ahead residuals
    resids <- OSAResids(objMARSS = peelDFA, fullDat = dfaDat, p = i,
                        colsRMSE = colsRMSE)
    
    # collect residuals for datum of interest
    peelRMSE <- resids %>% select(.rownames, t, resid.Naiv,
                                  resid.Inf, resid.Cont, resid.Proj) %>% 
                  filter(.rownames %in% colsRMSE, t == max(t)) %>%
                  mutate(peel = i) %>%
                  bind_rows(peelRMSE)
  } # end peel for-loop
    
  overPeel <- peelRMSE %>% 
                  pivot_longer(cols = c(resid.Naiv, resid.Inf, 
                                        resid.Cont, resid.Proj),
                               values_to = "resids", names_to = "resType") %>%
                  # Calc RMSE over all time points
                  group_by(.rownames, resType) %>%
                  summarize(sosRes = sum(resids^2, na.rm = TRUE),
                            nObs = sum(!is.na(resids)),
                            avgAbsRes = mean(abs(resids), na.rm = TRUE)) %>%
                  mutate(RMSE = sqrt(sosRes/nObs)) %>% arrange(RMSE)
  overPeelthenVar <- overPeel %>% group_by(resType) %>%
                  summarize(totRMSE = sum(RMSE),
                            meanRMSE = mean(RMSE))
  # Over variables then peels
  overVarthenPeel <- peelRMSE %>% 
                      pivot_longer(cols = c(resid.Naiv, resid.Inf,
                                            resid.Cont, resid.Proj),
                                   values_to = "resids", names_to = "resType") %>%
                      # Calc RMSE over all time points
                      group_by(peel, resType) %>%
                      summarize(sosRes = sum(resids^2, na.rm = TRUE),
                                nObs = sum(!is.na(resids)),
                                avgAbsRes = mean(abs(resids), na.rm = TRUE)) %>%
                      mutate(RMSE = sqrt(sosRes/nObs)) %>% arrange(RMSE) %>%
                      group_by(resType) %>%
                      summarize(totRMSE = sum(RMSE),
                                meanRMSE = mean(RMSE))
   # OR can calculate over all vars of interest AND peel
  overPeelandVar <- peelRMSE %>% 
                      pivot_longer(cols = c(resid.Naiv, resid.Inf,
                                            resid.Cont, resid.Proj),
                                   values_to = "resids", names_to = "resType") %>%
                      # Calc RMSE over all time points
                      group_by(resType) %>%
                      summarize(sosRes = sum(resids^2, na.rm = TRUE),
                                nObs = sum(!is.na(resids)),
                                avgAbsRes = mean(abs(resids), na.rm = TRUE)) %>%
                      mutate(RMSE = sqrt(sosRes/nObs)) %>% arrange(RMSE) 
  
  # set up ability to compare across calculation methods
  overPeelandVar <- overPeelandVar %>% mutate(compMethod = "peel & variable combined",
                                              variable = "total") %>% 
                        select(compMethod, variable, resType, RMSE)
  overPeel <- overPeel %>% rename(variable = .rownames) %>%
                mutate(compMethod = "individual") %>% 
                select(compMethod, variable, resType, RMSE)
  overPeelthenVar <- overPeelthenVar %>% mutate(compMethod = "peel then variable") %>%
                        pivot_longer(cols = c(totRMSE, meanRMSE), 
                                     names_to = "variable", values_to = "RMSE")
  
  LFOIC <- bind_rows(overPeel, overPeelandVar, overPeelthenVar)
            
  return(LFOIC)
}