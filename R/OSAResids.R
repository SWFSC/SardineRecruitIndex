# Function to calculate out of sample residuals given a prediction horizon
# where the DFA model used for prediction is only fit to 
# data available at the time of the forecast. Calculates naive and informed innovation 
# residuals, contemporaneous residuals, and projection residuals.
# Created 4/15/2024, Desiree Tommasi, Robert Wildermuth
#' @param objMARSS MARSS object fitted outside this function to data up to a specified peel
#' @param fullDat all available data
#' @param p peel used when fitting objMARSS
#' @param horizon integer for steps ahead in the prediction horizon
#' @param colsRMSE # column names to calculate RMSE for. Default: c("sardRec", "anchRec")
#' @return dataframe of residuals for the specified 
#' datasets, scaled observations fed into the model estimate (y.*),
#' original scaled observations (origDat), prediction (est.*), 
#' and timestep (t). Includes naive innovation residuals (Naiv) (persistence assumption),
#' informed innovation residuals with only 'colsRMSE' data missing from the forecast (Inf),
#' contemporaneous residuals with all data in the forecast (Cont), and projection residuals
#' with 'colsRMSE' data missing from full time series (Proj).
#' Note only the residuals for the last 'horizon' time steps where y=NA are 
#' out of sample and true innovation residuals

OSAResids <- function(objMARSS, fullDat, p, horizon = 1,
                          colsRMSE = c("sardRec", "anchRec") 
){
  
  
  origDat <- fullDat[,1:min(ncol(fullDat), # either full data set or ...
                            ncol(fullDat)-p+horizon)] # including 'horizon' steps ahead of peel
  # need to scale data for residual calcs and 'newdata' input
  origDat <- apply(origDat, 1, scale, simplify = TRUE) %>% t() 
  # !RW: scales slightly different from those used to fit the DFA model 'objMARSS'
  #      because of added data from 1 step ahead
  
  # create 'YStar' by removing last observation of 'colsRMSE' variables
  YStar <- origDat
  YStar[colsRMSE, (ncol(YStar)-p+1):ncol(YStar)] <- NA
  
  # naive innovations
  yExpctNaive <- predict(object = objMARSS, interval = "none", 
                    n.ahead = horizon, type = 'ytT')
  yExpctNaive <- yExpctNaive$pred %>% rename(y.Naiv = y,
                                             est.Naiv = estimate)
  
  # informed innovations
  yExpctInf <- predict(object = objMARSS, interval = "none", 
                    newdata = list(t = 1:ncol(origDat),
                                   y = YStar),
                    type = 'ytT', x0 = "use.model")
  yExpctInf <- yExpctInf$pred %>% rename(y.Inf = y,
                                         est.Inf = estimate)
  
  # contemporaneous
  yExpctCont <- predict(object = objMARSS, interval = "none", 
                       newdata = list(t = 1:ncol(origDat),
                                      y = origDat),
                       type = 'ytT', x0 = "use.model")
  yExpctCont <- yExpctCont$pred %>% rename(y.Cont = y,
                                           est.Cont = estimate)
  
  # projection residuals
  YStar[colsRMSE, ] <- NA
  yExpctProj <- predict(object = objMARSS, n.ahead = 0, interval = "none", 
                    newdata = list(t = 1:ncol(origDat),
                                   y = YStar),
                    type = 'ytT', x0 = "use.model")
  yExpctProj <- yExpctProj$pred %>% rename(y.Proj = y,
                                           est.Proj = estimate)
  
  # transpose for column addition
  trCol <- origDat %>% t() 
  trCol <- as.data.frame(trCol) %>% mutate(t = 1:nrow(trCol)) %>% 
    pivot_longer(cols = -t, names_to = ".rownames", 
                 values_to = "origDat")
  
  osaResids <- yExpctNaive %>% full_join(y = yExpctInf,
                                                   by = c("t", ".rownames")) %>% 
                      full_join(y = yExpctCont, by = c("t", ".rownames")) %>% 
                      full_join(y = yExpctProj, by = c("t", ".rownames")) %>% 
                      full_join(y = trCol, by = c("t", ".rownames"))
  
  osaResids <- osaResids %>% mutate(resid.Naiv = origDat - est.Naiv,
                                    resid.Inf = origDat - est.Inf,
                                    resid.Cont = origDat - est.Cont,
                                    resid.Proj = origDat - est.Proj)
  
  return(osaResids)
}