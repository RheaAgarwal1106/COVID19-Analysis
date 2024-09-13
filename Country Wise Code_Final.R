#LOADING LIBRARIES
#lag value is fixed in line 84
library(timetk)
library(tidyverse)
library(tvReg)
library(dplyr)
library(car)
library(lubridate)
library(broom)
library(knitr)
library(gt)
library(ggplot2)
library(segmented)
library(strucchange)
library(KernSmooth)
library(nprobust)
library(Ake)
library(plugdensity)
library(np)
library(ks)
library(iosmooth)
library(nprobust)
#------------------------------------------------------------------------------------------------------------------------------------
y_pred_lc_list <- list()
y_pred_ll_list <- list()
y_actual_list <- list()
MSPE_list_lc_bc <- list()
MSPE_list_lc_us <- list()
MSPE_list_ll_bc <- list()
MSPE_list_ll_us <- list()
MSPE <- list()
Lags <- list()
#-------------------------------------------------------------------------------------------------------------------------------------

#PREPARING DATASET

LocationsList <- list("Canada", "United States", "Israel", "South Africa", "Japan", "South Korea", "Brazil", "Pakistan", "Thailand")
Omicron_Start_Dates_Countries <- list("2021-12-12", "2021-12-11", "2021-12-30", "2021-11-27", "2022-01-10","2022-01-30", "2022-01-01", "2022-01-01", "2022-01-01")

#Reporting end dates for each country (found by visual inspection of the daily new cases plots)
Reporting_End_Dates_Countries <- list("2022-06-11","2022-10-21","2022-12-31","2022-07-30","2022-12-31","2022-12-31","2022-04-01","2022-10-01","2022-10-01")


for (region in LocationsList){
 
  y_pred_lc_list <- list()
  y_pred_ll_list <- list()
  y_actual_list <- list()
  MSPE_list_lc_bc <- list()
  MSPE_list_lc_us <- list()
  MSPE_list_ll_bc <- list()
  MSPE_list_ll_us <- list()
  y_pred_lc_bc_list <- list()
  y_pred_lc_us_list <- list()
  y_pred_ll_bc_list <- list()
  y_pred_lc_us_list <- list()
  daily_pred_lc_bc_list <- list()
  daily_pred_lc_us_list <- list()
  daily_pred_ll_bc_list <- list()
  daily_pred_ll_us_list <- list()
  best_lag_lc_bc <- list()
  best_lag_lc_us <- list()
  best_lag_ll_bc <- list()
  best_lag_ll_us <- list()
  
  urlnewOWIDdata <- 'https://covid.ourworldindata.org/data/owid-covid-data.csv'
  NewOWIDData <- read.csv(urlnewOWIDdata)
  NewOWIDData <- NewOWIDData %>% dplyr::select(date, location, total_cases, new_cases, total_deaths, new_deaths)
  
  NewOWIDData$Date <- as.Date(NewOWIDData$date, format="%Y-%m-%d") 
  
  cutoff_end <- "2022-12-31"
  NewOWIDData <- NewOWIDData %>% filter(Date <= cutoff_end)
  
  
  NewOWIDSingleRegionData <- NewOWIDData %>% filter(location %in% region)
  
  Omicron_Start_Date <-   Omicron_Start_Dates_Countries[which(LocationsList == region)]   #omicron start date
  Reporting_End_Date <- Reporting_End_Dates_Countries[which(LocationsList == region)]      #omicron end date
  
  #Filter the data to only use data after the start of the Omicron wave in the current region
  NewOWIDSingleRegionData <- NewOWIDSingleRegionData %>% filter(date >= Omicron_Start_Date)
  NewOWIDSingleRegionData <- NewOWIDSingleRegionData %>% filter(date <= Reporting_End_Date)
  
  #Select the Y variable
  death_long <- NewOWIDSingleRegionData %>% dplyr::select(date, new_deaths, total_deaths)  #total_deaths are cumulative_deaths
  names(death_long) <- c("Date","deaths", "cumulative_deaths")  #changing column names
  
  #Select the X variable
  cases_long <- NewOWIDSingleRegionData %>% dplyr::select(date, new_cases, total_cases)  #total cases are cumulative cases
  names(cases_long) <- c("Date","cases", "cumulative_cases")  #changing column names
  
  #Combine the Y and X variables into a single data frame
  df <-  merge (death_long[,c( "deaths", "cumulative_deaths","Date" )], 
                cases_long[,c( "cases", "cumulative_cases", "Date" ) ], 
                by = "Date") 
  
  #Reformat the dates in the data frame
  df$Date <- as.Date(df$Date, format="%Y-%m-%d")
  
  
  #Print the name of the current country to the console
  print(paste("For", region))
  
  names(df) <- c("Date", "y", "cumulative_y", "x", "cumulative_x")   #renaming column names
  
  #resetting x and y
  df$cumulative_y_reset <- df$cumulative_y - df$cumulative_y[1] 
  df$cumulative_x_reset <- df$cumulative_x - df$cumulative_x[1]
  
  #-------------------------------------------------------------------------------------------------------------------------------------
  max_lead <- 21
  for (lag in 5:21){
    y_pred_lc_bc_list <- list()
    y_pred_lc_us_list <- list()
    y_pred_ll_bc_list <- list()
    y_pred_ll_us_list <- list()
    
    y_actual_list <- list()
    
    D1 <- max(df$Date) -  2 * max_lead
    
    df_tv <- df %>%
      tk_augment_lags(c(cumulative_y_reset, y), .lags = -lag, .names = c("y_best_lead_tv", "daily_y_lead"))
    
    df_insample <- df_tv %>% filter (df$Date <= D1)   #insample data. if lag = 7, 2*7=14 dates prior to the last date.
    
    D2 <- max(df$Date) -  max_lead   #--------------------|--------|--------   D1 and D2 are dates
    #D1       D2
    
    #insample       outsample
    
    df_newsample <- df_tv %>% filter (df$Date > D1 & df$Date <= D2)#Between D1 and D2 predictions are to be made. After D2 we need only to calculate MPSE corresponding to values between D1 and D2. No predictions after D2.
    
    #Implementing lprobust for the first time on insample data
    
    X_insample <- df_tv$cumulative_x_reset[1: nrow(df_insample)]
    y_insample <- df_tv$cumulative_y_reset[(1+lag): (nrow(df_insample)+lag)]
    
    lc_result_new <- lprobust(y = y_insample, x = X_insample, eval = df_tv$cumulative_x_reset[nrow(df_insample)+1], p = 0, kernel ='tri')
    ll_result_new <- lprobust(y = y_insample, x = X_insample, eval = df_tv$cumulative_x_reset[nrow(df_insample)+1], p = 1, kernel = 'tri')
    
    y_pred_lc_bc_list <- append(y_insample, lc_result_new$Estimate[,6])
    y_pred_lc_us_list <- append(y_insample, lc_result_new$Estimate[,5])
    y_pred_ll_bc_list <- append(y_insample, ll_result_new$Estimate[,6])
    y_pred_ll_us_list <- append(y_insample, ll_result_new$Estimate[,5])
    
    y_actual <- df_tv$daily_y_lead[nrow(df_insample)+lag+1]
    y_actual_list <- append(y_actual_list, y_actual)
    
    #calling lprobust function one at a time
    for (i in 2:lag){
      
      X_new <- df_tv$cumulative_x_reset[1: (nrow(df_insample)+i-1)]  #i=1, original insample data, i=2, insample+1 datapoint, i=3 original insample +2 datapoints
      
      lc_bc_result_new <- lprobust(y = y_pred_lc_bc_list, x = X_new, eval = df_tv$cumulative_x_reset[nrow(df_insample)+i], p = 0, kernel = 'tri')
      
      lc_us_result_new <- lprobust(y = y_pred_lc_us_list, x = X_new, eval = df_tv$cumulative_x_reset[nrow(df_insample)+i], p = 0, kernel = 'tri')
      
      ll_bc_result_new <- lprobust(y = y_pred_ll_bc_list, x = X_new, eval = df_tv$cumulative_x_reset[nrow(df_insample)+i], p = 1, kernel = 'tri')
      
      ll_us_result_new <- lprobust(y = y_pred_ll_us_list, x = X_new, eval = df_tv$cumulative_x_reset[nrow(df_insample)+i], p = 1, kernel = 'tri')
      
      
      #store this prediction and the corresponding original y value. now take ith datapoint inside the sample and repeat the procedure.
      #After first iteration, we have 5 predictions corresponding to last 5 X values. Calculate MPSE. That will be the MPSE of lag =5.
      #Similarly for all lags upto 21. Find lowest MPSE. That will be the lag.
      
      y_pred_lc_bc <- lc_bc_result_new$Estimate[,6]
      y_pred_lc_us <- lc_us_result_new$Estimate[,5]
      y_pred_ll_bc <- ll_bc_result_new$Estimate[,6]
      y_pred_ll_us <- ll_us_result_new$Estimate[,5]
      
      y_actual <- df_tv$daily_y_lead[nrow(df_insample)+lag+i]   
      
      y_pred_lc_bc_list <- append(y_pred_lc_bc_list, y_pred_lc_bc)
      y_pred_lc_us_list <- append(y_pred_lc_us_list, y_pred_lc_us)
      y_pred_ll_bc_list <- append(y_pred_ll_bc_list, y_pred_ll_bc)
      y_pred_ll_us_list <- append(y_pred_ll_us_list , y_pred_ll_us)
      
      y_actual_list <- append(y_actual_list, y_actual)
    }
    
    y_actual_list <- rev(rev(df_tv$daily_y_lead)[(lag+1):(2*lag)])
    #calculating daily predictions
    daily_pred_lc_bc_list <- tail(y_pred_lc_bc_list, -1) - head(y_pred_lc_bc_list, -1)
    daily_pred_lc_us_list <- tail(y_pred_lc_us_list, -1) - head(y_pred_lc_us_list, -1)
    daily_pred_ll_bc_list <- tail(y_pred_ll_bc_list, -1) - head(y_pred_ll_bc_list, -1)
    daily_pred_ll_us_list <- tail(y_pred_ll_us_list, -1) - head(y_pred_ll_us_list, -1)
    
    #calculating MSPE
    MSPE_lc_bc <- mean((as.numeric(unlist(tail(daily_pred_lc_bc_list, lag)))-as.numeric(unlist(y_actual_list)))^2)
    MSPE_lc_us <- mean((as.numeric(unlist(tail(daily_pred_lc_us_list, lag)))-as.numeric(unlist(y_actual_list)))^2)
    MSPE_ll_bc <- mean((as.numeric(unlist(tail(daily_pred_ll_bc_list, lag)))-as.numeric(unlist(y_actual_list)))^2)
    MSPE_ll_us <- mean((as.numeric(unlist(tail(daily_pred_ll_us_list, lag)))-as.numeric(unlist(y_actual_list)))^2)
    
    #storing MSPE
    MSPE_list_lc_bc <- append(MSPE_list_lc_bc, MSPE_lc_bc)
    MSPE_list_lc_us <- append(MSPE_list_lc_us, MSPE_lc_us)
    MSPE_list_ll_bc <- append(MSPE_list_ll_bc, MSPE_ll_bc)
    MSPE_list_ll_us <- append(MSPE_list_ll_us, MSPE_ll_us)
  }
 
  #finding the index of the smallest lag
  best_lag_lc_bc <- which.min(MSPE_list_lc_bc)+4
  min(unlist(MSPE_list_lc_bc))
  best_lag_lc_us <- which.min(MSPE_list_lc_us)+4
  min(unlist(MSPE_list_lc_us))
  best_lag_ll_bc <- which.min(MSPE_list_ll_bc)+4
  min(unlist(MSPE_list_ll_bc))
  best_lag_ll_us <- which.min(MSPE_list_ll_us)+4
  min(unlist(MSPE_list_ll_us))
  
  
  #_______________________________________________________________________________________________________________________
  
  #Running the code with optimal lag (lc_bc)
  
  y_pred_lc_bc_list <- list()
  y_pred_lc_us_list <- list()
  y_pred_ll_bc_list <- list()
  y_pred_ll_us_list <- list()
  
  y_actual_list <- list()
  
  D1 <- max(df$Date) -  2 * best_lag_lc_bc
  
  df_tv <- df %>%
    tk_augment_lags(c(cumulative_y_reset, y), .lags = -best_lag_lc_bc, .names = c("y_best_lead_tv", "daily_y_lead"))
  
  df_insample <- df_tv %>% filter (df$Date <= D1)   #insample data. if lag = 7, 2*7=14 dates prior to the last date.
  
  D2 <- max(df$Date) -  best_lag_lc_bc   #--------------------|--------|--------   D1 and D2 are dates
  #D1       D2
  
  #insample       outsample
  
  df_newsample <- df_tv %>% filter (df$Date > D1 & df$Date <= D2)#Between D1 and D2 predictions are to be made. After D2 we need only to calculate MPSE corresponding to values between D1 and D2. No predictions after D2.
  
  #Implementing lprobust for the first time on insample data
  
  X_insample <- df_tv$cumulative_x_reset[1: nrow(df_insample)]
  y_insample <- df_tv$cumulative_y_reset[(1+best_lag_lc_bc): (nrow(df_insample)+best_lag_lc_bc)]
  
  lc_result_new <- lprobust(y = y_insample, x = X_insample, eval = df_tv$cumulative_x_reset[nrow(df_insample)+1], p = 0, kernel = 'tri')
  
  y_pred_lc_bc_list <- append(y_insample, lc_result_new$Estimate[,6])
  
  y_actual <- df_tv$daily_y_lead[(nrow(df_insample))+(best_lag_lc_bc+1)]
  y_actual_list <- append(y_actual_list, y_actual)
  
  #calling lprobust function one at a time
  for (i in 2:best_lag_lc_bc){
    
    X_new <- df_tv$cumulative_x_reset[1: (nrow(df_insample)+i-1)]  #i=1, original insample data, i=2, insample+1 datapoint, i=3 original insample +2 datapoints
    
    lc_result_new <- lprobust(y = y_pred_lc_bc_list, x = X_new, eval = df_tv$cumulative_x_reset[nrow(df_insample)+i], p = 0, kernel = 'tri')
    
    #store this prediction and the corresponding original y value. now take ith datapoint inside the sample and repeat the procedure.
    #After first iteration, we have 5 predictions corresponding to last 5 X values. Calculate MPSE. That will be the MPSE of lag =5.
    #Similarly for all lags upto 21. Find lowest MPSE. That will be the lag.
    
    y_pred_lc_bc <- lc_result_new$Estimate[,6]
    
    y_actual <- df_tv$daily_y_lead[(nrow(df_insample))+(best_lag_lc_bc+i)]
    
    y_pred_lc_bc_list <- append(y_pred_lc_bc_list, y_pred_lc_bc)
    
    y_actual_list <- append(y_actual_list, y_actual)
  }
  y_actual_list <- rev(rev(df_tv$daily_y_lead)[(best_lag_lc_bc+1):(2*best_lag_lc_bc)])
  
  daily_pred_lc_bc_list <- tail(y_pred_lc_bc_list, -1) - head(y_pred_lc_bc_list, -1)
  
  MSPE_lc_bc <- mean((as.numeric(unlist(tail(daily_pred_lc_bc_list, best_lag_lc_bc)))-as.numeric(unlist(y_actual_list)))^2)
  
  print(paste(MSPE_lc_bc, best_lag_lc_bc))
  
  # pdf(file = paste(region,"local_constant_unbiased.pdf"))
  # 
  # plot(df$cumulative_x_reset[(nrow(df_insample)+1) : (nrow(df_insample)+best_lag_lc_bc)], tail(daily_pred_lc_bc_list, best_lag_lc_bc), xlab ='Cumulative Cases', ylab='Cumulative Deaths', type = 'l', col = 'red')
  # lines(df$cumulative_x_reset[(nrow(df_insample)+1) : (nrow(df_insample)+best_lag_lc_bc)], y_actual_list, type='l', col='green', xlab='Cumulative Cases', ylab='Cumulative Deaths')
  # title(main='Local Constant Unbiased')
  # 
  # dev.off()
  #_______________________________________________________________________________________________________________________
  
  
  y_pred_lc_bc_list <- list()
  y_pred_lc_us_list <- list()
  y_pred_ll_bc_list <- list()
  y_pred_ll_us_list <- list()
  
  y_actual_list <- list()
  
  D1 <- max(df$Date) -  2 * best_lag_lc_us
  
  df_tv <- df %>%
    tk_augment_lags(c(cumulative_y_reset, y), .lags = -best_lag_lc_us, .names = c("y_best_lead_tv", "daily_y_lead"))
  
  
  df_insample <- df_tv %>% filter (df$Date <= D1)   #insample data. if lag = 7, 2*7=14 dates prior to the last date.
  
  D2 <- max(df$Date) -  best_lag_lc_us   #--------------------|--------|--------   D1 and D2 are dates
  #D1       D2
  
  #insample       outsample
  
  df_newsample <- df_tv %>% filter (df$Date > D1 & df$Date <= D2)#Between D1 and D2 predictions are to be made. After D2 we need only to calculate MPSE corresponding to values between D1 and D2. No predictions after D2.
  
  #Implementing lprobust for the first time on insample data
  
  X_insample <- df_tv$cumulative_x_reset[1: nrow(df_insample)]
  y_insample <- df_tv$cumulative_y_reset[(1+best_lag_lc_us): (nrow(df_insample)+best_lag_lc_us)]
  
  lc_result_new <- lprobust(y = y_insample, x = X_insample, eval = df_tv$cumulative_x_reset[nrow(df_insample)+1], p = 0, kernel = 'tri')
  
  y_pred_lc_us_list <- append(y_insample, lc_result_new$Estimate[,5])
  
  y_actual <- df_tv$daily_y_lead[nrow(df_insample)+best_lag_lc_us+1]
  y_actual_list <- append(y_actual_list, y_actual)
  
  #calling lprobust function one at a time
  for (i in 2:best_lag_lc_us){
    
    X_new <- df_tv$cumulative_x_reset[1: (nrow(df_insample)+i-1)]  #i=1, original insample data, i=2, insample+1 datapoint, i=3 original insample +2 datapoints
    
    lc_result_new <- lprobust(y = y_pred_lc_us_list, x = X_new, eval = df_tv$cumulative_x_reset[nrow(df_insample)+i], p = 0, kernel = 'tri')
    
    #store this prediction and the corresponding original y value. now take ith datapoint inside the sample and repeat the procedure.
    #After first iteration, we have 5 predictions corresponding to last 5 X values. Calculate MPSE. That will be the MPSE of lag =5.
    #Similarly for all lags upto 21. Find lowest MPSE. That will be the lag.
    
    y_pred_lc_us <- lc_result_new$Estimate[,5]
    
    y_actual <- df_tv$daily_y_lead[nrow(df_insample)+(best_lag_lc_us+i)]   
    
    y_pred_lc_us_list <- append(y_pred_lc_us_list, y_pred_lc_us)
    
    y_actual_list <- append(y_actual_list, y_actual)
  }
  y_actual_list <- rev(rev(df_tv$daily_y_lead)[(best_lag_lc_us+1):(2*best_lag_lc_us)])
  
  MSPE_lc_us <- mean((as.numeric(unlist(tail(y_pred_lc_us_list, best_lag_lc_us)))-as.numeric(unlist(y_actual_list)))^2)
  
  daily_pred_lc_us_list <- tail(y_pred_lc_us_list, -1) - head(y_pred_lc_us_list, -1)
  
  MSPE_lc_us <- mean((as.numeric(unlist(tail(daily_pred_lc_us_list, best_lag_lc_us)))-as.numeric(unlist(y_actual_list)))^2)
  
  print(paste(MSPE_lc_us, best_lag_lc_us))
  
  # pdf(file = paste(region,"local_constant_biased.pdf"))
  # 
  # plot(df$cumulative_x_reset[(nrow(df_insample)+1) : (nrow(df_insample)+best_lag_lc_us)], tail(daily_pred_lc_us_list, best_lag_lc_us), xlab ='Cumulative Cases', ylab='Cumulative Deaths', type = 'l', col = 'red')
  # lines(df$cumulative_x_reset[(nrow(df_insample)+1) : (nrow(df_insample)+best_lag_lc_us)], y_actual_list, type='l', col='green', xlab='Cumulative Cases', ylab='Cumulative Deaths')
  # title(main='Local Constant Biased')
  # 
  # dev.off()
  #_______________________________________________________________________________________________________________________
  
  y_pred_lc_bc_list <- list()
  y_pred_lc_us_list <- list()
  y_pred_ll_bc_list <- list()
  y_pred_ll_us_list <- list()
  
  y_actual_list <- list()
  
  D1 <- max(df$Date) -  2 * best_lag_ll_bc
  
  df_tv <- df %>%
    tk_augment_lags(c(cumulative_y_reset, y), .lags = -best_lag_ll_bc, .names = c("y_best_lead_tv", "daily_y_lead"))
  
  
  df_insample <- df_tv %>% filter (df$Date <= D1)   #insample data. if lag = 7, 2*7=14 dates prior to the last date.
  
  D2 <- max(df$Date) -  best_lag_ll_bc   #--------------------|--------|--------   D1 and D2 are dates
  #D1       D2
  
  #insample       outsample
  
  df_newsample <- df_tv %>% filter (df$Date > D1 & df$Date <= D2)#Between D1 and D2 predictions are to be made. After D2 we need only to calculate MPSE corresponding to values between D1 and D2. No predictions after D2.
  
  #Implementing lprobust for the first time on insample data
  
  X_insample <- df_tv$cumulative_x_reset[1: nrow(df_insample)]
  y_insample <- df_tv$cumulative_y_reset[(1+best_lag_ll_bc): (nrow(df_insample)+best_lag_ll_bc)]
  
  ll_result_new <- lprobust(y = y_insample, x = X_insample, eval = df_tv$cumulative_x_reset[nrow(df_insample)+1], p = 1, kernel = 'tri')
  
  y_pred_ll_bc_list <- append(y_insample, ll_result_new$Estimate[,6])
  
  y_actual <- df_tv$daily_y_lead[nrow(df_insample)+best_lag_ll_bc+1]
  y_actual_list <- append(y_actual_list, y_actual)
  
  #calling lprobust function one at a time
  for (i in 2:best_lag_ll_bc){
    
    X_new <- df_tv$cumulative_x_reset[1: (nrow(df_insample)+i-1)]  #i=1, original insample data, i=2, insample+1 datapoint, i=3 original insample +2 datapoints
    
    ll_result_new <- lprobust(y = y_pred_ll_bc_list, x = X_new, eval = df_tv$cumulative_x_reset[nrow(df_insample)+i], p = 1, kernel = 'tri')
    
    #store this prediction and the corresponding original y value. now take ith datapoint inside the sample and repeat the procedure.
    #After first iteration, we have 5 predictions corresponding to last 5 X values. Calculate MPSE. That will be the MPSE of lag =5.
    #Similarly for all lags upto 21. Find lowest MPSE. That will be the lag.
    
    y_pred_ll_bc <- ll_result_new$Estimate[,6]
    
    y_actual <- df_tv$cumulative_y_reset[nrow(df_insample)+best_lag_ll_bc+i]   
    
    y_pred_ll_bc_list <- append(y_pred_ll_bc_list, y_pred_ll_bc)
    
    y_actual_list <- append(y_actual_list, y_actual)
  }
  
  y_actual_list <- rev(rev(df_tv$daily_y_lead)[(best_lag_ll_bc+1):(2*best_lag_ll_bc)])
  
  daily_pred_lc_bc_list <- tail(y_pred_ll_bc_list, -1) - head(y_pred_ll_bc_list, -1)
  
  MSPE_lc_bc <- mean((as.numeric(unlist(tail(daily_pred_ll_bc_list, best_lag_ll_bc)))-as.numeric(unlist(y_actual_list)))^2)
  
  print(paste(MSPE_ll_bc, best_lag_ll_bc))
  
  # pdf(file = paste(region,"local_linear_unbiased.pdf"))
  # 
  # plot(df$cumulative_x_reset[(nrow(df_insample)+1) : (nrow(df_insample)+best_lag_ll_bc)], tail(daily_pred_ll_bc_list, best_lag_ll_bc), xlab ='Cumulative Cases', ylab='Cumulative Deaths', type = 'l', col = 'red')
  # lines(df$cumulative_x_reset[(nrow(df_insample)+1) : (nrow(df_insample)+best_lag_ll_bc)], y_actual_list, type='l', col='green', xlab='Cumulative Cases', ylab='Cumulative Deaths')
  # title(main='Local Linear Unbiased')
  # 
  # dev.off()
  #_______________________________________________________________________________________________________________________
  
  y_pred_lc_bc_list <- list()
  y_pred_lc_us_list <- list()
  y_pred_ll_bc_list <- list()
  y_pred_ll_us_list <- list()
  
  y_actual_list <- list()
  
  D1 <- max(df$Date) -  2 * best_lag_ll_us
  
  df_tv <- df %>%
    tk_augment_lags(c(cumulative_y_reset, y), .lags = -best_lag_ll_us, .names = c("y_best_lead_tv", "daily_y_lead"))
  
  df_insample <- df_tv %>% filter (df$Date <= D1)   #insample data. if lag = 7, 2*7=14 dates prior to the last date.
  
  D2 <- max(df$Date) -  best_lag_ll_us   #--------------------|--------|--------   D1 and D2 are dates
  #D1       D2
  
  #insample       outsample
  
  df_newsample <- df_tv %>% filter (df$Date > D1 & df$Date <= D2)#Between D1 and D2 predictions are to be made. After D2 we need only to calculate MPSE corresponding to values between D1 and D2. No predictions after D2.
  
  #Implementing lprobust for the first time on insample data
  
  X_insample <- df_tv$cumulative_x_reset[1: nrow(df_insample)]
  y_insample <- df_tv$cumulative_y_reset[(1+best_lag_ll_us): (nrow(df_insample)+best_lag_ll_us)]
  
  ll_result_new <- lprobust(y = y_insample, x = X_insample, eval = df_tv$cumulative_x_reset[nrow(df_insample)+1], p = 1, kernel = 'tri')
  
  y_pred_ll_us_list <- append(y_insample, ll_result_new$Estimate[,5])
  
  y_actual <- df_tv$daily_y_lead[(nrow(df_insample))+(best_lag_ll_us+1)]
  y_actual_list <- append(y_actual_list, y_actual)
  
  #calling lprobust function one at a time
  for (i in 2:best_lag_ll_us){
    
    X_new <- df_tv$cumulative_x_reset[1: (nrow(df_insample)+i-1)]  #i=1, original insample data, i=2, insample+1 datapoint, i=3 original insample +2 datapoints
    
    ll_result_new <- lprobust(y = y_pred_ll_us_list, x = X_new, eval = df_tv$cumulative_x_reset[nrow(df_insample)+i], p = 1, kernel = 'tri')
    
    #store this prediction and the corresponding original y value. now take ith datapoint inside the sample and repeat the procedure.
    #After first iteration, we have 5 predictions corresponding to last 5 X values. Calculate MPSE. That will be the MPSE of lag =5.
    #Similarly for all lags upto 21. Find lowest MPSE. That will be the lag.
    
    y_pred_ll_us <- ll_result_new$Estimate[,5]
    
    y_actual <- df_tv$daily_y_lead[(nrow(df_insample))+(best_lag_ll_us+i)]
    
    y_pred_ll_us_list <- append(y_pred_ll_us_list, y_pred_ll_us)
    
    y_actual_list <- append(y_actual_list, y_actual)
  }
  y_actual_list <- rev(rev(df_tv$daily_y_lead)[(best_lag_ll_us+1):(2*best_lag_ll_us)])
  
  daily_pred_ll_us_list <- tail(y_pred_ll_us_list, -1) - head(y_pred_ll_us_list, -1)
  
  MSPE_ll_us <- mean((as.numeric(unlist(tail(daily_pred_ll_us_list, best_lag_ll_us)))-as.numeric(unlist(y_actual_list)))^2)
  
  print(paste(MSPE_ll_us, best_lag_ll_us))
  
  # pdf(file = paste(region,"local_linear_biased.pdf"))
  # 
  # plot(df$cumulative_x_reset[(nrow(df_insample)+1) : (nrow(df_insample)+best_lag_ll_us)], tail(daily_pred_ll_us_list, best_lag_ll_us), xlab ='Cumulative Cases', ylab='Cumulative Deaths', type = 'l', col = 'red')
  # lines(df$cumulative_x_reset[(nrow(df_insample)+1) : (nrow(df_insample)+best_lag_ll_us)], y_actual_list, type='l', col='green', xlab='Cumulative Cases', ylab='Cumulative Deaths')
  # title(main='Local Linear Biased')
  # 
  # dev.off()
  #_______________________________________________________________________________________________________________________
  
  MSPE <- append(MSPE, c(MSPE_lc_bc, MSPE_lc_us, MSPE_ll_bc, MSPE_ll_us))
  Lags <- append(Lags, c(best_lag_lc_bc, best_lag_lc_us, best_lag_ll_bc, best_lag_ll_us))
}

library(knitr)
library(kableExtra)

Countries <- rep(LocationsList, each = 4)
MSPE_list <- unlist(MSPE)
options(scipen=9999)
Lags_list <- unlist(Lags)
Countries_list <- Countries

results <- data.frame('MSPE(lc_bc lc_us ll_bc ll_us)' = MSPE_list,'Lags' = Lags_list,'Countries' = unlist(Countries_list))

results %>%
  kable() %>%
  kable_styling()

print(results)
