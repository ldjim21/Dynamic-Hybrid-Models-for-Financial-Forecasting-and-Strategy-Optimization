library(forecast)
library(stats)
library(fracdiff)
library(lubridate)
library(dplyr)
library(lmtest)

btc_arfima <- read.csv("")

btc_arfima$Date <- as.Date(btc_arfima$Date)

rownames(btc_arfima) <- btc_arfima$Date

####################################################################################

CV_ARIMA <- function(df, train_months = 36, test_months = 6) {
  
  df$Date <- as.Date(df$Date)
  
  indices <- list()
  current_start_date <- min(df$Date)
  
  while (TRUE) {
    
    train_end_date <- seq(current_start_date, by = paste0(train_months, " months"), length.out = 2)[2] - 1
    test_start_date <- train_end_date + 1
    test_end_date <- seq(test_start_date, by = paste0(test_months, " months"), length.out = 2)[2] - 1
    
    if (test_end_date - 1 > max(df$Date)) {
      break
    }
    
    train_set <- df$Date[df$Date >= current_start_date & df$Date <= train_end_date]
    test_set <- df$Date[df$Date >= test_start_date & df$Date <= test_end_date]
    
    indices[[length(indices) + 1]] <- list(train_set, test_set)
    
    current_start_date <- seq(current_start_date, by = "6 months", length.out = 2)[2]
  }
  
  return(indices)
}

####################################################################################

indices <- CV_ARIMA(btc_arfima)

for (i in seq_along(indices)) {
  cat(sprintf("Window %d:\n", i))
  cat("Training set:\n")
  start_train = which(rownames(btc_arfima) == indices[[i]][[1]][1])
  end_train = which(rownames(btc_arfima) == indices[[i]][[1]][length(indices[[i]][[1]])])
  print(head(btc_arfima$Date[start_train:end_train], 3))
  print(tail(btc_arfima$Date[end_train:end_train], 3))
  cat("Testing set:\n")
  start_test = which(rownames(btc_arfima) == indices[[i]][[2]][1])
  end_test = which(rownames(btc_arfima) == indices[[i]][[2]][length(indices[[i]][[2]])])
  print(head(btc_arfima$Date[start_test:end_test], 3))
  print(tail(btc_arfima$Date[start_test:end_test], 3))
  cat("\n")
}

####################################################################################

fit_ARFIMA <- function(ts_data) {
  model <- arfima(ts_data, drange=c(0,0.5), start.p = 0, max.p = 5, start.q = 0, max.q = 5, max.order = 20, 
                  ic ="aic", seasonal = FALSE, stepwise = FALSE, allowdrift = FALSE)
  
  
  return(list(
    orders = model[c("ar", "d", "ma")]
  ))
}

####################################################################################

best_arfima_orders_btc <- list()

indices <- CV_ARIMA(btc_arfima)

for (i in seq_along(indices)) {
  
  print(i)
  
  start_train = which(rownames(btc_arfima) == indices[[i]][[1]][1])
  end_train = which(rownames(btc_arfima) == indices[[i]][[1]][length(indices[[i]][[1]])])
  train_set <- btc_arfima$Log_Returns[start_train:end_train]
  
  best_aic_model <- fit_ARFIMA(train_set)
  
  best_arfima_orders_btc[[i]] <- best_aic_model$orders
}

best_arfima_orders_btc

####################################################################################

best_arfima_orders_btc_ranks <- list()

for (i in seq_along(best_arfima_orders_btc)) {
  
  best_arfima_orders_btc_ranks[[i]] <- list(length(best_arfima_orders_btc[[i]][[1]]), length(best_arfima_orders_btc[[i]][[3]]))
  
}

best_arfima_orders_btc_ranks

####################################################################################

ARFIMA_forecasts <- function(df) {
  
  indices_1 <- CV_ARIMA(df, train_months = 42, test_months = 0)
  indices_2 <- CV_ARIMA(df, train_months = 36, test_months = 6)
  
  predictions <- c()
  
  for (i in seq_along(indices_1)) {
    
    start_train = which(rownames(df) == indices_2[[i]][[1]][1])
    end_train = which(rownames(df) == indices_2[[i]][[1]][length(indices_2[[i]][[1]])])
    
    start_test = which(rownames(df) == indices_2[[i]][[2]][1])
    end_test = which(rownames(df) == indices_2[[i]][[2]][length(indices_2[[i]][[2]])])
    
    len_train_series <- end_train - start_train + 1
    len_test_series <- end_test - start_test + 1
    
    start_series = which(rownames(df) == indices_1[[i]][[1]][1])
    end_series = which(rownames(df) == indices_1[[i]][[1]][length(indices_1[[i]][[1]])])
    
    series <- df$Log_Returns[start_series:end_series]
    
    best_aic_model <- fit_ARFIMA(df$Log_Returns[start_train:end_train])
    
    best_arfima_orders <- best_aic_model$orders
    
    best_arfima_orders_ranks <- list(length(best_arfima_orders[[1]]), length(best_arfima_orders[[3]]))
    
    for (i in 1:len_test_series) {
      train_series <- window(series, start=i, end=(len_train_series + i - 1))
      model <- fracdiff(na.omit(train_series), nar = best_arfima_orders_ranks[[1]], nma = best_arfima_orders_ranks[[2]])
      prediction <- forecast(model, h=1)$mean
      predictions <- c(predictions, as.numeric(prediction))
    }
    
  }
  
  return(predictions)
  
}

####################################################################################

best_arfima_orders_btc <- list()

indices <- CV_ARIMA(btc_arfima)

for (i in seq_along(indices)) {
  
  print(i)
  
  start_train = which(rownames(btc_arfima) == indices[[i]][[1]][1])
  end_train = which(rownames(btc_arfima) == indices[[i]][[1]][length(indices[[i]][[1]])])
  train_set <- btc_arfima$Log_Returns[start_train:end_train]
  
  best_aic_model <- fit_ARFIMA(train_set)
  
  best_arfima_orders_btc[[i]] <- best_aic_model$orders
}

best_arfima_orders_btc

#####################################################################################


predictions_btc = ARFIMA_forecasts(btc_arfima)
predictions_btc


####################################################################################

predictions_btc_df <- as.data.frame(predictions_btc)
write.csv(predictions_btc_df, file = "", row.names = FALSE)


####################################################################################

fit_ARIMA <- function(ts_data) {
  model <- arfima(ts_data, drange=c(0,0.5), ic ="aic")
  
  
  return(list(
    orders = model[c("ar", "d", "ma")]
  ))
}

results_df <- as.data.frame(predictions)

matplot(results_df, type = "b", pch = 1, col = 1:ncol(results_df),
        xlab = "Index", ylab = "Values", main = "Results of the List")
legend("topright", legend = names(results_df), col = 1:ncol(results_df), pch = 1)

###################################################################################

indices_1 <- CV_ARIMA(btc_arfima, train_years=6, test_years=0)
indices_2 <- CV_ARIMA(btc_arfima, train_years=5, test_years=1)

predictions <- c()

for (i in seq_along(indices_1)) {
  
  start_train = which(rownames(btc_arfima) == indices_2[[i]][[1]][1])
  end_train = which(rownames(btc_arfima) == indices_2[[i]][[1]][length(indices[[i]][[1]])])
  
  start_test = which(rownames(btc_arfima) == indices_2[[i]][[2]][1])
  end_test = which(rownames(btc_arfima) == indices_2[[i]][[2]][length(indices[[i]][[2]])])
  
  len_train_series <- end_train - start_train + 1
  len_test_series <- end_test - start_test + 1
  
  series <- btc_arfima$Log_Returns[start_test:end_test]
  
  print(c(len_train_series, len_test_series))
  print(c(start_train, end_train, start_test, end_test))
  #print(series)
  
}

########################################################################################

ARFIMA_forecasts <- function(df) {
  
  indices_1 <- CV_ARIMA(df, train_years=6, test_years=0)
  indices_2 <- CV_ARIMA(df, train_years=5, test_years=1)
  
  predictions <- c()
  
  for (i in seq_along(indices_1)) {
    
    start_train = which(rownames(df) == indices_2[[i]][[1]][1])
    end_train = which(rownames(df) == indices_2[[i]][[1]][length(indices_2[[i]][[1]])])
    
    start_test = which(rownames(df) == indices_2[[i]][[2]][1])
    end_test = which(rownames(df) == indices_2[[i]][[2]][length(indices_2[[i]][[2]])])
    
    len_train_series <- end_train - start_train + 1
    len_test_series <- end_test - start_test + 1
    
    start_series = which(rownames(df) == indices_1[[i]][[1]][1])
    end_series = which(rownames(df) == indices_1[[i]][[1]][length(indices_1[[i]][[1]])])
    
    series <- df$Log_Returns[start_series:end_series]
    
    for (i in 1:len_test_series) {
      train_series <- window(series, start=i, end=(len_train_series + i - 1))
      model <- auto.arima(na.omit(train_series), ic="aic")
      prediction <- forecast(model, h=1)$mean
      predictions <- c(predictions, as.numeric(prediction))
    }
    
  }
  
  return(predictions)
  
}

pred_sp = ARFIMA_forecasts(btc_arfima)

length(pred_sp)

results_df <- as.data.frame(pred_sp)

write.csv(results_df, file = "", row.names = FALSE)

matplot(results_df, type = "b", pch = 1, col = 1:ncol(results_df),
        xlab = "Index", ylab = "Values", main = "Results of the List")
legend("topright", legend = names(results_df), col = 1:ncol(results_df), pch = 1)

####################################################################################

fit_ARFIMA <- function(ts_data) {
  model <- arfima(ts_data, drange=c(0,1), ic ="aic")
  
  
  return(list(
    orders = model[c("ar", "d", "ma")]
  ))
}

######################################################################################

best_arfima_orders_btc <- list()

indices <- CV_ARIMA(btc_arfima)

for (i in seq_along(indices)) {
  
  print(i)
  
  start_train = which(rownames(btc_arfima) == indices[[i]][[1]][1])
  end_train = which(rownames(btc_arfima) == indices[[i]][[1]][length(indices[[i]][[1]])])
  train_set <- btc_arfima$Log_Close[start_train:end_train]
  
  best_aic_model <- fit_ARFIMA(train_set)
  
  best_arfima_orders_btc[[i]] <- best_aic_model$orders
}

best_arfima_orders_btc

####################################################################################

ARFIMA_forecasts <- function(df) {
  
  indices_1 <- CV_ARIMA(df, train_years=6, test_years=0)
  indices_2 <- CV_ARIMA(df, train_years=5, test_years=1)
  
  predictions <- c()
  
  for (i in seq_along(indices_1)) {
    
    start_train = which(rownames(df) == indices_2[[i]][[1]][1])
    end_train = which(rownames(df) == indices_2[[i]][[1]][length(indices_2[[i]][[1]])])
    
    start_test = which(rownames(df) == indices_2[[i]][[2]][1])
    end_test = which(rownames(df) == indices_2[[i]][[2]][length(indices_2[[i]][[2]])])
    
    len_train_series <- end_train - start_train + 1
    len_test_series <- end_test - start_test + 1
    
    start_series = which(rownames(df) == indices_1[[i]][[1]][1])
    end_series = which(rownames(df) == indices_1[[i]][[1]][length(indices_1[[i]][[1]])])
    
    series <- df$Log_Close[start_series:end_series]
    
    best_aic_model <- fit_ARFIMA(df$Log_Close[start_train:end_train])
    
    best_arfima_orders <- best_aic_model$orders
    
    best_arfima_orders_ranks <- list(length(best_arfima_orders[[1]]), length(best_arfima_orders[[3]]))
    
    for (i in 1:len_test_series) {
      train_series <- window(series, start=i, end=(len_train_series + i - 1))
      model <- fracdiff(na.omit(train_series), nar = best_arfima_orders_ranks[[1]], nma = best_arfima_orders_ranks[[2]])
      prediction <- forecast(model, h=1)$mean
      predictions <- c(predictions, as.numeric(prediction))
    }
    
  }
  
  return(predictions)
  
}

####################################################################################

predictions_btc = ARFIMA_forecasts(btc_arfima)
predictions_btc

predictions_btc_log_df <- as.data.frame(predictions_btc)
write.csv(predictions_btc_log_df, file = "", row.names = FALSE)

#####################################################################################
#
#
### Arfima predictions for the first train set for hybrid models
#
#
#####################################################################################

indices <- CV_ARIMA(btc_arfima)

start_train = which(rownames(btc_arfima) == indices[[1]][[1]][1])
end_train = which(rownames(btc_arfima) == indices[[1]][[1]][length(indices[[1]][[1]])])
train_set <- btc_arfima$Log_Returns[start_train:end_train]

best_aic_model <- arfima(train_set, drange=c(0,0.5), start.p = 0, max.p = 5, start.q = 0, max.q = 5, max.order = 20, ic ="aic", seasonal = FALSE, stepwise = FALSE, allowdrift = FALSE)
best_aic_model["fitted"]
best_aic_model

train_btc_arfima_for_hybrid <- as.data.frame(best_aic_model["fitted"])
write.csv(train_btc_arfima_for_hybrid, file = "", row.names = FALSE)
