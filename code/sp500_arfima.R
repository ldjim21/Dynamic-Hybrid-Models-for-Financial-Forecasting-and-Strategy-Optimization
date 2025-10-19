library(forecast)
library(stats)
library(fracdiff)
library(lubridate)
library(dplyr)

sp500_arfima <- read.csv("")

sp500_arfima$Date <- as.Date(sp500_arfima$Date)

rownames(sp500_arfima) <- sp500_arfima$Date

####################################################################################

CV_ARIMA <- function(df, train_years = 5, test_years = 1) {
  
  start_date <- as.numeric(format(min(df$Date), "%Y"))
  end_date <- as.numeric(format(max(df$Date), "%Y"))
  
  current_start_date <- start_date
  indices <- list()
  
  while (TRUE) {
    
    train_end_date <- current_start_date + train_years - 1
    test_end_date <- train_end_date + test_years
    
    if (test_end_date > end_date) {
      break
    }
    
    train_set <- df[(df$Date >= as.Date(paste0(current_start_date, "-01-01")) & df$Date <= as.Date(paste0(train_end_date, "-12-31"))), "Date"]
    test_set <- df[(df$Date >= as.Date(paste0(train_end_date + 1, "-01-01")) & df$Date <= as.Date(paste0(test_end_date, "-12-31"))), "Date"]
    
    indices[[length(indices) + 1]] <- list(train_set, test_set)
    
    current_start_date <- current_start_date + 1
  }
  
  return(indices)
}

####################################################################################

indices <- CV_ARIMA(sp500_arfima)

for (i in seq_along(indices)) {
  cat(sprintf("Window %d:\n", i))
  cat("Training set:\n")
  start_train = which(rownames(sp500_arfima) == indices[[i]][[1]][1])
  end_train = which(rownames(sp500_arfima) == indices[[i]][[1]][length(indices[[i]][[1]])])
  print(head(sp500_arfima$Date[start_train:end_train], 3))
  print(tail(sp500_arfima$Date[end_train:end_train], 3))
  cat("Testing set:\n")
  start_test = which(rownames(sp500_arfima) == indices[[i]][[2]][1])
  end_test = which(rownames(sp500_arfima) == indices[[i]][[2]][length(indices[[i]][[2]])])
  print(head(sp500_arfima$Date[start_test:end_test], 3))
  print(tail(sp500_arfima$Date[start_test:end_test], 3))
  cat("\n")
}

####################################################################################

fit_ARFIMA <- function(ts_data) {
  model <- arfima(ts_data, drange=c(0,0.5), ic ="aic")
  
  
  return(list(
    orders = model[c("ar", "d", "ma")]
  ))
}

####################################################################################

best_arfima_orders_sp500 <- list()

indices <- CV_ARIMA(sp500_arfima)

for (i in seq_along(indices)) {
  
  print(i)
  
  start_train = which(rownames(sp500_arfima) == indices[[i]][[1]][1])
  end_train = which(rownames(sp500_arfima) == indices[[i]][[1]][length(indices[[i]][[1]])])
  train_set <- sp500_arfima$Log_Returns[start_train:end_train]
  
  best_aic_model <- fit_ARFIMA(train_set)
  
  best_arfima_orders_sp500[[i]] <- best_aic_model$orders
}

best_arfima_orders_sp500

####################################################################################

best_arfima_orders_sp500_ranks <- list()

for (i in seq_along(best_arfima_orders_sp500)) {
  
  best_arfima_orders_sp500_ranks[[i]] <- list(length(best_arfima_orders_sp500[[i]][[1]]), length(best_arfima_orders_sp500[[i]][[3]]))

}

best_arfima_orders_sp500_ranks

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

best_arfima_orders_sp500 <- list()

indices <- CV_ARIMA(sp500_arfima)

for (i in seq_along(indices)) {
  
  print(i)
  
  start_train = which(rownames(sp500_arfima) == indices[[i]][[1]][1])
  end_train = which(rownames(sp500_arfima) == indices[[i]][[1]][length(indices[[i]][[1]])])
  train_set <- sp500_arfima$Log_Returns[start_train:end_train]
  
  best_aic_model <- fit_ARFIMA(train_set)
  
  best_arfima_orders_sp500[[i]] <- best_aic_model$orders
}

best_arfima_orders_sp500

####################################################################################

predictions_sp500 = ARFIMA_forecasts(sp500_arfima)
predictions_sp500

####################################################################################

predictions_sp500_df <- as.data.frame(predictions_sp500)
write.csv(predictions_sp500_df, file = "", row.names = FALSE)

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

indices_1 <- CV_ARIMA(sp500_arfima, train_years=6, test_years=0)
indices_2 <- CV_ARIMA(sp500_arfima, train_years=5, test_years=1)

predictions <- c()

for (i in seq_along(indices_1)) {
  
  start_train = which(rownames(sp500_arfima) == indices_2[[i]][[1]][1])
  end_train = which(rownames(sp500_arfima) == indices_2[[i]][[1]][length(indices[[i]][[1]])])
  
  start_test = which(rownames(sp500_arfima) == indices_2[[i]][[2]][1])
  end_test = which(rownames(sp500_arfima) == indices_2[[i]][[2]][length(indices[[i]][[2]])])
  
  len_train_series <- end_train - start_train + 1
  len_test_series <- end_test - start_test + 1
  
  series <- sp500_arfima$Log_Returns[start_test:end_test]
  
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

pred_sp = ARFIMA_forecasts(sp500_arfima)

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

best_arfima_orders_sp500 <- list()

indices <- CV_ARIMA(sp500_arfima)

for (i in seq_along(indices)) {
  
  print(i)
  
  start_train = which(rownames(sp500_arfima) == indices[[i]][[1]][1])
  end_train = which(rownames(sp500_arfima) == indices[[i]][[1]][length(indices[[i]][[1]])])
  train_set <- sp500_arfima$Log_Close[start_train:end_train]
  
  best_aic_model <- fit_ARFIMA(train_set)
  
  best_arfima_orders_sp500[[i]] <- best_aic_model$orders
}

best_arfima_orders_sp500

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

predictions_sp500 = ARFIMA_forecasts(sp500_arfima)
predictions_sp500

predictions_sp500_log_df <- as.data.frame(predictions_sp500)
write.csv(predictions_sp500_log_df, file = "", row.names = FALSE)

#####################################################################################
#
#
### Arfima predictions for the first train set for hybrid models
#
#
#####################################################################################

indices <- CV_ARIMA(sp500_arfima)

start_train = which(rownames(sp500_arfima) == indices[[1]][[1]][1])
end_train = which(rownames(sp500_arfima) == indices[[1]][[1]][length(indices[[1]][[1]])])
train_set <- sp500_arfima$Log_Returns[start_train:end_train]
  
best_aic_model <- arfima(train_set, drange=c(0,0.5), ic ="aic")
best_aic_model["fitted"]

train_sp500_arfima_for_hybrid <- as.data.frame(best_aic_model["fitted"])
write.csv(train_sp500_arfima_for_hybrid, file = "", row.names = FALSE)
