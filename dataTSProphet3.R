# Practicum #2 - Maggie Sleziak
# #1: Read data 

# Data: NCEI spatial average time series of monthly data
# Read the file:
# From: https://www.ncdc.noaa.gov/data-access/marineocean-data/noaa-global-surface-temperature-noaaglobaltemp
# Ascii Time Series: ftp://ftp.ncdc.noaa.gov/pub/data/noaaglobaltemp/operational/
# Monthly and annual land–ocean temperature time series 
# are available from 1880 to present for several zonal bands
# file name convention for areal average (aravg) time series:
# ann=annual average
# mon=monthly average
# land_ocean=merged land-ocean surface temperature
# land=land surface temperature for air
# ocean=ocean surface temperature for water
# latitudes=southern and northern limits of areal average
# v=version number
# yyyymm=date for the latest data

# Monthly data (aravg.mon.*) :
#  1st column = year
#  2nd column = month
#  3rd column = anomaly of temperature (K)
#  4th column = total error variance (K**2)
#  5th column = high-frequency error variance (K**2)
#  6th column = low-frequency error variance (K**2)
#  7th column = bias error variance (K**2)
#  8th column = diagnostic variable
#  9th column = diagnostic variable
#  10th column= diagnostic variable


rm(list = ls(all=TRUE))
library(ggplot2)
library(dplyr)
library(e1071)
library(lattice)
library(corrplot)
library(tseries) 
library(prophet)
library(vars)
library(astsa)
library(readr)
library(forecast)

setwd('/Users/maggie/Documents/dataScience/practicum2/data')
list.files()

# Read the file: 
monthDataRaw <- read_csv(file='monthlyClean3.csv')
head(monthDataRaw,10)
names(monthDataRaw)
monthDataRaw <- dplyr::rename(monthDataRaw, ds = 'YearMonthDay', y = 'TempAnomaly')
dataTS <- ts(monthDataRaw[,2], frequency=12, start=c(1880, 1))
plot(dataTS)

# In case we want to test with differenced data
# dataTsdiff <- diff(dataTS)
# plot(dataTsdiff)
# monthData  <- data.frame(ds=as.Date(as.yearmon(time(dataTsdiff))), y=as.matrix(dataTsdiff))

monthData  <- monthDataRaw
names(monthData)
head(monthData,10)
length(monthData$y)

# Save 80% for testing
monthData.end <- floor(0.80*length(monthData$y)) 
monthData.train <- monthData[1:monthData.end,] 
monthData.test <- monthData[(monthData.end+1):length(monthData$y),] 

monthData.end
length(monthData$y)
length(monthData.train$y)
length(monthData.test$y)

str(monthData.train)
head(monthData.train,10)
tail(monthData.train,10)

head(monthData.test,10)
tail(monthData.test,10)
str(monthData.test)

c(min(monthData.train$ds), max(monthData.train$ds))
c(min(monthData.test$ds), max(monthData.test$ds))

ntrain <- nrow(monthData.train)
ntrain
ntest <- nrow(monthData.test)
ntest

#This will be how far out we will need our training model to predict in order to compare with our observed values for evaluating the model accuracy.
forecast.horizon <- length(monthData.test$y) 
forecast.horizon

# RMSE
rmse <- function(actual,predicted)
          round((sum((actual-predicted)^2)/length(actual))^.5,2)

# Porphet train and predict
m <- prophet(monthData.train)
# Future includes this 20% (ntest+1) months
future <- make_future_dataframe(m, periods = ntest+1, freq = 'month', include_history = TRUE) 
head(future)
tail(future)

forecast <- predict(m, future)
tail(forecast[,c("ds", "yhat", 'yhat_lower', 'yhat_upper')])

monthData.train$yhat <- forecast$yhat[1:ntrain]
monthData.test$yhat <- forecast$yhat[(ntrain+1):(ntrain+ntest)]

# Trend and yearly seasonality
prophet_plot_components(m, forecast, uncertainty = TRUE, plot_cap = TRUE,
                        weekly_start = 0, yearly_start = 0, render_plot = TRUE)

dyplot.prophet(m, forecast, uncertainty = TRUE, plot_cap = TRUE,
               weekly_start = 0, yearly_start = 0, render_plot = TRUE)


# Plot train, plot test and plot fcst
monthTrainTS<- ts(monthData.train[,2], frequency=12, start=c(1880, 1))
plot.ts(monthTrainTS, main="Monthly Surface Air Temperature Anomalies ")
tail(monthData.test)

monthTestTS<- ts(monthData.test[,2], frequency=12, start=c(1991, 1))
plot.ts(monthTestTS, main="Monthly Surface Air Temperature Anomalies ")

monthlyProfForecast <- forecast[,c("ds", "yhat")]
tail(monthlyProfForecast)
str(monthlyProfForecast)

monthlyProfForecastEnd <- monthlyProfForecast[(monthData.end+1):length(monthlyProfForecast$y),] #assign the most recent 10% to the test set
head(monthlyProfForecastEnd)
str(monthlyProfForecastEnd)

monthForecastAllTS<- ts(monthlyProfForecast[,2], frequency=12, start=c(1880, 1))
monthForecastTS<- ts(monthlyProfForecastEnd[,2], frequency=12, start=c(1991, 1))
monthTestTSyhat<- ts(monthData.test[,3], frequency=12, start=c(1991, 1))
plot.ts(monthForecastTS, main="Monthly Surface Air Temperature Anomalies ")

plot(monthTrainTS, main="Monthly Surface Air Temperature Anomalies ")
lines(monthTestTS, col = "red")
lines(monthForecastAllTS, col = "green")
lines(monthForecastTS, col = "blue")
lines(monthTestTSyhat, col = "brown")

monthTestTS # This is the y for test data
monthForecastTS # This is yhat for test data
accuracy(monthForecastTS,monthTestTS)
accuracy(monthData.test$y,monthData.test$yhat)

# For train y and train yhat: 
accuracy(monthData.train$y,monthData.train$yhat)
c(rmse(monthData.train$y,monthData.train$yhat), rmse(monthData.test$y,monthData.test$yhat))

str(monthData.test)

library(scorer)
r2_score(monthData.test$y, monthData.test$yhat)
r2_score(monthData.train$y, monthData.train$yhat)
mean_squared_error(monthData.test$yhat, monthData.test$y)
sqrt(mean_squared_error(monthData.test$yhat, monthData.test$y))
sqrt(mean_squared_error(monthData.train$yhat, monthData.train$y))
