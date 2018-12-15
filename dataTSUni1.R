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
library(vars)
library(astsa)
library(readr)
library(forecast)

setwd('/Users/maggie/Documents/dataScience/practicum2/data')
list.files()

# Read the file: 
monthData <- read_csv(file='monthlyClean.csv')
head(monthData,10)
names(monthData)

# Time Series for monthly surface air temperature anomalies
anomalyTS<- ts(monthData[,3], frequency=12, start=c(1880, 1))
plot.ts(anomalyTS, main="Monthly Surface Air Temperature Anomalies ")

class(anomalyTS)
start(anomalyTS)
end(anomalyTS)
frequency(anomalyTS)

anomalyTSdec <- decompose(anomalyTS)
plot(anomalyTSdec)


# Showing outliers
# Box plot across months to explore seasonal effect
boxplot(anomalyTS~cycle(anomalyTS))

# Forecast Temperature Anomalies: 
# HoltWinters
# y default, HoltWinters() just makes forecasts for the same time period covered by 
# the original time series

## Seasonal Holt-Winters
anomalyTSmodelforecast <- HoltWinters(anomalyTS)

plot(anomalyTSmodelforecast, main="Temperture Anomaly TS against\n  HoltWintersForecast")
lines(fitted(anomalyTSmodelforecast)[,1], col = 3)

anomalyTSmodelforecast

# Summary Square of Errors
anomalyTSmodelforecast$SSE

acf(residuals(anomalyTSmodelforecast))
pacf(residuals(anomalyTSmodelforecast))

# https://stats.stackexchange.com/questions/110999/r-confused-on-residual-terminology
# # Mean squared error
# The mean squared error (MSE) is the mean of the square of the residuals:
tanomaly_mse <- mean(residuals(anomalyTSmodelforecast)^2)
tanomaly_mse

# Root mean squared error (RMSE) is then the square root of MSE:
tanomaly_rmse <- sqrt(tanomaly_mse)
tanomaly_rmse

# Residual sum of squares (RSS) is the sum of the squared residuals:
# The code below is the same as SSE: 
# tanomaly_rss <- sum(residuals(anomalyTSmodelforecast)^2)
# tanomaly_rss


anomalyTSExtended <- forecast(anomalyTSmodelforecast, h=24)
plot(anomalyTSExtended, main="Temperture Anomaly TS Extended 2 Years \n  HoltWintersForecast")
# anomalyTSExtended

plot.ts(anomalyTSExtended$residuals, main="Temperture Anomaly TS Extended 5 Years \n Holt Winters Forecast Residuals") 
# anomalyTSExtended$residuals
dim(residuals(anomalyTSExtended))
dim(residuals(anomalyTSmodelforecast))
dim(anomalyTS)


Box.test(window(anomalyTSExtended$residuals,start = c(1880,1)),lag = 1)

# Plot provided by Kevin Maher in lieu of plotForecastErrors
# http://t-redactyl.io/blog/2016/02/creating-plots-in-r-using-ggplot2-part-7-histograms.html # except args was arg in the example link above, so the link's code did not run as-is
library("ggplot2")
forecastResiduals <- as.data.frame(anomalyTSExtended$residuals)
ggplot(forecastResiduals, aes(x = x)) + geom_histogram(aes(y = ..density..), col='grey', fill='#6699ff') + stat_function(fun = dnorm, colour = "red", args = list(mean = mean(forecastResiduals$x, na.rm = TRUE), sd = sd(forecastResiduals$x, na.rm = TRUE)))

# Stationarity

# Check time series and if stationary: 
# Time series is stationary if its mean level and variance stay steady over time
# Check for stationarity
# http://www.statosphere.com.au/check-time-series-stationary-r/
# This data is stationary. 

# The Ljung-Box test examines whether there is significant evidence 
# for non-zero correlations at lags 1-20. 
# Small p-values (i.e., less than 0.05) suggest that the series is stationary.
# If P less than 0.05 then stationary: 
Box.test(anomalyTS,lag = 20, type = "Ljung-Box")
# Result p-value < 2.2e-16 = stationary
Box.test(diff(anomalyTS),lag = 20, type = "Ljung-Box")
# Result p-value < 2.2e-16 = stationary

# The Augmented Dickey–Fuller (ADF) t-statistic test: 
# small p-values suggest the data is stationary and doesn’t need to be differenced stationarity.
# If P less than 0.05 then it is stationary: 
# For the warnings: 
suppressWarnings(adf.test(anomalyTS, alternative = "stationary"))
# With differencing also
suppressWarnings(adf.test(diff(anomalyTS), alternative="stationary", k=0))


# The Kwiatkowski-Phillips-Schmidt-Shin (KPSS) test; 
# here accepting the null hypothesis means that the series is stationarity, 
# and small p-values suggest that the series is not stationary and a differencing is required.
# If P greater than 0.05 then it is stationary: 
suppressWarnings(kpss.test(anomalyTS))
# This gave a p-value = 0.01 thus not stationary, thus differencing is required: 
suppressWarnings(kpss.test(diff(anomalyTS)))
# Now we have a p-value of 0.1 = So we have to diff the dataset for it to be stationary
# Important Inferences

require(graphics)
# Autocorrelation
# Plot the autocorrelation function
# Plots lags on the horizontal and the correlations on vertical axis.
acf(diff(anomalyTS), lag.max=20)

# We can see with differencing there are have a few significant lags but these 
# die out quickly

# Partial Autocorrelation Function - PACF: 
# Gives the partial correlation of a time series with its own lagged values, 
pacf(diff(anomalyTS), lag.max=20)
anomalyTSdiff <- diff(anomalyTS)
plot(anomalyTS)
plot(anomalyTSdiff)

# Arima with 1 lag: 
anomalyTSsarima1 <- arima(anomalyTS, order=c(1,1,0)) # fit an ARIMA(0,1,1) model
tsdiag(anomalyTSsarima1)

# Arima with 3 lags: 
anomalyTSsarima3 <- arima(anomalyTS, order=c(3,1,0)) # fit an ARIMA(0,1,1) model
tsdiag(anomalyTSsarima3)

# Arima with 5 lags: 
anomalyTSsarima5 <- arima(anomalyTS, order=c(5,1,0)) # fit an ARIMA(0,1,1) model
tsdiag(anomalyTSsarima5)

# Comparison of different Arima models
BIC(anomalyTSsarima1, anomalyTSsarima3, anomalyTSsarima10)
## compare a whole set of models; BIC() would choose the smallest

# Comparison of different Arima models
AIC(anomalyTSsarima1, arima(anomalyTS, c(2,1,0)),
    arima(anomalyTS, c(2,1,1)), # <- chosen (barely) by AIC
    anomalyTSsarima3, arima(anomalyTS, c(3,1,1)))

anomalyTSsarima1

acf(residuals(anomalyTSsarima1))
pacf(residuals(anomalyTSsarima1))

# https://stats.stackexchange.com/questions/110999/r-confused-on-residual-terminology
# # Mean squared error
# The mean squared error (MSE) is the mean of the square of the residuals:
tarima_mse <- mean(residuals(anomalyTSsarima1)^2)
tarima_mse

# Root mean squared error (RMSE) is then the square root of MSE:
tarima_rmse <- sqrt(tarima_mse)
tarima_rmse

# Residual sum of squares (RSS) is the sum of the squared residuals:
# Same as SSE: 
tarima_rss <- sum(residuals(anomalyTSsarima1)^2)
tarima_rss



