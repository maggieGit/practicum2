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
monthData <- read_csv(file='monthlyClean3.csv')
head(monthData,10)
names(monthData)
monthData <- dplyr::rename(monthData, ds = 'YearMonthDay', y = 'TempAnomaly')
monthData$yorig = monthData$y
names(monthData)
head(monthData,10)

m <- prophet(monthData)
future <- make_future_dataframe(m, periods = 24, freq = 'month')
tail(future)
fcst <- predict(m, future)
fcstOrig <- fcst
tail(fcst[,c("ds", "yhat", 'yhat_lower', 'yhat_upper')])
plot(m, fcst, cex = .2)

# Trend and yearly seasonality
prophet_plot_components(m, fcst, uncertainty = TRUE, plot_cap = TRUE,
                        weekly_start = 0, yearly_start = 0, render_plot = TRUE)

dyplot.prophet(m, fcst, uncertainty = TRUE, plot_cap = TRUE,
               weekly_start = 0, yearly_start = 0, render_plot = TRUE)


# Try: Out-of-sample forecast
#df.cv <- cross_validation(m, horizon = (365/12), units = "days", period = (362/12),
#                               initial = 132*12*(365/12))

df.cv <- cross_validation(m, horizon = 2*(365/12), units = "days", period = 2*(362/12),
                          initial = 132*12*(365/12))
head(df.cv)
tail(df.cv)
df.p <- performance_metrics(df.cv)
head(df.p)
tail(df.p)

# Mean Absolute Percentage Error
plot_cross_validation_metric(df.cv, metric = 'mape')

# Cross validation performance metrics can be visualized with plot_cross_validation_metric, 
# here shown for MAPE. Dots show the absolute percent error for each prediction in df_cv. 
# The blue line shows the MAPE, where the mean is taken over a rolling window of the dots. 
# We see for this forecast that errors around 5% are typical for predictions one month into 
# the future, and that errors increase up to around 11% for predictions that are a year out.

plot_cross_validation_metric(df.cv, metric = 'rmse')

