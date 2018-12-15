# Practicum #2 - Maggie Sleziak
# Regression Excercises

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

# Some code of adapted from the tutorial: 
# R PROGRAMMING FOR CLIMATE DATA ANALYSIS AND VISUALIZATION
# Samuel S.P. Shen
# San Diego State University and Scripps Institution of Oceanography
# http://scrippsscholars.ucsd.edu/s4shen/files/r-textbysamshenjune2017.pdf


rm(list = ls(all=TRUE))
library(ggplot2)
library(tseries) 
library(pscl)


setwd('/Users/maggie/Documents/dataScience/practicum2/data')
list.files()

# Read the file: 
monthData <- read_csv(file='monthlyClean.csv')
head(monthData,10)
names(monthData)

# Time Series for monthly surface air temperature anomalies
monthOnlyTS<- ts(monthData[,3], frequency=12, start=c(1880, 1))
monthOnlyTS
str(monthOnlyTS)

##
#linear model
# timemo=seq(1880,2018,length=1665)
# Can also be done with timemo instead of time(monthOnlyTS), but decided to use time(monthOnlyTS) 
# str(timemo)

plot.new()
str(time(monthOnlyTS))
plot(time(monthOnlyTS),monthOnlyTS,type="l", cex.lab=1.4,
     xlab="Year", ylab="Temperature Anomaly [oC]",
     main="Global Avg of Monthly Surface Air Temp Anomalies: Jan 1880 - Sept 2018")
abline(lm(monthOnlyTS ~ time(monthOnlyTS)),col="blue",lwd=2)
text(1930,0.7, paste("Linear Trend per Century [oC]: = ",
                     round(digits=3,
                           (100*coefficients(lm(monthOnlyTS ~ time(monthOnlyTS)))[2]))), 
     col="red")



#Linear model 
monthOnlyTS.lm = lm(monthOnlyTS~time(monthOnlyTS))
summary(monthOnlyTS.lm) #r square 0.6832 not so good

#time series plot of the standardized residuals
plot(y=rstandard(monthOnlyTS.lm), x=as.vector(time(monthOnlyTS)), type = 'l')

# All plots:  
# plot(monthOnlyTS.lm, main="Temperture Anomaly TS Linear Regression")

acf(residuals(monthOnlyTS.lm))
pacf(residuals(monthOnlyTS.lm))

lmR_mse <- mean(residuals(monthOnlyTS.lm)^2)
lmR_mse

# Root mean squared error (RMSE) is then the square root of MSE:
lmR_rmse <- sqrt(lmR_mse)
lmR_rmse

plot.new()
plot(time(monthOnlyTS),monthOnlyTS,type="s",
     cex.lab=1.4, lwd=1,
     xlab="Year", ylab="Temperature anomaly [oC]",
     main="Monthly SAT Time Series - Linear Regression Line Fit")
abline(lm(monthOnlyTS ~ time(monthOnlyTS)),col="blue",lwd=2)
text(1940,0.7, paste("RMSE: = ",
                     round(digits=3,
                           (lmR_rmse))), 
     col="red")




# Change to anual so we can do orthogonal polynomial regression: 
monthlyMatrix = matrix(monthOnlyTS[1:1656], ncol=12, byrow=TRUE)
#compute annual average
yearData=rowMeans(monthlyMatrix)
#Plot the annual mean global average temp
timeyr<-seq(1880, 2017)
str(timeyr)
str(yearData)


#Quadratic model trend

yearlyOnlyTS.qm = lm(yearData ~ timeyr + I(timeyr^2))
summary(yearlyOnlyTS.qm) # r square: 0.8717 - better

#time series plot of the standardized residuals
plot(y=rstandard(yearlyOnlyTS.qm), x=as.vector(timeyr), type = 'l')

# All plots:  
# plot(yearlyOnlyTS.qm, main="Temperture Anomaly TS Quadratic Regression")

acf(residuals(yearlyOnlyTS.qm))
pacf(residuals(yearlyOnlyTS.qm))

qmR_mse <- mean(residuals(yearlyOnlyTS.qm)^2)
qmR_mse

# Root mean squared error (RMSE) is then the square root of MSE:
qmR_rmse <- sqrt(qmR_mse)
qmR_rmse

plot.new()
plot(timeyr,yearData,type="s",
     cex.lab=1.4, lwd=1,
     xlab="Year", ylab="Temperature anomaly [oC]",
     main="Yearly SAT Time Series - Quadratic Regression Line Fit")
lines(timeyr,predict(yearlyOnlyTS.qm),col="blue", lwd=2)
text(1940,0.7, paste("RMSE: = ",
                     round(digits=3,
                           (qmR_rmse))), 
     col="red")
  


#polynomial model trend

yearlyOnlyPol.qm = lm(yearData ~ poly(timeyr,9, raw= TRUE))
summary(yearlyOnlyPol.qm) # r square: 0.8717 - better

#time series plot of the standardized residuals
plot(y=rstandard(yearlyOnlyPol.qm), x=as.vector(timeyr), type = 'l')

# All plots:  
# plot(yearlyOnlyPol.qm, main="Temperture Anomaly TS Polynomial Regression")

acf(residuals(yearlyOnlyPol.qm))
pacf(residuals(yearlyOnlyPol.qm))

polR_mse <- mean(residuals(yearlyOnlyPol.qm)^2)
polR_mse

# Root mean squared error (RMSE) is then the square root of MSE:
polR_rmse <- sqrt(polR_mse)
polR_rmse

plot.new()
plot(timeyr,yearData,type="s",
     cex.lab=1.4, lwd=1,
     xlab="Year", ylab="Temperature anomaly [oC]",
     main="Yearly SAT Time Series - Polynomial Orthogonal Regression Line Fit")
lines(timeyr,predict(yearlyOnlyPol.qm),col="blue", lwd=2)
text(1940,0.7, paste("RMSE: = ",
                     round(digits=3,
                           (polR_rmse))), 
     col="red")



yearlyOnlyPolOrtho.qm = lm(yearData ~ poly(timeyr,9, raw= FALSE))
summary(yearlyOnlyPolOrtho.qm) # r square: 0.8717 - better

#time series plot of the standardized residuals
plot(y=rstandard(yearlyOnlyPolOrtho.qm), x=as.vector(timeyr), type = 'l')

# All plots:  
# plot(yearlyOnlyPolOrtho.qm, main="Temperture Anomaly TS Polynomial Regression")

acf(residuals(yearlyOnlyPolOrtho.qm))
pacf(residuals(yearlyOnlyPolOrtho.qm))

polQR_mse <- mean(residuals(yearlyOnlyPolOrtho.qm)^2)
polQR_mse

# Root mean squared error (RMSE) is then the square root of MSE:
polQR_rmse <- sqrt(polQR_mse)
polQR_rmse

plot.new()
plot(timeyr,yearData,type="s",
     cex.lab=1.4, lwd=1,
     xlab="Year", ylab="Temperature anomaly [oC]",
     main="Yearly SAT Time Series - Polynomial Non-Orthogonal Regression Line Fit")
lines(timeyr,predict(yearlyOnlyPolOrtho.qm),col="blue", lwd=2)
text(1940,0.7, paste("RMSE: = ",
                     round(digits=3,
                           (polQR_rmse))), 
     col="red")

