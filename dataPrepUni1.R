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
library(dplyr)
library(mice)
library(lubridate)
library(ggplot2)
library(readr)
library(stringr)
setwd('/Users/maggie/Documents/dataScience/practicum2/data')
list.files()

monthTempAve=read.table("aravg.mon.land_ocean.90S.90N.v4.0.1.201809.asc")
class(monthTempAve)
str(monthTempAve)
summary(monthTempAve)
head(monthTempAve)
head(monthTempAve[,1]) # Year
head(monthTempAve[,2]) # Month
head(monthTempAve[,3]) # Temperature anomalies
dim(monthTempAve)
monthTempAve <- monthTempAve[, c('V1', 'V2', 'V3')]
head(monthTempAve)
monthTempAve <- dplyr::rename(monthTempAve, Year = 'V1', Month = 'V2', TempAnomaly = 'V3')
monthTempAve$MonthLong <- str_pad(monthTempAve$Month, 2, pad = "0")
head(monthTempAve,10)
monthTempAve$YearMonth <- paste(monthTempAve$Year, monthTempAve$Month, sep="-" ) 
monthTempAve$YearMonthDay <- paste(monthTempAve$Year, monthTempAve$MonthLong, "01", sep="-" ) 
head(monthTempAve,10)

# Check for NA values: No NA values in data
apply(monthTempAve, 2, function(x) any(is.na(x)))

monthTempAve2<- monthTempAve[, c('YearMonth', 'TempAnomaly')]
head(monthTempAve2,10)

monthTempAve3<- monthTempAve[, c('YearMonthDay', 'TempAnomaly')]
head(monthTempAve3,10)

# Write to file for data prep for later use
write_csv(monthTempAve, 'monthlyClean.csv')
write_csv(monthTempAve2, 'monthlyClean2.csv')
write_csv(monthTempAve3, 'monthlyClean3.csv')



monthTempAve2$TempAnomaly1 <- c(NA, monthTempAve2$TempAnomaly[1:{nrow(monthTempAve2) - 1}])
monthTempAve2$TempAnomaly2 <- c(NA, NA, monthTempAve2$TempAnomaly[1:{nrow(monthTempAve2) - 2}])
monthTempAve2$TempAnomaly3 <- c(NA, NA, NA, monthTempAve2$TempAnomaly[1:{nrow(monthTempAve2) - 3}])
head(monthTempAve2,10)


# Lag one day: 
monthTempLag1 <- monthTempAve2[, c('TempAnomaly1', 'TempAnomaly' )]
head(monthTempLag1,10)
# Remove the 1st row of NAs: 
monthTempLag1<- monthTempLag1[-c(1), ]
head(monthTempLag1,10)

# Lag 3 days: 
monthTempLag3 <- monthTempAve2[, c('TempAnomaly3', 'TempAnomaly2', 'TempAnomaly1', 'TempAnomaly' )]
head(monthTempLag3,10)
# Remove the 1st row of NAs: 
monthTempLag3<- monthTempLag3[-c(1,2,3), ]
head(monthTempLag3,10)

# Write data to files
write_csv(monthTempLag1, 'monthTempLag1.csv')
write_csv(monthTempLag3, 'monthTempLag3.csv')

par(mar=c(4,4,3,1))
x<-timemo # These are the mmonths
y<-monthTempAve[,3] # This is the temperature anomalies

# Y is Temperature Anomalies
n1<-which(y>=0) # Assign where average temp is greater than zero
x1<-x[n1] # Years for temp greater than zero
y1<-y[n1] # Temps greater than zero

n2<-which(y<0) # vector elements where temps less than zerp
x2<-x[n2] # years less than zero temps
y2<-y[n2] # less than zero temps

plot(x1,y1,type="h",xlim=c(1880,2018),lwd=3,
     tck=0.02, ylim=c(-1.0,1.0), #tck>0 makes ticks inside the plot
     ylab="Temperature [deg C]",
     xlab="Year",col="red",
     main="NOAA Global Annual Temperature Anomalies")
lines(x2,y2,type="h",
      lwd=3, tck=-0.02,  col="blue")

hist(monthTempAve$TempAnomaly, freq=NULL, density=NULL, breaks=20, prob = T, xlim = c(-1 , 1), col = "orange", main = "Histogram: Surface Air Temperature Montly Anomalies ", xlab = "Celcius")


# Looking at data patterns throughout the years: 
scatter <- ggplot(data=monthTempAve, aes(x = Year, y = TempAnomaly))
scatter +  geom_point(aes(color=TempAnomaly)) + xlab("Year") +  ylab("%") + ggtitle("Yearly Patterns: SAT Anomalies")


# Looking at data patterns throughout the months: 
scatter <- ggplot(data=monthTempAve, aes(x = Month, y = TempAnomaly))
scatter +  geom_point(aes(color=TempAnomaly)) + xlab("Month") +  ylab("%") + ggtitle("Monthly Patterns: SAT Anomalies")

