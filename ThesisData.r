## Thesis Data management -------------------------------------------------
##                                                                       ##
## File name: ThesisData.R                                               ##  
## Author:    Mathias Concha Modin                                       ##
## Email:     mathias.modin@ne.su.se                                     ##
##                                                                       ##




# Definitions -------------------------------------------------------------

# The following abbreviations  are used in "rawdata.csv":
# DS == TRUE implies that the series have been deseasonalised pre downlaod    

# "L"     = Aggregate hours worked (in 10^6 SEK) => DS == TRUE 
# "nC"    = Nominal Consumption (in 10^6 SEK) => DS == TRUE
# "rC"    = Real Consumption in 2014 prices (in 10^6 SEK) => DS == TRUE
# "frI"   = Fixed Real Investments in 2014 prices (in 10^6 SEK) => DS == TRUE
# "nY"    = Nominal GDP (in 10^6 SEK) => DS == FALSE
# "rY"    = Real GDP in 2014 prices (in 10^6 SEK) => DS == FALSE
# "rYds"  = Real GDP in 2014 prices (in 10^6 SEK) => DS == TRUE
# "repo"  = Riksbank Repo Rate in percent
# "pop"   = Population age 16-64 in 1000s
# "wc"    = Firms' nominal wage costs (in 10^6 SEK)
# "sf"    = Firms' nominal social fees (in 10^6 SEK)




# Importing Data ----------------------------------------------------------

# Clearing Workspace
rm(list=ls())

#setwd("Z:/NEK MASTEr/Thesis/Data")
data <- ts(read.csv2("rawdata.csv", header = TRUE)[,-1], start = c(1980, 1), freq = 4)
# cpi <- ts(read.csv2("CPI.csv", header=F), start=1980, freq=12 )
# Subsetting the dataset to start in 1994Q2
data <- window(data, start = c(1994,1))




# Loading packages --------------------------------------------------------

library(seasonal) # arima x13 seasonal estimation
library(dynlm)    # used for OLS estimation of trend components
library(mFilter)  # used for filtering trend components
library(urca)     # used for unit root testing




# Data Manipulations ------------------------------------------------------

# "rawdata.csv" contains the combined raw data downloaded directly from
# Statistics Sweden and the Riksbank. All transformations of data that are 
# later used in the estimated model is documented in this script.

# Storing variables as separate objects instead
Y <- data[,"rYds"]/(data[,"pop"]*1000)
C <- data[,"rC"]/(data[,"pop"]*1000)
I <- data[,"frI"]/(data[,"pop"]*1000)
AL <- data[,"L"]
R <- data[,"repo"]
pop <- data[,"pop"]

# Creating the GDP deflator from nominal and real GDP both of which are not
# deseasonalised prior to construction
Ydef <-  data[,"nY"]/data[,"rY"]*100

# Creating the proxy for the hourly wage
# nominal wage bill
nw <- data[,"wc"] + data[,"sf"] #sum of wage costs + social fees
W <- (nw*100)/(AL*Ydef) # real wage


# Aggregating CPI
cpi<- aggregate(cpi, FUN=mean, nfrequency=4)


# Subsetting and normalising reference to 2014
dtime <- time(cpi)
cpi <- window(cpi/cpi[which(dtime==2014)]*100, start=c(1994, 1), end=c(2015,4))




# Plotting raw model variables --------------------------------------------

# Storing Legend Names
var_names <- c("GDP", "Consumption", "Investments", "Labour", 
             "Repo Rate", "Hourly Wage Bill", "GDP-deflator")
# Storing unit names
unit <- c("M SEK", "M SEK", "M SEK", "M Hours", "Percent", "SEK", "2014Q4 = 100")

# Creating a color palette
red <- "#E53935"
purple <- "#5E35B1"
teal <- "#00897B"
indigo <- "#3949AB"
orange <- "#FB8C00"
bgrey <- "#546E7A"


# Creating a matrix of only relevant variables
new_data <- cbind(Y, C, I, AL, R , W, Ydef)




# Seasonal adjustments ----------------------------------------------------

# Because nominal GDP is not available in de-seasonalsied form at SCB, 
# the GDP deflator has to be deseasonalsied. Equivalently the firms' wage
# also has a seasonal pattern that has to be removed.

# Estimation of ARIMA x13 models
w_season <- seas(W, transform.function = "none")
ydef_season <- seas(Ydef, transform.function = "none")

# Extracting the seasonally adjusted series
ds.W <- w_season$data[,"seasonaladj"]
ds.ydef <- ydef_season$data[,"seasonaladj"]

# Plotting results
plot(w_season)
plot(ydef_season)




# Trend Adjustments -------------------------------------------------------
# Real variables are detrended using a linear trend
Y.model <- dynlm(log(Y) ~ trend(Y))
C.model <- dynlm(log(C) ~ trend(C))
I.model <- dynlm(log(I) ~ trend(I))
W.model <- dynlm(log(ds.W) ~ trend(ds.W))
R.model <- dynlm(R~ trend(R)) # not estimated using logs


# Saving trends for plotting
t.Y <- Y.model$fitted.values
t.C <- C.model$fitted.values
t.I <- I.model$fitted.values
t.W <- W.model$fitted.values
t.R <- R.model$fitted.values



# Detrended data / residuals
dt.Y <- Y.model$residuals
dt.C <- C.model$residuals
dt.I <- I.model$residuals
dt.W <- W.model$residuals
dt.R <- R.model$residuals

# Removing population growth from growth in hours worked
L.model <- dynlm(log(AL) ~ 1 + log(pop))
dt.L <-L.model$residuals 
t.L <- L.model$fitted.values




# Log linearised steady state data ----------------------------------------

# The log of the detrended series are used to proxy the log-deviation from
# steady state in the theoretical model. The deviations are also demeaned.

lnY <- dt.Y-mean(dt.Y)
lnC <- dt.C-mean(dt.C)
lnI <- dt.I-mean(dt.I)
lnW <- dt.W-mean(dt.W)
lnRR <- log(1 + (R+1)/400) - log(1+(t.R)/400)
lnR <- lnRR - mean(lnRR)
lnL <- dt.L-mean(dt.L)
PIE <- (log(ds.ydef)-log(lag(ds.ydef,-1)))
PIE <- PIE-mean(PIE)
inf=(log(cpi)-log(lag(cpi,-1)))
inf=inf-mean(inf)




# Creating final data -----------------------------------------------------

# log lin data used for model estimation
data <- cbind(lnY, lnC, lnI, lnL, lnR, lnW, PIE)
data <- window(data, start=c(1994, 2)) # subsetting again to get rid of first NA

datainf <- cbind(lnY, lnC, lnI, lnL, lnR, lnW, inf)
datainf <- window(datainf, start=c(1994, 2)) 

new_data <- cbind(log(Y), log(C), log(I), dt.L, R , log(ds.W), log(ds.ydef), log(cpi))
# trends saved for plotting
#t.data <- cbind(t.Y, t.C, t.I, t.W, t.R, t.L)
#t.data <- window(data, start=c(1994, 2))

# raw data saved for plotting
#new_data <- log(cbind(Y, C, I, AL, exp(R), ds.W))
#new_data <- window(new_data, start=c(1994, 2))

# naming columns
datanames <- c("Y", "C", "I", "L", "R", "W", "PIE")
colnames(data) <- datanames
colnames(datainf) <- datanames
colnames(new_data) <- c(datanames, "PIEinf")
# saving data
if (getwd() == "Z:/NEK MASTEr/Thesis/Data"){
  write.csv2(data, "data.csv")
  write.csv2(datainf, "datainf.csv")
  write.csv2(new_data, "ndata.csv")   
}




# Plotting final data -----------------------------------------------------


# GDP plots
par(mfrow=c(1,2), font.lab=2, font=2)
ts.plot(log(Y), t.Y, col = c(teal, orange), lwd = c(2,2))
legend("topleft", c("Log GDP", "Trend"), col = c(teal, orange), lwd=c(2,2), bty="n")
axis(1, font=2)
axis(2, font=2)

ts.plot(lnY, col = teal, lwd=2, ylab="")
legend("topleft", "Y", col = teal, lwd=2, bty="n")
abline(h=0, col=orange, lwd=2)
axis(1, font=2)
axis(2, font=2)

# Consumption
par(mfrow=c(1,2), font.lab=2, font=2)
ts.plot(log(C), t.C, col = c(teal, orange), lwd = c(2,2), ylab="")
legend("topleft", c("Log Consumption   ", "Trend"), col = c(teal, orange), lwd=c(2,2), bty="n")
axis(1, font=2)
axis(2, font=2)

ts.plot(lnC, col = teal, lwd=2, ylab="")
legend("topleft", "C", col = teal, lwd=2, bty="n")
abline(h=0, col=orange, lwd=2)
axis(1, font=2)
axis(2, font=2)

# Investments
ts.plot(log(I), t.I, col = c(teal, orange), lwd = c(2,2), ylab="")
legend("topleft", c("Log Investments   ", "Trend"), col = c(teal, orange), lwd=c(2,2), bty="n")
axis(1, font=2)
axis(2, font=2)

ts.plot(lnI, col = teal, lwd=2, ylab="")
legend("topleft", "I", col = teal, lwd=2, bty="n")
abline(h=0, col=orange, lwd=2)
axis(1, font=2)
axis(2, font=2)


# Wages
ts.plot(log(ds.W), t.W, col = c(teal, orange), lwd = c(2,2), ylab="")
legend("topleft", c("Log Hourly Wage", "Trend"), col = c(teal, orange), lwd=c(2,2), bty="n")
axis(1, font=2)
axis(2, font=2)

ts.plot(lnW, col = teal, lwd=2, ylab="")
legend("topleft", "W", col = teal, lwd=2, bty="n")
abline(h=0, col=orange, lwd=2)
axis(1, font=2)
axis(2, font=2)


# Repo rate
ts.plot(R, t.R, col = c(teal, orange), lwd = c(2,2), ylab="")
legend("topright", c("Repo Rate        ", "Trend   "), col = c(teal, orange), lwd=c(2,2), bty="n")
axis(1, font=2)
axis(2, font=2)

ts.plot(lnR, col = teal, lwd=2, ylab="")
legend("topright", "R", col = teal, lwd=2, bty="n")
abline(h=0, col=orange, lwd=2)
axis(1, font=2)
axis(2, font=2)

# GDP deflator/Inflation
ts.plot(Ydef, ds.ydef, col = c(teal, orange), lwd = c(2,2), ylab="")
legend("topleft", c("GDP deflator", "Season adj.\nGDP deflator"), col = c(teal, orange), lwd=c(2,2), bty="n")
axis(1, font=2)
axis(2, font=2)

ts.plot(PIE, col = teal, lwd=2, ylab="")
legend("bottomright", "GDP inflation            ", col = teal, lwd=2, bty="n")
abline(h=0, col=orange, lwd=2)
axis(1, font=2)
axis(2, font=2)

# Aggregate Labour Supply
par(mfrow=c(1,2), font.lab=2, font=2)
ts.plot(log(AL), t.L, col = c(teal, orange), lwd = c(2,2), ylab="")
legend("topleft", c("Log Hours Worked", "Log Population\nTrend"), col = c(teal, orange), lwd=c(2,2), bty="n")
axis(1, font=2)
axis(2, font=2)

ts.plot(lnL, col = teal, lwd=2, ylab="")
legend("topleft", "L", col = teal, lwd=2, bty="n")
abline(h=0, col=orange, lwd=2)
axis(1, font=2)
axis(2, font=2)

# CPI inflation
ts.plot(cpi, col = c(teal), lwd = c(2), ylab="")
legend("topleft", c("CPI"), col = c(teal), lwd=c(2), bty="n")
axis(1, font=2)
axis(2, font=2)

ts.plot(inf, col = teal, lwd=2, ylab="")
legend("topright", "CPI inflation ", col = teal, lwd=2, bty="n")
abline(h=0, col=orange, lwd=2)
axis(1, font=2)
axis(2, font=2)




