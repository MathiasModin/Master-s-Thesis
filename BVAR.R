## Thesis: Bayesian Estimation of VAR -------------------------------------
##                                                                       ##
## File name: BVAR.R                                                     ##  
## Author:    Mathias Concha Modin                                       ##
## Email:     mathias.modin@ne.su.se                                     ##
##                                                                       ##




# Loading data ------------------------------------------------------------

#setwd("Z:/NEK MASTEr/Thesis/Data")
rm(list=ls())
ndata <- ts(read.csv2("ndata.csv", header=T)[,-1], start=c(1994,1), freq=4)
x <- diff(ndata[,-c(4,5)])


n <- colnames(ndata)
data <- cbind(x[,"Y"], x[,"W"], x[,"PIE"], ndata[,"L"], 
                 x[,"C"], x[,"I"], ndata[,"R"])

colnames(data) <- cbind("Output", "Wages", "Inflation", "Labour", 
                        "Consumption", "Investment", "Repo Rate")

data <- window(data, start = c(1994, 2))
# Estimation --------------------------------------------------------------

library(BMR)

obj <- BVARM(data, coefprior = c(1, 1, 1, 1, 1, 1, 1), p=2, constant=T, 
             irf.periods = 20, VType=2, HP1=0.2, HP2=0.5, HP3=10^5, HP4=2
             , keep=15000, burnin=5000)

source("plot.IRF.R")
IRFs(obj)
