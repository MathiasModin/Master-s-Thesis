## Function for plotting IRFs separetely ----------------------------------

IRFs<- function (obj, percentiles = c(0.05, 0.5, 0.95), save = FALSE, 
                     height = 13, width = 13) 
{
  irf.periods <- dim(obj$IRFs)[1]
  M <- dim(obj$IRFs)[4]
  mydata <- obj$data
  IRFs <- obj$IRFs
  plots <- list()
  if (class(dev.list()) != "NULL") {
    dev.off()
  }
  
  total <- dim(IRFs)[3]
  IRFUpper <- round(percentiles[3] * total)
  IRFMid <- round(percentiles[2] * total)
  IRFLower <- round(percentiles[1] * total)
  IRFPlot <- array(NA, dim = c(irf.periods, 4, M, M))
  for (i in 1:M) {
    for (k in 1:M) {
      IRFPData <- data.frame(IRFs[, k, IRFLower, i], 
                             IRFs[, k, IRFMid, i], IRFs[, k, IRFUpper, i], 
                             1:(irf.periods))
      IRFPData <- as.matrix(IRFPData)
      IRFPlot[, , k, i] <- IRFPData
    }
  }
  if (class(dev.list()) != "NULL") {
    dev.off()
  }
  if (save == TRUE) {
    cairo_ps(filename = "IRFs.eps", height = height, width = width)
  }
  
  Plist <- list()
  IRFL <- IRFM <- IRFU <- Time <- NULL
  for (i in 1:M) {
    for (k in 1:M) {
      NameImpulse <- colnames(mydata)[k]
      NameRespone <- colnames(mydata)[i]
      IRFDF <- IRFPlot[, , k, i]
      IRFDF <- data.frame(IRFDF)
      colnames(IRFDF) <- c("IRFL", "IRFM", "IRFU", "Time")
      gg1 <- ggplot(data = (IRFDF), aes(x = Time)) + xlab("") + 
          ggtitle(paste(NameImpulse)) + 
          theme(plot.title = element_text(lineheight=.8, face="bold"))
      gg2 <- gg1 + geom_ribbon(aes(ymin = IRFL, ymax = IRFU), 
                               color = "black", lty = 1, fill = "black", alpha = 0.4, 
                               size = 0.1) + geom_hline(yintercept = 0) + geom_line(aes(y = IRFM), 
                                                                                    color = "black", size = 2) + ylab(" ")
      gg3 <- gg2 + theme(panel.background = element_rect(fill = "white", 
                                                         colour = "grey5")) + theme(panel.grid.major = element_line(colour = "grey89"))
      namm <- paste("irf", i, k, sep="_")
      assign(namm, gg3)
      Sys.sleep(0.3)
      
    }
  }
  require(gridExtra)
  require(grid)
  require(ggplot2)
  require(lattice)
  
  
  # TThis is a bad way of plotting, figure out something smarter....pls
  grid.arrange(irf_1_1, irf_1_2,irf_1_3,irf_1_4,irf_1_5,irf_1_6,irf_1_7, ncol=3)
  grid.arrange(irf_2_1, irf_2_2,irf_2_3,irf_2_4,irf_2_5,irf_2_6,irf_2_7, ncol=3)
  grid.arrange(irf_3_1, irf_3_2,irf_3_3,irf_3_4,irf_3_5,irf_3_6,irf_3_7, ncol=3)
  grid.arrange(irf_4_1, irf_4_2,irf_4_3,irf_4_4,irf_4_5,irf_4_6,irf_4_7, ncol=3)
  grid.arrange(irf_5_1, irf_5_2,irf_5_3,irf_5_4,irf_5_5,irf_5_6,irf_5_7, ncol=3)
  grid.arrange(irf_6_1, irf_6_2,irf_6_3,irf_6_4,irf_6_5,irf_6_6,irf_6_7, ncol=3)
  grid.arrange(irf_7_1, irf_7_2,irf_7_3,irf_7_4,irf_7_5,irf_7_6,irf_7_7, ncol=3)
  
  return = list(IRFs = IRFPlot)
}




