library(coda)
library(gdata)
library(boot)

version <- "i-009f9c75a9206a0ac"

#######################################################################################################################################
# load data

load(file = "../Library/SalixBorerWide.rda")

load(file = paste0("../ModelBuild/GlobalModel/", version, "/MCMCsamplesCore.rda"))

chains <- rbind(as.matrix(MCMCsamplesCore[[1]]), as.matrix(MCMCsamplesCore[[2]]), as.matrix(MCMCsamplesCore[[3]]))

load(file = paste0("../ModelBuild/GlobalModel/", version, "/MCMCsamplesBDGamma.rda"))

chains2 <- rbind(as.matrix(MCMCsamplesBDGamma[[1]]), as.matrix(MCMCsamplesBDGamma[[2]]), as.matrix(MCMCsamplesBDGamma[[3]]))
                 
#######################################################################################################################################
# determine rgr to use by habitat

growthUpland <- c(SalixBorerWide$gr09.sd[which(SalixBorerWide$habitat==1 & SalixBorerWide$site10==1)], 
            
                 SalixBorerWide$gr10.sd[which(SalixBorerWide$habitat==1 & SalixBorerWide$site11==1)])

growthRiparian <- c(SalixBorerWide$gr09.sd[which(SalixBorerWide$habitat==0 & SalixBorerWide$site10==1)], 
              
                  SalixBorerWide$gr10.sd[which(SalixBorerWide$habitat==0 & SalixBorerWide$site11==1)])

(growthUpland <- median(growthUpland))

(growthRiparian <- median(growthRiparian))

#######################################################################################################################################
# get data to plot

dPlotData <- function (intercept, rgr, d, repro, reproVal, growthVal, minX, maxX, mesh, intervalType, y1, y2){
  
  tbd2.mean <- unique(SalixBorerWide$tbd2.mean)

  tbd2.sd <- unique(SalixBorerWide$tbd2.sd)
  
  tbd <- seq(minX, maxX, (maxX - minX) / mesh)
  
  tbd.sd <- (tbd - tbd2.mean)/tbd2.sd
  
  intGroup1 <- chains[, matchcols(chains, with = intercept[1])]

  intGroup2 <- chains[, matchcols(chains, with = intercept[2])]

  slopeGroup1Repro <- chains[, matchcols(chains, with = repro[1])]

  slopeGroup2Repro <- chains[, matchcols(chains, with = repro[2])]
  
  slopeGroup1rgr <- chains[, matchcols(chains, with = rgr[1])]
  
  slopeGroup2rgr <- chains[, matchcols(chains, with = rgr[2])]

  slopeGroup1d <- chains[, matchcols(chains, with = d[1])]
  
  slopeGroup2d <- chains[, matchcols(chains, with = d[2])]
  
  samplesGroup1 <- samplesGroup2 <- array(NA, dim = c(4500, length(tbd.sd)))
  
  for (i in 1:4500){
    
    samplesGroup1[i,] <- inv.logit(intGroup1[i] + reproVal[1] * slopeGroup1Repro[i] + growthVal[1] * slopeGroup1rgr[i] + tbd.sd * slopeGroup1d[i])
  
    samplesGroup2[i,] <- inv.logit(intGroup2[i] + reproVal[2] * slopeGroup2Repro[i] + growthVal[2] * slopeGroup2rgr[i] + tbd.sd * slopeGroup2d[i])

  }
  
  y <- array(NA, dim = c(length(tbd.sd), 6))
  
  for (i in 1:length(tbd.sd)){
    
    y[i, 1] <- median(samplesGroup1[,i])
    
    y[i, 4] <- median(samplesGroup2[,i])
    
    if (intervalType == "HPD") {
    
      y[i, 2] <- HPDinterval(as.mcmc(samplesGroup1[,i]), prob = .95)[1]
      
      y[i, 3] <- HPDinterval(as.mcmc(samplesGroup1[,i]), prob = .95)[2]
      
      y[i, 5] <- HPDinterval(as.mcmc(samplesGroup2[,i]), prob = .95)[1]
      
      y[i, 6] <- HPDinterval(as.mcmc(samplesGroup2[,i]), prob = .95)[2]
    
    }
    
    else {
      
      y[i, 2] <- quantile(samplesGroup1[,i], c(.025, .975))[1]
      
      y[i, 3] <- quantile(samplesGroup1[,i], c(.025, .975))[2]
    
      y[i, 5] <- quantile(samplesGroup2[,i], c(.025, .975))[1]
      
      y[i, 6] <- quantile(samplesGroup2[,i], c(.025, .975))[2]
    
    } 
    
  }

  g1 <- (quantile(chains2[, matchcols(chains2, with = y1)], prob = c(.125, .25, .5, .75, .875)) - tbd2.mean) / tbd2.sd
  
  g2 <- (quantile(chains2[, matchcols(chains2, with = y2)], prob = c(.125, .25, .5, .75, .875)) - tbd2.mean) / tbd2.sd

  out <- list(y, tbd.sd, tbd2.mean, tbd2.sd, g1, g2, samplesGroup1, samplesGroup2)

  return(out)

}

#######################################################################################################################################
# plot it

dPlot <- function (out, minX2, maxX2, cexAll = .79, cexP = 1.5){
  
  tbd2.mean <- out[[3]]
  
  tbd2.sd <- out[[4]]
  
  minX <- round((minX2 - out[[3]]) / out[[4]], 2)
  
  maxX <- round((maxX2 - out[[3]]) / out[[4]], 2)
  
  xLabels <- seq(minX2, maxX2, 10)
  
  xTicks <- (xLabels - out[[3]]) / out[[4]]
  
  plot(0, bty="n", type = "n" ,ylim = c(-.3, 1), xlim = c(minX, maxX), ylab = NA, xlab = NA, axes = FALSE, cex.lab = cexAll)

  g1H <- "gray68"
  
  g1L <- "gray88"
  
  g2H <- "gray68"
  
  g2L <- "gray88"

  axis(side = 1, line = -.1, at = xTicks, tick = TRUE, labels = xLabels, cex.axis = cexAll)

  axis(side = 2, at = c(0, .2, .4, .6, .8, 1), tick = TRUE, cex.axis = cexAll)

  yp <- 0

  yt <- -.15

  ys <- -.015

  ybox <- c(yp + yt, yp, yp, yp + yt)

  polygon(c(out[[5]][1], out[[5]][1], out[[5]][5], out[[5]][5]), ybox, col = g2L, border = NA)

  polygon(c(out[[6]][1], out[[6]][1], out[[6]][5], out[[6]][5]), ybox + yt + ys, col = g2L, border = NA)

  polygon(c(out[[5]][2], out[[5]][2], out[[5]][4], out[[5]][4]), ybox, col = g2H, border = NA)

  polygon(c(out[[6]][2], out[[6]][2], out[[6]][4], out[[6]][4]), ybox + yt + ys, col = g2H, border = NA)

  points(out[[5]][3], (yt - yp)/2, cex = cexP, col = "black", pch = 15)

  points(out[[6]][3], (yp + yt + ys) + (yt - yp)/2, cex = cexP, col = "black", pch = 16)

  polygon(c(out[[2]],rev(out[[2]])), c(out[[1]][,3], rev(out[[1]][,2])), col = g1L, border = NA)

  polygon(c(out[[2]],rev(out[[2]])), c(out[[1]][,6], rev(out[[1]][,5])), col = g2L, border = NA)

  lines(out[[2]], out[[1]][,4], type = "l" , cex = 0.4 , lwd = 1.5 , lty = 2, col = "black")

  lines(out[[2]], out[[1]][,1], type = "l" , cex = 0.4 , lwd = 1.5 , col = "black")

}

#######################################################################################################################################

pdf(file = "fig5.pdf", width = 8, height = 8)

par(mfrow = c(2, 2))

out <- dPlotData(c("mu.gamma\\[1", "mu.gamma\\[1"), c("beta.rgr.gamma\\[1", "beta.rgr.gamma\\[1"), c("beta.d.gamma\\[1", "beta.d.gamma\\[1"), 
             
             c("beta.repro.gamma\\[1", "beta.repro.gamma\\[1"), c(1,0), c(growthRiparian, growthRiparian), 10, 80, 1000, "quantile",
             
             y1 = "tbdNew\\[1,1,1\\]", y2 = "tbdNew\\[1,1,2\\]")

par(pin = c(3, 3), mai = c(.5, 1, 1, .25))

dPlot(out, 10, 80, 1, 2)

mtext("(a)", side = 3, line = 0, at = -.7, cex = 1, font = 2)

out <- dPlotData(c("mu.gamma\\[3", "mu.gamma\\[3"), c("beta.rgr.gamma\\[2", "beta.rgr.gamma\\[2"), c("beta.d.gamma\\[2", "beta.d.gamma\\[2"), 
                 
                 c("beta.repro.gamma\\[3", "beta.repro.gamma\\[3"), c(1,0), c(growthUpland, growthUpland), 10, 50, 1000, "quantile",
                 
                 y1 = "tbdNew\\[2,1,1\\]", y2 = "tbdNew\\[2,1,2\\]")

par(pin = c(3, 3), mai = c(.5, .25, 1, 1))

dPlot(out, 10, 50, 1, 2)

mtext("(b)", side = 3, line = 0, at = -.7, cex = 1, font = 2)

out <- dPlotData(c("mu.gamma\\[2", "mu.gamma\\[2"), c("beta.rgr.gamma\\[1", "beta.rgr.gamma\\[1"), c("beta.d.gamma\\[1", "beta.d.gamma\\[1"), 
                 
                 c("beta.repro.gamma\\[2", "beta.repro.gamma\\[2"), c(1,0), c(growthRiparian, growthRiparian), 10, 80, 1000, "quantile",
                 
                 y1 = "tbdNew\\[1,2,1\\]", y2 = "tbdNew\\[1,2,2\\]")

par(pin = c(3, 3), mai = c(1, 1, .5, .25))

dPlot(out, 10, 80, 1, 2)

mtext("(c)", side = 3, line = 0, at = -.7, cex = 1, font = 2)

out <- dPlotData(c("mu.gamma\\[4", "mu.gamma\\[4"), c("beta.rgr.gamma\\[2", "beta.rgr.gamma\\[2"), c("beta.d.gamma\\[2", "beta.d.gamma\\[2"), 
                 
                 c("beta.repro.gamma\\[4", "beta.repro.gamma\\[4"), c(1,0), c(growthUpland, growthUpland), 10, 50, 1000, "quantile",
                 
                 y1 = "tbdNew\\[2,2,1\\]", y2 = "tbdNew\\[2,2,2\\]")

par(pin = c(3, 3), mai = c(1, .5, .5, 1))

dPlot(out, 10, 50, 1, 2)

legend(-.4, .55, inset = .05, cex = 1, title = "", c("flowering", "non-reproductive"), lty = c(1, 2), lwd = 1.5, bty = "n")

mtext("(d)", side = 3, line = 0, at = -.7, cex = 1, font = 2)

par(oma = c(2, 2, 0, 0))

mtext("Stem diameter (mm)", outer = TRUE, cex = 1, side = 1)
       
mtext("Colonization rate", outer = TRUE, cex = 1, side = 2)

dev.off()

#######################################################################################################################################
# load data

load(file = "../Library/SalixBorerWide.rda")

load(file = paste0("../ModelBuild/GlobalModel/", version, "/MCMCsamplesCore.rda"))

chains <- rbind(as.matrix(MCMCsamplesCore[[1]]), as.matrix(MCMCsamplesCore[[2]]), as.matrix(MCMCsamplesCore[[3]]))

load(file = paste0("../ModelBuild/GlobalModel/", version, "/MCMCsamplesBDPhi.rda"))

chains2 <- rbind(as.matrix(MCMCsamplesBDPhi[[1]]), as.matrix(MCMCsamplesBDPhi[[2]]), as.matrix(MCMCsamplesBDPhi[[3]]))

#######################################################################################################################################

pdf(file = "fig6.pdf", width = 8, height = 8)

par(mfrow = c(2, 2))

out <- dPlotData(c("mu.phi\\[1", "mu.phi\\[1"), c("beta.rgr.phi\\[1", "beta.rgr.phi\\[1"), c("beta.d.phi\\[1", "beta.d.phi\\[1"), 
                 
                 c("beta.repro.phi\\[1", "beta.repro.phi\\[1"), c(1,0), c(growthRiparian, growthRiparian), 10, 90, 1000, "quantile",
                 
                 y1 = "tbdNew\\[1,1,1\\]", y2 = "tbdNew\\[1,1,2\\]")

par(pin = c(3, 3), mai = c(.5, 1, 1, .25))

dPlot(out, 10, 100, 1, 2)

mtext("(a)", side = 3, line = 0, at = -.7, cex = 1, font = 2)

out <- dPlotData(c("mu.phi\\[3", "mu.phi\\[3"), c("beta.rgr.phi\\[2", "beta.rgr.phi\\[2"), c("beta.d.phi\\[2", "beta.d.phi\\[2"), 
                 
                 c("beta.repro.phi\\[3", "beta.repro.phi\\[3"), c(1,0), c(growthUpland, growthUpland), 10, 80, 1000, "quantile",
                 
                 y1 = "tbdNew\\[2,1,1\\]", y2 = "tbdNew\\[2,1,2\\]")

par(pin = c(3, 3), mai = c(.5, .25, 1, 1))

dPlot(out, 10, 80, 1, 2)

mtext("(b)", side = 3, line = 0, at = -.7, cex = 1, font = 2)

out <- dPlotData(c("mu.phi\\[2", "mu.phi\\[2"), c("beta.rgr.phi\\[1", "beta.rgr.phi\\[1"), c("beta.d.phi\\[1", "beta.d.phi\\[1"), 
                 
                 c("beta.repro.phi\\[2", "beta.repro.phi\\[2"), c(1,0), c(growthRiparian, growthRiparian), 10, 90, 1000, "quantile",
                 
                 y1 = "tbdNew\\[1,2,1\\]", y2 = "tbdNew\\[1,2,2\\]")

par(pin = c(3, 3), mai = c(1, 1, .5, .25))

dPlot(out, 10, 100, 1, 2)

mtext("(c)", side = 3, line = 0, at = -.7, cex = 1, font = 2)

out <- dPlotData(c("mu.phi\\[4", "mu.phi\\[4"), c("beta.rgr.phi\\[2", "beta.rgr.phi\\[2"), c("beta.d.phi\\[2", "beta.d.phi\\[2"), 
                 
                 c("beta.repro.phi\\[4", "beta.repro.phi\\[4"), c(1,0), c(growthUpland, growthUpland), 10, 80, 1000, "quantile",
                 
                 y1 = "tbdNew\\[2,2,1\\]", y2 = "tbdNew\\[2,2,2\\]")

par(pin = c(3, 3), mai = c(1, .5, .5, 1))

dPlot(out, 10, 80, 1, 2)

legend(-.4, .55, inset = .05, cex = 1, title = "", c("flowering", "non-reproductive"), lty = c(1, 2), lwd = 1.5, bty = "n")

mtext("(d)", side = 3, line = 0, at = -.7, cex = 1, font = 2)

par(oma = c(2, 2, 0, 0))

mtext("Stem diameter (mm)", outer = TRUE, cex = 1, side = 1)

mtext("Re-attack rate", outer = TRUE, cex = 1, side = 2)

dev.off()
