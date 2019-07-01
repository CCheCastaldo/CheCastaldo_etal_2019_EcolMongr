library(coda)
library(gdata)
library(boot)

version <- "i-009f9c75a9206a0ac" 

#######################################################################################################################################
# load data

load(file = "../Library/SalixBorerWide.rda")

load(file = paste0("../ModelBuild/GlobalModel/", version, "/MCMCsamplesCore.rda"))

chains <- rbind(as.matrix(MCMCsamplesCore[[1]]), as.matrix(MCMCsamplesCore[[2]]), as.matrix(MCMCsamplesCore[[3]]))

#######################################################################################################################################
# determine x axis range
# 10 - 70 mm

upland <- c(SalixBorerWide$tbd10[which(SalixBorerWide$habitat==1 & SalixBorerWide$site10==1)], 
            
            SalixBorerWide$tbd11[which(SalixBorerWide$habitat==1 & SalixBorerWide$site11==1)])

riparian <- c(SalixBorerWide$tbd10[which(SalixBorerWide$habitat==0 & SalixBorerWide$site10==1)], 
              
              SalixBorerWide$tbd11[which(SalixBorerWide$habitat==0 & SalixBorerWide$site11==1)])

quantile(upland, prob = c(.125,.875))

quantile(riparian, prob = c(.125,.875))

#######################################################################################################################################
# get data to plot

dPlot <- function (intercept, slope, minX, maxX, mesh, intervalType){
  
  tbd2.mean <- unique(SalixBorerWide$tbd2.mean)

  tbd2.sd <- unique(SalixBorerWide$tbd2.sd)
  
  tbd <- seq(minX, maxX, (maxX - minX) / mesh)
  
  tbd.sd <- (tbd - tbd2.mean)/tbd2.sd
  
  int.upland <- chains[, matchcols(chains, with = intercept[1])]
  
  slope.upland <- chains[, matchcols(chains, with = slope[1])]
  
  int.riparian <- chains[, matchcols(chains, with = intercept[2])]
  
  slope.riparian <- chains[, matchcols(chains, with = slope[2])]
  
  samples.upland <- samples.riparian <- array(NA, dim = c(4500, length(tbd.sd)))
  
  for (i in 1:4500){
    
    samples.upland[i,] <- inv.logit(int.upland[i] + tbd.sd * slope.upland[i])
  
    samples.riparian[i,] <- inv.logit(int.riparian[i] + tbd.sd * slope.riparian[i])
    
  }
  
  y <- array(NA, dim = c(length(tbd.sd), 6))
  
  for (i in 1:length(tbd.sd)){
    
    y[i, 1] <- median(samples.upland[,i])
    
    y[i, 4] <- median(samples.riparian[,i])
    
    if (intervalType == "HPD") {
    
      y[i, 2] <- HPDinterval(as.mcmc(samples.upland[,i]), prob = .95)[1]
      
      y[i, 3] <- HPDinterval(as.mcmc(samples.upland[,i]), prob = .95)[2]
      
      y[i, 5] <- HPDinterval(as.mcmc(samples.riparian[,i]), prob = .95)[1]
      
      y[i, 6] <- HPDinterval(as.mcmc(samples.riparian[,i]), prob = .95)[2]
    
    }
    
    else {
      
      y[i, 2] <- quantile(samples.upland[,i], c(.025, .975))[1]
      
      y[i, 3] <- quantile(samples.upland[,i], c(.025, .975))[2]
      
      y[i, 5] <- quantile(samples.riparian[,i], c(.025, .975))[1]
      
      y[i, 6] <- quantile(samples.riparian[,i], c(.025, .975))[2]
    
    } 
    
  }
  
  out <- list(y, tbd.sd, tbd2.mean, tbd2.sd)
  
  return(out)

}

out <- dPlot(c("alpha.p\\[2", "alpha.p\\[3"), c("beta.d.p\\[2", "beta.d.p\\[3"), 10, 70, 1000, "quantile")

#######################################################################################################################################
# set axes limits

tbd2.mean <- out[[3]]

tbd2.sd <- out[[4]]

minX <- round((10 - out[[3]]) / out[[4]], 2)

maxX <- round((70 - out[[3]]) / out[[4]], 2)

xLabels <- c(10, 20, 30, 40, 50, 60, 70)

xTicks <- (xLabels - out[[3]]) / out[[4]]

#######################################################################################################################################
# make rugs

upland <- c(SalixBorerWide$tbd10.sd[which(SalixBorerWide$habitat==1 & SalixBorerWide$site10==1)], 
            
            SalixBorerWide$tbd11.sd[which(SalixBorerWide$habitat==1 & SalixBorerWide$site11==1)])

riparian <- c(SalixBorerWide$tbd10.sd[which(SalixBorerWide$habitat==0 & SalixBorerWide$site10==1)], 
              
              SalixBorerWide$tbd11.sd[which(SalixBorerWide$habitat==0 & SalixBorerWide$site11==1)])

upland <- quantile(upland, prob = c(.125, .25, .5, .75, .875))

riparian <- quantile(riparian, prob = c(.125, .25, .5, .75, .875))

#######################################################################################################################################
# plot it

pdf(file = "fig3.pdf", width = 4, height = 4)

plot(0, bty="n", type = "n" ,ylim = c(.5, 1), xlim = c(minX, maxX), ylab = "         Detection rate", xlab = "Stem diameter (mm)", axes = FALSE, cex.lab = .8)

g1H <- "gray68"

g1L <- "gray88"

g2H <- "gray68"

g2L <- "gray88"

axis(side = 1, line = -0.05 , at = xTicks, tick = TRUE, labels = xLabels, cex.axis = .8, cex.lab = .4)

axis(side = 2, at = c(.6, .7, .8, .9, 1), tick = TRUE, cex.axis = .8)

yp <- .595

yt <- -.05

ys <- -.005

ybox <- c(yp + yt, yp, yp, yp + yt)

polygon(c(upland[1], upland[1], upland[5], upland[5]), ybox, col = g2L, border = NA)

polygon(c(riparian[1], riparian[1], riparian[5], riparian[5]), ybox + yt + ys, col = g2L, border = NA)

polygon(c(upland[2], upland[2], upland[4], upland[4]), ybox, col = g2H, border = NA)

polygon(c(riparian[2], riparian[2], riparian[4], riparian[4]), ybox + yt + ys, col = g2H, border = NA)

points(upland[3], (yp + (yt + yp))/2, cex = 2, col = "black", pch = 16)

points(riparian[3], ((yp + yt + ys) + (yp + 2 * yt + ys))/2, cex = 2, col = "black", pch = 18)

polygon(c(out[[2]],rev(out[[2]])), c(out[[1]][,3], rev(out[[1]][,2])), col = g1L, border = NA)

polygon(c(out[[2]],rev(out[[2]])), c(out[[1]][,6], rev(out[[1]][,5])), col = g1L, border = NA)

lines(out[[2]], out[[1]][,1], type = "l" , cex = 0.4 , lwd = 1.5 , col = "black")

lines(out[[2]], out[[1]][,4], type = "l" , cex = 0.4 , lwd = 1.5 , lty = 2, col = "black")

legend(.58, .8, inset = .05, cex = .8, pt.cex = 1.5, title = "", c("upland", "riparian"), lty = c(1, 2), pch = c(16, 18), lwd = 1.5, bty = "n")

dev.off()
