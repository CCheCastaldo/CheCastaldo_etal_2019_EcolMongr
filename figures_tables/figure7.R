library(coda)
library(gdata)
library(boot)
library(ggplot2)

version <- "i-009f9c75a9206a0ac"

#######################################################################################################################################
# load data

load(file = "../library/salix_borer_wide.rda")

load(file = paste0("../model_build/global_model/", version, "/mcmc_samples_core.rda"))

chains <- rbind(as.matrix(MCMCsamplesCore[[1]]), as.matrix(MCMCsamplesCore[[2]]), as.matrix(MCMCsamplesCore[[3]]))

load(file = "../stem_biomass/mcmc_stem_biomass.rda")

chains2 <- rbind(as.matrix(stemBiomass[[1]]), as.matrix(stemBiomass[[2]]), as.matrix(stemBiomass[[3]]))

#######################################################################################################################################
# determine rgr to use by habitat

growthUpland <- c(SalixBorerWide$gr09.sd[which(SalixBorerWide$habitat==1 & SalixBorerWide$site10==1)], 
            
                 SalixBorerWide$gr10.sd[which(SalixBorerWide$habitat==1 & SalixBorerWide$site11==1)])

(growthUpland <- median(growthUpland))

#######################################################################################################################################
# get data to plot

dPlotData <- function (growthVal, minX, maxX, mesh){
  
  tbd2.mean <- unique(SalixBorerWide$tbd2.mean)

  tbd2.sd <- unique(SalixBorerWide$tbd2.sd)
  
  tbd <- seq(minX, maxX, (maxX - minX) / mesh)
  
  tbd.sd <- (tbd - tbd2.mean)/tbd2.sd
  
  intGroup1 <- chains[, matchcols(chains, with = "mu.gamma\\[3")]

  intGroup2 <- chains[, matchcols(chains, with = "mu.gamma\\[4")]

  slopeGroup1Repro <- chains[, matchcols(chains, with = "beta.repro.gamma\\[3")]

  slopeGroup2Repro <- chains[, matchcols(chains, with = "beta.repro.gamma\\[4")]
  
  slopeGroup1rgr <- chains[, matchcols(chains, with = "beta.rgr.gamma\\[2")]
  
  slopeGroup2rgr <- chains[, matchcols(chains, with = "beta.rgr.gamma\\[2")]

  slopeGroup1d <- chains[, matchcols(chains, with = "beta.d.gamma\\[2")]
  
  slopeGroup2d <- chains[, matchcols(chains, with = "beta.d.gamma\\[2")]
  
  samplesGroup1 <- samplesGroup2 <- samplesGroup3 <- samplesGroup4 <- array(NA, dim = c(4500, length(tbd.sd)))
  
  for (i in 1:4500){
    
    samplesGroup1[i,] <- inv.logit(intGroup1[i] + 0 * slopeGroup1Repro[i] + growthVal * slopeGroup1rgr[i] + tbd.sd * slopeGroup1d[i])
  
    samplesGroup2[i,] <- inv.logit(intGroup2[i] + 1 * slopeGroup1Repro[i] + growthVal * slopeGroup1rgr[i] + tbd.sd * slopeGroup1d[i])

    samplesGroup3[i,] <- inv.logit(intGroup1[i] + 0 * slopeGroup1Repro[i] + growthVal * slopeGroup1rgr[i] + tbd.sd * slopeGroup1d[i])
    
    samplesGroup4[i,] <- inv.logit(intGroup2[i] + 1 * slopeGroup1Repro[i] + growthVal * slopeGroup1rgr[i] + tbd.sd * slopeGroup1d[i])
  
  }
  
  y <- 1 - (samplesGroup1 + samplesGroup2 + samplesGroup3 + samplesGroup4)/4
      
  out <- list(y, tbd)

  return(out)

}

#######################################################################################################################################
# estimate n2

out <- dPlotData(growthUpland, 10, 40, 1000)

beta <- fit <- NA

for (i in 1:4500){
  
  y <- log(-log(out[[1]][i,]))
  
  x <- log(out[[2]])

  beta[i] <- lm(y ~ x)[[1]][2]
  
  fit[i] <- summary(lm(y ~ x))[[8]]

}

quantile(fit, c(.025, .5, .975))

quantile(beta, c(.025, .5, .975))

length(beta[which(beta>=1)])/4500

#######################################################################################################################################
# n1

n1 <- chains2[, matchcols(chains2, with = "n1")]

quantile(n1, c(.025, .5, .975))

#######################################################################################################################################
# estimate n3 as ratio of n2/flux

n3 <- beta / n1

quantile(n3, c(.025, .5, .975))

length(n3[which(n3>=1)])/4500

#######################################################################################################################################
# plot these

data_summary <- function(x) {
  
  m <- mean(x)
  
  ymin <- as.numeric(quantile(x, prob = .025))
  
  ymax <- as.numeric(quantile(x, prob = .975))
  
  return(c(y=m,ymin=ymin,ymax=ymax))
  
}

setEPS()

postscript("fig9n1.eps")

work1 <- as.data.frame(n1)

work1$type <- "n1"

names(work1) <- c("samples", "type")

work2 <- work1

work2$samples <- 2^work1$samples

work2$type <- "delta1"

work3 <- rbind(work1, work2)

work3$type <- as.factor(work3$type)

p <- ggplot(work3, aes(x = type, y = samples))+ geom_violin(trim = FALSE, fill = '#A4A4A4', color = "black")

p <- p + stat_summary(fun.data=data_summary) + theme_minimal() 

p

dev.off()

setEPS()

postscript("fig9n2.eps")

work1 <- as.data.frame(beta)

work1$type <- "n2"

names(work1) <- c("samples", "type")

work2 <- work1

work2$samples <- 2^work1$samples

work2$type <- "delta2"

work3 <- rbind(work1, work2)

work3$type <- as.factor(work3$type)

p <- ggplot(work3, aes(x = type, y = samples))+ geom_violin(trim = FALSE, fill = '#A4A4A4', color = "black")

p <- p + stat_summary(fun.data=data_summary) + theme_minimal() 

p

dev.off()

setEPS()

postscript("fig9n3.eps")

work1 <- as.data.frame(n3)

work1$type <- "n3"

names(work1) <- c("samples", "type")

work2 <- work1

work2$samples <- 2^work1$samples

work2$type <- "delta3"

work3 <- rbind(work1, work2)

work3$type <- as.factor(work3$type)

p <- ggplot(work3, aes(x = type, y = samples))+ geom_violin(trim = FALSE, fill = '#A4A4A4', color = "black")

p <- p + stat_summary(fun.data=data_summary) + theme_minimal() 

p

dev.off()
