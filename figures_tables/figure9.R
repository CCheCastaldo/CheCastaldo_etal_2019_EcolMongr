library(coda)
library(gdata)
library(boot)
library(coda)

version <- "i-009f9c75a9206a0ac"

#######################################################################################################################################
# load data

load(file = "../library/salix_borer_wide.rda")

site <- subset(SalixBorerWide, select = c(stem, site09, site10, site11, habitat, sex))

i <- which(site$site09==0); site$site09[i] <- NA

i <- which(site$site10==0); site$site10[i] <- NA

i <- which(site$site11==0); site$site11[i] <- NA

load(file = paste0("../model_build/global_model/", version, "/mcmc_samples_z.rda"))

load(file = paste0("../model_build/global_model/", version, "/mcmc_samples_psi_1.rda"))

load(file = paste0("../model_build/global_model/", version, "/mcmc_samples_psi_2.rda"))

load(file = paste0("../model_build/global_model/", version, "/mcmc_samples_psi_3.rda"))

load(file = paste0("../model_build/global_model/", version, "/mcmc_samples_gamma_1.rda"))

load(file = paste0("../model_build/global_model/", version, "/mcmc_samples_gamma_2.rda"))

load(file = paste0("../model_build/global_model/", version, "/mcmc_samples_phi_1.rda"))

load(file = paste0("../model_build/global_model/", version, "/mcmc_samples_phi_2.rda"))

psi1 <- t(rbind(as.matrix(MCMCsamplesPsi1[[1]]), as.matrix(MCMCsamplesPsi1[[2]]), as.matrix(MCMCsamplesPsi1[[3]])))

psi2 <- t(rbind(as.matrix(MCMCsamplesPsi2[[1]]), as.matrix(MCMCsamplesPsi2[[2]]), as.matrix(MCMCsamplesPsi2[[3]])))

psi3 <- t(rbind(as.matrix(MCMCsamplesPsi3[[1]]), as.matrix(MCMCsamplesPsi3[[2]]), as.matrix(MCMCsamplesPsi3[[3]])))

gamma1 <- t(rbind(as.matrix(MCMCsamplesGamma1[[1]]), as.matrix(MCMCsamplesGamma1[[2]]), as.matrix(MCMCsamplesGamma1[[3]])))

gamma2 <- t(rbind(as.matrix(MCMCsamplesGamma2[[1]]), as.matrix(MCMCsamplesGamma2[[2]]), as.matrix(MCMCsamplesGamma2[[3]])))

phi1 <- t(rbind(as.matrix(MCMCsamplesPhi1[[1]]), as.matrix(MCMCsamplesPhi1[[2]]), as.matrix(MCMCsamplesPhi1[[3]])))

phi2 <- t(rbind(as.matrix(MCMCsamplesPhi2[[1]]), as.matrix(MCMCsamplesPhi2[[2]]), as.matrix(MCMCsamplesPhi2[[3]])))

z1 <- t(rbind(as.matrix(MCMCsamplesZ[, unlist(matchcols(MCMCsamplesZ[[1]], with = c("z\\[", "1\\]")))])))

z2 <- t(rbind(as.matrix(MCMCsamplesZ[, unlist(matchcols(MCMCsamplesZ[[1]], with = c("z\\[", "2\\]")))])))

#######################################################################################################################################
# load data

makeChains <- function (input1, input2, rate, type){
  
  probs1 <- probs2 <-  array(NA, dim = c(2087, 4500))
  
  w1 <- w2 <- array(1, dim = c(2087, 4500))
  
  for (i in 1:dim(z1)[1]) {
    
    for (j in 1:dim(z1)[2]) {
      
      if (z1[i, j]==rate) w1[i, j] <- NA
      
      if (z2[i, j]==rate) w2[i, j] <- NA

      probs1[i, j] <- site$site10[i] * input1[i, j] * w1[i, j]

      probs2[i, j] <- site$site11[i] * input2[i, j] * w2[i, j]
    
    }
    
  }

  j <- which(site$habitat==0 & site$sex==0)

  k <- which(site$habitat==0 & site$sex==1)

  l <- which(site$habitat==1 & site$sex==0)

  m <- which(site$habitat==1 & site$sex==1)

  chains <- array(NA, dim = c(4500, 8))
  
  for (i in 1:dim(probs2)[2]) {
  
    chains[i, 1] <- type(probs1[j, i], na.rm=TRUE)
    chains[i, 2] <- type(probs2[j, i], na.rm=TRUE)
    chains[i, 3] <- type(probs1[k, i], na.rm=TRUE)
    chains[i, 4] <- type(probs2[k, i], na.rm=TRUE)
    chains[i, 5] <- type(probs1[l, i], na.rm=TRUE)
    chains[i, 6] <- type(probs2[l, i], na.rm=TRUE)
    chains[i, 7] <- type(probs1[m, i], na.rm=TRUE)
    chains[i, 8] <- type(probs2[m, i], na.rm=TRUE)
  }

  return(chains)

}

#####################################################################################################################################################
# gamma

chains <- makeChains(gamma1, gamma2, 1, median)

dif <- chains[,1] - chains[,3]; length(dif[which(dif>0)])/4500
dif <- chains[,2] - chains[,4]; length(dif[which(dif>0)])/4500
dif <- chains[,5] - chains[,7]; length(dif[which(dif>0)])/4500
dif <- chains[,6] - chains[,8]; length(dif[which(dif>0)])/4500

dif <- chains[,1] - chains[,3]; length(dif[which(dif>0)])/4500
dif <- chains[,2] - chains[,4]; length(dif[which(dif>0)])/4500
dif <- chains[,5] - chains[,7]; length(dif[which(dif>0)])/4500
dif <- chains[,6] - chains[,8]; length(dif[which(dif>0)])/4500

dif <- chains[,1] - chains[,2]; length(dif[which(dif>0)])/4500
dif <- chains[,3] - chains[,4]; length(dif[which(dif>0)])/4500
dif <- chains[,5] - chains[,6]; length(dif[which(dif>0)])/4500
dif <- chains[,7] - chains[,8]; length(dif[which(dif>0)])/4500

gamma.m <- gamma.LCB <- gamma.UCB <- matrix(NA, nrow = 2, ncol = 4)

gamma.m[1, 1] <- median(chains[, 1])
gamma.m[2, 1] <- median(chains[, 2])
gamma.m[1, 2] <- median(chains[, 3])
gamma.m[2, 2] <- median(chains[, 4])
gamma.m[1, 3] <- median(chains[, 5])
gamma.m[2, 3] <- median(chains[, 6])
gamma.m[1, 4] <- median(chains[, 7])
gamma.m[2, 4] <- median(chains[, 8])

gamma.LCB[1, 1] <- quantile(chains[, 1], prob = .025)
gamma.LCB[2, 1] <- quantile(chains[, 2], prob = .025)
gamma.LCB[1, 2] <- quantile(chains[, 3], prob = .025)
gamma.LCB[2, 2] <- quantile(chains[, 4], prob = .025)
gamma.LCB[1, 3] <- quantile(chains[, 5], prob = .025)
gamma.LCB[2, 3] <- quantile(chains[, 6], prob = .025)
gamma.LCB[1, 4] <- quantile(chains[, 7], prob = .025)
gamma.LCB[2, 4] <- quantile(chains[, 8], prob = .025)

gamma.UCB[1, 1] <- quantile(chains[, 1], prob = .975)
gamma.UCB[2, 1] <- quantile(chains[, 2], prob = .975)
gamma.UCB[1, 2] <- quantile(chains[, 3], prob = .975)
gamma.UCB[2, 2] <- quantile(chains[, 4], prob = .975)
gamma.UCB[1, 3] <- quantile(chains[, 5], prob = .975)
gamma.UCB[2, 3] <- quantile(chains[, 6], prob = .975)
gamma.UCB[1, 4] <- quantile(chains[, 7], prob = .975)
gamma.UCB[2, 4] <- quantile(chains[, 8], prob = .975)

#####################################################################################################################################################
# phi

chains <- makeChains(phi1, phi2, 0, median)

dif <- chains[,1] - chains[,3]; length(dif[which(dif>0)])/4500
dif <- chains[,2] - chains[,4]; length(dif[which(dif>0)])/4500
dif <- chains[,5] - chains[,7]; length(dif[which(dif>0)])/4500
dif <- chains[,6] - chains[,8]; length(dif[which(dif>0)])/4500

dif <- chains[,1] - chains[,3]; length(dif[which(dif>0)])/4500
dif <- chains[,2] - chains[,4]; length(dif[which(dif>0)])/4500
dif <- chains[,5] - chains[,7]; length(dif[which(dif>0)])/4500
dif <- chains[,6] - chains[,8]; length(dif[which(dif>0)])/4500

dif <- chains[,1] - chains[,2]; length(dif[which(dif>0)])/4500
dif <- chains[,3] - chains[,4]; length(dif[which(dif>0)])/4500
dif <- chains[,5] - chains[,6]; length(dif[which(dif>0)])/4500
dif <- chains[,7] - chains[,8]; length(dif[which(dif>0)])/4500

phi.m <- phi.LCB <- phi.UCB <- matrix(NA, nrow = 2, ncol = 4)

phi.m[1, 1] <- median(chains[, 1])
phi.m[2, 1] <- median(chains[, 2])
phi.m[1, 2] <- median(chains[, 3])
phi.m[2, 2] <- median(chains[, 4])
phi.m[1, 3] <- median(chains[, 5])
phi.m[2, 3] <- median(chains[, 6])
phi.m[1, 4] <- median(chains[, 7])
phi.m[2, 4] <- median(chains[, 8])

phi.LCB[1, 1] <- quantile(chains[, 1], prob = .025)
phi.LCB[2, 1] <- quantile(chains[, 2], prob = .025)
phi.LCB[1, 2] <- quantile(chains[, 3], prob = .025)
phi.LCB[2, 2] <- quantile(chains[, 4], prob = .025)
phi.LCB[1, 3] <- quantile(chains[, 5], prob = .025)
phi.LCB[2, 3] <- quantile(chains[, 6], prob = .025)
phi.LCB[1, 4] <- quantile(chains[, 7], prob = .025)
phi.LCB[2, 4] <- quantile(chains[, 8], prob = .025)

phi.UCB[1, 1] <- quantile(chains[, 1], prob = .975)
phi.UCB[2, 1] <- quantile(chains[, 2], prob = .975)
phi.UCB[1, 2] <- quantile(chains[, 3], prob = .975)
phi.UCB[2, 2] <- quantile(chains[, 4], prob = .975)
phi.UCB[1, 3] <- quantile(chains[, 5], prob = .975)
phi.UCB[2, 3] <- quantile(chains[, 6], prob = .975)
phi.UCB[1, 4] <- quantile(chains[, 7], prob = .975)
phi.UCB[2, 4] <- quantile(chains[, 8], prob = .975)

#####################################################################################################################################################
# psi
  
probs1 <- site$site09 * psi1
  
probs2 <- site$site10 * psi2
  
probs3 <- site$site11 * psi3

j <- which(site$habitat==0 & site$sex==0)
  
k <- which(site$habitat==0 & site$sex==1)
  
l <- which(site$habitat==1 & site$sex==0)
  
m <- which(site$habitat==1 & site$sex==1)
  
chains <- array(NA, dim = c(4500, 12))
  
for (i in 1:dim(probs2)[2]) {
    
  chains[i, 1] <- median(probs1[j, i], na.rm=TRUE)
  chains[i, 2] <- median(probs2[j, i], na.rm=TRUE)
  chains[i, 3] <- median(probs3[j, i], na.rm=TRUE)
  chains[i, 4] <- median(probs1[k, i], na.rm=TRUE)
  chains[i, 5] <- median(probs2[k, i], na.rm=TRUE)
  chains[i, 6] <- median(probs3[k, i], na.rm=TRUE)
  chains[i, 7] <- median(probs1[l, i], na.rm=TRUE)
  chains[i, 8] <- median(probs2[l, i], na.rm=TRUE)
  chains[i, 9] <- median(probs3[l, i], na.rm=TRUE)
  chains[i, 10] <- median(probs1[m, i], na.rm=TRUE)
  chains[i, 11] <- median(probs2[m, i], na.rm=TRUE)
  chains[i, 12] <- median(probs3[m, i], na.rm=TRUE)
  
}
  
dif <- chains[,4] - chains[,1]; length(dif[which(dif>0)])/4500
dif <- chains[,5] - chains[,2]; length(dif[which(dif>0)])/4500
dif <- chains[,6] - chains[,3]; length(dif[which(dif>0)])/4500

dif <- chains[,10] - chains[,7]; length(dif[which(dif>0)])/4500
dif <- chains[,11] - chains[,8]; length(dif[which(dif>0)])/4500
dif <- chains[,12] - chains[,9]; length(dif[which(dif>0)])/4500


psi.m <- psi.LCB <- psi.UCB <- matrix(NA, nrow = 3, ncol = 4)

psi.m[1, 1] <- median(chains[, 1])
psi.m[2, 1] <- median(chains[, 2])
psi.m[3, 1] <- median(chains[, 3])
psi.m[1, 2] <- median(chains[, 4])
psi.m[2, 2] <- median(chains[, 5])
psi.m[3, 2] <- median(chains[, 6])
psi.m[1, 3] <- median(chains[, 7])
psi.m[2, 3] <- median(chains[, 8])
psi.m[3, 3] <- median(chains[, 9])
psi.m[1, 4] <- median(chains[, 10])
psi.m[2, 4] <- median(chains[, 11])
psi.m[3, 4] <- median(chains[, 12])

psi.LCB[1, 1] <- quantile(chains[, 1], prob = .025)
psi.LCB[2, 1] <- quantile(chains[, 2], prob = .025)
psi.LCB[3, 1] <- quantile(chains[, 3], prob = .025)
psi.LCB[1, 2] <- quantile(chains[, 4], prob = .025)
psi.LCB[2, 2] <- quantile(chains[, 5], prob = .025)
psi.LCB[3, 2] <- quantile(chains[, 6], prob = .025)
psi.LCB[1, 3] <- quantile(chains[, 7], prob = .025)
psi.LCB[2, 3] <- quantile(chains[, 8], prob = .025)
psi.LCB[3, 3] <- quantile(chains[, 9], prob = .025)
psi.LCB[1, 4] <- quantile(chains[, 10], prob = .025)
psi.LCB[2, 4] <- quantile(chains[, 11], prob = .025)
psi.LCB[3, 4] <- quantile(chains[, 12], prob = .025)

psi.UCB[1, 1] <- quantile(chains[, 1], prob = .975)
psi.UCB[2, 1] <- quantile(chains[, 2], prob = .975)
psi.UCB[3, 1] <- quantile(chains[, 3], prob = .975)
psi.UCB[1, 2] <- quantile(chains[, 4], prob = .975)
psi.UCB[2, 2] <- quantile(chains[, 5], prob = .975)
psi.UCB[3, 2] <- quantile(chains[, 6], prob = .975)
psi.UCB[1, 3] <- quantile(chains[, 7], prob = .975)
psi.UCB[2, 3] <- quantile(chains[, 8], prob = .975)
psi.UCB[3, 3] <- quantile(chains[, 9], prob = .975)
psi.UCB[1, 4] <- quantile(chains[, 10], prob = .975)
psi.UCB[2, 4] <- quantile(chains[, 11], prob = .975)
psi.UCB[3, 4] <- quantile(chains[, 12], prob = .975)


#####################################################################################################################################################
# plot it

pdf(file = "fig9.pdf", width = 8, height = 8)

par(mfrow = c(3, 2))

cexAll <- 1.1

d <- .04
d2 <- .03
d3 <- .01
d4 <- .03

f1 <- c(.1, .8); f1s <- c(.1, .3, .5, .7, .9)
f2 <- c(.5, .9); f2s <- c(.5, .6, .7, .8, .9)
f3 <- c(0,  .4); f3s <- c(0, .1, .2, .3, .4)
f4 <- c(.4, .8); f4s <- c(.4, .5, .6, .7, .8)
f5 <- c(.5, .9); f5s <- c(.5, .6, .7, .8, .9)
f6 <- c(.6,  1); f6s <- c(.6, .7, .8, .9, 1)

#### Riparian, Occupancy

par(pin = c(2, 2), mai = c(.75, 1, .25, .25))

plot(0, bty = "n", type = "n" ,ylim = f1, xlim = c(.9, 3.1), ylab = "Occupancy", xlab = "Years", axes = FALSE, cex.lab = cexAll)

axis(side = 1, at = c(1, 2, 3), tick = TRUE, labels = c(2009, 2010, 2011), cex.axis = cexAll)

axis(side = 2, at = f1s, tick = TRUE, cex.axis = cexAll)

lines(c(1 - d, 2 - d, 3 - d), psi.m[,1], type = "b" , cex = cexAll, lwd = 1.5 , lty = 1, col = "black", pch = 16)

lines(c(1 + d , 2 + d, 3 + d), psi.m[,2], type = "b" , cex = cexAll, lwd = 1.5 , lty = 2, col = "black", pch = 16)

arrows(c(1 - d, 2 - d, 3 - d), psi.LCB[,1], c(1 - d, 2 - d, 3 - d), psi.UCB[,1], code = 3, angle = 90, length = 0.05)

arrows(c(1 + d, 2 + d, 3 + d), psi.LCB[,2], c(1 + d, 2 + d, 3 + d), psi.UCB[,2], code = 3, angle = 90, length = 0.05)

mtext("(a)", side = 3, line = 0, at = 1, cex = cexAll, font = 1)

#### Upland, Occupancy

par(pin = c(2, 2), mai = c(.75, .25, .25, 1))

plot(0, bty="n", type = "n", ylim = f2, xlim = c(.9, 3.1), ylab = NA, xlab = "Years", axes = FALSE, cex.lab = cexAll)

axis(side = 1, at = c(1, 2, 3), tick = TRUE, labels = c(2009, 2010, 2011), cex.axis = cexAll)

axis(side = 2, at = f2s, tick = TRUE, cex.axis = cexAll)

lines(c(1 - d, 2 - d, 3 - d), psi.m[,3], type = "b" , cex = cexAll, lwd = 1.5 , lty = 1, col = "black", pch = 16)

lines(c(1 + d , 2 + d, 3 + d), psi.m[,4], type = "b" , cex = cexAll, lwd = 1.5 , lty = 2, col = "black", pch = 16)

arrows(c(1 - d, 2 - d, 3 - d), psi.LCB[,3], c(1 - d, 2 - d, 3 - d), psi.UCB[,3], code = 3, angle = 90, length = 0.05)

arrows(c(1 + d, 2 + d, 3 + d), psi.LCB[,4], c(1 + d, 2 + d, 3 + d), psi.UCB[,4], code = 3, angle = 90, length = 0.05)

mtext("(b)", side = 3, line = 0, at = 1, cex = cexAll, font = 1)

#### Riparian, Colonization

par(pin = c(2, 2), mai = c(.75, 1, .25, .25))

plot(0, bty = "n", type = "n", ylim = f3, xlim = c(.8, 2.2), ylab = "Colonization rate", xlab = "Between years", axes = FALSE, cex.lab = cexAll)

axis(side = 1, at = c(1, 2), tick = TRUE, labels = c("2009-2010","2010-2011"), cex.axis = cexAll)

axis(side = 2, at = f3s, tick = TRUE, cex.axis = cexAll)

lines(c(1 - d2, 2 - d2), gamma.m[,1], type = "b" , cex = cexAll, lwd = 1.5 , lty = 1, col = "black", pch = 16)

lines(c(1 + d2 , 2 + d2), gamma.m[,2], type = "b" , cex = cexAll, lwd = 1.5 , lty = 2, col = "black", pch = 16)

arrows(c(1 - d2, 2 - d2, 3 - d2), gamma.LCB[,1], c(1 - d2, 2 - d2, 3 - d2), gamma.UCB[,1], code = 3, angle = 90, length = 0.05)

arrows(c(1 + d2, 2 + d2, 3 + d2), gamma.LCB[,2], c(1 + d2, 2 + d2, 3 + d2), gamma.UCB[,2], code = 3, angle = 90, length = 0.05)

mtext("(c)", side = 3, line = 0, at = 1, cex = cexAll, font = 1)

#### Upland, Colonization

par(pin = c(2, 2), mai = c(.75, .25, .25, 1))

plot(0, bty = "n", type = "n", ylim = f4, xlim = c(.8, 2.2), ylab = NA, xlab = "Between years", axes = FALSE, cex.lab = cexAll)

axis(side = 1, at = c(1, 2), tick = TRUE, labels = c("2009-2010","2010-2011"), cex.axis = cexAll)

axis(side = 2, at = f4s, tick = TRUE, cex.axis = cexAll)

lines(c(1 - d2, 2 - d2), gamma.m[,3], type = "b" , cex = cexAll, lwd = 1.5 , lty = 1, col = "black", pch = 16)

lines(c(1 + d2, 2 + d2), gamma.m[,4], type = "b" , cex = cexAll, lwd = 1.5 , lty = 2, col = "black", pch = 16)

arrows(c(1 - d2, 2 - d2, 3 - d2), gamma.LCB[,3], c(1 - d2, 2 - d2, 3 - d2), gamma.UCB[,3], code = 3, angle = 90, length = 0.05)

arrows(c(1 + d2, 2 + d2, 3 + d2), gamma.LCB[,4], c(1 + d2, 2 + d2, 3 + d2), gamma.UCB[,4], code = 3, angle = 90, length = 0.05)

mtext("(d)", side = 3, line = 0, at = 1, cex = cexAll, font = 1)

#### Riparian, Persistence

par(pin = c(2, 2), mai = c(.75, 1, .25, .25))

plot(0, bty = "n", type = "n", ylim = f5, xlim = c(.8, 2.2), ylab = "Re-attack rate", xlab = "Between years", axes = FALSE, cex.lab = cexAll)

axis(side = 1, at = c(1, 2), tick = TRUE, labels = c("2009-2010","2010-2011"), cex.axis = cexAll)

axis(side = 2, at = f5s, tick = TRUE, cex.axis = cexAll)

lines(c(1 - d2, 2 - d2), phi.m[,1], type = "b" , cex = cexAll, lwd = 1.5 , lty = 1, col = "black", pch = 16)

lines(c(1 + d2, 2 + d2), phi.m[,2], type = "b" , cex = cexAll, lwd = 1.5 , lty = 2, col = "black", pch = 16)

arrows(c(1 - d2, 2 - d2, 3 - d2), phi.LCB[,1], c(1 - d2, 2 - d2, 3 - d2), phi.UCB[,1], code = 3, angle = 90, length = 0.05)

arrows(c(1 + d2, 2 + d2, 3 + d2), phi.LCB[,2], c(1 + d2, 2 + d2, 3 + d2), phi.UCB[,2], code = 3, angle = 90, length = 0.05)

mtext("(e)", side = 3, line = 0, at = 1, cex = cexAll, font = 1)


#### Upland, Persistence

par(pin = c(2, 2), mai = c(.75, .25, .25, 1))

plot(0, bty = "n", type = "n", ylim = f6, xlim = c(.8, 2.2), ylab = NA, xlab = "Between years", axes = FALSE, cex.lab = cexAll)

axis(side = 1, at = c(1, 2), tick = TRUE, labels = c("2009-2010","2010-2011"), cex.axis = cexAll)

axis(side = 2, at = f6s, tick = TRUE, cex.axis = cexAll)

lines(c(1 - d2, 2 - d2), phi.m[,3], type = "b" , cex = cexAll, lwd = 1.5 , lty = 1, col = "black", pch = 16)

lines(c(1 + d2 , 2 + d2), phi.m[,4], type = "b" , cex = cexAll, lwd = 1.5 , lty = 2, col = "black", pch = 16)

arrows(c(1 - d2, 2 - d2, 3 - d2), phi.LCB[,3], c(1 - d2, 2 - d2, 3 - d2), phi.UCB[,3], code = 3, angle = 90, length = 0.05)

arrows(c(1 + d2, 2 + d2, 3 + d2), phi.LCB[,4], c(1 + d2, 2 + d2, 3 + d2), phi.UCB[,4], code = 3, angle = 90, length = 0.05)

mtext("(f)", side = 3, line = 0, at = 1, cex = cexAll, font = 1)

dev.off()
