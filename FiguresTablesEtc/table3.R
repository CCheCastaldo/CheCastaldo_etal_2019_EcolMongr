library(coda)
library(gdata)
library(rjags)
library(broman)
library(plyr)
library(reshape)
library(xtable)

version <- "i-009f9c75a9206a0ac"

#######################################################################################################################################
# load data

load(file = "../Library/SalixBorerWide.rda")

load(file = paste0("../ModelBuild/GlobalModel/", version, "/MCMCsamplesZ.rda"))

z1 <- t(rbind(as.matrix(MCMCsamplesZ[, unlist(matchcols(MCMCsamplesZ[[1]], with = c("z\\[", "1\\]")))])))

z2 <- t(rbind(as.matrix(MCMCsamplesZ[, unlist(matchcols(MCMCsamplesZ[[1]], with = c("z\\[", "2\\]")))])))

z3 <- t(rbind(as.matrix(MCMCsamplesZ[, unlist(matchcols(MCMCsamplesZ[[1]], with = c("z\\[", "3\\]")))])))

site <- subset(SalixBorerWide, select = c(stem, site09, site10, site11, habitat, sex, repro10, repro11))

#######################################################################################################################################

w1 <- site$site10 * z1

w2 <- site$site11 * z2

j1 <- which(site$habitat==0 & site$sex==0)

j2 <- which(site$habitat==0 & site$sex==1)

j3 <- which(site$habitat==1 & site$sex==0)

j4 <- which(site$habitat==1 & site$sex==1)

chains <- chains2 <- array(NA, dim = c(4500, 8))

for (i in 1:dim(w1)[2]) {
  
  chains[i, 1] <- sum(w1[j1, i])/sum(site$site10[j1])
  chains[i, 2] <- sum(w1[j2, i])/sum(site$site10[j2])
  chains[i, 3] <- sum(w1[j3, i])/sum(site$site10[j3])
  chains[i, 4] <- sum(w1[j4, i])/sum(site$site10[j4])
  
  chains[i, 5] <- sum(w2[j1, i])/sum(site$site11[j1])
  chains[i, 6] <- sum(w2[j2, i])/sum(site$site11[j2])
  chains[i, 7] <- sum(w2[j3, i])/sum(site$site11[j3])
  chains[i, 8] <- sum(w2[j4, i])/sum(site$site11[j4])

  chains2[i, 1] <- 1-chains[i, 1]
  chains2[i, 2] <- 1-chains[i, 2]
  chains2[i, 3] <- 1-chains[i, 3]
  chains2[i, 4] <- 1-chains[i, 4]
  
  chains2[i, 5] <- 1-chains[i, 5]
  chains2[i, 6] <- 1-chains[i, 6]
  chains2[i, 7] <- 1-chains[i, 7]
  chains2[i, 8] <- 1-chains[i, 8]
  
}

#######################################################################################################################################

z.m <- z.LCB <- z.UCB <- matrix(NA, nrow = 4, ncol = 2)

z.m[1, 1] <- median(chains[, 1])
z.m[2, 1] <- median(chains[, 2])
z.m[3, 1] <- median(chains[, 3])
z.m[4, 1] <- median(chains[, 4])
z.m[1, 2] <- median(chains[, 5])
z.m[2, 2] <- median(chains[, 6])
z.m[3, 2] <- median(chains[, 7])
z.m[4, 2] <- median(chains[, 8])

z.LCB[1, 1] <- quantile(chains[, 1], prob = .025)
z.LCB[2, 1] <- quantile(chains[, 2], prob = .025)
z.LCB[3, 1] <- quantile(chains[, 3], prob = .025)
z.LCB[4, 1] <- quantile(chains[, 4], prob = .025)
z.LCB[1, 2] <- quantile(chains[, 5], prob = .025)
z.LCB[2, 2] <- quantile(chains[, 6], prob = .025)
z.LCB[3, 2] <- quantile(chains[, 7], prob = .025)
z.LCB[4, 2] <- quantile(chains[, 8], prob = .025)

z.UCB[1, 1] <- quantile(chains[, 1], prob = .975)
z.UCB[2, 1] <- quantile(chains[, 2], prob = .975)
z.UCB[3, 1] <- quantile(chains[, 3], prob = .975)
z.UCB[4, 1] <- quantile(chains[, 4], prob = .975)
z.UCB[1, 2] <- quantile(chains[, 5], prob = .975)
z.UCB[2, 2] <- quantile(chains[, 6], prob = .975)
z.UCB[3, 2] <- quantile(chains[, 7], prob = .975)
z.UCB[4, 2] <- quantile(chains[, 8], prob = .975)

habitat <- c("Riparian", "Riparian", "Upland", "Upland")

sex <- c("Male", "Female", "Male", "Female")

CI2010 <- paste0("(", round(z.LCB[,1], 2), " -- ", round(z.UCB[,1], 2), ")")

CI2011 <- paste0("(", round(z.LCB[,2], 2), " -- ", round(z.UCB[,2], 2), ")")

Dtable2 <- as.data.frame(cbind(habitat, sex, round(z.m[,1],2), CI2010, round(z.m[,2],2), CI2011))

xtable(Dtable2, digits = 1)

