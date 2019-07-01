library(coda)
library(gdata)
library(rjags)
library(plyr)
library(reshape)
library(xtable)
library(truncdist)

version <- "i-009f9c75a9206a0ac"
  
#######################################################################################################################################
# load data
  
load(file = "../Library/SalixBorerWide.rda")
  
load(file = paste0("../ModelBuild/GlobalModel/", version, "/MCMCsamplesZ.rda"))
  
chains <- apply(rbind(as.matrix(MCMCsamplesZ[[1]]), as.matrix(MCMCsamplesZ[[2]]), as.matrix(MCMCsamplesZ[[3]])), 2, mean)
  
zWeights <- as.data.frame(cbind(chains[1:2087], chains[2088:4174], chains[4175:6261]))
  
SalixBorerWide <- cbind(SalixBorerWide, zWeights)
  
#######################################################################################################################################
# define weights and create data object for jags
  
work1 <- subset(SalixBorerWide, site10==1, select = c(habitat, sex, repro10, tbd10, V1))
  
work2 <- subset(SalixBorerWide, site11==1, select = c(habitat, sex, repro11, tbd11, V2))
  
names(work1) <- c("habitat", "sex", "repro", "tbd", "w")                          
  
names(work2) <- c("habitat", "sex", "repro", "tbd", "w")                          
  
work1$year <- 1
  
work2$year <- 2
  
work3 <- rbind(work1, work2)
  
work3$check <- 0; i <- which(work3$w <= .1 | work3$w >= .9); work3$check[i] <- 1

table(work3$check)/dim(work3)[1]

work4 <- subset(work3, w <= .1)

DataQuery <- list(
  
  tbd = work4$tbd,
  habitat = work4$habitat + 1,
  sex = work4$sex + 1,
  year = work4$year,
  n = dim(work4)[1])

(work5 <- ddply(work4, c("habitat", "sex", "year"), function(x) c(
  
  meanTBD = mean(x$tbd, na.rm = TRUE),
  
  sdeTBD = sd(x$tbd, na.rm = TRUE))))

#######################################################################################################################################
# do weighted regression on diameter
  
sink("StemDiameterJAGS.R")
  
cat("
      
model {
      
# priors    
      
for (i in 1:2) {
  for (j in 1:2) {
    for (k in 1:2) {
      
      sigma[i, j, k] ~ dunif(0, 5)
      tau[i, j, k] <- pow(sigma[i, j, k], -2) 
      mu[i, j, k] ~ dunif(0, 5)    

    }
  }
}
      
# likelihood
      
for (i in 1:n) {
      
  tbd[i] ~ dlnorm(mu[habitat[i], sex[i], year[i]], tau[habitat[i], sex[i], year[i]])T(10, )
  tbdSim[i] ~ dlnorm(mu[habitat[i], sex[i], year[i]], tau[habitat[i], sex[i], year[i]])T(10, )
  errorRaw[i] <- (tbd[i] - exp(mu[habitat[i], sex[i], year[i]]))^2
  errorSim[i] <- (tbdSim[i] - exp(mu[habitat[i], sex[i], year[i]]))^2
  
}
      
# derived quantites

for (i in 1:2) {
  for (j in 1:2) {
    for (k in 1:2) {

      tbdNew[i, j, k] ~ dlnorm(mu[i, j, k], tau[i, j, k])T(10, )

    }
  }
}

# ppc
      
tbdRawSum <- sum(errorRaw[])
tbdSimSum <- sum(errorSim[])
tbdRawM <- mean(tbd[])
tbdSimM <- mean(tbdSim[])
bpVar <- step(tbdRawSum - tbdSimSum)
bpMean <- step(tbdRawM - tbdSimM)
 
}    
",fill = TRUE)
sink()

  
# Initial values
  
initsFinal <- list(
  
  list(
    
    mu = array (c(3.07, 0.02794,	3.086, 0.01808, 3.148, 0.09498, 3.19, 0.395), dim = c(2,2,2)),
    sigma = array(c(0.43, 0.70, 0.48, 0.92, 0.48, 0.47, 0.55, 0.49), dim = c(2,2,2))),
  
  list(
    
    mu = array (c(3.17, 0.6282, 3.206, 0.4663, 3.266, 1.453, 3.317, 1.966), dim = c(2,2,2)),
    sigma = array(c(0.49, 0.99, 0.55, 1.16, 0.55, 0.76, 0.62, 0.69), dim = c(2,2,2))),
  
  list(
    
    mu = array (c(3.25, 1.81, 3.30, 1.47, 3.36, 2.35, 3.42, 2.45), dim = c(2,2,2)),
    sigma = array(c(0.574, 1.182, 0.6508, 1.321, 0.6479, 1.086, 0.728, 1.093), dim = c(2,2,2))))
  
 
# Parameters monitored
  
params <- c(
  "mu", 
  "sigma", 
  "bpVar",
  "bpMean",
  "tbdRawSum", 
  "tbdSimSum", 
  "tbdRawM", 
  "tbdSimM",
  "tbdNew")

# MCMC settings
  
n.adapt = 5000
n.update = 20000
n.iter = 15000
  
# JAGS run
  
jm = jags.model("StemDiameterJAGS.R", data = DataQuery, inits = initsFinal, n.chains = 3, n.adapt = n.adapt)
  
update(jm, n.iter = n.update)
  
zm = coda.samples(jm, variable.names = params, n.iter = n.iter, thin = 5)

# summarize JAGS run

summary(zm)

gelman.diag(zm, multivariate = FALSE)

#######################################################################################################################################
# ppc plot

chains2 <- rbind(as.matrix(zm[[1]]), as.matrix(zm[[2]]), as.matrix(zm[[3]]))

fitActual <- chains2[, matchcols(chains2, with = "tbdRawSum")]

fitNew <- chains2[, matchcols(chains2, with = "tbdSimSum")]

bppc <- array(0, dim = dim(chains2)[1])

i <- which(fitNew > fitActual)

bppc[i] <- 1

minValue <- min(c(fitActual, fitNew))

maxValue <- max(c(fitActual, fitNew))

mv1 <- round_any(minValue, 1, f = floor)

mv2 <- round_any(maxValue, 1, f = ceiling)

par(mar = c(5, 5, 1, 1))

plot(fitActual, fitNew, main = NULL, xlab = expression(T^{obs}), ylab = "", xlim = c(mv1, mv2), ylim = c(mv1, mv2), cex.lab = 1, cex.axis = 1, las = 1)

mtext(expression(T^{rep}), side = 2, line = 4, cex = 1)

abline(0, 1, lwd = 2)

#######################################################################################################################################
# calculate adjusted means and do mean comparisons

chains3 <- chains2[, 3:18]

mu <- array(NA, dim = c(dim(chains3)[1], 8))

for (i in 1:dim(chains3)[1]){
  for (j in 1:8){
  
    mu[i, j] <- extrunc("lnorm", chains3[i, j], chains3[i, (j + 8)], a = 10, b = Inf)
    
  }
}

# year comparison, male, riparian
dif <- mu[,5] - mu[,1]
length(dif[which(dif>0)])/dim(chains3)[1]

# year comparison, upland, male
dif <- mu[,6] - mu[,2]
length(dif[which(dif>0)])/dim(chains3)[1]

# year comparison, female, riparian
dif <- mu[,7] - mu[,3]
length(dif[which(dif>0)])/dim(chains3)[1]

# year comparison, female, upland
dif <- mu[,8] - mu[,4]
length(dif[which(dif>0)])/dim(chains3)[1]


# sex comparison, 2010, riparian
dif <- mu[,3] - mu[,1]
length(dif[which(dif>0)])/dim(chains3)[1]

# sex comparison, 2010, upland
dif <- mu[,4] - mu[,2]
length(dif[which(dif>0)])/dim(chains3)[1]

# sex comparison, 2011, riparian
dif <- mu[,7] - mu[,5]
length(dif[which(dif>0)])/dim(chains3)[1]

# sex comparison, 2011, upland
dif <- mu[,8] - mu[,6]
length(dif[which(dif>0)])/dim(chains3)[1]

#######################################################################################################################################
# make table

Dtable <- array(NA, dim = c(4, 6))

Dtable[1, 1] <- median(mu[, 1])
Dtable[2, 1] <- median(mu[, 3])
Dtable[3, 1] <- median(mu[, 2])
Dtable[4, 1] <- median(mu[, 4])

Dtable[1, 2] <- quantile(mu[, 1], prob = .025)
Dtable[2, 2] <- quantile(mu[, 3], prob = .025)
Dtable[3, 2] <- quantile(mu[, 2], prob = .025)
Dtable[4, 2] <- quantile(mu[, 4], prob = .025)

Dtable[1, 3] <- quantile(mu[, 1], prob = .975)
Dtable[2, 3] <- quantile(mu[, 3], prob = .975)
Dtable[3, 3] <- quantile(mu[, 2], prob = .975)
Dtable[4, 3] <- quantile(mu[, 4], prob = .975)

Dtable[1, 4] <- median(mu[, 5])
Dtable[2, 4] <- median(mu[, 7])
Dtable[3, 4] <- median(mu[, 6])
Dtable[4, 4] <- median(mu[, 8])

Dtable[1, 5] <- quantile(mu[, 5], prob = .025)
Dtable[2, 5] <- quantile(mu[, 7], prob = .025)
Dtable[3, 5] <- quantile(mu[, 6], prob = .025)
Dtable[4, 5] <- quantile(mu[, 8], prob = .025)

Dtable[1, 6] <- quantile(mu[, 5], prob = .975)
Dtable[2, 6] <- quantile(mu[, 7], prob = .975)
Dtable[3, 6] <- quantile(mu[, 6], prob = .975)
Dtable[4, 6] <- quantile(mu[, 8], prob = .975)

habitat <- c("Riparian", "Riparian", "Upland", "Upland")

sex <- c("Male", "Female", "Male", "Female")

CI2010 <- paste0("(", round(Dtable[,2], 1), " -- ", round(Dtable[,3], 1), ")")

CI2011 <- paste0("(", round(Dtable[,5], 1), " -- ", round(Dtable[,6], 1), ")")

Dtable2 <- as.data.frame(cbind(habitat, sex, round(Dtable[,1],1), CI2010, round(Dtable[,4],1), CI2011))

xtable(Dtable2, digits = 1)

#######################################################################################################################################
# posterior predictive distribution for tbd

params <- c("tbdNew")

zm = coda.samples(jm, variable.names = params, n.iter = n.iter, thin = 5)

summary(zm)

MCMCsamplesBDGamma <- zm

save(MCMCsamplesBDGamma, file = paste0("GlobalModel/", version, "/MCMCsamplesBDGamma.rda"))

unlist("StemDiamaterJAGS.R")
