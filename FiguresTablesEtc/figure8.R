library(rjags)
library(gdata)
library(plyr)

load(file = "../Library/weevilPhenology.rda")

###############################################################################################################################################################################

sink("model.txt")

cat("

model {
	
# priors	

for (i in 1:2){

  for (j in 1:2){

    alpha[i, j] ~ dnorm(0,0.0001)

  }

}

tau <- 1/(sigma*sigma)	

sigma~dunif(0,10)
	
for (i in 1:n){
	
	C[i] ~ dpois(lambda[i])T(1,)
	
  log(lambda[i]) <- alpha[habitat[i], cohort[i]] + eps[i]
		
	eps[i] ~ dnorm(0, tau)

	Presi[i] <- (C[i]-lambda[i])/sqrt(lambda[i])
	
	C.new[i] ~ dpois(lambda[i])T(1,)
	
	Presi.new[i] <- (C.new[i]-lambda[i])/sqrt(lambda[i])
	
	D[i] <- pow(Presi[i],2)

	D.new[i] <- pow(Presi.new[i],2)

}

	zfit.actual <- sum(D[])
	
	zfit.new <- sum(D.new[])


for (i in 1:2){

  for (j in 1:2){

    weevils[i, j] <- exp(alpha[i, j])/(1 - exp(-exp(alpha[i, j])))

  }

}

  dif[1] <- weevils[1,1] - weevils[2,1]
  dif[2] <- weevils[1,2] - weevils[2,2]

}

",fill=TRUE)

sink()

inits <- function(){list(alpha = array(runif(4,-1,1), dim=c(2,2)), sigma=runif(1,3,7))}

params <- c("alpha","sigma","zfit.actual","zfit.new","weevils", "dif")


# MCMC settings
n.adapt = 5000
n.update = 100000
n.iter = 15000

# Call JAGS from R

###############################################################################################################################################################################
# early instars

set.seed(1)

work7 <- subset(weevilPhenology, stage=="early instar")

work8 <- subset(work7, month %in% c(7,9,10) & count > 0)

work8$cohort <- 1

work8$cohort[which(work8$month %in% c(9,10))] <- 2

DataQuery <- list(C=work8$count, habitat=as.numeric(as.factor(work8$habitat)), cohort=work8$cohort, n = length(work8$habitat))

jm = jags.model("model.txt", data = DataQuery, inits = inits(), n.chains = 3, n.adapt = n.adapt)

update(jm, n.iter = n.update)

zm = coda.samples(jm, variable.names = params, n.iter = n.iter, thin = 10)

# summarize JAGS run

summary(zm)

gelman.diag(zm, multivariate = FALSE)

chains2 <- rbind(as.matrix(zm[[1]]), as.matrix(zm[[2]]), as.matrix(zm[[3]]))

fitActual <- chains2[, matchcols(chains2, with = "zfit.act")]

fitNew <- chains2[, matchcols(chains2, with = "zfit.new")]

bppc <- array(0, dim = dim(chains2)[1])

i <- which(fitNew > fitActual)

bppc[i] <- 1

mean(bppc)

minValue <- min(c(fitActual, fitNew))

maxValue <- max(c(fitActual, fitNew))

mv1 <- round_any(minValue, 1, f = floor)

mv2 <- round_any(maxValue, 1, f = ceiling)

par(mar = c(5, 5, 1, 1))

plot(fitActual, fitNew, main = NULL, xlab = expression(T^{obs}), ylab = "", xlim = c(mv1, mv2), ylim = c(mv1, mv2), cex.lab = 1, cex.axis = 1, las = 1)

mtext(expression(T^{rep}), side = 2, line = 4, cex = 1)

abline(0, 1, lwd = 2)

# different based on 95% CI

instar.flush.riparian <- chains2[, 8]
instar.flush.upland <- chains2[, 9]
instar.senescent.riparian <- chains2[, 10]
instar.senescent.upland <- chains2[, 11]

effectHabitat <- instar.flush.riparian - instar.flush.upland
effectCohort <-  instar.flush.riparian - instar.senescent.riparian
interactionHabitatCohort <-  instar.flush.riparian - instar.flush.upland - instar.senescent.riparian + instar.senescent.upland

length(effectHabitat[which(effectHabitat>=1)])/4500
length(effectCohort[which(effectCohort>=1)])/4500
1 - length(interactionHabitatCohort[which(interactionHabitatCohort>=1)])/4500

a <- instar.senescent.upland/instar.senescent.riparian
length(a[which(a>=1)])/4500

###############################################################################################################################################################################
# eggs

set.seed(1)

work7 <- subset(weevilPhenology, stage=="egg")

work8 <- subset(work7, month %in% c(6, 7, 9, 10) & count > 0)

work8$cohort <- 1

work8$cohort[which(work8$month %in% c(9,10))] <- 2

DataQuery <- list(C=work8$count, habitat=as.numeric(as.factor(work8$habitat)), cohort=work8$cohort, n = length(work8$habitat))

jm = jags.model("model.txt", data = DataQuery, inits = inits(), n.chains = 3, n.adapt = n.adapt)

update(jm, n.iter = n.update)

zm = coda.samples(jm, variable.names = params, n.iter = n.iter, thin = 10)

# summarize JAGS run

summary(zm)

gelman.diag(zm, multivariate = FALSE)

chains2 <- rbind(as.matrix(zm[[1]]), as.matrix(zm[[2]]), as.matrix(zm[[3]]))

fitActual <- chains2[, matchcols(chains2, with = "zfit.act")]

fitNew <- chains2[, matchcols(chains2, with = "zfit.new")]

bppc <- array(0, dim = dim(chains2)[1])

i <- which(fitNew > fitActual)

bppc[i] <- 1

mean(bppc)

minValue <- min(c(fitActual, fitNew))

maxValue <- max(c(fitActual, fitNew))

mv1 <- round_any(minValue, 1, f = floor)

mv2 <- round_any(maxValue, 1, f = ceiling)

par(mar = c(5, 5, 1, 1))

plot(fitActual, fitNew, main = NULL, xlab = expression(T^{obs}), ylab = "", xlim = c(mv1, mv2), ylim = c(mv1, mv2), cex.lab = 1, cex.axis = 1, las = 1)

mtext(expression(T^{rep}), side = 2, line = 4, cex = 1)

abline(0, 1, lwd = 2)

# different based on 95% CI

egg.flush.riparian <- chains2[,8]
egg.flush.upland <- chains2[,9]
egg.senescent.riparian <- chains2[,10]
egg.senescent.upland <- chains2[,11]

effectHabitat <- egg.flush.riparian - egg.flush.upland
effectCohort <-  egg.flush.riparian - egg.senescent.riparian
interactionHabitatCohort <-  egg.flush.riparian - egg.flush.upland - egg.senescent.riparian + egg.senescent.upland

length(effectHabitat[which(effectHabitat>=1)])/4500
length(effectCohort[which(effectCohort>=1)])/4500
length(interactionHabitatCohort[which(interactionHabitatCohort>=1)])/4500

b <- egg.senescent.upland/egg.senescent.riparian
length(b[which(b>=1)])/4500

###############################################################################################################################################################################

ME.instar <- c(mean(instar.flush.riparian), mean(instar.flush.upland), 0, mean(instar.senescent.riparian), mean(instar.senescent.upland))

ME.eggs <- c(mean(egg.flush.riparian), mean(egg.flush.upland), 0, mean(egg.senescent.riparian), mean(egg.senescent.upland))

LCB.instar <- c(quantile(instar.flush.riparian, prob = .025), quantile(instar.flush.upland, prob = .025), 0, quantile(instar.senescent.riparian, prob = .025), quantile(instar.senescent.upland, prob = .025))

LCB.eggs <- c(quantile(egg.flush.riparian, prob = .025), quantile(egg.flush.upland, prob = .025), 0, quantile(egg.senescent.riparian, prob = .025), quantile(egg.senescent.upland, prob = .025))

UCB.instar <- c(quantile(instar.flush.riparian, prob = .975), quantile(instar.flush.upland, prob = .975), 0, quantile(instar.senescent.riparian, prob = .975), quantile(instar.senescent.upland, prob = .975))

UCB.eggs <- c(quantile(egg.flush.riparian, prob = .975), quantile(egg.flush.upland, prob = .975), 0, quantile(egg.senescent.riparian, prob = .975), quantile(egg.senescent.upland, prob = .975))

###############################################################################################################################################################################

pdf(file = "fig8.pdf", width = 6, height = 6)

par(mfrow=c(2,1))

cexAll <- .79

par(pin = c(3, 3), mai = c(.25, 1, 1, .25))

barplot(ME.instar, axis.lty=0, xpd=FALSE, xaxt="n", ylim=c(0,30), ylab=expression('Early instars per stem'), cex.names= cexAll,  cex.lab=cexAll, cex.axis=cexAll, space = 0,

        col = c("white","white","white"))

abline(h=0,lty=1,lwd=1)

arrows(c(1:5)*1-0.5,LCB.instar, c(1:5)*1-0.5,UCB.instar, code=3, angle=90, length=0.05)

#text(x=4, y = UCB.instar[5]+1, "*", cex = 2)

mtext("(a)", side = 3, line = 0, at = .25, cex = cexAll, font = 2)

par(pin = c(3, 3), mai = c(1, 1, .25, .25))

barplot(ME.eggs, axis.lty=0, xpd=FALSE, xaxt="n", ylim=c(0,14), ylab=expression('Eggs per stem'), cex.names= cexAll,  cex.lab=cexAll, cex.axis=cexAll, space = 0,
        
        col = c("white","white","white"))

abline(h=0,lty=1,lwd=1)

arrows(c(1:5)*1-0.5,LCB.eggs, c(1:5)*1-0.5,UCB.eggs, code=3, angle=90, length=0.05)

axis(side=1, at=.5:4.5, labels=c("Riparian", "Upland", "", "Riparian", "Upland"), tick=FALSE, cex.axis = cexAll)

mtext("Flush-feeders", side = 1, line = 3 , cex = cexAll, at = 1)
mtext("Senescent-feeders", side = 1, line = 3 , cex = cexAll, at = 4)

#text(x=4, y = UCB.eggs[5]+1, "*", cex = 2)

mtext("(b)", side = 3, line = 0, at = .2, cex = cexAll, font = 2)

dev.off()

###############################################################################################################################################################################
