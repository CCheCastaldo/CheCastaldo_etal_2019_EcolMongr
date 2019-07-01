library(gdata)
library(reshape)
library(plyr)
library(rjags)

##########################################################################################################
# plot and bundle data

load(file = "../Library/pairedStemMeasurements.rda")

plot(pairedStemMeasurements[, 1], pairedStemMeasurements[, 2], pch = 16, cex = 1.5, col = "black")

data <- list(y = pairedStemMeasurements, n = dim(pairedStemMeasurements)[1])

##########################################################################################################
# model the error as truncated normal

sink("TruncatedNormalError.R")

cat("
    
model {
    
# priors    
  sigma ~ dunif(0, 100)
  tau <- 1/pow(sigma, 2)  
    
  for (i in 1:n) {  
    z[i] ~ dunif(0, 200)
  }
    
# likelihood
  for (i in 1:n) {
    y.pred[i] ~ dnorm(z[i], tau)T(0,)
    for (j in 1:2) {
      y[i,j] ~ dnorm(z[i], tau)T(0,)
      y.new[i,j] ~ dnorm(z[i], tau)T(0,)
      sq.error.dat[i,j] <- (y[i,j] - z[i])^2
      sq.error.sim[i,j] <- (y.new[i,j] - z[i])^2
      }
    }
    
# derived quantities
  p.discrep <- step(sum(sq.error.sim) - sum(sq.error.dat))

}    
",fill = TRUE)
sink()


# Initial values
inits <- function() {list(z = pairedStemMeasurements[,1], sigma = runif(1,0,2))} 

# Parameters monitored
params1 <- c("sigma", "p.discrep")
params2 <- c("z")
params3 <- c("y.pred")

# MCMC settings
n.adapt = 5000
n.update = 20000
n.iter = 20000

jm1 = jags.model("TruncatedNormalError.R", data = data, inits = inits(), n.chains = 3, n.adapt = n.adapt)

update(jm1, n.iter = n.update)

zm1 = coda.samples(jm1, variable.names = params1, n.iter = n.iter, thin = 5)

zm2 = coda.samples(jm1, variable.names = params3, n.iter = n.iter, thin = 5)

summary(zm1)

gelman.diag(zm1, multivariate=FALSE)

densplot(zm1)

##########################################################################################################
# plot it

cols <- matchcols(zm2[[1]], with = "y.pred\\[")

a <- as.data.frame(summary(zm2[, cols])[[2]])

b <- as.data.frame(seq(1, 106))

plotData <- cbind(a, b, pairedStemMeasurements)

colnames(plotData) <- c("q25", "q25", "q50", "q75", "q975","stem","bdObs1","bdObs2")

plotData <- plotData[order(plotData$q50),]

xrange = c(0, 60)

stems <- dim(plotData)[1]

col = "black"

cex = .8

par(mar = c(5, 5, 1, 1))

plot(y = xrange, x = c(1, stems), type = "n", xlab = "stem", ylab = "BD", axes = FALSE, cex.lab = cex, cex.axis = cex)

axis(2, cex.lab = cex, cex.axis = cex)

axis(1, at = 1:stems, labels = plotData[,6], padj = 0.5, las = 1, cex.lab = cex, cex.axis = cex)

segments(y0 = plotData[, 1], y1 = plotData[, 5], x0 = 1:length(cols), x1 = 1:length(cols), col = col)

segments(y0 = plotData[, 2], y1 = plotData[, 4], x0 = 1:length(cols), x1 = 1:length(cols), lwd = 3, col = col)

points((1:stems), plotData[, 3], col = col)

points((1:stems), plotData[, 7], col = "red")

points((1:stems), plotData[, 8], col = "red")

