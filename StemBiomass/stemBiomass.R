library(gdata)
library(reshape)
library(plyr)
library(rjags)

##########################################################################################################################################################################
# input data

load(file = "../Library/stemBiomass.rda")

DataQuery <- list(
  
  bd = log(stemBiomass$bd),
  soma = log(stemBiomass$soma),
  bd.new = log(c(12, 15, 20, 30, 35)),
  n = dim(stemBiomass)[1],
  n.new = 5)

##########################################################################################################################################################################
# model biomass

sink("StemBiomassJAGS.R")

cat("
    
model {

# priors

a ~ dnorm(0, .0001)    	
b ~ dnorm(0, .0001)    	
sigma ~ dunif(0, 10)
tau <- 1/(sigma * sigma)
    
# likelihood

for (i in 1:n) {
  soma[i] ~ dnorm(mu[i], tau) 
  y.new[i] ~ dnorm(mu[i], tau)
  sq.error.dat[i] <- (soma[i] - mu[i])^2
  sq.error.sim[i] <- (y.new[i] - mu[i])^2
  mu[i] <- a + b * bd[i]
}

for (i in 1:n.new) {
  soma.new[i] ~ dnorm(mu.new[i], tau)

  mu.new[i] <- a + b * bd.new[i]
  soma.pred[i] <- exp(soma.new[i])
}

# derived quantities

p.discrep <- step(sum(sq.error.sim) - sum(sq.error.dat))
n1 <- b - 1

} 
    
",fill = TRUE)
sink()

# Initial values

inits <- function() {list(

  a = runif(1,-5,5), 
  b = runif(1,-5,5), 
  sigma = runif(1,1,10))} 

# Parameters monitored

params <- c(
  "a", 
  "b", 
  "sigma", 
  "n1", 
  "soma.pred",
  "p.discrep")

# MCMC settings

n.adapt = 5000
n.update = 10000
n.iter = 15000

# JAGS run

jm = jags.model("StemBiomassJAGS.R", data = DataQuery, inits = inits(), n.chains = 3, n.adapt = n.adapt)

update(jm, n.iter = n.update)

zm = coda.samples(jm, variable.names = params, n.iter = n.iter, thin = 10)

# summarize JAGS run

summary(zm)

gelman.diag(zm, multivariate = FALSE)

stemBiomass <- zm

save(stemBiomass, file = "MCMCstemBiomass.rda")

unlist("StemBiomassJAGS.R")
