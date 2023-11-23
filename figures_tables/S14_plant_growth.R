library(coda)
library(gdata)
library(rjags)
library(plyr)
library(reshape)
library(xtable)

#######################################################################################################################################
# load data

load(file = "../library/salix_borer_wide.rda")

work1 <- ddply(SalixBorerWide, c("plant", "habitat", "sex"), function(x) c(tbd09 = sum(x$tbd09, na.rm = TRUE), tbd10 = sum(x$tbd10, na.rm = TRUE), tbd11 = sum(x$tbd11, na.rm = TRUE)))

work1$ratio <- work1$tbd11/work1$tbd09

i <- which(work1$ratio==0); work1$ratio[i] <- NA

work2 <- subset(work1, is.na(ratio)==FALSE)

DataQuery <- list(
  
  ratio = work2$ratio^.5,
  n = dim(work2)[1],
  habitat = work2$habitat + 1,
  sex = work2$sex + 1)

#######################################################################################################################################

sink("PlantGrowthJAGS.R")

cat("
    
model {
    
# priors    

for (i in 1:2) {
  for (j in 1:2) {
  alpha[i, j] ~ dunif(0, 100)
  beta[i, j] ~ dunif(0, 100)
  mu[i, j] <- alpha[i, j]/beta[i, j]
  }
}

# likelihood

for (i in 1:n) {
  ratio[i] ~ dgamma(alpha[habitat[i], sex[i]], beta[habitat[i], sex[i]])
  ratio.new[i] ~ dgamma(alpha[habitat[i], sex[i]], beta[habitat[i], sex[i]])
  errorRaw[i] <- (ratio[i] - mu[habitat[i], sex[i]])^2
  errorSim[i] <- (ratio.new[i] - mu[habitat[i], sex[i]])^2
}

ratioRawSum <- sum(errorRaw[])
ratioSimSum <- sum(errorSim[])
bpVar <- step(ratioRawSum - ratioSimSum)

sex.riparian <- mu[1, 1] - mu[1, 2]
sex.upland <- mu[2, 1] - mu[2, 2]
habitat.male <- mu[1, 1] - mu[2, 1]
habitat.female <- mu[1, 2] - mu[2, 2]

}    
",fill = TRUE)
sink()

# Initial values

inits <- function() {list(
  alpha = array(runif(2, 0, 5), dim = c(2, 2)),
  beta = array(runif(2, 0, 5), dim = c(2, 2)))}
  
# Parameters monitored

params <- c(
  "mu",
  "alpha",
  "beta",
  "bpVar",
  "sex.riparian",
  "sex.upland",
  "habitat.male",
  "habitat.female")
  
# MCMC settings

n.adapt = 5000
n.update = 10000
n.iter = 15000

# JAGS run

jm = jags.model("PlantGrowthJAGS.R", data = DataQuery, n.chains = 3, inits = inits(), n.adapt = n.adapt)

update(jm, n.iter = n.update)

zm = coda.samples(jm, variable.names = params, n.iter = n.iter, thin = 10)

# summarize JAGS run

summary(zm)

gelman.diag(zm, multivariate = FALSE)

