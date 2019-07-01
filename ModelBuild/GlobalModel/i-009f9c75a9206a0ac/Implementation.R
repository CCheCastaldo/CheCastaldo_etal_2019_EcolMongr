library(rjags)
library(gdata)
library(coda)

load(file = "../../../Library/SalixBorerWide.rda")

work1 <- SalixBorerWide

# habitat = 0: riparian
# habitat = 1: upland
# sex = 0: male
# sex = 1: female

#####################################################################################################################################################
# group variables

work1$group.stem <- NA

i <- which(work1$habitat==0 & work1$sex==0); work1$group.stem[i] <- 1

i <- which(work1$habitat==0 & work1$sex==1); work1$group.stem[i] <- 2

i <- which(work1$habitat==1 & work1$sex==0); work1$group.stem[i] <- 3

i <- which(work1$habitat==1 & work1$sex==1); work1$group.stem[i] <- 4

# make unique list of plant indicator in same order as main data frame

work2 <- subset(work1, select = c(plant, habitat, group.stem))

work3 <- unique(work2)

#####################################################################################################################################################
# covariates

group.stem <- as.numeric(work1$group.stem)

group.plant <- as.numeric(work3$group.stem)

habitat.plant <- as.numeric(work3$habitat) + 1

habitat <- as.numeric(work1$habitat) + 1

sex <- as.numeric(work1$sex) + 1

plant <- as.numeric(work1$plant)

site <- as.matrix(subset(work1, select = c(site09, site10, site11)))

d <- as.matrix(subset(work1, select = c(tbd09.sd, tbd10.sd, tbd11.sd)))

repro <- as.matrix(subset(work1, select = c(repro10, repro11)))

rgr <- as.matrix(subset(work1, select = c(gr09.sd, gr10.sd)))

#####################################################################################################################################################
# convert NA to 0 so JAGS doesn't choke on site = 0 cases

i <- which(site[, 1]==0) 
  
d[i, 1] <- 0
  
i <- which(site[, 2]==0) 
  
d[i, 2] <- repro[i, 1] <- rgr[i, 1] <- 0
  
i <- which(site[, 3]==0) 
  
d[i, 3] <- repro[i, 2] <- rgr[i, 2] <- 0
  
#####################################################################################################################################################
# dimensions

j <- dim(work1)[1]
s <- 2 
t <- 3 
nplant <- length(unique(work1$plant))
U <- dim(subset(work1, habitat==0))[1] + 1

#####################################################################################################################################################
# response

y <- array(NA, dim = c(j, s, t))	

y[,1,1] <- as.numeric(work1$borer1.09)
y[,2,1] <- as.numeric(work1$borer2.09)
y[,1,2] <- as.numeric(work1$borer1.10)
y[,2,2] <- as.numeric(work1$borer2.10)
y[,1,3] <- as.numeric(work1$borer1.11)
y[,2,3] <- as.numeric(work1$borer2.11)

#####################################################################################################################################################
# data

jags.data <- list(

y = y,
nstem = j, 
nsurvey = s, 
nyear = t, 
nplant = nplant, 
U = U,
site = site,
plant = plant,
repro = repro, 
d = d, 
rgr = rgr,
group.stem = group.stem,
group.plant = group.plant,
habitat.plant = habitat.plant,
habitat = habitat,
sex = sex)

#####################################################################################################################################################
# Specify model in BUGS language

sink("weeviloccupancy.jags")

cat("	

model {

# priors

for (i in 1:3) { alpha.p[i] ~ dunif(-10, 10) }

for (i in 1:3) { beta.d.p[i] ~ dnorm(0, .386) }

for (i in 1:4) { beta.repro.phi[i] ~ dnorm(0, .386) }

for (i in 1:4) { beta.repro.gamma[i] ~ dnorm(0, .386) }

for (i in 1:2) { beta.d.psi[i] ~ dunif(-10, 10) }

for (i in 1:2) { beta.d.phi[i] ~ dunif(-10, 10) }

for (i in 1:2) { beta.d.gamma[i] ~ dunif(-10, 10) }

for (i in 1:2) { beta.rgr.phi[i] ~ dnorm(0, .386) }

for (i in 1:2) { beta.rgr.gamma[i] ~ dnorm(0, .386) }

for (i in 1:4) { mu.psi[i] ~ dunif(-10, 10) }
  
for (i in 1:4) { mu.phi[i] ~ dunif(-10, 10) }

for (i in 1:4) { mu.gamma[i] ~ dunif(-10, 10) }
  
for (i in 1:2) { sigma.psi[i] ~ dunif(0, 10) }
  
sigma.phi ~ dunif(0, 10)

sigma.gamma ~ dunif(0, 10)

for (i in 1:2) { tau.psi[i] <- 1/(sigma.psi[i] * sigma.psi[i]) }

tau.phi <- 1/(sigma.phi * sigma.phi) 

tau.gamma <- 1/(sigma.gamma * sigma.gamma)


# process model

for (i in 1:nplant) { 

  alpha.psi[i] ~ dnorm(mu.psi[group.plant[i]], tau.psi[habitat.plant[i]])           

  alpha.phi[i] ~ dnorm(mu.phi[group.plant[i]], tau.phi)           

  alpha.gamma[i] ~ dnorm(mu.gamma[group.plant[i]], tau.gamma)  

} #i

for (j in 1:(U-1)) {
  	
	z[j, 1] ~ dbern(psi[j, 1] * site[j, 1])

	logit(psi[j, 1]) <- alpha.psi[plant[j]] + beta.d.psi[1] * d[j, 1]

} #j

for (j in U:nstem) {
  	
  z[j, 1] ~ dbern(psi[j, 1] * site[j, 1])

  logit(psi[j, 1]) <- alpha.psi[plant[j]] + beta.d.psi[2] * d[j, 1]

} #j

for (j in 1:nstem) {
 
 	for (t in 2:nyear) {
  
    	z[j, t] ~ dbern(psi[j, t] * site[j, t])

    	psi[j, t] <- z[j, t - 1] * phi[j, t - 1] + (1 - z[j, t - 1]) * gamma[j, t - 1]

    	logit(gamma[j, t - 1]) <- alpha.gamma[plant[j]] + beta.repro.gamma[group.stem[j]] * repro[j, t - 1] + beta.rgr.gamma[habitat[j]] * rgr[j, t - 1] + beta.d.gamma[habitat[j]] * d[j, t]        

    	logit(phi[j, t - 1]) <- alpha.phi[plant[j]] + beta.repro.phi[group.stem[j]] * repro[j, t - 1] + beta.rgr.phi[habitat[j]] * rgr[j, t - 1] + beta.d.phi[habitat[j]] * d[j, t]

  	} #t

} #j

 
# observation model

# all the riparian stems

for (j in 1:(U - 1)) {

  for (s in 1:nsurvey) {

    for (t in 1:nyear) {

      y[j, s, t] ~ dbern(muy[j, s, t] * site[j, t])

      muy[j, s, t] <- z[j, t] * p[j, s, t]

      logit(p[j, s, t]) <- alpha.p[3] + beta.d.p[3] * d[j, t]

      y.new[j, s, t] ~ dbern(muy[j, s, t] * site[j, t])

    } #t

  } #s

} #j

# the upland stems in 2009

for (j in U:nstem) {

	for (s in 1:nsurvey) {
   
    	for (t in 1:1) {
     
        y[j, s, t] ~ dbern(muy[j, s, t] * site[j, t])

 		    muy[j, s, t] <- z[j, t] * p[j, s, t]
       
       	logit(p[j, s, t]) <- alpha.p[1] + beta.d.p[1] * d[j, t]
         
       	y.new[j, s, t] ~ dbern(muy[j, s, t] * site[j, t])

	    } #t
    
      #  the upland stems in 2010 and 2011

	  	for (t in 2:nyear) {
     
        y[j, s, t] ~ dbern(muy[j, s, t] * site[j, t])

 		    muy[j, s, t] <- z[j, t] * p[j, s, t]

        logit(p[j, s, t]) <- alpha.p[2] + beta.d.p[2] * d[j, t]

        y.new[j, s, t] ~ dbern(muy[j, s, t] * site[j, t])

      } #t

  } #s

} #j


# ppc

for (j in 1:nstem) {

  for (t in 1:nyear) {

    sum.y[j, t] <- max(0.01, sum(y[j,, t]))
    
	  eval[j, t] <- max(0.01, sum(muy[j,, t]))
    
	  E[j, t] <- pow((sum.y[j, t] - eval[j, t]), 2) / (eval[j, t] + 0.01) 
    
    sum.y.new[j, t] <- max(0.01, sum(y.new[j,, t]))
    
    E.new[j, t] <- pow((sum.y.new[j, t] - eval[j, t]), 2) / (eval[j, t] + 0.01)

  } #t

} #j

zActual <- sum(E[,])

zNew <- sum(E.new[,])


# metapopulation derived quantities

for (j in 1:nstem) {
  
  m.r.psi[j, 1] <- z[j, 1] * ifelse(habitat[j]==1, 1, 0) * ifelse(sex[j]==1, 1, 0) * site[j, 1]

  m.r.psi.sites[j, 1] <- ifelse(habitat[j]==1, 1, 0) * ifelse(sex[j]==1, 1, 0) * site[j, 1]

  f.r.psi[j, 1] <- z[j, 1] * ifelse(habitat[j]==1, 1, 0) * ifelse(sex[j]==2, 1, 0) * site[j, 1]

  f.r.psi.sites[j, 1] <- ifelse(habitat[j]==1, 1, 0) * ifelse(sex[j]==2, 1, 0) * site[j, 1]

  m.u.psi[j, 1] <- z[j, 1] * ifelse(habitat[j]==2, 1, 0) * ifelse(sex[j]==1, 1, 0) * site[j, 1]

  m.u.psi.sites[j, 1] <- ifelse(habitat[j]==2, 1, 0) * ifelse(sex[j]==1, 1, 0) * site[j, 1]

  f.u.psi[j, 1] <- z[j, 1] * ifelse(habitat[j]==2, 1, 0) * ifelse(sex[j]==2, 1, 0) * site[j, 1]

  f.u.psi.sites[j, 1] <- ifelse(habitat[j]==2, 1, 0) * ifelse(sex[j]==2, 1, 0) * site[j, 1]

  for (t in 2:nyear) {

    m.r.psi[j, t] <- z[j, t] * ifelse(habitat[j]==1, 1, 0) * ifelse(sex[j]==1, 1, 0) * site[j, t]

    m.r.psi.sites[j, t] <- ifelse(habitat[j]==1, 1, 0) * ifelse(sex[j]==1, 1, 0) * site[j, t]

    f.r.psi[j, t] <- z[j, t] * ifelse(habitat[j]==1, 1, 0) * ifelse(sex[j]==2, 1, 0) * site[j, t]

    f.r.psi.sites[j, t] <- ifelse(habitat[j]==1, 1, 0) * ifelse(sex[j]==2, 1, 0) * site[j, t]

    m.u.psi[j, t] <- z[j, t] * ifelse(habitat[j]==2, 1, 0) * ifelse(sex[j]==1, 1, 0) * site[j, t]

    m.u.psi.sites[j, t] <- ifelse(habitat[j]==2, 1, 0) * ifelse(sex[j]==1, 1, 0) * site[j, t]

    f.u.psi[j, t] <- z[j, t] * ifelse(habitat[j]==2, 1, 0) * ifelse(sex[j]==2, 1, 0) * site[j, t]

    f.u.psi.sites[j, t] <- ifelse(habitat[j]==2, 1, 0) * ifelse(sex[j]==2, 1, 0) * site[j, t]

    m.r.gamma[j, t - 1] <- z[j, t] * (1 - z[j, t - 1]) * ifelse(habitat[j]==1, 1, 0) * ifelse(sex[j]==1, 1, 0) * site[j, t]

    m.r.gamma.sites[j, t - 1] <- (1 - z[j, t - 1]) * ifelse(habitat[j]==1, 1, 0) * ifelse(sex[j]==1, 1, 0) * site[j, t]

    f.r.gamma[j, t - 1] <- z[j, t] * (1 - z[j, t - 1]) * ifelse(habitat[j]==1, 1, 0) * ifelse(sex[j]==2, 1, 0) * site[j, t]

    f.r.gamma.sites[j, t - 1] <- (1 - z[j, t - 1]) * ifelse(habitat[j]==1, 1, 0) * ifelse(sex[j]==2, 1, 0) * site[j, t]

    m.u.gamma[j, t - 1] <- z[j, t] * (1 - z[j, t - 1]) * ifelse(habitat[j]==2, 1, 0) * ifelse(sex[j]==1, 1, 0) * site[j, t]

    m.u.gamma.sites[j, t - 1] <- (1 - z[j, t - 1]) * ifelse(habitat[j]==2, 1, 0) * ifelse(sex[j]==1, 1, 0) * site[j, t]

    f.u.gamma[j, t - 1] <- z[j, t] * (1 - z[j, t - 1]) * ifelse(habitat[j]==2, 1, 0) * ifelse(sex[j]==2, 1, 0) * site[j, t]

    f.u.gamma.sites[j, t - 1] <- (1 - z[j, t - 1]) * ifelse(habitat[j]==2, 1, 0) * ifelse(sex[j]==2, 1, 0) * site[j, t]

    m.r.phi[j, t - 1] <- z[j, t] * z[j, t - 1] * ifelse(habitat[j]==1, 1, 0) * ifelse(sex[j]==1, 1, 0) * site[j, t]

    m.r.phi.sites[j, t - 1] <- z[j, t - 1] * ifelse(habitat[j]==1, 1, 0) * ifelse(sex[j]==1, 1, 0) * site[j, t]

    f.r.phi[j, t - 1] <- z[j, t] * z[j, t - 1] * ifelse(habitat[j]==1, 1, 0) * ifelse(sex[j]==2, 1, 0) * site[j, t]

    f.r.phi.sites[j, t - 1] <- z[j, t - 1] * ifelse(habitat[j]==1, 1, 0) * ifelse(sex[j]==2, 1, 0) * site[j, t]

    m.u.phi[j, t - 1] <- z[j, t] * z[j, t - 1] * ifelse(habitat[j]==2, 1, 0) * ifelse(sex[j]==1, 1, 0) * site[j, t]

    m.u.phi.sites[j, t - 1] <- z[j, t - 1] * ifelse(habitat[j]==2, 1, 0) * ifelse(sex[j]==1, 1, 0) * site[j, t]

    f.u.phi[j, t - 1] <- z[j, t] * z[j, t - 1] * ifelse(habitat[j]==2, 1, 0) * ifelse(sex[j]==2, 1, 0) * site[j, t]

    f.u.phi.sites[j, t - 1] <- z[j, t - 1] * ifelse(habitat[j]==2, 1, 0) * ifelse(sex[j]==2, 1, 0) * site[j, t]

  }  #t

}  #j

psi.male.riparian[1] <- sum(m.r.psi[, 1]) / sum(m.r.psi.sites[, 1])

psi.female.riparian[1] <- sum(f.r.psi[, 1]) / sum(f.r.psi.sites[, 1])

psi.male.upland[1] <- sum(m.u.psi[, 1]) / sum(m.u.psi.sites[, 1])

psi.female.upland[1] <- sum(f.u.psi[, 1]) / sum(f.u.psi.sites[, 1])

for (t in 2:nyear) {

  psi.male.riparian[t] <- sum(m.r.psi[, t]) / sum(m.r.psi.sites[, t])

  psi.female.riparian[t] <- sum(f.r.psi[, t]) / sum(f.r.psi.sites[, t])

  psi.male.upland[t] <- sum(m.u.psi[, t]) / sum(m.u.psi.sites[, t])

  psi.female.upland[t] <- sum(f.u.psi[, t]) / sum(f.u.psi.sites[, t])

  gamma.male.riparian[t - 1] <- sum(m.r.gamma[, t - 1]) /  sum(m.r.gamma.sites[, t - 1])

  gamma.female.riparian[t - 1] <- sum(f.r.gamma[, t - 1]) /  sum(f.r.gamma.sites[, t - 1])

  gamma.male.upland[t - 1] <- sum(m.u.gamma[, t - 1]) / sum(m.u.gamma.sites[, t - 1])

  gamma.female.upland[t - 1] <- sum(f.u.gamma[, t - 1]) / sum(f.u.gamma.sites[, t - 1])

  phi.male.riparian[t - 1] <- sum(m.r.phi[, t - 1]) /  sum(m.r.phi.sites[, t - 1])

  phi.female.riparian[t - 1] <- sum(f.r.phi[, t - 1]) /  sum(f.r.phi.sites[, t - 1])

  phi.male.upland[t - 1] <- sum(m.u.phi[, t - 1]) /  sum(m.u.phi.sites[, t - 1])

  phi.female.upland[t - 1] <- sum(f.u.phi[, t - 1]) /  sum(f.u.phi.sites[, t - 1])

} 

}

",fill = TRUE)

sink()

#####################################################################################################################################################
# Initial values

initialInits <- function() {list(

# start at constant values near zero for covariates

z = site, 
alpha.p = runif(3,-.25, .25),
alpha.psi = runif(nplant, -.25, .25),
alpha.phi = runif(nplant, -.25, .25),
alpha.gamma = runif(nplant, -.25, .25),
beta.repro.phi = runif(4, -.25, .25),
beta.repro.gamma = runif(4, -.25, .25),
beta.rgr.phi = runif(2, -.25, .25),
beta.rgr.gamma = runif(2, -.25, .25),
beta.d.p = runif(3, -.25, .25),
beta.d.psi = runif(2, -.25, .25),
beta.d.phi = runif(2, -.25, .25),
beta.d.gamma = runif(2, -.25, .25),
mu.alpha.psi = runif(4, -.25, .25),
mu.alpha.phi = runif(4, -.25, .25),
mu.alpha.gamma = runif(4, -.25, .25),
sigma.psi = runif(2, .5, .7),
sigma.phi = runif(1, .5, .7),
sigma.gamma = runif(1, .5, .7))}

# Parameters monitored

# parameter values and fit statistics only

params <- c(
"alpha.gamma",
"alpha.phi",
"alpha.psi",
"alpha.p",
"beta.repro.phi",
"beta.repro.gamma",
"beta.rgr.phi",
"beta.rgr.gamma",
"beta.d.p",
"beta.d.psi",
"beta.d.phi",
"beta.d.gamma",
"mu.psi",
"mu.phi",
"mu.gamma",
"sigma.psi",
"sigma.phi",
"sigma.gamma",
"zActual",
"zNew",
"z",
"psi",
"gamma",
"phi",
"gamma.male.riparian",
"gamma.female.riparian",
"gamma.male.upland",
"gamma.female.upland",
"phi.male.riparian",
"phi.female.riparian",
"phi.male.upland",
"phi.female.upland",
"psi.male.riparian",
"psi.female.riparian",
"psi.male.upland",
"psi.female.upland")

n.adapt <- 5000
n.update <- 100000
n.iter <- 15000

#####################################################################################################################################################
# Call JAGS from R

# multiple cores

library(parallel)

cl <- makeCluster(3) # Request 3 cores

clusterExport(cl, c("jags.data", "initialInits", "params", "nplant", "site", "n.adapt", "n.update", "n.iter")) 

system.time({

out <- clusterEvalQ(cl, {

  library(rjags)
  
  inits <- initialInits()  
    
  jm <- jags.model("weeviloccupancy.jags", jags.data, inits = inits, n.adapt = n.adapt, n.chains = 1)

  update(jm, n.iter = n.update)
  
  samples <- coda.samples(jm, n.iter = n.iter, variable.names = params, thin = 10)

  return(as.mcmc(samples))
    
  })

}) 

stopCluster(cl)
 
#####################################################################################################################################################

MCMCsamples <- mcmc.list(out)

MCMCsamplesCore <- MCMCsamples[, unlist(matchcols(MCMCsamples[[1]], with = c("alpha.p\\[", "beta.repro", "beta.rgr", "beta.d", "mu.", "sigma."), method = "or"))]

MCMCsamplesPlant <- MCMCsamples[, unlist(matchcols(MCMCsamples[[1]], with = c("alpha.psi\\[", "alpha.gamma\\[", "alpha.phi\\["), method = "or"))]

MCMCsamplesZ <- MCMCsamples[, unlist(matchcols(MCMCsamples[[1]], with = c("z\\[")))]

MCMCsamplesPsi1 <- MCMCsamples[, unlist(matchcols(MCMCsamples[[1]], with = c("psi\\[", "1\\]"), without = (".psi")))]

MCMCsamplesPsi2 <- MCMCsamples[, unlist(matchcols(MCMCsamples[[1]], with = c("psi\\[", "2\\]"), without = (".psi")))]

MCMCsamplesPsi3 <- MCMCsamples[, unlist(matchcols(MCMCsamples[[1]], with = c("psi\\[", "3\\]"), without = (".psi")))]

MCMCsamplesPhi1 <- MCMCsamples[, unlist(matchcols(MCMCsamples[[1]], with = c("phi\\[", "1\\]"), without = (".phi")))]

MCMCsamplesPhi2 <- MCMCsamples[, unlist(matchcols(MCMCsamples[[1]], with = c("phi\\[", "2\\]"), without = (".phi")))]

MCMCsamplesGamma1 <- MCMCsamples[, unlist(matchcols(MCMCsamples[[1]], with = c("gamma\\[", "1\\]"), without = (".gamma")))]

MCMCsamplesGamma2 <- MCMCsamples[, unlist(matchcols(MCMCsamples[[1]], with = c("gamma\\[", "2\\]"), without = (".gamma")))]

MCMCsamplesMeta <- MCMCsamples[, unlist(matchcols(MCMCsamples[[1]], with = c("gamma.male.riparian", "gamma.female.riparian",
                                                                             
                      "gamma.male.upland", "gamma.female.upland", "phi.male.riparian", "phi.female.riparian", "phi.male.upland",
                      
                      "phi.female.upland", "psi.male.riparian", "psi.female.riparian", "psi.male.upland","psi.female.upland"), method = "or"))]

MCMCsamplesPPC <- MCMCsamples[, unlist(matchcols(MCMCsamples[[1]], with = c("zActual", "zNew"), method = "or"))]

#####################################################################################################################################################

# setwd("/Users/coldwater/Documents/Projects/SalixStressVigor/ModelBuild")

save(MCMCsamplesCore, file = "MCMCsamplesCore.rda")

save(MCMCsamplesPlant, file = "MCMCsamplesPlant.rda")

save(MCMCsamplesZ, file = "MCMCsamplesZ.rda")

save(MCMCsamplesPsi1, file = "MCMCsamplesPsi1.rda")

save(MCMCsamplesPsi2, file = "MCMCsamplesPsi2.rda")

save(MCMCsamplesPsi3, file = "MCMCsamplesPsi3.rda")

save(MCMCsamplesPhi1, file = "MCMCsamplesPhi1.rda")

save(MCMCsamplesPhi2, file = "MCMCsamplesPhi2.rda")

save(MCMCsamplesGamma1, file = "MCMCsamplesGamma1.rda")

save(MCMCsamplesGamma2, file = "MCMCsamplesGamma2.rda")

save(MCMCsamplesMeta, file = "MCMCsamplesMeta.rda")

save(MCMCsamplesPPC, file = "MCMCsamplesPPC.rda")

unlink("weeviloccupancy.jags")