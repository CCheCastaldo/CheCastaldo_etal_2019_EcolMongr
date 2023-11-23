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

load(file = "../library/salix_borer_wide.rda")

load(file = paste0("../model_build/global_model/", version, "/mcmc_samples_z.rda"))

chains <- apply(rbind(as.matrix(MCMCsamplesZ[[1]]), as.matrix(MCMCsamplesZ[[2]]), as.matrix(MCMCsamplesZ[[3]])), 2, mean)

zWeights <- as.data.frame(cbind(chains[1:2087], chains[2088:4174], chains[4175:6261]))

SalixBorerWide <- cbind(SalixBorerWide, zWeights)

#######################################################################################################################################
# define weights and create data object for jags

work1 <- subset(SalixBorerWide, site09==1, select = c(habitat, sex, site09, site10, V1))

work2 <- subset(SalixBorerWide, site10==1, select = c(habitat, sex, site10, site11, V2))

names(work1) <- c("habitat", "sex", "denom", "num", "w")                          

names(work2) <- c("habitat", "sex", "denom", "num", "w")                          

work1$year <- 1

work2$year <- 2

work3 <- rbind(work1, work2)

work3$check <- 0; i <- which(work3$w <= .1 | work3$w >= .9); work3$check[i] <- 1

table(work3$check)/dim(work3)[1]

work3a <- subset(work3, w <= .1)

work4 <- ddply(work3a, c("habitat", "sex", "year"), function(x) c(
  
  alive = sum(x$num, na.rm = TRUE),
  
  total = sum(x$denom, na.rm = TRUE)))

work4$proportion <- work4$alive/work4$total

DataQuery <- list(
  
  alive = work4$alive,
  n = work4$total,
  habitat = work4$habitat + 1,
  sex = work4$sex + 1,
  year = work4$year)

#######################################################################################################################################
# do weighted regression on diameter

sink("SurvivalStems.R")

cat("
    
model {
    
# priors    

for (i in 1:2) {
  for (j in 1:2) {
    for (k in 1:2) {
      p[i, j, k] ~ dunif (0, 1)
    }
  }
}

# likelihood

for (i in 1:8) {

  alive[i] ~ dbinom(p[habitat[i], sex[i], year[i]], n[i])

}

# mean comparisons

sex.riparian.2010 <- p[1, 1, 1] - p[1, 2, 1]
sex.riparian.2011 <- p[1, 1, 2] - p[1, 2, 2]
sex.upland.2010 <- p[2, 1, 1] - p[2, 2, 1]
sex.upland.2011 <- p[2, 1, 2] - p[2, 2, 2]

year.riparian.male <- p[1, 1, 1] - p[1, 1, 2]
year.riparian.female <- p[1, 2, 1] - p[1, 2, 2]
year.upland.male <- p[2, 1, 1] - p[2, 1, 2]
year.upland.female <- p[2, 2, 1] - p[2, 2, 2]

habitat.male.2010 <- p[1, 1, 1] - p[2, 1, 1]
habitat.female.2010 <- p[1, 2, 1] - p[2, 2, 1]
habitat.male.2011 <- p[1, 1, 2] - p[2, 1, 2]
habitat.female.2011 <- p[1, 2, 2] - p[2, 2, 2]

}    
",fill = TRUE)
sink()

# Initial values

inits <- function() {list(
  
  p = array(runif(8,0,1), dim = c(2,2,2)))} 

# Parameters monitored

params <- c(
  "p", 
  "sex.riparian.2010", 
  "sex.upland.2010",
  "sex.riparian.2011", 
  "sex.upland.2011",
  "year.riparian.male",
  "year.riparian.female",
  "year.upland.male",
  "year.upland.female",
  "habitat.male.2010",
  "habitat.female.2010",
  "habitat.male.2011",
  "habitat.female.2011")
  
# MCMC settings

n.adapt = 5000
n.update = 10000
n.iter = 15000

# JAGS run

jm = jags.model("SurvivalStems.R", data = DataQuery, inits = inits(), n.chains = 3, n.adapt = n.adapt)

update(jm, n.iter = n.update)

zm = coda.samples(jm, variable.names = params, n.iter = n.iter, thin = 10)

# summarize JAGS run

summary(zm)

gelman.diag(zm, multivariate = FALSE)

MCMCsamplesSurvivalA <- zm

unlist("SurvivalStems.R")

#######################################################################################################################################
# define weights and create data object for jags

work3a <- subset(work3, w >= .9)

work4 <- ddply(work3a, c("habitat", "sex", "year"), function(x) c(
  
  alive = sum(x$num, na.rm = TRUE),
  
  total = sum(x$denom, na.rm = TRUE)))

work4$proportion <- work4$alive/work4$total

DataQuery <- list(
  
  alive = work4$alive,
  n = work4$total,
  habitat = work4$habitat + 1,
  sex = work4$sex + 1,
  year = work4$year)

#######################################################################################################################################
# do weighted regression on diameter

sink("SurvivalStems.R")

cat("
    
    model {
    
    # priors    
    
    for (i in 1:2) {
    for (j in 1:2) {
    for (k in 1:2) {
    p[i, j, k] ~ dunif (0, 1)
    }
    }
    }
    
    # likelihood
    
    for (i in 1:8) {
    
    alive[i] ~ dbinom(p[habitat[i], sex[i], year[i]], n[i])
    
    }
    
    # mean comparisons
    
    sex.riparian.2010 <- p[1, 1, 1] - p[1, 2, 1]
    sex.riparian.2011 <- p[1, 1, 2] - p[1, 2, 2]
    sex.upland.2010 <- p[2, 1, 1] - p[2, 2, 1]
    sex.upland.2011 <- p[2, 1, 2] - p[2, 2, 2]
    
    year.riparian.male <- p[1, 1, 1] - p[1, 1, 2]
    year.riparian.female <- p[1, 2, 1] - p[1, 2, 2]
    year.upland.male <- p[2, 1, 1] - p[2, 1, 2]
    year.upland.female <- p[2, 2, 1] - p[2, 2, 2]
 
  habitat.male.2010 <- p[1, 1, 1] - p[2, 1, 1]
  habitat.female.2010 <- p[1, 2, 1] - p[2, 2, 1]
  habitat.male.2011 <- p[1, 1, 2] - p[2, 1, 2]
  habitat.female.2011 <- p[1, 2, 2] - p[2, 2, 2]
    
}    
",fill = TRUE)
sink()

# Initial values

inits <- function() {list(
  
  p = array(runif(8,0,1), dim = c(2,2,2)))} 

# Parameters monitored

params <- c(
  "p", 
  "sex.riparian.2010", 
  "sex.upland.2010",
  "sex.riparian.2011", 
  "sex.upland.2011",
  "year.riparian.male",
  "year.riparian.female",
  "year.upland.male",
  "year.upland.female",
  "habitat.male.2010",
  "habitat.female.2010",
  "habitat.male.2011",
  "habitat.female.2011")

# MCMC settings

n.adapt = 5000
n.update = 10000
n.iter = 15000

# JAGS run

jm = jags.model("SurvivalStems.R", data = DataQuery, inits = inits(), n.chains = 3, n.adapt = n.adapt)

update(jm, n.iter = n.update)

zm = coda.samples(jm, variable.names = params, n.iter = n.iter, thin = 10)

# summarize JAGS run

summary(zm)

gelman.diag(zm, multivariate = FALSE)

MCMCsamplesSurvivalNA <- zm

unlist("SurvivalStems.R")

#######################################################################################################################################

chains2 <- rbind(as.matrix(MCMCsamplesSurvivalA[[1]]), as.matrix(MCMCsamplesSurvivalA[[2]]), as.matrix(MCMCsamplesSurvivalA[[3]]))

Dtable <- array(NA, dim = c(4, 6))

Dtable[1, 1] <- median(chains2[, matchcols(chains2, with = "p\\[1,1,1")])
Dtable[2, 1] <- median(chains2[, matchcols(chains2, with = "p\\[1,2,1")])
Dtable[3, 1] <- median(chains2[, matchcols(chains2, with = "p\\[2,1,1")])
Dtable[4, 1] <- median(chains2[, matchcols(chains2, with = "p\\[2,2,1")])

Dtable[1, 2] <- quantile(chains2[, matchcols(chains2, with = "p\\[1,1,1")], prob = .025)
Dtable[2, 2] <- quantile(chains2[, matchcols(chains2, with = "p\\[1,2,1")], prob = .025)
Dtable[3, 2] <- quantile(chains2[, matchcols(chains2, with = "p\\[2,1,1")], prob = .025)
Dtable[4, 2] <- quantile(chains2[, matchcols(chains2, with = "p\\[2,2,1")], prob = .025)

Dtable[1, 3] <- quantile(chains2[, matchcols(chains2, with = "p\\[1,1,1")], prob = .975)
Dtable[2, 3] <- quantile(chains2[, matchcols(chains2, with = "p\\[1,2,1")], prob = .975)
Dtable[3, 3] <- quantile(chains2[, matchcols(chains2, with = "p\\[2,1,1")], prob = .975)
Dtable[4, 3] <- quantile(chains2[, matchcols(chains2, with = "p\\[2,2,1")], prob = .975)

Dtable[1, 4] <- median(chains2[, matchcols(chains2, with = "p\\[1,1,2")])
Dtable[2, 4] <- median(chains2[, matchcols(chains2, with = "p\\[1,2,2")])
Dtable[3, 4] <- median(chains2[, matchcols(chains2, with = "p\\[2,1,2")])
Dtable[4, 4] <- median(chains2[, matchcols(chains2, with = "p\\[2,2,2")])

Dtable[1, 5] <- quantile(chains2[, matchcols(chains2, with = "p\\[1,1,2")], prob = .025)
Dtable[2, 5] <- quantile(chains2[, matchcols(chains2, with = "p\\[1,2,2")], prob = .025)
Dtable[3, 5] <- quantile(chains2[, matchcols(chains2, with = "p\\[2,1,2")], prob = .025)
Dtable[4, 5] <- quantile(chains2[, matchcols(chains2, with = "p\\[2,2,2")], prob = .025)

Dtable[1, 6] <- quantile(chains2[, matchcols(chains2, with = "p\\[1,1,2")], prob = .975)
Dtable[2, 6] <- quantile(chains2[, matchcols(chains2, with = "p\\[1,2,2")], prob = .975)
Dtable[3, 6] <- quantile(chains2[, matchcols(chains2, with = "p\\[2,1,2")], prob = .975)
Dtable[4, 6] <- quantile(chains2[, matchcols(chains2, with = "p\\[2,2,2")], prob = .975)

#######################################################################################################################################

habitat <- c("Riparian", "Riparian", "Upland", "Upland")

sex <- c("Male", "Female", "Male", "Female")

S2010 <- paste0("(", round(Dtable[,2], 2), " -- ", round(Dtable[,3], 2), ")")

S2011 <- paste0("(", round(Dtable[,5], 2), " -- ", round(Dtable[,6], 2), ")")

Dtable2 <- as.data.frame(cbind(habitat, sex, round(Dtable[,1],2), S2010, round(Dtable[,4],2), S2011))

xtable(Dtable2, digits = 2)

#######################################################################################################################################

chains2 <- rbind(as.matrix(MCMCsamplesSurvivalNA[[1]]), as.matrix(MCMCsamplesSurvivalNA[[2]]), as.matrix(MCMCsamplesSurvivalNA[[3]]))

Dtable <- array(NA, dim = c(4, 6))

Dtable[1, 1] <- median(chains2[, matchcols(chains2, with = "p\\[1,1,1")])
Dtable[2, 1] <- median(chains2[, matchcols(chains2, with = "p\\[1,2,1")])
Dtable[3, 1] <- median(chains2[, matchcols(chains2, with = "p\\[2,1,1")])
Dtable[4, 1] <- median(chains2[, matchcols(chains2, with = "p\\[2,2,1")])

Dtable[1, 2] <- quantile(chains2[, matchcols(chains2, with = "p\\[1,1,1")], prob = .025)
Dtable[2, 2] <- quantile(chains2[, matchcols(chains2, with = "p\\[1,2,1")], prob = .025)
Dtable[3, 2] <- quantile(chains2[, matchcols(chains2, with = "p\\[2,1,1")], prob = .025)
Dtable[4, 2] <- quantile(chains2[, matchcols(chains2, with = "p\\[2,2,1")], prob = .025)

Dtable[1, 3] <- quantile(chains2[, matchcols(chains2, with = "p\\[1,1,1")], prob = .975)
Dtable[2, 3] <- quantile(chains2[, matchcols(chains2, with = "p\\[1,2,1")], prob = .975)
Dtable[3, 3] <- quantile(chains2[, matchcols(chains2, with = "p\\[2,1,1")], prob = .975)
Dtable[4, 3] <- quantile(chains2[, matchcols(chains2, with = "p\\[2,2,1")], prob = .975)

Dtable[1, 4] <- median(chains2[, matchcols(chains2, with = "p\\[1,1,2")])
Dtable[2, 4] <- median(chains2[, matchcols(chains2, with = "p\\[1,2,2")])
Dtable[3, 4] <- median(chains2[, matchcols(chains2, with = "p\\[2,1,2")])
Dtable[4, 4] <- median(chains2[, matchcols(chains2, with = "p\\[2,2,2")])

Dtable[1, 5] <- quantile(chains2[, matchcols(chains2, with = "p\\[1,1,2")], prob = .025)
Dtable[2, 5] <- quantile(chains2[, matchcols(chains2, with = "p\\[1,2,2")], prob = .025)
Dtable[3, 5] <- quantile(chains2[, matchcols(chains2, with = "p\\[2,1,2")], prob = .025)
Dtable[4, 5] <- quantile(chains2[, matchcols(chains2, with = "p\\[2,2,2")], prob = .025)

Dtable[1, 6] <- quantile(chains2[, matchcols(chains2, with = "p\\[1,1,2")], prob = .975)
Dtable[2, 6] <- quantile(chains2[, matchcols(chains2, with = "p\\[1,2,2")], prob = .975)
Dtable[3, 6] <- quantile(chains2[, matchcols(chains2, with = "p\\[2,1,2")], prob = .975)
Dtable[4, 6] <- quantile(chains2[, matchcols(chains2, with = "p\\[2,2,2")], prob = .975)

#######################################################################################################################################

habitat <- c("Riparian", "Riparian", "Upland", "Upland")

sex <- c("Male", "Female", "Male", "Female")

S2010 <- paste0("(", round(Dtable[,2], 2), " -- ", round(Dtable[,3], 2), ")")

S2011 <- paste0("(", round(Dtable[,5], 2), " -- ", round(Dtable[,6], 2), ")")

Dtable2 <- as.data.frame(cbind(habitat, sex, round(Dtable[,1],2), S2010, round(Dtable[,4],2), S2011))

xtable(Dtable2, digits = 2)

unlist("SurvivalStems.R")