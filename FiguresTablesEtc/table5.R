library(coda)
library(gdata)
library(rjags)
library(broman)
library(plyr)
library(reshape)
library(xtable)

#######################################################################################################################################
# load data

load(file = "../Library/SalixBorerWide.rda")

#######################################################################################################################################
# define weights and create data object for jags

work1 <- subset(SalixBorerWide, site10==1, select = c(habitat, sex, site09, site10))

work2 <- subset(SalixBorerWide, site11==1, select = c(habitat, sex, site10, site11))

work1$num <- 1 - work1$site09

work2$num <- 1 - work2$site10

work1$denom <- 1

work2$denom <- 1

work1$year <- 1

work2$year <- 2

work1 <- subset(work1, select = c("habitat", "sex", "year", "denom", "num"))                         

work2 <- subset(work2, select = c("habitat", "sex", "year", "denom", "num"))                         

work3 <- rbind(work1, work2)

work4 <- ddply(work3, c("habitat", "sex", "year"), function(x) c(
  
  recruited = sum(x$num, na.rm = TRUE),
  
  total = sum(x$denom, na.rm = TRUE)))

work4$proportion <- work4$recruited/work4$total

DataQuery <- list(
  
  recruited = work4$recruited,
  n = work4$total,
  habitat = work4$habitat + 1,
  sex = work4$sex + 1,
  year = work4$year)

#######################################################################################################################################
# do weighted regression on diameter

sink("RecruitmentStems.R")

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

  recruited[i] ~ dbinom(p[habitat[i], sex[i], year[i]], n[i])

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

jm = jags.model("RecruitmentStems.R", data = DataQuery, inits = inits(), n.chains = 3, n.adapt = n.adapt)

update(jm, n.iter = n.update)

zm = coda.samples(jm, variable.names = params, n.iter = n.iter, thin = 10)

# summarize JAGS run

summary(zm)

gelman.diag(zm, multivariate = FALSE)

unlist("RecruitmentStems.R")

#######################################################################################################################################

chains2 <- rbind(as.matrix(zm[[1]]), as.matrix(zm[[2]]), as.matrix(zm[[3]]))

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

