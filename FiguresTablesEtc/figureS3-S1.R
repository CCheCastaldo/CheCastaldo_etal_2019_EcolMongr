# ________________________________________________________________________________
# load packages

if('pacman' %in% rownames(installed.packages()) == FALSE) {
  
  install.packages('pacman', repos = "http://cran.case.edu")
  
}

pacman::p_load(gdata, dplyr, rjags, MCMCvis)

# ________________________________________________________________________________
# load data

load(file = "../Library/leafGasExchange.rda")

load(file = "../Library/leafP.rda")

work7a <- leafGasExchange %>%
  
  filter(year==2011) %>%
  
  group_by(plant_id, section) %>%
  
  summarise(PHOTO = mean(photosynthesis), TRANS = mean(transpiration), COND = mean(conductance),  WUE = mean(wue))

work7 <- leafP %>%
  
  filter(year==2011) %>%
  
  select(plant_id, section, P) %>%
  
  left_join(work7a, by = c("plant_id", "section")) %>%
  
  arrange(plant_id, section)

work8a <- leafGasExchange %>%
  
  filter(year==2010) %>%
  
  group_by(plant_id, section) %>%
  
  summarise(PHOTO = mean(photosynthesis), TRANS = mean(transpiration), COND = mean(conductance),  WUE = mean(wue))

work8 <- leafP %>%
  
  filter(year==2010) %>%
  
  select(plant_id, section, P) %>%
  
  left_join(work8a, by = c("plant_id", "section")) %>%
  
  arrange(plant_id, section)

  # ________________________________________________________________________________
# jags models

sink("jagsModel.jags")
  
cat("

model {
      
# priors

for (i in 1:3) {

  mu[i] ~ dnorm(0, 0.001)
  a[i] <- mu[i]^2/sigma^2
  b[i] <- mu[i]/sigma^2

}

sigma ~ dunif(0,100)

# likelihood

for (i in 1:n){
  y[i] ~ dgamma(a[section[i]], b[section[i]])
  rnew[i] ~ dgamma(a[section[i]], b[section[i]])
  rsqActual[i] <- pow(y[i] - mu[section[i]], 2)
  rsqNew[i] <- pow(rnew[i] - mu[section[i]], 2)
}

delta <- mu[1] - mu[3]

#  posterior predictive check
zActual <- sum(rsqActual[])   
zNew <- sum(rsqNew[]) 
test <- step(zNew - zActual) 	
bpvalue <- mean(test) 	

}", fill = TRUE)
sink()

sink("jagsModel1.jags")

cat("
    
    model {
    
    # priors
    
    for (i in 1:3) {
    
      alpha[i] ~ dnorm(0, .0001)
      beta[i] ~ dnorm(0, .0001)
    
    }
    
    sigma ~ dunif(0,100)
    
    # likelihood
    
    for (i in 1:n){
    mu[i] <- alpha[section[i]] + beta[section[i]] * g[i]
    a[i] <- mu[i]^2/sigma^2
    b[i] <- mu[i]/sigma^2
    photo[i] ~ dgamma(a[section[i]], b[section[i]])
    rnew[i] ~ dgamma(a[section[i]], b[section[i]])
    rsqActual[i] <- pow(photo[i] - mu[section[i]], 2)
    rsqNew[i] <- pow(rnew[i] - mu[section[i]], 2)
    }
    
    delta <- beta[1] - beta[3]
    
    #  posterior predictive check
    zActual <- sum(rsqActual[])   
    zNew <- sum(rsqNew[]) 
    test <- step(zNew - zActual) 	
    bpvalue <- mean(test) 	
    
    }", fill = TRUE)
sink()

sink("jagsModel2.jags")

cat("
    
model {
    
# priors
    
for (i in 1:3) {
    
mu[i] ~ dnorm(0, 0.001)

}
    
sigma ~ dunif(0,100)
tau <- pow(sigma, -2)

# likelihood
    
for (i in 1:n){
  y[i] ~ dnorm(mu[section[i]], tau)
  rnew[i] ~ dnorm(mu[section[i]], tau)
  rsqActual[i] <- pow(y[i] - mu[section[i]], 2)
  rsqNew[i] <- pow(rnew[i] - mu[section[i]], 2)
}
    
delta <- mu[1] - mu[2]
    
#  posterior predictive check
zActual <- sum(rsqActual[])   
zNew <- sum(rsqNew[]) 
test <- step(zNew - zActual) 	
bpvalue <- mean(test) 	
    
}", fill = TRUE)
sink()

sink("jagsModel3.jags")

cat("
    
model {
    
# priors
    
for (i in 1:3) {
  mu[i] ~ dnorm(0, 0.001)
  a[i] <- (mu[i]^2 - mu[i]^3 - mu[i] * sigma^2)  / sigma^2
  b[i] <- (mu[i] - 2 * mu[i]^2 + mu[i]^3 - sigma^2 + mu[i] * sigma^2)  / sigma^2
}
    
sigma ~ dunif(0,100)

# likelihood
    
for (i in 1:n){
  y[i] ~ dbeta(a[section[i]], b[section[i]])
  rnew[i] ~ dbeta(a[section[i]], b[section[i]])
  rsqActual[i] <- pow(y[i] - mu[section[i]], 2)
  rsqNew[i] <- pow(rnew[i] - mu[section[i]], 2)
}
    
delta <- mu[1] - mu[3]
    
#  posterior predictive check
zActual <- sum(rsqActual[])   
zNew <- sum(rsqNew[]) 
test <- step(zNew - zActual) 	
bpvalue <- mean(test) 	
    
}", fill = TRUE)
sink()

sink("jagsModel4.jags")

cat("
    
    model {
    
    # priors
    
    for (i in 1:2) {
    
    mu[i] ~ dnorm(0, 0.001)
    a[i] <- mu[i]^2/sigma^2
    b[i] <- mu[i]/sigma^2
    
    }
    
    sigma ~ dunif(0,100)
    
    # likelihood
    
    for (i in 1:n){
    y[i] ~ dgamma(a[section[i]], b[section[i]])
    rnew[i] ~ dgamma(a[section[i]], b[section[i]])
    rsqActual[i] <- pow(y[i] - mu[section[i]], 2)
    rsqNew[i] <- pow(rnew[i] - mu[section[i]], 2)
    }
    
    delta <- mu[1] - mu[2]
    
    #  posterior predictive check
    zActual <- sum(rsqActual[])   
    zNew <- sum(rsqNew[]) 
    test <- step(zNew - zActual) 	
    bpvalue <- mean(test) 	
    
    }", fill = TRUE)
sink()

sink("jagsModel5.jags")

cat("
    
    model {
    
    # priors
    
    for (i in 1:2) {
    
    mu[i] ~ dnorm(0, 0.001)
    
    }
    
    sigma ~ dunif(0,100)
    tau <- pow(sigma, -2)
    
    # likelihood
    
    for (i in 1:n){
    y[i] ~ dnorm(mu[section[i]], tau)
    rnew[i] ~ dnorm(mu[section[i]], tau)
    rsqActual[i] <- pow(y[i] - mu[section[i]], 2)
    rsqNew[i] <- pow(rnew[i] - mu[section[i]], 2)
    }
    
    delta <- mu[1] - mu[2]
    
    #  posterior predictive check
    zActual <- sum(rsqActual[])   
    zNew <- sum(rsqNew[]) 
    test <- step(zNew - zActual) 	
    bpvalue <- mean(test) 	
    
    }", fill = TRUE)
sink()

sink("jagsModel6.jags")

cat("
    
    model {
    
    # priors
    
    for (i in 1:2) {
    mu[i] ~ dnorm(0, 0.001)
    a[i] <- (mu[i]^2 - mu[i]^3 - mu[i] * sigma^2)  / sigma^2
    b[i] <- (mu[i] - 2 * mu[i]^2 + mu[i]^3 - sigma^2 + mu[i] * sigma^2)  / sigma^2
    }
    
    sigma ~ dunif(0,100)
    
    # likelihood
    
    for (i in 1:n){
    y[i] ~ dbeta(a[section[i]], b[section[i]])
    rnew[i] ~ dbeta(a[section[i]], b[section[i]])
    rsqActual[i] <- pow(y[i] - mu[section[i]], 2)
    rsqNew[i] <- pow(rnew[i] - mu[section[i]], 2)
    }
    
    delta <- mu[1] - mu[2]
    
    #  posterior predictive check
    zActual <- sum(rsqActual[])   
    zNew <- sum(rsqNew[]) 
    test <- step(zNew - zActual) 	
    bpvalue <- mean(test) 	
    
    }", fill = TRUE)
sink()


params <- c(
  "mu", 
  "sigma", 
  "delta", 
  "zActual", 
  "zNew",
  "bpvalue")

# ________________________________________________________________________________
# photo rate, 2010

pdf(file = "figureS3-S1.pdf", width = 8, height = 11)

par(mfrow = c(5, 2))

DataQuery <- list(section = as.numeric(as.factor(work8$section)), y = work8$PHOTO, n = dim(work8)[1])  

randomInits <- function() {mu <- runif(3, 5, 10); sigma <- runif(1, 5, 10); return(list(mu = mu, sigma = sigma))}

jm = rjags::jags.model(data = DataQuery, file = "jagsModel.jags", inits = randomInits(), n.chains = 3, n.adapt = 1000)

update(jm, n.iter = 5000)

out = coda.samples(jm, n.iter = 20000, variable.names = params, thin = 10)

#MCMCvis::MCMCtrace(out, params = params, pdf = FALSE)

MCMCvis::MCMCsummary(out, params = params, n.eff = TRUE, Rhat = TRUE, digits = 2)

MCMCvis::MCMCplot(out, params = c("mu"), labels = c("wet\n(260% f.c.)", "intermediate\n(105% f.c)", "intermediate\n(125% f.c)"), main = "",
                  
                  mar = c(5, 5, 4, 2), xlab = "Photosynthetic rate, 2010", axis_text_sz = 1, tick_text_sz = 1, main_text_sz = 1, labels_sz = 1, xlim= c(4, 22), tick_pos = seq(4, 22, by = 3))

# ________________________________________________________________________________
# photo rate, 2011

DataQuery <- list(section = as.numeric(as.factor(work7$section)), y = work7$PHOTO, n = dim(work7)[1])  

randomInits <- function() {mu <- runif(3, 5, 10); sigma <- runif(1, 5, 10); return(list(mu = mu, sigma = sigma))}

jm = rjags::jags.model(data = DataQuery, file = "jagsModel.jags", inits = randomInits(), n.chains = 3, n.adapt = 1000)

update(jm, n.iter = 5000)

out = coda.samples(jm, n.iter = 20000, variable.names = params, thin = 10)

#MCMCvis::MCMCtrace(out, params = params, pdf = FALSE)

MCMCvis::MCMCsummary(out, params = params, n.eff = TRUE, Rhat = TRUE, digits = 2)

MCMCvis::MCMCplot(out, params = c("mu"), labels = c("wet\n(220% f.c.)", "intermediate\n(95% f.c)", "dry\n(35% f.c)"), main = "",
                  
                  mar = c(5, 5, 4, 2), xlab = "Photosynthetic rate, 2011", axis_text_sz = 1, tick_text_sz = 1, main_text_sz = 1, labels_sz = 1, xlim= c(4, 22), tick_pos = seq(4, 22, by = 3))

# ________________________________________________________________________________
# conductance rate, 2010

DataQuery <- list(section = as.numeric(as.factor(work8$section)), y = work8$COND, n = dim(work8)[1])  

randomInits <- function() {mu <- runif(3, 5, 10); sigma <- runif(1, 5, 10); return(list(mu = mu, sigma = sigma))}

jm = rjags::jags.model(data = DataQuery, file = "jagsModel.jags", inits = randomInits(), n.chains = 3, n.adapt = 1000)

update(jm, n.iter = 5000)

out = coda.samples(jm, n.iter = 20000, variable.names = params, thin = 10)

#MCMCvis::MCMCtrace(out, params = params, pdf = FALSE)

MCMCvis::MCMCsummary(out, params = params, n.eff = TRUE, Rhat = TRUE, digits = 2)

MCMCvis::MCMCplot(out, params = c("mu"), labels = c("wet\n(260% f.c.)", "intermediate\n(105% f.c)", "intermediate\n(125% f.c)"), main = "",
                  
                  mar = c(5, 5, 4, 2), ref = NULL, xlab = "Conductance, 2010", axis_text_sz = 1, tick_text_sz = 1, main_text_sz = 1, labels_sz = 1, xlim= c(0, .7), tick_pos = seq(0, .7, by = .1))

# ________________________________________________________________________________
# conductance rate, 2011

DataQuery <- list(section = as.numeric(as.factor(work7$section)), y = work7$COND, n = dim(work7)[1])  

randomInits <- function() {mu <- runif(3, 5, 10); sigma <- runif(1, 5, 10); return(list(mu = mu, sigma = sigma))}

jm = rjags::jags.model(data = DataQuery, file = "jagsModel.jags", inits = randomInits(), n.chains = 3, n.adapt = 1000)

update(jm, n.iter = 5000)

out = coda.samples(jm, n.iter = 20000, variable.names = params, thin = 10)

#MCMCvis::MCMCtrace(out, params = params, pdf = FALSE)

MCMCvis::MCMCsummary(out, params = params, n.eff = TRUE, Rhat = TRUE, digits = 2)

MCMCvis::MCMCplot(out, params = c("mu"), labels = c("wet\n(220% f.c.)", "intermediate\n(95% f.c)", "dry\n(35% f.c)"), main = "",
                  
                  mar = c(5, 5, 4, 2), ref = NULL, xlab = "Conductance, 2011", axis_text_sz = 1, tick_text_sz = 1, main_text_sz = 1, labels_sz = 1, xlim= c(0, .7), tick_pos = seq(0, .7, by = .1))

# ________________________________________________________________________________
# trans, 2010

DataQuery <- list(section = as.numeric(as.factor(work8$section)), y = work8$TRANS, n = dim(work8)[1])  

randomInits <- function() {mu <- runif(3, 5, 10); sigma <- runif(1, 5, 10); return(list(mu = mu, sigma = sigma))}

jm = rjags::jags.model(data = DataQuery, file = "jagsModel.jags", inits = randomInits(), n.chains = 3, n.adapt = 1000)

update(jm, n.iter = 5000)

out = coda.samples(jm, n.iter = 20000, variable.names = params, thin = 10)

#MCMCvis::MCMCtrace(out, params = params, pdf = FALSE)

MCMCvis::MCMCsummary(out, params = params, n.eff = TRUE, Rhat = TRUE, digits = 2)

MCMCvis::MCMCplot(out, params = c("mu"), labels = c("wet\n(260% f.c.)", "intermediate\n(105% f.c)", "intermediate\n(125% f.c)"), main = "",
                  
                  mar = c(5, 5, 4, 2), ref = NULL, xlab = "Transpiration, 2010", axis_text_sz = 1, tick_text_sz = 1, main_text_sz = 1, labels_sz = 1, xlim= c(1, 6), tick_pos = seq(1, 6, by = 1))

# ________________________________________________________________________________
# trans, 2011

DataQuery <- list(section = as.numeric(as.factor(work7$section)), y = work7$TRANS, n = dim(work7)[1])  

randomInits <- function() {mu <- runif(3, 5, 10); sigma <- runif(1, 5, 10); return(list(mu = mu, sigma = sigma))}

jm = rjags::jags.model(data = DataQuery, file = "jagsModel.jags", inits = randomInits(), n.chains = 3, n.adapt = 1000)

update(jm, n.iter = 5000)

out = coda.samples(jm, n.iter = 20000, variable.names = params, thin = 10)

#MCMCvis::MCMCtrace(out, params = params, pdf = FALSE)

MCMCvis::MCMCsummary(out, params = params, n.eff = TRUE, Rhat = TRUE, digits = 2)

MCMCvis::MCMCplot(out, params = c("mu"), labels = c("wet\n(220% f.c.)", "intermediate\n(95% f.c)", "dry\n(35% f.c)"), main = "",
                  
                  mar = c(5, 5, 4, 2), ref = NULL, xlab = "Transpiration, 2011", axis_text_sz = 1, tick_text_sz = 1, main_text_sz = 1, labels_sz = 1, xlim= c(1, 6), tick_pos = seq(1, 6, by = 1))

# ________________________________________________________________________________
# wue, 2010

DataQuery <- list(section = as.numeric(as.factor(work8$section)), y = work8$WUE, n = dim(work8)[1])  

randomInits <- function() {mu <- runif(3, .5, .75); sigma <- runif(1, .25, .5); return(list(mu = mu, sigma = sigma))}

jm = rjags::jags.model(data = DataQuery, file = "jagsModel.jags", inits = randomInits(), n.chains = 3, n.adapt = 1000)

update(jm, n.iter = 5000)

out = coda.samples(jm, n.iter = 20000, variable.names = params, thin = 10)

#MCMCvis::MCMCtrace(out, params = params, pdf = FALSE)

MCMCvis::MCMCsummary(out, params = params, n.eff = TRUE, Rhat = TRUE, digits = 2)

MCMCvis::MCMCplot(out, params = c("mu"), labels = c("wet\n(260% f.c.)", "intermediate\n(105% f.c)", "intermediate\n(125% f.c)"), main = "",
                  
                  mar = c(5, 5, 4, 2), ref = NULL, xlab = "WUE, 2010", axis_text_sz = 1, tick_text_sz = 1, main_text_sz = 1, labels_sz = 1, xlim= c(25, 105), tick_pos = seq(25, 105, by = 20))

# ________________________________________________________________________________
# wue, 2011

DataQuery <- list(section = as.numeric(as.factor(work7$section)), y = work7$WUE, n = dim(work7)[1])  

randomInits <- function() {mu <- runif(3, .5, .75); sigma <- runif(1, .25, .5); return(list(mu = mu, sigma = sigma))}

jm = rjags::jags.model(data = DataQuery, file = "jagsModel.jags", inits = randomInits(), n.chains = 3, n.adapt = 1000)

update(jm, n.iter = 5000)

out = coda.samples(jm, n.iter = 20000, variable.names = params, thin = 10)

#MCMCvis::MCMCtrace(out, params = params, pdf = FALSE)

MCMCvis::MCMCsummary(out, params = params, n.eff = TRUE, Rhat = TRUE, digits = 2)

MCMCvis::MCMCplot(out, params = c("mu"), labels = c("wet\n(220% f.c.)", "intermediate\n(95% f.c)", "dry\n(35% f.c)"), main = "",
                  
                  mar = c(5, 5, 4, 2), ref = NULL, xlab = "WUE, 2011", axis_text_sz = 1, tick_text_sz = 1, main_text_sz = 1, labels_sz = 1, xlim= c(25, 105), tick_pos = seq(25, 105, by = 20))

# ________________________________________________________________________________
# P, 2010

DataQuery <- list(section = as.numeric(as.factor(work8$section)), y = work8$P, n = dim(work8)[1])  

randomInits <- function() {mu <- runif(3, .002, .004); sigma <- runif(1, .001, .002); return(list(mu = mu, sigma = sigma))}

jm = rjags::jags.model(data = DataQuery, file = "jagsModel3.jags", inits = randomInits(), n.chains = 3, n.adapt = 1000)

update(jm, n.iter = 5000)

out = coda.samples(jm, n.iter = 20000, variable.names = params, thin = 10)

#MCMCvis::MCMCtrace(out, params = params, pdf = FALSE)

MCMCvis::MCMCsummary(out, params = params, n.eff = TRUE, Rhat = TRUE, digits = 10)

MCMCvis::MCMCplot(out, params = c("mu"), labels = c("wet\n(260% f.c.)", "intermediate\n(105% f.c)", "intermediate\n(125% f.c)"), main = "",
                  
                  mar = c(5, 5, 4, 2), ref = NULL, xlab = "Leaf percent P, 2010", axis_text_sz = 1, tick_text_sz = 1, main_text_sz = 1, labels_sz = 1, xlim= c(.2, .6), tick_pos = seq(.2, .6, by = .1))

# ________________________________________________________________________________
# P, 2011

DataQuery <- list(section = as.numeric(as.factor(work7$section)), y = work7$P, n = dim(work7)[1])  

randomInits <- function() {mu <- runif(3, .002, .004); sigma <- runif(1, .001, .002); return(list(mu = mu, sigma = sigma))}

jm = rjags::jags.model(data = DataQuery, file = "jagsModel3.jags", inits = randomInits(), n.chains = 3, n.adapt = 1000)

update(jm, n.iter = 5000)

out = coda.samples(jm, n.iter = 20000, variable.names = params, thin = 10)

#MCMCvis::MCMCtrace(out, params = params, pdf = FALSE)

MCMCvis::MCMCsummary(out, params = params, n.eff = TRUE, Rhat = TRUE, digits = 10)

MCMCvis::MCMCplot(out, params = c("mu"), labels = c("wet\n(220% f.c.)", "intermediate\n(95% f.c)", "dry\n(35% f.c)"), main = "",
                  
                  mar = c(5, 5, 4, 2), ref = NULL, xlab = "Leaf percent P, 2011", axis_text_sz = 1, tick_text_sz = 1, main_text_sz = 1, labels_sz = 1, xlim= c(.2, .6), tick_pos = seq(.2, .6, by = .1))

dev.off()
