####### Single season site-occupancy model on weasel Data from Dieren
rm(list=ls())

## Open library
library(jagsUI)

###### Model for 2017 ######

# Define model in BUGS language
cat(file="occ-Weasel17month.txt", 
    "
    model{
    
    # Priors and model for params

    int.psi ~ dunif(0,1)        #Intercept for occupancy probability 

    int.p ~ dunif(0,1)          #Intercept for detection probability 

    # prior for random intercept per location
    sigma_r ~ dunif(0,100)
    tau_r <- 1/(sigma_r*sigma_r)
    for(j in 1:n.loc){
    r.psi[j] ~ dnorm(0,tau_r)
    }

    # Priors for random intercepts per month
    sigma_r ~ dunif(0,100)
    tau_r <- 1/(sigma_r*sigma_r)
    sigma_r2 ~ dunif(0,100)
    tau_r2 <- 1/(sigma_r2*sigma_r2)
    for(t in 1:n.month){
      r.psi[t] ~ dnorm(0,tau_r)
      r.p[t] ~ dnorm(0,tau_r2)
    }
    
    #Slope of p covariate :: diameter
    
    beta.lp ~ dnorm(0,0.001)
    
    
     #Likelihood (or basic model structure)
    for(i in 1:n.site){
    
    #Occurrence in site i
    z[i] ~ dbern(psi[i])
    logit(psi[i]) <- logit(int.psi) + r.psi[loc[i]] + r.psi2[month[i]]
    
    
    for(k in 1:n.days){
    
    # detection probability on day k
    y[i,k] ~ dbern(mu.y[i,k])
    mu.y[i,k] <- z[i] * p[i,k]
    logit(p[i,k]) <- logit(int.p) + r.p[month[i]] +
    beta.lp * diameter[i]
    }
    }

    # model for missing covariates
    for(i in 1:n.site){
    diameter[i] ~ dbern(q) #prior model for wheter diameter is 0 (8cm) or 1 (10cm)
    }
    # vague prior for prob of diameter being 1
    q ~ dbeta(1,1)

    # derived parameters

    logit(psi.mar) <- logit(int.psi) + r.psi[1]
    logit(psi.apr) <- logit(int.psi) + r.psi[2]
    logit(psi.may) <- logit(int.psi) + r.psi[3]
    logit(psi.jun) <- logit(int.psi) + r.psi[4]
    logit(psi.jul) <- logit(int.psi) + r.psi[5]
    logit(psi.aug) <- logit(int.psi) + r.psi[6]
    logit(psi.sep) <- logit(int.psi) + r.psi[7]
    logit(psi.oct) <- logit(int.psi) + r.psi[8]

    logit(p.mar8) <- logit(int.p) + r.p[1]
    logit(p.apr8) <- logit(int.p) + r.p[2]
    logit(p.may8) <- logit(int.p) + r.p[3]
    logit(p.jun8) <- logit(int.p) + r.p[4]
    logit(p.jul8) <- logit(int.p) + r.p[5]
    logit(p.aug8) <- logit(int.p) + r.p[6]
    logit(p.sep8) <- logit(int.p) + r.p[7]
    logit(p.oct8) <- logit(int.p) + r.p[8]

    logit(p.mar10) <- logit(int.p) + beta.lp + r.p[1]
    logit(p.apr10) <- logit(int.p) + beta.lp + r.p[2]
    logit(p.may10) <- logit(int.p) + beta.lp + r.p[3]
    logit(p.jun10) <- logit(int.p) + beta.lp + r.p[4]
    logit(p.jul10) <- logit(int.p) + beta.lp + r.p[5]
    logit(p.aug10) <- logit(int.p) + beta.lp + r.p[6]
    logit(p.sep10) <- logit(int.p) + beta.lp + r.p[7]
    logit(p.oct10) <- logit(int.p) + beta.lp + r.p[8]

    det.2wk.mar10 <- 1-(1-p.mar10)^14
    det.2wk.apr10 <- 1-(1-p.apr10)^14
    det.2wk.may10 <- 1-(1-p.may10)^14
    det.2wk.jun10 <- 1-(1-p.jun10)^14
    det.2wk.jul10 <- 1-(1-p.jul10)^14
    det.2wk.aug10 <- 1-(1-p.aug10)^14
    det.2wk.sep10 <- 1-(1-p.sep10)^14
    det.2wk.oct10 <- 1-(1-p.oct10)^14

    det.2wk.mar8 <- 1-(1-p.mar8)^14
    det.2wk.apr8 <- 1-(1-p.apr8)^14
    det.2wk.may8 <- 1-(1-p.may8)^14
    det.2wk.jun8 <- 1-(1-p.jun8)^14
    det.2wk.jul8 <- 1-(1-p.jul8)^14
    det.2wk.aug8 <- 1-(1-p.aug8)^14
    det.2wk.sep8 <- 1-(1-p.sep8)^14
    det.2wk.oct8 <- 1-(1-p.oct8)^14

        } #End Model
    ")

# Parameters monitored
params <- c("int.psi","int.p","beta.lp","psi.mar","psi.apr","psi.may","psi.jun","psi.jul",
            "psi.aug","psi.sep","psi.oct","p.mar8","p.apr8","p.may8","p.jun8","p.jul8",
            "p.aug8","p.sep8","p.oct8","p.mar10","p.apr10","p.may10","p.jun10","p.jul10",
            "p.aug10","p.sep10","p.oct10","det.2wk.mar8","det.2wk.apr8","det.2wk.may8",
            "det.2wk.jun8","det.2wk.jul8","det.2wk.aug8","det.2wk.sep8","det.2wk.oct8",
            "det.2wk.mar10","det.2wk.apr10","det.2wk.may10",
            "det.2wk.jun10","det.2wk.jul10","det.2wk.aug10","det.2wk.sep10","det.2wk.oct10")

# MCMC settings
ni <- 120000; nt <- 30; nb <- 30000; nc <- 3

## load data
load("occ-weasel2017.RData")
attach(data)

win.data <- list(y=y, n.site=dim(y)[1],n.days=dim(y)[2],diameter=diameter,
                 month=month, n.month=n.month, loc=loc, n.loc=n.loc)
# Initial values
zst <- apply(y,1,max,na.rm=T)      #inits for presence (z)
inits <- function()list(z=zst, int.psi=0.5)
#outJ.w17month <- jags(win.data, inits, params, "occ-Weasel17month.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, parallel=T)
#save(outJ.w17month,file="outJ-w17monthv20200111.RData")
detach(data); rm("data")

###### Model for 2018 ######


# Define model in BUGS language
cat(file="occ-Weasel18month.txt", 
    "
    model{
    
    # Priors and model for params
    int.psi ~ dunif(0,1)        #Intercept for occupancy probability
    
    int.p ~ dunif(0,1)          #Intercept for detection probability
    
    # prior for random intercept
    sigma_r ~ dunif(0,100)
    tau_r <- 1/(sigma_r*sigma_r)
    for(j in 1:n.loc){
    r.psi[j] ~ dnorm(0,tau_r)
    }

    sigma_r2 ~ dunif(0,100)
    tau_r2 <- 1/(sigma_r2*sigma_r2)
    sigma_r3 ~ dunif(0,100)
    tau_r3 <- 1/(sigma_r3*sigma_r3)
    for(j in 1:n.month){
    r.psi2[j] ~ dnorm(0,tau_r2)
    r.p[j] ~ dnorm(0,tau_r3)
    }

    #Slope of p covariate ::  diameter
    beta.lp ~ dnorm(0,0.001)
    
    #Likelihood (or basic model structure)
    for(i in 1:n.site){
    
    #Occurrence in site i
    z[i] ~ dbern(psi[i])
    logit(psi[i]) <- logit(int.psi) + r.psi[loc[i]] + r.psi2[month[i]]
    
    for(k in 1:n.days){
    
    # detection probability on day k
    y[i,k] ~ dbern(mu.y[i,k])
    mu.y[i,k] <- z[i] * p[i,k]
    logit(p[i,k]) <- logit(int.p) + r.p[month[i]] +
    beta.lp * diameter[i]
    }
    }
    
    # derived parameters
    logit(psi.feb) <- logit(int.psi) + r.psi2[1]
    logit(psi.mar) <- logit(int.psi) + r.psi2[2]
    logit(psi.apr) <- logit(int.psi) + r.psi2[3]
    logit(psi.may) <- logit(int.psi) + r.psi2[4]
    logit(psi.jun) <- logit(int.psi) + r.psi2[5]
    logit(psi.jul) <- logit(int.psi) + r.psi2[6]
    logit(psi.aug) <- logit(int.psi) + r.psi2[7]
    logit(psi.sep) <- logit(int.psi) + r.psi2[8]
    
    logit(p.feb8) <- logit(int.p) + r.p[1]
    logit(p.mar8) <- logit(int.p) + r.p[2]
    logit(p.apr8) <- logit(int.p) + r.p[3]
    logit(p.may8) <- logit(int.p) + r.p[4]
    logit(p.jun8) <- logit(int.p) + r.p[5]
    logit(p.jul8) <- logit(int.p) + r.p[6]
    logit(p.aug8) <- logit(int.p) + r.p[7]
    logit(p.sep8) <- logit(int.p) + r.p[8]
    
    logit(p.feb10) <- logit(int.p) + beta.lp[1] + r.p[1]
    logit(p.mar10) <- logit(int.p) + beta.lp[1] + r.p[2]
    logit(p.apr10) <- logit(int.p) + beta.lp[1] + r.p[3]
    logit(p.may10) <- logit(int.p) + beta.lp[1] + r.p[4]
    logit(p.jun10) <- logit(int.p) + beta.lp[1] + r.p[5]
    logit(p.jul10) <- logit(int.p) + beta.lp[1] + r.p[6]
    logit(p.aug10) <- logit(int.p) + beta.lp[1] + r.p[7]
    logit(p.sep10) <- logit(int.p) + beta.lp[1] + r.p[8]
    
    det.2wk.feb10 <- 1-(1-p.feb10)^14
    det.2wk.mar10 <- 1-(1-p.mar10)^14
    det.2wk.apr10 <- 1-(1-p.apr10)^14
    det.2wk.may10 <- 1-(1-p.may10)^14
    det.2wk.jun10 <- 1-(1-p.jun10)^14
    det.2wk.jul10 <- 1-(1-p.jul10)^14
    det.2wk.aug10 <- 1-(1-p.aug10)^14
    det.2wk.sep10 <- 1-(1-p.sep10)^14
    
    det.2wk.feb8 <- 1-(1-p.feb8)^14
    det.2wk.mar8 <- 1-(1-p.mar8)^14
    det.2wk.apr8 <- 1-(1-p.apr8)^14
    det.2wk.may8 <- 1-(1-p.may8)^14
    det.2wk.jun8 <- 1-(1-p.jun8)^14
    det.2wk.jul8 <- 1-(1-p.jul8)^14
    det.2wk.aug8 <- 1-(1-p.aug8)^14
    det.2wk.sep8 <- 1-(1-p.sep8)^14
    
    } #End Model
    ")

# Parameters monitored
params <- c("int.psi","int.p","beta.lp","psi.feb","psi.mar","psi.apr","psi.may","psi.jun","psi.jul",
            "psi.aug","psi.sep","p.feb8","p.mar8","p.apr8","p.may8","p.jun8","p.jul8",
            "p.aug8","p.sep8","p.feb10","p.mar10","p.apr10","p.may10","p.jun10","p.jul10",
            "p.aug10","p.sep10","det.2wk.mar8","det.2wk.apr8","det.2wk.may8",
            "det.2wk.jun8","det.2wk.jul8","det.2wk.aug8","det.2wk.sep8","det.2wk.oct8",
            "det.2wk.mar10","det.2wk.apr10","det.2wk.may10",
            "det.2wk.jun10","det.2wk.jul10","det.2wk.aug10","det.2wk.sep10","det.2wk.oct10")

# MCMC settings
ni <- 120000; nt <- 30; nb <- 30000; nc <- 3

## load data
load("occ-weasel2018.RData")
attach(data)

win.data <- list(y=y, n.site=dim(y)[1],n.days=dim(y)[2],diameter=diameter,loc=loc,n.loc=n.loc,
                 month=month,n.month=n.month)
# Initial values
zst <- apply(y,1,max,na.rm=T)      #inits for presence (z)
inits <- function()list(z=zst, int.psi=0.5, beta.lpsi=0)
#outJ.w18month.2 <- jags(win.data, inits, params, "occ-Weasel18month.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, parallel=T)
#save(outJ.w18month.2,file="outJ-w18month-20200111.RData")

detach(data); rm("data")
