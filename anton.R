library(mvtBinaryEP)
set.seed(1001)

setwd("/users/dustin/desktop/") 
load("Data.Sim.50.RData") 
mydata <- Data.Sim[[1]][[1]]

## Estimated from previous model
mu.pm <- 47.71
beta.pm <- 8.96
sigma.y <- 22.25
sigma.alpha <- 11.45
sigma.gamma <- 1.44
n <- 

## Generate some data
x.l <- ep(mu=0.5, rho=0.95, n=1386, nRep=90)[[1]]
delta <- c(10,9,8,7,6,5,4,3,2,1,rep(0,1376))
gamma.l.mean <- vector(length=90)
gamma.l <- vector(length=90)
for (i in 1:90) gamma.l.mean[i] <- x.l[i,]%*%delta
for (i in 1:90) gamma.l[i] <- rnorm(1,gamma.l[i],sigma.gamma)
