

library(MASS)
library(mclust)
library(gdata)
library(phyclust)

N = 200
tau = 0.5
mu1 = c(0,0)
mu2 = c(8,0)
sig1 = matrix(c(64,0,0,64),nrow=2)
sig2 = matrix(c(16,0,0,16),nrow=2)
k = 9
p = 0.1
nseed = 100

#check.seed(N,tau,mu1,mu2,sig1,sig2,k,p)

setwd("Z://Research Project/simulations/Sim/example 3/p=0.1")

ptm <- proc.time()
sim.driver(N,tau,mu1,mu2,sig1,sig2,k,p,nseed)
t01 = proc.time() - ptm






N = 200
tau = 0.5
mu1 = c(0,0)
mu2 = c(8,0)
sig1 = matrix(c(64,0,0,64),nrow=2)
sig2 = matrix(c(16,0,0,16),nrow=2)
k = 9
p = 0.3
nseed = 100

#check.seed(N,tau,mu1,mu2,sig1,sig2,k,p)
  
setwd("Z://Research Project/simulations/Sim/example 3/p=0.3")

ptm <- proc.time()
sim.driver(N,tau,mu1,mu2,sig1,sig2,k,p,nseed)
t03 = proc.time() - ptm






N = 200
tau = 0.5
mu1 = c(0,0)
mu2 = c(8,0)
sig1 = matrix(c(64,0,0,64),nrow=2)
sig2 = matrix(c(16,0,0,16),nrow=2)
k = 9
p = 0.7
nseed = 100

#check.seed(N,tau,mu1,mu2,sig1,sig2,k,p)

setwd("Z://Research Project/simulations/Sim/example 3/p=0.7")

ptm <- proc.time()
sim.driver(N,tau,mu1,mu2,sig1,sig2,k,p,nseed)
t07 = proc.time() - ptm






N = 200
tau = 0.5
mu1 = c(0,0)
mu2 = c(8,0)
sig1 = matrix(c(64,0,0,64),nrow=2)
sig2 = matrix(c(16,0,0,16),nrow=2)
k = 9
p = 0.9
nseed = 100

#check.seed(N,tau,mu1,mu2,sig1,sig2,k,p)

setwd("Z://Research Project/simulations/Sim/example 3/p=0.9")

ptm <- proc.time()
sim.driver(N,tau,mu1,mu2,sig1,sig2,k,p,nseed)
t09 = proc.time() - ptm



