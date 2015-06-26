
N = 1000             

#Sample N random uniforms U
U =runif(N)

#Variable to store the samples from the mixture distribution                                             
rand.samples = matrix(rep(NA,2*N),nrow=2)

mu1 = c(0,0)
mu2 = c(1,1)
mu3 = c(-1,-1)

set.seed(0)
sig1 = matrix(c(2,0,0,2),nrow=2)
sig2 = matrix(c(4,-2,-2,4),nrow=2)
sig3 = matrix(c(4,2,2,4),nrow=2)

library(MASS)
#Sampling from the mixture
for(i in 1:N){
    if(U[i]<.3){
        rand.samples[,i] = mvrnorm(1,mu1,sig1)
    }else if(U[i]<.8){
        rand.samples[,i] = mvrnorm(1,mu2,sig2)
    }else{
        rand.samples[,i] = mvrnorm(1,mu3,sig3)
    }
}

plot(rand.samples[1,],rand.samples[2,],cex=0.5,xlab="",ylab="")


par(mfrow=c(1,3))

