


### Checking whether meVVV is using loglikelihood increment as stopping criterion ###


library(mclust)
library(gdata)
library(MASS)

## Generate two clusters
  mu1 = c(-12,-3)
  mu2 = c(12,3)
  sig1 = matrix(c(25,0,0,4),nrow=2)
  sig2 = matrix(c(4,0,0,25),nrow=2)
  
  # Sample size:
  n1 = 500
  n2 = 500
  s = n1 + n2

  # Construct initial classification. Use true membership:
  z.ini = matrix(c(rep(c(1,0),n1),rep(c(0,1),n2)), byrow=T, nrow=s)

  # Generate random errors:
  g = 15

  # Generate samples:
  mean1 = mvrnorm(n1, mu1, sig1)
  mean2 = mvrnorm(n2, mu2, sig2)
  
  # Combine two groups:
  truemean = rbind(mean1, mean2)
  
  
  # Observations with same errors:
  DAT = matrix(0,s,2)
  ID = matrix(c(1,0,0,1),nrow=2)
  for(i in 1:s){
    sig = g*ID
    DAT[i,] = mvrnorm(1,truemean[i,],sig)
  }
  


  #------------------------------------------------------------------------------------#



ONE = meVVV(DAT,z.ini,control=emControl(itmax=1))
#> ONE$loglik
#[1] -6859.297

TWO = meVVV(DAT,z.ini,control=emControl(itmax=2))
#> TWO$loglik
#[1] -6859.227

THREE = meVVV(DAT,z.ini,control=emControl(itmax=3))
#> THREE$loglik
#[1] -6859.218


 res = meVVV(DAT,z.ini)   
 attributes(res)

#$info
# iterations       error 
#3.00000e+00 1.32887e-06 # This number is clearly not = loglik[3]-loglik[2] from above.

#$returnCode
#[1] 0

#> res$loglik
#[1] -6859.218

#> res$loglik==THREE$loglik
#[1] TRUE
> 