library(MASS)
library(gdata)
library(mclust)

files = list.files(path="../../R functions for VVV, with errors", pattern="*.R");
files = paste0("../../R functions for VVV, with errors/", files);
for (file in files) {
  source(file)
}


### Simulation 3: varying error structures

## Structure 1: identical errors across observations
simulate.data.1 = function(e, seed=0) {
                                        # set sample size
  nvec = c(3,3,4) * 5
  n = sum(nvec)
  p = 2
  G = 3

  tau = c(0.3,0.3,0.4)

  set.seed(seed)
  mu1 = c(1,1)
  mu2 = c(10,-10)
  mu3 = c(-25,25)
  mu = cbind(mu1,mu2,mu3)

                                        # set errors to be zero
  err = array(0, dim=c(p,p,n))  

  for(i in 1:n){
    err[,,i] = matrix(0,p,p)
    diag(err[,,i]) = e
  }

  sigma1 = matrix(c(15,-2,-2,15),nrow=2)
  sigma2 = matrix(c(23,3,3,23),nrow=2)
  sigma3 = matrix(c(31,-4,-4,31),nrow=2)

  set.seed(seed)
  s1 = mvrnorm(nvec[1], mu1, (sigma1+err[,,1]))
  s2 = mvrnorm(nvec[2], mu2, (sigma2+err[,,1]))
  s3 = mvrnorm(nvec[3], mu3, (sigma3+err[,,1]))

  # true membership matrix:
  temp = c(rep(c(1,0,0),nvec[1]),rep(c(0,1,0),nvec[2]),rep(c(0,0,1),nvec[3]))
  z.ini = matrix(temp, nrow=n, byrow=TRUE)

  samp = rbind(s1,s2,s3);
  
  ini.par = ini.par.iderr(samp,z.ini,err)

  list(z.ini=z.ini, samp = samp, err = err, ini.par=ini.par);
}





## Structure 2: errors the same within clusters
simulate.data.2 = function(seed=0) {
  nvec = c(3,3,4) * 5
  n = sum(nvec)
  p = 2
  G = 3

  tau = c(0.3,0.3,0.4)

  set.seed(seed)
  mu1 = c(1,1)
  mu2 = c(10,-10)
  mu3 = c(-25,25)
  mu = cbind(mu1,mu2,mu3)

  err = array(0, dim=c(p,p,n))
  err1 = matrix(c(0.1,0,0,0.1),nrow=2)
  err2 = matrix(c(0.3,0,0,0.3),nrow=2)
  err3 = matrix(c(0.5,0,0,0.5),nrow=2)

  for(i in 1:nvec[1]){
    err[,,i] = err1
  }

  for(i in (nvec[1]+1):(n-nvec[3])){
    err[,,i] = err2
  }

  for(i in (n-nvec[3]+1):n){
    err[,,i] = err3
  }
  
  sigma1 = matrix(c(15,-2,-2,15),nrow=2)
  sigma2 = matrix(c(23,3,3,23),nrow=2)
  sigma3 = matrix(c(31,-4,-4,31),nrow=2)
  sigma = array(0,dim=c(p,p,G))
  sigma[,,1] = sigma1
  sigma[,,2] = sigma2
  sigma[,,3] = sigma3
  
  cholsig = array(0,dim=c(p,p,G))
  cholsig[,,1] = chol(sigma1+err1)
  cholsig[,,2] = chol(sigma2+err2)
  cholsig[,,3] = chol(sigma3+err3)

  set.seed(seed)
  s1 = mvrnorm(nvec[1], mu1, (sigma1+err1))
  s2 = mvrnorm(nvec[2], mu2, (sigma2+err2))
  s3 = mvrnorm(nvec[3], mu3, (sigma3+err3))

  # true membership matrix:
  temp = c(rep(c(1,0,0),nvec[1]),rep(c(0,1,0),nvec[2]),rep(c(0,0,1),nvec[3]))
  z.ini = matrix(temp, nrow=n, byrow=TRUE)

  samp = rbind(s1,s2,s3)


  ini.par = ini.par.clust(samp,z.ini,err)
  
  list(z.ini=z.ini, samp = samp, err = err, ini.par=ini.par, tau=tau, mu=mu, sigma=sigma, cholsig=cholsig);
}





## Verify when the error matrices are identical, mclust and mclustme
## should give equivalent results
test.1 = function() {

  ## Simulate data with no error
  d = simulate.data.1(e=0.1); 
  
  z.me = z.mc = d$z.ini
  samp = d$samp
  err = d$err
  mu = d$mu
  sigma = d$sigma
  cholsig = d$cholsig
  tau = d$tau
  
  param.me = list(mu=mu,phat=tau,sigma=sigma)
  param.mc = list(pro=tau,mean=mu,variance=list(modelName="VVV",d=2,G=3,cholsigma=cholsig))
  
  
  ini.par = d$ini.par
  lik.me = lik.mc = numeric()
  
  k = 1
  
  ## step-by-step iteration for both methods:
  while(k < 51){
    
    ## Run mstep of mclustme
    res.me = MstepVVV.err(z.me, samp, err, ini.par);
    
    ## Run mstep of mclust
    res = mstepVVV(samp, z.mc);
    
    ## Compare means
    print(res$parameters$mean - res.me$muhat)
    
    ## Compare varaince convariance
    print(res$parameter$variance$sigma - res.me$sigma)
    
    
    param.me = res.me
    param.mc = res$parameters
    
 
    ## Run E-step of mclustme
    es.me = EstepVVV.err(res.me,samp,err)
    
    ## Run E-step of mclust
    es = estepVVV(samp,res$parameters)
    
    #es.me[[1]] - es$z
    ## Compare results
    
    lik.me[k] = es.me[[3]]
    lik.mc[k] = es$loglik
    
    print(c(lik.me[k],lik.mc[k]))
    
    z.me = es.me[[1]]
    z.mc = es$z
    
    k = k + 1
  }
  
  ## Plot the error for two methods:
  me.error = -log(diff(lik.me))
  mc.error = -log(diff(lik.mc))
  
  par(lwd=2.5)
  plot(me.error,type="l",lty="dashed",col="red",xlab="iteration",ylab="-log(error)",main="Convergence: Same Errors \n (step-by-step iteration)")
  lines(mc.error,col="green",lty="dotted")
  abline(h=-log(1e-5),lty="dashed",col="blue")
  legend(30,6, c("new method","mclust","-log(tolerance)"),lty=c("dashed","dotted","dashed"),lwd=c(3,3,3),col=c("red","green","blue"))
  
  ## Plot difference between errors:
  par(mfrow=c(1,1))
  plot(me.error-mc.error, type="l", lty="dashed",xlab="iteration",ylab="difference in error",col="green",main="Difference in Errors \n (log ratios mc/me)")
  abline(h=0,col="blue",lty="dashed") 
 
  
  ## Plot error for our method:
  me.result = ME.VVV.err(samp,z.me,err,"identical")
  likvec = me.result$likvec
  plot(-log(diff(likvec)),type="l",ylab="-log(error)",main="Convergence of EM Wrapper",col="red")
  abline(h=-log(1e-5),lty="dashed",col="blue")
  legend(22,4,"-log(tolerance)",lty="dashed",lwd=3,col="blue")
 
  
  mc.result = meVVV(samp,z.ini)
  
}






## Verify when the error matrices are identical within each cluster,
## mclust and mclustme should give equivalent results
test.2 = function() {
  
  ## Simulate data with no error
  d = simulate.data.2(); 
  
  z.me = z.mc = d$z.ini
  samp = d$samp
  err = d$err
  sigma = d$sigma
  cholsig = d$cholsig
  tau = d$tau
  mu = d$mu
  ini.par = d$ini.par
  
  param.me = list(mu=mu,phat=tau,sigma=sigma)
  param.mc = list(pro=tau,mean=mu,variance=list(modelName="VVV",d=2,G=3,cholsigma=cholsig))
    
  
  
  
  ## See if loglikelihood is same at each step:
  mye = EstepVVV.err(param.me,samp,err) # -337.6775
  mce = estepVVV(samp,param.mc) # -337.6746
  # not exactly the same (I wonder why...) Could it be that mclust computes the "complete" loglike?
  
  sig = mye[[2]]$sigmahat
  param = numeric()
  for(i in 1:3){
    param = c(param,upperTriangle(chol(sig[,,i]),diag=TRUE))  
  }  
  z = mye[[1]] 
  -1*obj.fun.VVV.err(param,z,samp,err) # -339.0867
  # So no, mclust is NOT using the "complete" loglike.
  
  
  
  
  ## step-by-step iteration for both methods:
  lik.me = lik.mc = numeric()
  k = 1  
  while(k < 51){
    ## Run mstep of mclustme
    res.me = MstepVVV.err(z.me, samp, err, ini.par);
    
    ## Run mstep of mclust
    res = mstepVVV(samp, z.mc);
    
    ## Compare means
    print(res$parameters$mean - res.me$muhat)
    
    ## Compare varaince convariance
    print(res$parameter$variance$sigma - res.me$sigma)
    
    
    
    ## Run E-step of mclustme
    es.me = EstepVVV.err(res.me,samp,err)
    
    ## Run E-step of mclust
    es = estepVVV(samp,res$parameters)
    
    #es.me[[1]] - es$z
    ## Compare results
    
    lik.me[k] = es.me[[3]]
    lik.mc[k] = es$loglik
    
    print(c(lik.me[k],lik.mc[k]))
    
    z.me = es.me[[1]]
    z.mc = es$z
    
    k = k + 1
  }
  
  ## Plot the error for two methods:
  me.error = -log(diff(lik.me))
  mc.error = -log(diff(lik.mc))
  
  par(mfrow=c(1,1))
  plot(me.error,type="l",lty="dashed",col="red",xlab="iteration",ylab="-log(error)",main="Convergence: Same Error within Clusters \n (step-by-step iteration)")
  lines(mc.error,col="green",lty="dotted")
  abline(h=-log(1e-5),lty="dashed",col="blue")
  legend(26,6, c("new method","mclust","-log(tolerance)"),lty=c("dashed","dotted","dashed"),lwd=c(3,3,3),col=c("red","green","blue"))
  
  ## Difference between errors:
  plot(me.error-mc.error, type="l", lty="dashed",xlab="iteration",ylab="log(mc/me)",col="green",main="Difference in Errors \n (log ratios mc/me)")
  abline(h=0,lty="dashed",col="blue")
    
  ## Plot error for our method:
  me.result = ME.VVV.err(samp,z.me,err,"cluster")
  likvec = me.result$likvec
  plot(-log(diff(likvec)),type="l",ylab="-log(error)",main="Convergence of EM Wrapper",col="red")
  abline(h=-log(1e-5),lty="dashed",col="blue")
  legend(22,4,"-log(tolerance)",lty="dashed",lwd=3,col="blue")
  
  par(mfrow=c(1,2))
  mc.result = meVVV(samp,z.ini)
  
}



