

      #---------------------Cluster Simulation 4--------------------------#
## When measurement errors are generated from a random distribution,
## will the results be the same as when errors equal to their expectation?

library(mclust)
library(gdata)
library(MASS)
library(phyclust)

col.blue = rgb(30,144,255,max=255)

## Simulate two clusters far away from each other.
  mu1 = c(-7,0)
  mu2 = c(7,0)
  sig1 = matrix(c(36,0,0,36),nrow=2)
  sig2 = matrix(c(36,0,0,36),nrow=2)
  
  # Sample size:
  n1 = 100
  n2 = 100
  s = n1 + n2

  # Construct initial classification. Use true membership:
  z.ini = matrix(c(rep(c(1,0),n1),rep(c(0,1),n2)), byrow=T, nrow=s)

  # Generate random errors:
  k = runif(s, 0, 25) # expectation = 12.5
  g = 12.5

  # Generate samples:
  set.seed(1)
  mean1 = mvrnorm(n1, mu1, sig1)
  mean2 = mvrnorm(n2, mu2, sig2)
  
  # Combine two groups:
  truemean = rbind(mean1, mean2)
  
  
  # Observations with random errors:
  dat1 = matrix(0,n1,2)
  id = matrix(c(1,0,0,1),nrow=2)
  for(i in 1:n1){
    sig = k[i]*id
    dat1[i,] = mvrnorm(1,truemean[i,],sig)
  }
  
  dat2 = mean2
  dat = rbind(dat1,dat2)
  

  # Observations with same errors:
  DAT1 = matrix(0,n1,2)
  ID = matrix(c(1,0,0,1),nrow=2)
  for(i in 1:n1){
    sig = g*ID
    DAT1[i,] = mvrnorm(1,truemean[i,],sig)
  }
  
  DAT2 = mean2
  DAT = rbind(DAT1,DAT2)

  
  
  xlow = min(dat[,1],DAT[,1])-0.5
  xup = max(dat[,1],DAT[,1])+0.5
  ylow = min(dat[,2],DAT[,2])-0.2
  yup = max(dat[,2],DAT[,2])+0.2
  
  par(mfrow=c(2,2))
  plot(dat1,xlim = c(xlow,xup), ylim = c(ylow,yup),col=col.blue,pch=19
    ,xlab="",ylab="",main="Uniform Errors")
  points(dat2,col="red",pch=0)
  
  plot(DAT1,xlim = c(xlow,xup), ylim = c(ylow,yup),col=col.blue,pch=19
    ,xlab="",ylab="",main="Same Errors")
  points(DAT2,col="red",pch=0)


 #--------------------------------------------------------------------------#


  
  ## When errors are randomly distributed:
  # Construct the true error matrices:
  errmat = array(0,dim=c(2,2,s))
  
  for(i in 1:n1){
    errmat[,,i] = k[i]*id
  }
  
  # Perform our clustering method:
  tmp = proc.time()
  my.result1 = mcmeVVV(dat, z.ini, errmat)
  mytime1 = proc.time() - tmp

  names(my.result1)
  gr1 = my.result1$z
  
  t = nrow(gr1)
  G = ncol(gr1)
  
  for(i in 1:t){
    rowmax = max(gr1[i,])
    for(k in 1:G){
      gr1[i,k] = ifelse(gr1[i,k]==rowmax,2,1)
    }  
  }
    
  group1.new0 = cbind(dat1,gr1[1:n1,])
  group2.new0 = cbind(dat2,gr1[(n1+1):s,])
  
  plot(dat1, xlim = c(xlow,xup), ylim = c(ylow,yup), pch=ifelse(group1.new0[,3]==2,19,0),
    col=ifelse(group1.new0[,3]==2,col.blue,"red"),xlab="",ylab="",main="MCME, Uniform")
  points(dat2, pch=ifelse(group2.new0[,3]==2,19,0),
    col=ifelse(group2.new0[,3]==2,col.blue,"red"))
  
  

  #------------------------------------------------------------------------#


  
  ## When errors equal to their expectation:
  ERRMAT = array(0,dim=c(2,2,s))
  
  for(i in 1:n1){
    ERRMAT[,,i] = g*id
  }
  
  # Perform our clustering method:
  tmp = proc.time()
  my.result2 = mcmeVVV(DAT, z.ini, ERRMAT)
  mytime2 = proc.time() - tmp

  names(my.result2)
  gr2 = my.result2$z
  
  t = nrow(gr2)
  G = ncol(gr2)
  
  for(i in 1:t){
    rowmax = max(gr2[i,])
    for(k in 1:G){
      gr2[i,k] = ifelse(gr2[i,k]==rowmax,2,1)
    }  
  }
    
  group1.new = cbind(DAT1,gr2[1:n1,])
  group2.new = cbind(DAT2,gr2[(n1+1):s,])
  
  plot(DAT1, xlim = c(xlow,xup), ylim = c(ylow,yup), pch=ifelse(group1.new[,3]==2,19,0),
    col=ifelse(group1.new[,3]==2,col.blue,"red") ,xlab="", ylab="",main="MCME, Same")
  points(DAT2, pch=ifelse(group2.new[,3]==2,19,0),
    col=ifelse(group2.new[,3]==2,col.blue,"red"))
 
  par(mfrow=c(1,1))
  
  
  
  #------------------------------------------------------------------------#
  
  
  ## Compute Rand index
  z.ini = z.ini + 1
  r1 = RRand(z.ini,gr1);r1
  r2 = RRand(z.ini,gr2);r2
  z.ini = z.ini - 1
  
  #------------------------------------------------------------------#
  

  ## How different are the clustering results in terms of membership?
  mbdiff = (my.result1$z-my.result2$z)[,1]
  max(abs(mbdiff))

  lbdiff = (gr1-gr2)[,1]
  mean(lbdiff!=0)


