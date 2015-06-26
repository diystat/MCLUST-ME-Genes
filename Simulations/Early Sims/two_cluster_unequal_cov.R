

library(mclust)
library(gdata)
library(MASS)
library(phyclust)

col.blue = rgb(30,144,255,max=255)


sim8 = function(seed,n){

## Simulate two clusters far away from each other.
  mu1 = c(-5,0)
  mu2 = c(5,0)
  sig1 = matrix(c(16,0,0,16),nrow=2)
  sig2 = matrix(c(36,0,0,36),nrow=2)  
  I = matrix(c(1,0,0,1),nrow=2)

  # Sample size:
  n1 = n2 = n = 200
  s = n1 + n2
  
  seed = 92
  
  # Generate measurement errors:
  set.seed(seed); k = runif(s,0,3.9)

  # Generate means:
  set.seed(seed); mean1 = mvrnorm(n1, mu1, sig1)
  set.seed(seed); mean2 = mvrnorm(n2, mu2, sig2)
  
  # Combine two groups:
  truemean = rbind(mean1, mean2)
  dat = matrix(0,s,2)
  for(i in 1:s){
    sig = (k[i])^2*I
    set.seed(seed); dat[i,] = mvrnorm(1,truemean[i,],sig)
  }
  
  dat1 = dat[1:n1,]
  dat2 = dat[(n1+1):s,]
  
  # Plot generated clusters:
  xlow = min(dat[,1])-0.5
  xup = max(dat[,1])+0.5
  ylow = min(dat[,2])-0.2
  yup = max(dat[,2])+0.2
  
  par(mfrow=c(1,3))
  
  plot(dat1,xlim = c(xlow,xup), ylim = c(ylow,yup),col=col.blue,pch=19
    ,xlab="",ylab="",main="Generated Clusters",cex=1)
  points(dat2,col="red",pch=0,cex=1)
  
  

  #-------------------------------------------------------------#




  # Construct initial classification. Use true membership:
  z.ini = matrix(c(rep(c(1,0),n1),rep(c(0,1),n2)), byrow=T, nrow=s)

  # Construct the true error matrices:
  errmat = array(0,dim=c(2,2,s))
  
  for(i in 1:s){
    errmat[,,i] = (k[i])^2*I
  }


  
  #-------------------------------------------------------------#



  # Use meVVV
  res.me = meVVV(dat,z.ini)
  q1 = res.me$z

  t = nrow(q1)
  G = ncol(q1)
  
  z1 = matrix(0,t,G)
  for(i in 1:t){
    rowmax = max(q1[i,])
    for(k in 1:G){
      z1[i,k] = ifelse(q1[i,k]==rowmax,2,1)
    }  
  }
    
  group1.new1 = cbind(dat1,z1[1:n1,])
  group2.new1 = cbind(dat2,z1[(n1+1):s,])
  
  plot(dat1, xlim = c(xlow,xup), ylim = c(ylow,yup), pch=ifelse(group1.new1[,3]==2,19,0),
    col=ifelse(group1.new1[,3]==2,col.blue,"red") ,xlab="", ylab="",main="meVVV",cex=1)
  points(dat2, pch=ifelse(group2.new1[,3]==2,19,0),
    col=ifelse(group2.new1[,3]==2,col.blue,"red"),cex=1)
 
   
  
    #-------------------------------------------------------------#

  
  
  
  ## Use MCME:
  my.result = mcmeVVV(dat, z.ini, errmat)
  q2 = my.result$z
  
  t = nrow(q2)
  G = ncol(q2)
  
  z2 = matrix(0,t,G)
  for(i in 1:t){
    rowmax = max(q2[i,])
    for(k in 1:G){
      z2[i,k] = ifelse(q2[i,k]==rowmax,2,1)
    }  
  }
   
  group1.new2 = cbind(dat1,z2[1:n1,])
  group2.new2 = cbind(dat2,z2[(n1+1):s,])
  
  plot(dat1, xlim = c(xlow,xup), ylim = c(ylow,yup), pch=ifelse(group1.new2[,3]==2,19,0),
    col=ifelse(group1.new2[,3]==2,col.blue,"red") ,xlab="", ylab="",main="MCME",cex=1)
  points(dat2, pch=ifelse(group2.new2[,3]==2,19,0),
    col=ifelse(group2.new2[,3]==2,col.blue,"red"),cex=1)
  
  par(mfrow=c(1,1))
  
  
  
    #-----------------------------------------------------------------------#


  
  
  
  ### Plot both clustering results
  ### Highlight points that are grouped differently,
  ### make all other points gray.
  
  z.ini1 = z.ini + 1  
  group1.new = cbind(dat1,z1[1:n1,1],z2[1:n1,1],z.ini1[1:n1,1])
  group2.new = cbind(dat2,z1[(n1+1):s,1],z2[(n1+1):s,1],z.ini1[(n1+1):s,1])

  
  
  par(mfrow=c(3,1))
  
  plot(dat1,xlim = c(xlow,xup), ylim = c(ylow,yup),col=col.blue,pch=19
    ,xlab="",ylab="",main="Generated Clusters")
  points(dat2,col="red",pch=0)
  
  
  plot(dat1, xlim = c(xlow,xup), ylim = c(ylow,yup), pch=ifelse(group1.new[,3]==2,19,0),
    col=ifelse(group1.new[,3]==group1.new[,4],"gray",
      ifelse((group1.new[,3]==2 & group1.new[,5]==2),col.blue,
        ifelse((group1.new[,3]==2 & group1.new[,5]==1),"red",
          ifelse((group1.new[,3]==1 & group1.new[,5]==1),"red",col.blue)))),
    xlab="", ylab="",main="meVVV")  
  points(dat2, pch=ifelse(group2.new[,3]==2,19,0),
    col=ifelse(group2.new[,3]==group2.new[,4],"gray",
      ifelse((group2.new[,3]==2 & group2.new[,5]==2),col.blue,
        ifelse((group2.new[,3]==2 & group2.new[,5]==1),"red",
          ifelse((group2.new[,3]==1 & group2.new[,5]==1),"red",col.blue)))))
  
  
  plot(dat1, xlim = c(xlow,xup), ylim = c(ylow,yup), pch=ifelse(group1.new[,4]==2,19,0),
    col=ifelse(group1.new[,3]==group1.new[,4],"gray",
      ifelse((group1.new[,4]==2 & group1.new[,5]==2),col.blue,
        ifelse((group1.new[,4]==2 & group1.new[,5]==1),"red",
          ifelse((group1.new[,4]==1 & group1.new[,5]==1),"red",col.blue)))),
    xlab="", ylab="",main="MCME")
  points(dat2, pch=ifelse(group2.new[,4]==2,19,0),
    col=ifelse(group2.new[,3]==group2.new[,4],"gray",
      ifelse((group2.new[,4]==2 & group2.new[,5]==2),col.blue,
        ifelse((group2.new[,4]==2 & group2.new[,5]==1),"red",
          ifelse((group2.new[,4]==1 & group2.new[,5]==1),"red",col.blue)))))
  
   par(mfrow=c(1,1))

  #-----------------------------------------------------------------------#


  
  ## Compute Rand index
  z.ini = z.ini + 1
  rand1.me = RRand(z.ini[1:n1,],z1[1:n1,])
  rand1.nn = RRand(z.ini[(n1+1):s,],z1[(n1+1):s,])
  
  rand2.me = RRand(z.ini[1:n1,],z2[1:n1,])
  rand2.nn = RRand(z.ini[(n1+1):s,],z2[(n1+1):s,])
  
  rand1 = RRand(z.ini,z1)
  rand2 = RRand(z.ini,z2)
  z.ini = z.ini - 1
  
  

  #-----------------------------------------------------------------------#



  ## Misclassification rate
  z.ini = z.ini + 1
  tab1 = table(z.ini,z1)
  tab2 = table(z.ini,z2)
  mis1 = 1-sum(diag(tab1))/sum(tab1)
  mis2 = 1-sum(diag(tab2))/sum(tab2)
  z.ini = z.ini - 1

  
  out = list(seed=seed,n=n,rand1=rand1$Rand,rand1.me=rand1.me$Rand,rand1.nn=rand1.nn$Rand,
    rand2=rand2$Rand,rand2.me=rand2.me$Rand,rand2.nn=rand2.nn$Rand,mis1=mis1,mis2=mis2)
  return(out)
  
  
  
}

