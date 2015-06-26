


          #---------------------Cluster Simulation 5--------------------------#

          ### One cluster has random measurement error, the other does not.
          ### Errors are uniformly distributed.


library(mclust)
library(gdata)
library(MASS)
library(phyclust)

col.blue = rgb(30,144,255,max=255)


sim5 = function(seed,n){

  n = n
  set.seed(seed)
  
## Simulate two clusters far away from each other.
  mu1 = c(-7,0)
  mu2 = c(7,0)
  sig = matrix(c(36,0,0,36),nrow=2) 
  I = matrix(c(1,0,0,1),nrow=2)

  # Sample size:
  n1 = n2 = n
  s = n1 + n2

  # Generate samples:
  mean1 = mvrnorm(n1, mu1, sig)
  mean2 = mvrnorm(n2, mu2, sig)
   
  # Generate random errors:
  k = runif(s, 0, 25)
  dat1 = matrix(0,n1,2)
  for(i in 1:n1){
    sigma = k[i]*I
    dat1[i,] = mvrnorm(1,mean1[i,],sigma)
  }
  
  dat2 = mean2
  
  # Combine two groups:
  dat = rbind(dat1, dat2)
  
  
  xlow = min(dat[,1])-0.5
  xup = max(dat[,1])+0.5
  ylow = min(dat[,2])-0.2
  yup = max(dat[,2])+0.2
  
  #par(mfrow=c(2,2))
  
  plot(dat1,xlim = c(xlow,xup), ylim = c(ylow,yup),col=col.blue,pch=19
    ,xlab="",ylab="",main="Generated Clusters")
  points(dat2,col="red",pch=0)
  # Observe that the two clusters overlap.
  
  

  #-------------------------------------------------------------#




  # Construct initial classification. Use true membership:
  z.ini = matrix(c(rep(c(1,0),n1),rep(c(0,1),n2)), byrow=T, nrow=s)

  # Construct the true error matrices:
  errmat = array(0,dim=c(2,2,s))
  
  for(i in 1:n1){
    errmat[,,i] = k[i]*I
  }


  
  #-------------------------------------------------------------#



  # Use meVVV
  res.me = meVVV(dat,z.ini)
  #names(res.me)
  m1 = res.me$z

  t1 = nrow(m1)
  G1 = ncol(m1)
  
  z1 = matrix(0,t1,G1)
  for(i in 1:t1){
    rowmax = max(m1[i,])
    for(k in 1:G1){
      z1[i,k] = ifelse(m1[i,k]==rowmax,2,1)
    }  
  }
    
  group1.new1 = cbind(dat1,z1[1:n1,])
  group2.new1 = cbind(dat2,z1[(n1+1):s,])
  
  plot(dat1, xlim = c(xlow,xup), ylim = c(ylow,yup), pch=ifelse(group1.new1[,3]==2,19,0),
    col=ifelse(group1.new1[,3]==2,col.blue,"red") ,xlab="", ylab="",main="meVVV")
  points(dat2, pch=ifelse(group2.new1[,3]==2,19,0),
    col=ifelse(group2.new1[,3]==2,col.blue,"red"))
  

  
  
    #-------------------------------------------------------------#

  
  
  
  ## Use MCME:
  my.result = mcmeVVV(dat, z.ini, errmat, "none")
  m2 = my.result$z
  
  t2 = nrow(m2)
  G2 = ncol(m2)
  
  for(i in 1:t2){
    rowmax = max(m2[i,])
    for(k in 1:G2){
      z2[i,k] = ifelse(m2[i,k]==rowmax,2,1)
    }  
  }
    
  group1.new2 = cbind(dat1,z2[1:n1,])
  group2.new2 = cbind(dat2,z2[(n1+1):s,])
  
  plot(dat1, xlim = c(xlow,xup), ylim = c(ylow,yup), pch=ifelse(group1.new2[,3]==2,19,0),
    col=ifelse(group1.new2[,3]==2,col.blue,"red") ,xlab="", ylab="",main="MCME")
  points(dat2, pch=ifelse(group2.new2[,3]==2,19,0),
    col=ifelse(group2.new2[,3]==2,col.blue,"red"))

  #par(mfrow=c(1,1))



  #-----------------------------------------------------------------------#

  
  
  z.ini1 = z.ini + 1  
  group1.new = cbind(dat1,z1[1:n1,1],z2[1:n1,1],z.ini1[1:n1,1])
  group2.new = cbind(dat2,z1[(n1+1):s,1],z2[(n1+1):s,1],z.ini1[(n1+1):s,1])
  
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
  

  
    #-----------------------------------------------------------------------#

  
  
  ## Compute Rand index
  rand1.me = RRand(z.ini1[1:n1,],z1[1:n1,])
  rand1.nn = RRand(z.ini1[(n1+1):s,],z1[(n1+1):s,])
  
  rand2.me = RRand(z.ini1[1:n1,],z2[1:n1,])
  rand2.nn = RRand(z.ini1[(n1+1):s,],z2[(n1+1):s,])
  
  rand1 = RRand(z.ini1,z1)
  rand2 = RRand(z.ini1,z2)
  
  

  #-----------------------------------------------------------------------#



  ## Misclassification rate
  tab1 = table(z.ini,z1)
  tab2 = table(z.ini,z2)
  mis1 = 1-sum(diag(tab1))/sum(tab1)
  mis2 = 1-sum(diag(tab2))/sum(tab2)

  
  out = list(seed=seed,n=n,rand1=rand1$Rand,rand1.me=rand1.me$Rand,rand1.nn=rand1.nn$Rand,
    rand2=rand2$Rand,rand2.me=rand2.me$Rand,rand2.nn=rand2.nn$Rand,mis1=mis1,mis2=mis2,
    m1=m1,m2=m2)
  return(out)
  
  
  
}




  #-----------------------------------------------------------------------#




seedvec = c(91,92,93,94,95)

png("50-R.png",900,1500)
par(mfrow=c(5,3))
res50_r = sapply(seedvec,sim5,n=50)
par(mfrow=c(1,1))
dev.off()


png("100-R.png",900,1500)
par(mfrow=c(5,3))
res100_r = sapply(seedvec,sim5,n=100)
par(mfrow=c(1,1))
dev.off()


png("200-R.png",900,1500)
par(mfrow=c(5,3))
res200_r = sapply(seedvec,sim5,n=200)
par(mfrow=c(1,1))
dev.off()


png("500-R.png",900,1500)
par(mfrow=c(5,3))
res500_r = sapply(seedvec,sim5,n=500)
par(mfrow=c(1,1))
dev.off()


png("800-R.png",900,1500)
par(mfrow=c(5,3))
res800_r = sapply(seedvec,sim5,n=800)
par(mfrow=c(1,1))
dev.off()


out = list(res50=res50_r,res100=res100_r,res200=res200_r,res500=res500_r,res800=res800_r)
save(out,file="sim5res")


