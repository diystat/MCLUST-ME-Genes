

      #---------------------Cluster Simulation 6--------------------------#
      
      ### Asymmetric clusters with random measurement errors.


library(mclust)
library(gdata)
library(MASS)
library(phyclust)

col.blue = rgb(30,144,255,max=255)


sim6 = function(seed,n1,n2,l1,l2){

  set.seed(seed)
  n1 = n1
  n2 = n2
  I = matrix(c(1,0,0,1),nrow=2)
  sig1 = (l1^2)*I
  sig2 = (l2^2)*I

  ## Simulate two clusters:
  mu1 = c(-7,0)
  mu2 = c(7,0)
  
  # Sample size:
  s = n1 + n2

  # Generate random errors:
  k = runif(s, 0, 25)

  # Generate samples:
  mean1 = mvrnorm(n1, mu1, sig1)
  mean2 = mvrnorm(n2, mu2, sig2)
  
  # Combine two groups:
  truemean = rbind(mean1, mean2)
  
  dat = matrix(0,s,2)
  id = matrix(c(1,0,0,1),nrow=2)
  for(i in 1:s){
    sig = k[i]*id
    dat[i,] = mvrnorm(1,truemean[i,],sig)
  }
  
  dat1 = dat[1:n1,]
  dat2 = dat[-(1:n1),]
  
  
  xlow = min(dat[,1])-0.5
  xup = max(dat[,1])+0.5
  ylow = min(dat[,2])-0.2
  yup = max(dat[,2])+0.2
  
  #par(mfrow=c(3,1))
  plot(dat1,xlim = c(xlow,xup), ylim = c(ylow,yup),col=col.blue,pch=19
    ,xlab="",ylab="",main="Generated Clusters",xaxt="n",yaxt="n")
  points(dat2,col="red",pch=0)
  


  #-----------------------------------------------------------------------------#
  



  # Construct initial classification. Use true membership:
  z.ini = matrix(c(rep(c(1,0),n1),rep(c(0,1),n2)), byrow=T, nrow=s)

  # Construct the true error matrices:
  errmat = array(0,dim=c(2,2,s))
  for(i in 1:s){
    errmat[,,i] = k[i]*id
  }
  
  


  #-----------------------------------------------------------------------------#


  
  ## meVVV clustering by MCLUST
  
  res1 = meVVV(dat, z.ini)
  m1 = res1$z
  
  t1 = nrow(m1)
  G1 = ncol(m1)
  
  z1 = matrix(0,t1,G1)
  for(i in 1:t1){
    rowmax = max(m1[i,])
    for(j in 1:G1){
      z1[i,j] = ifelse(m1[i,j]==rowmax,2,1)
    }  
  }
    
  
  
  #-------------------------------------------------------------#


  
  ## Perform MCME:
  my.result = mcmeVVV(dat, z.ini, errmat)
  m2 = my.result$z
  
  t2 = nrow(m2)
  G2 = ncol(m2)
  
  z2 = matrix(0,t2,G2)
  for(i in 1:t2){
    rowmax = max(m2[i,])
    for(k in 1:G2){
      z2[i,k] = ifelse(m2[i,k]==rowmax,2,1)
    }  
  }
    
  
  
  
    #------------------------------------------------------------------------#


  
  
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
  
  
  
  #------------------------------------------------------------------------#
  
  
  
  ## Compute Rand index
  rand1.me = RRand(z.ini1[1:n1,],z1[1:n1,])
  rand1.nn = RRand(z.ini1[(n1+1):s,],z1[(n1+1):s,])
  
  rand2.me = RRand(z.ini1[1:n1,],z2[1:n1,])
  rand2.nn = RRand(z.ini1[(n1+1):s,],z2[(n1+1):s,])
  
  rand1 = RRand(z.ini1,z1)
  rand2 = RRand(z.ini1,z2)

  
  
  #------------------------------------------------------------------#
  
  
  
  ## Misclassification rate
  tab1 = table(z.ini1,z1)
  tab2 = table(z.ini1,z2)
  mis1 = 1-sum(diag(tab1))/sum(tab1)
  mis2 = 1-sum(diag(tab2))/sum(tab2)

  
  out = list(seed=seed,n=n,rand1=rand1$Rand,rand1.me=rand1.me$Rand,rand1.nn=rand1.nn$Rand,
    rand2=rand2$Rand,rand2.me=rand2.me$Rand,rand2.nn=rand2.nn$Rand,mis1=mis1,mis2=mis2,
    m1=m1,m2=m2)
  return(out)
  
  
  
}  
  
  
  


  #-----------------------------------------------------------------------------#



ff = function(){
  
seedvec = c(91,92,93,94,95)



## Case 1: Same covariance, unequal sample sizes.

png("50-size.png",900,1500)
par(mfrow=c(5,3))
res50_rr = sapply(seedvec,sim6,n1=100,n2=50,l1=4,l2=4)
par(mfrow=c(1,1))
dev.off()  

png("100-size.png",900,1500)
par(mfrow=c(5,3))
res100_rr = sapply(seedvec,sim6,n1=200,n2=100,l1=4,l2=4)
par(mfrow=c(1,1))
dev.off()


png("200-size.png",900,1500)
par(mfrow=c(5,3))
res200_rr = sapply(seedvec,sim6,n1=400,n2=200,l1=4,l2=4)
par(mfrow=c(1,1))
dev.off()


png("500-size.png",900,1500)
par(mfrow=c(5,3))
res500_rr = sapply(seedvec,sim6,n1=1000,n2=500,l1=4,l2=4)
par(mfrow=c(1,1))
dev.off()


png("800-size.png",900,1500)
par(mfrow=c(5,3))
res800_rr = sapply(seedvec,sim6,n1=1600,n2=800,l1=4,l2=4)
par(mfrow=c(1,1))
dev.off()


out1 = list(res50=res50_rr,res100=res100_rr,res200=res200_rr,res500=res500_rr,res800=res800_rr)
save(out1,file="sim6res_unequal_size")




  #-----------------------------------------------------------------------------#




## Case 2: Same size, unequal covariances.

png("50-cov.png",900,1500)
par(mfrow=c(5,3))
res50_rr = sapply(seedvec,sim6,n1=50,n2=50,l1=6,l2=3)
par(mfrow=c(1,1))
dev.off()  

png("100-cov.png",900,1500)
par(mfrow=c(5,3))
res100_rr = sapply(seedvec,sim6,n1=100,n2=100,l1=6,l2=3)
par(mfrow=c(1,1))
dev.off()


png("200-cov.png",900,1500)
par(mfrow=c(5,3))
res200_rr = sapply(seedvec,sim6,n1=200,n2=200,l1=6,l2=3)
par(mfrow=c(1,1))
dev.off()


png("500-cov.png",900,1500)
par(mfrow=c(5,3))
res500_rr = sapply(seedvec,sim6,n1=500,n2=500,l1=6,l2=3)
par(mfrow=c(1,1))
dev.off()


png("800-cov.png",900,1500)
par(mfrow=c(5,3))
res800_rr = sapply(seedvec,sim6,n1=800,n2=800,l1=6,l2=3)
par(mfrow=c(1,1))
dev.off()


out2 = list(res50=res50_rr,res100=res100_rr,res200=res200_rr,res500=res500_rr,res800=res800_rr)
save(out2,file="sim6res_unequal_cov")


}




