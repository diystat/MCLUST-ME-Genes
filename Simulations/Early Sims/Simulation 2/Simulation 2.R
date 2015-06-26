

          #---------------------Cluster Simulation 2--------------------------#

          ### One cluster has measurement error, the other does not.


library(mclust)
library(gdata)
library(MASS)
library(phyclust)

col.blue = rgb(30,144,255,max=255)


sim2 = function(seed,n,lambda){

## Simulate two clusters far away from each other.
  mu1 = c(-3,0)
  mu2 = c(3,0)
  sig1 = matrix(c(4,0,0,4),nrow=2)
  sig2 = matrix(c(4,0,0,4),nrow=2)
  
  I = matrix(c(1,0,0,1),nrow=2)
  E = (lambda^2)*I

  # Sample size:
  n1 = n2 = n
  s = n1 + n2

  # Generate samples:
  set.seed(seed); dat1 = mvrnorm(n1, mu1, (sig1+E))
  set.seed(seed); dat2 = mvrnorm(n2, mu2, sig2)
  
  # Combine two groups:
  dat = rbind(dat1, dat2)
  
  xlow = min(dat[,1])-0.5
  xup = max(dat[,1])+0.5
  ylow = min(dat[,2])-0.2
  yup = max(dat[,2])+0.2
  
  #par(mfrow=c(2,2))
  
  #plot(dat1,xlim = c(xlow,xup), ylim = c(ylow,yup),col=col.blue,pch=19
  #  ,xlab="",ylab="",main="Generated Clusters")
  #points(dat2,col="red",pch=0)
  #abline(v=0,lty="dashed")
  # Observe that the two clusters overlap.
  
  

  #-------------------------------------------------------------#




  # Construct initial classification. Use true membership:
  z.ini = matrix(c(rep(c(1,0),n1),rep(c(0,1),n2)), byrow=T, nrow=s)

  # Construct the true error matrices:
  errmat = array(0,dim=c(2,2,s))
  
  for(i in 1:n1){
    errmat[,,i] = E
  }


  
  #-------------------------------------------------------------#



  # Use meVVV
  res.me = meVVV(dat,z.ini)
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
    
 
   
  
    #-------------------------------------------------------------#

  
  
  
  ## Use MCME:
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
   
  
  
  
  
    #-----------------------------------------------------------------------#


  
  
  
  ### Plot both clustering results
  ### Highlight points that are grouped differently,
  ### make all other points gray.
  
  z.ini1 = z.ini + 1  
  group1.new = cbind(dat1,z1[1:n1,1],z2[1:n1,1],z.ini1[1:n1,1])
  group2.new = cbind(dat2,z1[(n1+1):s,1],z2[(n1+1):s,1],z.ini1[(n1+1):s,1])

  
  
  #par(mfrow=c(3,1))
  
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
  
   #par(mfrow=c(1,1))

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
  tab1 = table(z.ini1,z1)
  tab2 = table(z.ini1,z2)
  mis1 = 1-sum(diag(tab1))/sum(tab1)
  mis2 = 1-sum(diag(tab2))/sum(tab2)

  
  out = list(seed=seed,n=n,rand1=rand1$Rand,rand1.me=rand1.me$Rand,rand1.nn=rand1.nn$Rand,
    rand2=rand2$Rand,rand2.me=rand2.me$Rand,rand2.nn=rand2.nn$Rand,mis1=mis1,mis2=mis2,
    m1=m1,m2=m2)
  return(out)
  
  
  
}


  #-----------------------------------------------------------------------#



rrun = function(){
  
seedvec = c(91,92,93,94,95)



## lambda = 0.1
png("50-1.png",900,1500)
par(mfrow=c(5,3))
res50 = sapply(seedvec,sim2,n=50,lambda=0.1)
par(mfrow=c(1,1))
dev.off()


png("100-1.png",900,1500)
par(mfrow=c(5,3))
res100 = sapply(seedvec,sim2,n=100,lambda=0.1)
par(mfrow=c(1,1))
dev.off()


png("200-1.png",900,1500)
par(mfrow=c(5,3))
res200 = sapply(seedvec,sim2,n=200,lambda=0.1)
par(mfrow=c(1,1))
dev.off()


png("500-1.png",900,1500)
par(mfrow=c(5,3))
res500 = sapply(seedvec,sim2,n=500,lambda=0.1)
par(mfrow=c(1,1))
dev.off()


png("800-1.png",900,1500)
par(mfrow=c(5,3))
res800 = sapply(seedvec,sim2,n=800,lambda=0.1)
par(mfrow=c(1,1))
dev.off()



out = list(res50=res50,res100=res100,res200=res200,res500=res500,res800=res800)
save(out,file="sim2res-0.1")







## lambda = 1
## large overlap between clusters

png("50-11.png",900,1500)
par(mfrow=c(5,3))
res50 = sapply(seedvec,sim2,n=50,lambda=1)
par(mfrow=c(1,1))
dev.off()


png("100-11.png",900,1500)
par(mfrow=c(5,3))
res100 = sapply(seedvec,sim2,n=100,lambda=1)
par(mfrow=c(1,1))
dev.off()


png("200-11.png",900,1500)
par(mfrow=c(5,3))
res200 = sapply(seedvec,sim2,n=200,lambda=1)
par(mfrow=c(1,1))
dev.off()


png("500-11.png",900,1500)
par(mfrow=c(5,3))
res500 = sapply(seedvec,sim2,n=500,lambda=1)
par(mfrow=c(1,1))
dev.off()


png("800-11.png",900,1500)
par(mfrow=c(5,3))
res800 = sapply(seedvec,sim2,n=800,lambda=1)
par(mfrow=c(1,1))
dev.off()



out = list(res50=res50,res100=res100,res200=res200,res500=res500,res800=res800)
save(out,file="sim2res-1")
  
  
  
  
  
  
  
## lambda = 2
## large overlap between clusters

png("50-2.png",900,1500)
par(mfrow=c(5,3))
res50 = sapply(seedvec,sim2,n=50,lambda=2)
par(mfrow=c(1,1))
dev.off()


png("100-2.png",900,1500)
par(mfrow=c(5,3))
res100 = sapply(seedvec,sim2,n=100,lambda=2)
par(mfrow=c(1,1))
dev.off()


png("200-2.png",900,1500)
par(mfrow=c(5,3))
res200 = sapply(seedvec,sim2,n=200,lambda=2)
par(mfrow=c(1,1))
dev.off()


png("500-2.png",900,1500)
par(mfrow=c(5,3))
res500 = sapply(seedvec,sim2,n=500,lambda=2)
par(mfrow=c(1,1))
dev.off()


png("800-2.png",900,1500)
par(mfrow=c(5,3))
res800 = sapply(seedvec,sim2,n=800,lambda=2)
par(mfrow=c(1,1))
dev.off()



out = list(res50=res50,res100=res100,res200=res200,res500=res500,res800=res800)
save(out,file="sim2res-2")


  
  
  
## lambda = 3
## large overlap between clusters

png("50-3.png",900,1500)
par(mfrow=c(5,3))
res50 = sapply(seedvec,sim2,n=50,lambda=3)
par(mfrow=c(1,1))
dev.off()


png("100-3.png",900,1500)
par(mfrow=c(5,3))
res100 = sapply(seedvec,sim2,n=100,lambda=3)
par(mfrow=c(1,1))
dev.off()


png("200-3.png",900,1500)
par(mfrow=c(5,3))
res200 = sapply(seedvec,sim2,n=200,lambda=3)
par(mfrow=c(1,1))
dev.off()


png("500-3.png",900,1500)
par(mfrow=c(5,3))
res500 = sapply(seedvec,sim2,n=500,lambda=3)
par(mfrow=c(1,1))
dev.off()


png("800-3.png",900,1500)
par(mfrow=c(5,3))
res800 = sapply(seedvec,sim2,n=800,lambda=3)
par(mfrow=c(1,1))
dev.off()



out = list(res50=res50,res100=res100,res200=res200,res500=res500,res800=res800)
save(out,file="sim2res-3")


}
