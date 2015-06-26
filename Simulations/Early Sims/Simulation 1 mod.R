
      #----------------------Cluster Simulation 1-mod-------------------------#

      ### Error outside; no error inside


library(mclust)
library(gdata)
library(MASS)
library(phyclust)

col.blue = rgb(30,144,255,max=255)

sim1.mod = function(seed,n){

  set.seed(seed)
  n = n

## Simulate two clusters that are clearly separated but close to each other.
  mu1 = c(-2,0)
  mu2 = c(2,0)
  # 3.2 is chosen as square root of chi(0.01) with df=2, so that 99% of data are expected
  # to fall into the region
  sig = matrix(c(2,0,0,2),nrow=2)
 
  dat1 = mvrnorm(n, mu1, sig)
  dat2 = mvrnorm(n, mu2, sig)

  #plot(dat1, xlim = c(-7,7), ylim = c(-5,5), col="blue")
  #points(dat2, col="red")
  # Observe that the two clusters are separated but close.
  
  
  # Now we truncate each group and keep the observations with x-coord <=-3 and >=3:
  dat1.tr = dat1[which(dat1[,1]>=-2),]
  dat2.tr = dat2[which(dat2[,1]<=2),]
  
  
  # Define an error matrix:
  err = matrix(c(4,0,0,4),nrow=2)
  
  sig.new = sig + err
  dat1.new = mvrnorm(n, mu1, sig.new)
  dat2.new = mvrnorm(n, mu2, sig.new)
  
  dat1.new.tr = dat1.new[which(dat1.new[,1]<(-2)),]
  dat2.new.tr = dat2.new[which(dat2.new[,1]>2),]
  
  group1 = rbind(dat1.tr, dat1.new.tr)
  group2 = rbind(dat2.tr, dat2.new.tr)
  
  n1 = nrow(group1)
  n2 = nrow(group2)
  
  # Combine two groups:
  dat = rbind(group1, group2)
  s = nrow(dat)


  xlow = min(dat[,1])-0.5
  xup = max(dat[,1])+0.5
  ylow = min(dat[,2])-0.2
  yup = max(dat[,2])+0.2


  plot(group1, xlim=c(xlow,xup), ylim=c(ylow,yup), pch=19, main="Generated clusters"
    ,xlab="",ylab="",col=col.blue,xaxt="n",yaxt="n")
  points(group2, pch=0, col="red")
  abline(v=2, lty="dashed")
  abline(v=-2, lty="dashed")
  # The combined dataset has points in the middle are commingled.
  
  


  #-----------------------------------------------------------------------------#




  # Construct initial classification:
  z.ini = matrix(c(rep(c(1,0),n1),rep(c(0,1),n2)),byrow=T, nrow=s)
  
  # Construct the true error matrices:
  temp = array(0,dim=c(2,3,s))
  
  for(i in 1:s){
    temp[,1,i] = dat[i,]
  }
  
  for(i in 1:s){
    if(temp[1,1,i]<(-2) || temp[1,1,i]>2){
      temp[1,2,i] = temp[2,3,i] = 4
    }
  }
  
  errmat = array(0,dim=c(2,2,s))
  for(i in 1:s){
    errmat[,,i] = temp[,2:3,i]
  }
  
  
  #-----------------------------------------------------------------------------#
  
  
  
  ## Use meVVV:
 
  res1 = meVVV(dat, z.ini)
  z1 = res1$z
  
  t = nrow(z1)
  G = ncol(z1)
  
  for(i in 1:t){
    rowmax = max(z1[i,])
    for(k in 1:G){
      z1[i,k] = ifelse(z1[i,k]==rowmax,2,1)
    }  
  }
   

  
  
  #--------------------------------------------------------------------------#
  
  
  
  
  ## Use MCME:
  my.result = mcmeVVV(dat, z.ini, errmat)
  names(my.result)
  z2 = my.result$z
  
  t = nrow(z2)
  G = ncol(z2)
  
  for(i in 1:t){
    rowmax = max(z2[i,])
    for(k in 1:G){
      z2[i,k] = ifelse(z2[i,k]==rowmax,2,1)
    }  
  }
  

  
  
   
      #-----------------------------------------------------------------------#


  
  
  
  ### Plot both clustering results
  ### Highlight points that are grouped differently,
  ### make all other points gray.
  
    
  group1.new = cbind(group1,z1[1:n1,1],z2[1:n1,1])
  group2.new = cbind(group2,z1[(n1+1):s,1],z2[(n1+1):s,1])

  
  
  
  par(mfrow=c(3,1))
  
  plot(group1, xlim=c(xlow,xup), ylim=c(ylow,yup), pch=19, main="Generated clusters"
    ,xlab="",ylab="",col=col.blue,xaxt="n",yaxt="n")
  points(group2, pch=0, col="red")
  abline(v=1, lty="dashed")
  abline(v=-1, lty="dashed")
  
  plot(group1, xlim = c(xlow,xup), ylim = c(ylow,yup), pch=ifelse(group1.new[,3]==2,19,0),
    col=ifelse(group1.new[,3]==group1.new[,4],"gray",ifelse(group1.new[,3]==2,col.blue,"red")),
    xlab="", ylab="",main="meVVV")
  points(group2, pch=ifelse(group2.new[,3]==2,19,0),
    col=ifelse(group2.new[,3]==group2.new[,4],"gray",ifelse(group2.new[,3]==2,col.blue,"red")))
  abline(v=1, lty="dashed")
  abline(v=-1, lty="dashed")
  
  plot(group1, xlim = c(xlow,xup), ylim = c(ylow,yup), pch=ifelse(group1.new[,4]==2,19,0),
    col=ifelse(group1.new[,3]==group1.new[,4],"gray",ifelse(group1.new[,4]==2,col.blue,"red")),
    xlab="", ylab="",main="MCME")
  points(group2, pch=ifelse(group2.new[,4]==2,19,0),
    col=ifelse(group2.new[,3]==group2.new[,4],"gray",ifelse(group2.new[,4]==2,col.blue,"red")))
  abline(v=1, lty="dashed")
  abline(v=-1, lty="dashed")
   
  par(mfrow=c(1,1))


  #-------------------------------------------------------------------------------------#


  
  ## Compute Rand index
  z.ini = z.ini + 1
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

  
    
  out = list(seed=seed,n=n,rand1=rand1$Rand,rand2=rand2$Rand,mis1=mis1,mis2=mis2)
  return(out)
  
}




  #-----------------------------------------------------------------------#




seedvec = c(91,92,93,94,95)

png("50.png",900,1500)
par(mfrow=c(5,3))
res50 = sapply(seedvec,sim1.mod,n=50)
par(mfrow=c(1,1))
dev.off()


png("100.png",900,1500)
par(mfrow=c(5,3))
res100 = sapply(seedvec,sim1.mod,n=100)
par(mfrow=c(1,1))
dev.off()


png("200.png",900,1500)
par(mfrow=c(5,3))
res200 = sapply(seedvec,sim1.mod,n=200)
par(mfrow=c(1,1))
dev.off()


png("500.png",900,1500)
par(mfrow=c(5,3))
res500 = sapply(seedvec,sim1.mod,n=500)
par(mfrow=c(1,1))
dev.off()


png("800.png",900,1500)
par(mfrow=c(5,3))
res800 = sapply(seedvec,sim1.mod,n=800)
par(mfrow=c(1,1))
dev.off()



out = list(res50=res50,res100=res100,res200=res200,res500=res500,res800=res800)
save(out,file="sim1modres")

