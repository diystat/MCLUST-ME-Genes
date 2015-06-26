
          #-----------------Cluster Simulation 2-fix center-------------------#

          ### One cluster has measurement error, the other does not.
          ### At each iteration, fix center estimate at the true value,
          ### and examine the clustering results.


library(mclust)
library(gdata)
library(MASS)
library(phyclust)

col.blue = rgb(30,144,255,max=255)


sim2.fc = function(seed,n,lambda){

  n = n
  set.seed(seed)

## Simulate two clusters far away from each other.
  mu1 = c(-3,0)
  mu2 = c(3,0)
  sig = matrix(c(4,0,0,4),nrow=2)
  
  I = matrix(c(1,0,0,1),nrow=2)
  E = (lambda^2)*I

  # Sample size:
  n1 = n2 = n
  s = n1 + n2

  # Generate samples:
  dat1 = mvrnorm(n1, mu1, (sig+E))
  dat2 = mvrnorm(n2, mu2, sig)
  
  # Combine two groups:
  dat = rbind(dat1, dat2)
  
  xlow = min(dat[,1])-0.5
  xup = max(dat[,1])+0.5
  ylow = min(dat[,2])-0.2
  yup = max(dat[,2])+0.2
  
  
  plot(dat1,xlim = c(xlow,xup), ylim = c(ylow,yup),col=col.blue,pch=19
    ,xlab="",ylab="",main="Generated Clusters")
  points(dat2,col="red",pch=0)
  
  

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
  z1 = res.me$z

  t = nrow(z1)
  G = ncol(z1)
  
  for(i in 1:t){
    rowmax = max(z1[i,])
    for(k in 1:G){
      z1[i,k] = ifelse(z1[i,k]==rowmax,2,1)
    }  
  }
    
  group1.new1 = cbind(dat1,z1[1:n1,])
  group2.new1 = cbind(dat2,z1[(n1+1):s,])
  
  plot(dat1, xlim = c(xlow,xup), ylim = c(ylow,yup), pch=ifelse(group1.new1[,3]==2,19,0),
    col=ifelse(group1.new1[,3]==2,col.blue,"red") ,xlab="", ylab="",main="meVVV")
  points(dat2, pch=ifelse(group2.new1[,3]==2,19,0),
    col=ifelse(group2.new1[,3]==2,col.blue,"red"))
  

  
  

  #-------------------------------------------------------------#



  # meVVV with centers fixed
  res.me.fc = meVVV.fc(dat,z.ini,mu1=mu1,mu2=mu2)
  z2 = res.me.fc$z

  t = nrow(z2)
  G = ncol(z2)
  
  for(i in 1:t){
    rowmax = max(z2[i,])
    for(k in 1:G){
      z2[i,k] = ifelse(z2[i,k]==rowmax,2,1)
    }  
  }
    
  group1.new2 = cbind(dat1,z2[1:n1,])
  group2.new2 = cbind(dat2,z2[(n1+1):s,])
  
  plot(dat1, xlim = c(xlow,xup), ylim = c(ylow,yup), pch=ifelse(group1.new2[,3]==2,19,0),
    col=ifelse(group1.new2[,3]==2,col.blue,"red") ,xlab="", ylab="",main="meVVV-fixed center")
  points(dat2, pch=ifelse(group2.new2[,3]==2,19,0),
    col=ifelse(group2.new2[,3]==2,col.blue,"red"))







    #-------------------------------------------------------------#





  ## Use MCME:
  my.result = mcmeVVV(dat, z.ini, errmat)
  z3 = my.result$z
  
  t = nrow(z3)
  G = ncol(z3)
  
  for(i in 1:t){
    rowmax = max(z3[i,])
    for(k in 1:G){
      z3[i,k] = ifelse(z3[i,k]==rowmax,2,1)
    }  
  }
    
  group1.new3 = cbind(dat1,z3[1:n1,])
  group2.new3 = cbind(dat2,z3[(n1+1):s,])
  
  plot(dat1, xlim = c(xlow,xup), ylim = c(ylow,yup), pch=ifelse(group1.new3[,3]==2,19,0),
    col=ifelse(group1.new3[,3]==2,col.blue,"red") ,xlab="", ylab="",main="MCME")
  points(dat2, pch=ifelse(group2.new3[,3]==2,19,0),
    col=ifelse(group2.new3[,3]==2,col.blue,"red"))
  
  


    #-------------------------------------------------------------#



  
  ## MCME with fixed centers:
    
  my.result.fc = mcmeVVV.fc(dat, z.ini, errmat, mu1=mu1, mu2=mu2)
  z4 = my.result.fc$z
  
  t = nrow(z4)
  G = ncol(z4)
  
  for(i in 1:t){
    rowmax = max(z4[i,])
    for(k in 1:G){
      z4[i,k] = ifelse(z4[i,k]==rowmax,2,1)
    }  
  }
    
  group1.new4 = cbind(dat1,z4[1:n1,])
  group2.new4 = cbind(dat2,z4[(n1+1):s,])
  
  plot(dat1, xlim = c(xlow,xup), ylim = c(ylow,yup), pch=ifelse(group1.new4[,3]==2,19,0),
    col=ifelse(group1.new4[,3]==2,col.blue,"red") ,xlab="", ylab="",main="MCME-fixed center")
  points(dat2, pch=ifelse(group2.new4[,3]==2,19,0),
    col=ifelse(group2.new4[,3]==2,col.blue,"red"))




  #-----------------------------------------------------------------------#


  
  ## Compute Rand index
  z.ini = z.ini + 1
  rand1 = RRand(z.ini,z1)
  rand2 = RRand(z.ini,z2)
  rand3 = RRand(z.ini,z3)
  rand4 = RRand(z.ini,z4)
  z.ini = z.ini - 1
  
  

  #-----------------------------------------------------------------------#



  ## Misclassification rate
  z.ini = z.ini + 1
  tab1 = table(z.ini,z1)
  tab2 = table(z.ini,z2)
  tab3 = table(z.ini,z3)
  tab4 = table(z.ini,z4)
  mis1 = 1-sum(diag(tab1))/sum(tab1)
  mis2 = 1-sum(diag(tab2))/sum(tab2)
  mis3 = 1-sum(diag(tab3))/sum(tab3)
  mis4 = 1-sum(diag(tab4))/sum(tab4)
  z.ini = z.ini - 1

  
  out = list(seed=seed,n=n,rand1=rand1$Rand,rand2=rand2$Rand,rand3=rand3$Rand,
    rand4=rand4$Rand,mis1=mis1,mis2=mis2,mis3=mis3,mis4=mis4)
  return(out)
  
  
  
}




  #-----------------------------------------------------------------------#




seedvec = c(91,92,93,94,95)


## lambda = 1
png("50-2.png",1500,1500)
par(mfrow=c(5,5))
res50 = sapply(seedvec,sim2.fc,n=50,lambda=1)
par(mfrow=c(1,1))
dev.off()


png("100-2.png",1500,1500)
par(mfrow=c(5,5))
res100 = sapply(seedvec,sim2.fc,n=100,lambda=1)
par(mfrow=c(1,1))
dev.off()


png("200-2.png",1500,1500)
par(mfrow=c(5,5))
res200 = sapply(seedvec,sim2.fc,n=200,lambda=1)
par(mfrow=c(1,1))
dev.off()


png("500-2.png",1500,1500)
par(mfrow=c(5,5))
res500 = sapply(seedvec,sim2.fc,n=500,lambda=1)
par(mfrow=c(1,1))
dev.off()


png("800-2.png",1500,1500)
par(mfrow=c(5,5))
res800 = sapply(seedvec,sim2.fc,n=800,lambda=1)
par(mfrow=c(1,1))
dev.off()






## lambda = 2
png("50-3.png",1500,1500)
par(mfrow=c(5,5))
res50_2 = sapply(seedvec,sim2.fc,n=50,lambda=2)
par(mfrow=c(1,1))
dev.off()


png("100-3.png",1500,1500)
par(mfrow=c(5,5))
res100_2 = sapply(seedvec,sim2.fc,n=100,lambda=2)
par(mfrow=c(1,1))
dev.off()


png("200-3.png",1500,1500)
par(mfrow=c(5,5))
res200_2 = sapply(seedvec,sim2.fc,n=200,lambda=2)
par(mfrow=c(1,1))
dev.off()


png("500-3.png",1500,1500)
par(mfrow=c(5,5))
res500_2 = sapply(seedvec,sim2.fc,n=500,lambda=2)
par(mfrow=c(1,1))
dev.off()


png("800-3.png",1500,1500)
par(mfrow=c(5,5))
res800_2 = sapply(seedvec,sim2.fc,n=800,lambda=2)
par(mfrow=c(1,1))
dev.off()



out = list(res50=res50,res100=res100,res200=res200,res500=res500,res800=res800,
  res50_2=res50_2,res100_2=res100_2,res200_2=res200_2,res500_2=res500_2,res800_2=res800_2)
save(out,file="sim2res.fc")







