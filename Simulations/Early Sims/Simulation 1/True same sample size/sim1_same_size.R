
library(mclust)
library(gdata)
library(MASS)
library(phyclust)
library(ggplot2)

col.blue = rgb(30,144,255,max=255)

sim1.ss = function(seed,n,k){

  mu1 = c(-10,0)
  mu2 = c(10,0)
  sig = matrix(c(36,0,0,36),nrow=2)
 
  set.seed(seed); dat1 = mvrnorm(n, mu1, sig)
  set.seed(seed); dat2 = mvrnorm(n, mu2, sig)
  
  # Define an error matrix:
  id = matrix(c(1,0,0,1),nrow=2)
  err = 25*id
  err1 = k*id
  
  set.seed(seed); dat1.err = mvrnorm(n, mu1, (sig+err))
  set.seed(seed); dat2.err = mvrnorm(n, mu2, (sig+err1))
  
  
  # Error-free obs:
  group1.ef = dat1[which(dat1[,1]<=-10),]
  group2.ef = dat2[which(dat2[,1]>=10),]
  
  # Erroneous obs:
  group1.err = dat1.err[which(dat1.err[,1]>-10 & dat1.err[,1]<10),]
  group2.err = dat2.err[which(dat2.err[,1]<10 & dat2.err[,1]>-10),]
  
  # Equalize the number of obs in two groups:
  n1.ef = nrow(group1.ef)
  n1.err = nrow(group1.err)
  n2.ef = nrow(group2.ef)
  n2.err = nrow(group2.err)
  
  ndiff = (n1.ef+n1.err)-(n2.ef+n2.err)
  
  if(ndiff > 0){
    count = 1
    while(count <= ndiff){
      a = mvrnorm(1, mu2, sig)
      if(a[1] >= 10){
        group2.ef = rbind(group2.ef,a)
        count = count + 1
      }
    }
  } else if(ndiff < 0){
    count = 1
    while(count <= abs(ndiff)){
      a = mvrnorm(1, mu1, sig)
      if(a[1] <= -10){
        group1.ef = rbind(group1.ef,a)
        count = count + 1
      }
    }
  }
  
  n1.ef = nrow(group1.ef)
  n2.ef = nrow(group2.ef)
  
   
  # Form two new clusters:
  group1 = rbind(group1.ef, group1.err)
  group2 = rbind(group2.ef, group2.err)
  
  n1 = nrow(group1)
  n2 = nrow(group2)
  
  # Combine two groups:
  dat = rbind(group1, group2)
  s = nrow(dat)


  xlow = min(dat[,1])-0.5
  xup = max(dat[,1])+0.5
  ylow = min(dat[,2])-0.2
  yup = max(dat[,2])+0.2


  plot(group1, xlim=c(xlow,xup), ylim=c(ylow,yup), pch=19, main="Generated clusters",
    xlab="",ylab="",col=col.blue)
  points(group2, pch=0, col="red")
  abline(v=10, lty="dashed")
  abline(v=-10, lty="dashed")
  
  


  #-----------------------------------------------------------------------------#




  # Construct initial classification:
  z.ini = matrix(c(rep(c(1,0),n1),rep(c(0,1),n2)),byrow=T, nrow=s)
  
  # Construct the true error matrices:
  id = matrix(c(1,0,0,1),nrow=2)
  
  errmat = array(0,dim=c(2,2,s))
  for(i in 1:n1.ef){
    errmat[,,i] = 0*id
  }
  for(i in (n1.ef+1):n1){
    errmat[,,i] = err
  }
  for(i in (n1+1):(n1+n2.ef)){
    errmat[,,i] = 0*id
  }
  for(i in (n1+n2.ef+1):s){
    errmat[,,i] = err1
  }
  
  
  #-----------------------------------------------------------------------------#
  
  
  
  ## Use meVVV:
 
  res1 = meVVV(dat, z.ini)
  m1 = res1$z
  
  t1 = nrow(m1)
  G1 = ncol(m1)
  
  z1 = matrix(0,t1,G1)
  for(i in 1:t1){
    rowmax = max(m1[i,])
    for(k in 1:G1){
      z1[i,k] = ifelse(m1[i,k]==rowmax,2,1)
    }  
  }
  
  dat.new = cbind(dat,z1)
  group1.new1 = cbind(group1,z1[1:n1,])
  group2.new1 = cbind(group2,z1[(n1+1):s,])
  
  plot(group1, xlim=c(xlow,xup), ylim=c(ylow,yup), pch=ifelse(group1.new1[,3]==2,19,0)
    ,col=ifelse(group1.new1[,3]==2,col.blue,"red"), xlab="", ylab="",main="meVVV",xaxt="n",yaxt="n")
  points(group2, pch=ifelse(group2.new1[,3]==2,19,0), col=ifelse(group2.new1[,3]==2,col.blue,"red"))
  abline(v=10, lty="dashed")
  abline(v=-10, lty="dashed")
  

  
  #--------------------------------------------------------------------------#
  
  
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
  
  dat.new = cbind(dat,z2)
  n1 = nrow(group1)
  group1.new2 = cbind(group1,z2[1:n1,])
  group2.new2 = cbind(group2,z2[(n1+1):s,])
  
  plot(group1, xlim=c(xlow,xup), ylim=c(ylow,yup), pch=ifelse(group1.new2[,3]==2,19,0)
    ,col=ifelse(group1.new2[,3]==2,col.blue,"red"), xlab="", ylab="",main="MCME",xaxt="n",yaxt="n")
  points(group2, pch=ifelse(group2.new2[,3]==2,19,0), col=ifelse(group2.new2[,3]==2,col.blue,"red"))
  abline(v=10, lty="dashed")
  abline(v=-10, lty="dashed")
  
  
  

  #-------------------------------------------------------------------------------------#


  
  ## Compute Rand index
  z.ini1 = z.ini + 1
  if(k>0){
    z1.me = rbind(z1[(n1.ef+1):n1,],z1[(n1+n2.ef+1):s,])
    z1.nn = rbind(z1[1:n1.ef,],z1[(n1+1):(n1+n2.ef),])
    
    z2.me = rbind(z2[(n1.ef+1):n1,],z2[(n1+n2.ef+1):s,])
    z2.nn = rbind(z2[1:n1.ef,],z2[(n1+1):(n1+n2.ef),])
    
    z.ini.me = rbind(z.ini1[(n1.ef+1):n1,],z.ini1[(n1+n2.ef+1):s,])
    z.ini.nn = rbind(z.ini1[1:n1.ef,],z.ini1[(n1+1):(n1+n2.ef),])
    
    rand1.me = RRand(z.ini.me,z1.me)
    rand1.nn = RRand(z.ini.nn,z1.nn)
  
    rand2.me = RRand(z.ini.me,z2.me)
    rand2.nn = RRand(z.ini.me,z2.me)
  
    rand1 = RRand(z.ini1,z1)
    rand2 = RRand(z.ini1,z2)
  } else if(k==0){
    z1.me = z1[(n1.ef+1):n1,]
    z1.nn = rbind(z1[1:n1.ef,],z1[(n1+1):s,])
    
    z2.me = z2[(n1.ef+1):n1,]
    z2.nn = rbind(z2[1:n1.ef,],z2[(n1+1):s,])
    
    z.ini.me = z.ini1[(n1.ef+1):n1,]
    z.ini.nn = rbind(z.ini1[1:n1.ef,],z.ini1[(n1+1):s,])
    
    rand1.me = RRand(z.ini.me,z1.me)
    rand1.nn = RRand(z.ini.nn,z1.nn)
  
    rand2.me = RRand(z.ini.me,z2.me)
    rand2.nn = RRand(z.ini.me,z2.me)
  
    rand1 = RRand(z.ini1,z1)
    rand2 = RRand(z.ini1,z2)
  }


  #-----------------------------------------------------------------------#



  ## Misclassification rate
  tab1 = table(z.ini1,z1)
  tab2 = table(z.ini1,z2)
  mis1 = 1-sum(diag(tab1))/sum(tab1)
  mis2 = 1-sum(diag(tab2))/sum(tab2)

  evaluation = list(seed=seed,n=n,rand1=rand1$Rand,rand1.me=rand1.me$Rand,rand1.nn=rand1.nn$Rand,
    rand2=rand2$Rand,rand2.me=rand2.me$Rand,rand2.nn=rand2.nn$Rand,mis1=mis1,mis2=mis2)
  
  out = list(evaluation=evaluation,m1=m1,m2=m2,z1=z1,z2=z2,z.ini=z.ini,group1=group1,group2=group2)
  return(out)

}



## Comparison plots:
sim1.complot = function(z.ini,z1,z2,dat1,dat2){
  
  z.ini1 = z.ini + 1 
  n1 = nrow(dat1)
  n2 = nrow(dat2)
  s = n1 + n2
  dat = rbind(dat1,dat2)
  
  xlow = min(dat[,1])-0.5
  xup = max(dat[,1])+0.5
  ylow = min(dat[,2])-0.2
  yup = max(dat[,2])+0.2
  
  group1.new = cbind(dat1,z1[1:n1,1],z2[1:n1,1],z.ini1[1:n1,1])
  group2.new = cbind(dat2,z1[(n1+1):s,1],z2[(n1+1):s,1],z.ini1[(n1+1):s,1])

  plot(dat1,xlim = c(xlow,xup), ylim = c(ylow,yup),col=col.blue,pch=19
    ,xlab="",ylab="",main="Generated Clusters")
  points(dat2,col="red",pch=0)
  abline(v=10, lty="dashed")
  abline(v=-10, lty="dashed")
    
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
  abline(v=10, lty="dashed")
  abline(v=-10, lty="dashed")
    
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
  abline(v=10, lty="dashed")
  abline(v=-10, lty="dashed")
  
}











seedvec = c(91,92,93,94,95)

png("50_ss.png",900,1500)
par(mfrow=c(5,3))
res50.ss = sapply(seedvec,sim1.ss,n=50,k=25)
par(mfrow=c(1,1))
dev.off()

png("50_ss_com.png",900,1500)
par(mfrow=c(5,3))
for(i in 1:5){
  z.ini = res50.ss[,i]$z.ini
  z1 = res50.ss[,i]$z1
  z2 = res50.ss[,i]$z2
  group1 = res50.ss[,i]$group1
  group2 = res50.ss[,i]$group2
  cp = sim1.complot(z.ini,z1,z2,group1,group2)
}
par(mfrow=c(1,1))
dev.off()



png("100_ss.png",900,1500)
par(mfrow=c(5,3))
res100.ss = sapply(seedvec,sim1.ss,n=100,k=25)
par(mfrow=c(1,1))
dev.off()

png("100_ss_com.png",900,1500)
par(mfrow=c(5,3))
for(i in 1:5){
  z.ini = res100.ss[,i]$z.ini
  z1 = res100.ss[,i]$z1
  z2 = res100.ss[,i]$z2
  group1 = res100.ss[,i]$group1
  group2 = res100.ss[,i]$group2
  cp = sim1.complot(z.ini,z1,z2,group1,group2)
}
par(mfrow=c(1,1))
dev.off()



png("200_ss.png",900,1500)
par(mfrow=c(5,3))
res200.ss = sapply(seedvec,sim1.ss,n=200,k=25)
par(mfrow=c(1,1))
dev.off()

png("200_ss_com.png",900,1500)
par(mfrow=c(5,3))
for(i in 1:5){
  z.ini = res200.ss[,i]$z.ini
  z1 = res200.ss[,i]$z1
  z2 = res200.ss[,i]$z2
  group1 = res200.ss[,i]$group1
  group2 = res200.ss[,i]$group2
  cp = sim1.complot(z.ini,z1,z2,group1,group2)
}
par(mfrow=c(1,1))
dev.off()



png("500_ss.png",900,1500)
par(mfrow=c(5,3))
res500.ss = sapply(seedvec,sim1.ss,n=500,k=25)
par(mfrow=c(1,1))
dev.off()

png("500_ss_com.png",900,1500)
par(mfrow=c(5,3))
for(i in 1:5){
  z.ini = res500[,i]$z.ini
  z1 = res500.ss[,i]$z1
  z2 = res500.ss[,i]$z2
  group1 = res500.ss[,i]$group1
  group2 = res500.ss[,i]$group2
  cp = sim1.complot(z.ini,z1,z2,group1,group2)
}
par(mfrow=c(1,1))
dev.off()



png("800_ss.png",900,1500)
par(mfrow=c(5,3))
res800.ss = sapply(seedvec,sim1.ss,n=800,k=25)
par(mfrow=c(1,1))
dev.off()
  
png("800_ss_com.png",900,1500)
par(mfrow=c(5,3))
for(i in 1:5){
  z.ini = res800.ss[,i]$z.ini
  z1 = res800.ss[,i]$z1
  z2 = res800.ss[,i]$z2
  group1 = res800.ss[,i]$group1
  group2 = res800.ss[,i]$group2
  cp = sim1.complot(z.ini,z1,z2,group1,group2)
}
par(mfrow=c(1,1))
dev.off()


out = list(res50=res50.ss,res100=res100.ss,res200=res200.ss,res500=res500.ss,res800=res800.ss)
save(out,file="sim1_ss_res")


