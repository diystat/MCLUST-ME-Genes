

library(mclust)
library(gdata)
library(MASS)
library(phyclust)
col.blue = rgb(30,144,255,max=255)


sim.new = function(p,N){
  
  set.seed(91)
  
  tau = 0.5
  mu1 = c(0,0)
  mu2 = c(0,0)
  sig1 = matrix(c(64,0,0,64),nrow=2)
  sig2 = matrix(c(16,0,0,16),nrow=2)

  ## Randomly assign constant measurement error to half of observations:
  index = rbinom(N,1,p)
  errmat = array(0,dim=c(2,2,N))
  E = matrix(c(16,0,0,16),nrow=2)
    for(i in 1:N){
      errmat[,,i] = index[i]*E
    }

  ## Generate sample from mixture distribution:
  U = runif(N)
  rand.samples = matrix(0,N,2)
  z.ini = matrix(0,N,2)
  for(i in 1:N){
    if(U[i]<tau){
      rand.samples[i,] = mvrnorm(1,mu1,(sig1+errmat[,,i]))
      z.ini[i,] = c(1,0)
    } else{
      rand.samples[i,] = mvrnorm(1,mu2,(sig2+errmat[,,i]))
      z.ini[i,] = c(0,1)
    }
  }

  xlow = min(rand.samples[,1])-0.5
  xup = max(rand.samples[,1])+0.5
  ylow = min(rand.samples[,2])-0.2
  yup = max(rand.samples[,2])+0.2

  rs.index = cbind(rand.samples,index)
  plot(rand.samples[,1],rand.samples[,2],xlim=c(xlow,xup),ylim=c(ylow,yup),xlab="",ylab="",
    pch=ifelse(index==1,2,1),cex=1.5,main="Generated Sample")



  rs.id = cbind(rand.samples,z.ini,index)  
  plot(rand.samples, xlim=c(xlow,xup), ylim=c(ylow,yup), 
    pch=ifelse((rs.id[,3]==1 & rs.id[,5]==1),17,
      ifelse((rs.id[,3]==1 & rs.id[,5]==0),19,
        ifelse((rs.id[,3]==0 & rs.id[,5]==1),2,0)))
    ,col=ifelse(rs.id[,3]==1,col.blue,"red"), xlab="", ylab="",main="True Grouping",cex=1.5)
  
  
  ## Run MCME:
  my.result = mcmeVVV(rand.samples, z.ini, errmat)
  m = my.result$z
  
  t = nrow(m)
  G = ncol(m)
  
  z = matrix(0,t,G)
  for(i in 1:t){
    rowmax = max(m[i,])
    for(k in 1:G){
      z[i,k] = ifelse(m[i,k]==rowmax,1,0)
    }  
  }

  rs.new = cbind(rand.samples,z,index)  
  #par(mfrow=c(1,2))
  plot(rand.samples, xlim=c(xlow,xup), ylim=c(ylow,yup), 
    pch=ifelse((rs.new[,3]==1 & rs.new[,5]==1),17,
      ifelse((rs.new[,3]==1 & rs.new[,5]==0),19,
        ifelse((rs.new[,3]==0 & rs.new[,5]==1),2,0)))
    ,col=ifelse(rs.new[,3]==1,col.blue,"red"), xlab="", ylab="",main="MCME",cex=1.5)
  
  

  ## Run mevvv:
  res.me = meVVV(rand.samples,z.ini)
  mm = res.me$z

  tt = nrow(mm)
  GG = ncol(mm)
  
  zz = matrix(0,tt,GG)
  for(i in 1:tt){
    rowmax = max(mm[i,])
    for(k in 1:GG){
      zz[i,k] = ifelse(mm[i,k]==rowmax,1,0)
    }  
  }
    
  rs.new1 = cbind(rand.samples,zz,index)
  
  plot(rand.samples, xlim=c(xlow,xup), ylim=c(ylow,yup), 
    pch=ifelse((rs.new1[,3]==1 & rs.new1[,5]==1),17,
      ifelse((rs.new1[,3]==1 & rs.new1[,5]==0),19,
        ifelse((rs.new1[,3]==0 & rs.new1[,5]==1),2,0)))
    ,col=ifelse(rs.new1[,3]==1,col.blue,"red"), xlab="", ylab="",main="meVVV",cex=1.5)
  #par(mfrow=c(1,1))


  ## Calculate Rand index:
  rand.mcme = round(RRand(z.ini+1,z+1)$Rand,4)
  rand.mcme.me = round(RRand(z.ini[which(index==1),]+1,z[which(index==1),]+1)$Rand,4)
  rand.mcme.nn = round(RRand(z.ini[which(index==0),]+1,z[which(index==0),]+1)$Rand,4)
  rand.mevvv = round(RRand(z.ini+1,zz+1)$Rand,4)
  rand.mevvv.me = round(RRand(z.ini[which(index==1),]+1,zz[which(index==1),]+1)$Rand,4)
  rand.mevvv.nn = round(RRand(z.ini[which(index==0),]+1,zz[which(index==0),]+1)$Rand,4)
  
  out = list(rmcme=rand.mcme,rmcme.me=rand.mcme.me,rmcme.nn=rand.mcme.nn,
    rmevvv=rand.mevvv,rmevvv.me=rand.mevvv.me,rmevvv.nn=rand.mevvv.nn)
  return(out)
}






p = c(0.1,0.3,0.5,0.7,0.9)

png("50.png",1200,1500)
par(mfrow=c(5,4))
res50 = sapply(p,sim.new,N=50)
par(mfrow=c(1,1))
dev.off()

png("100.png",1200,1500)
par(mfrow=c(5,4))
res100 = sapply(p,sim.new,N=100)
par(mfrow=c(1,1))
dev.off()

png("200.png",1200,1500)
par(mfrow=c(5,4))
res200 = sapply(p,sim.new,N=200)
par(mfrow=c(1,1))
dev.off()

png("300.png",1200,1500)
par(mfrow=c(5,4))
res300 = sapply(p,sim.new,N=300)
par(mfrow=c(1,1))
dev.off()


out = list(res100=res100,res200=res200,res300=res300)
save(out,file="simnewres")






