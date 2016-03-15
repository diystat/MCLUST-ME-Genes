
## Obtain coordinates for two boundaries
get.boundary = function(y,simrun,k){
  E = k*matrix(c(1,0,0,1),nrow=2)
  param = simrun$res.mcme$param
  errmat = simrun$errmat
  tau1 = param$tauhat[1]
  tau2 = param$tauhat[2]
  sig1 = param$sigmahat[,,1]
  sig2 = param$sigmahat[,,2]
  mu1 = as.numeric(param$muhat[,1])
  mu2 = as.numeric(param$muhat[,2])
  
  part1 = tau1/sqrt(det(sig1+E))*exp((-1/2)*t(c(y[1],y[2])-mu1)%*%solve(sig1+E)%*%(c(y[1],y[2])-mu1))
  part2 = tau2/sqrt(det(sig2+E))*exp((-1/2)*t(c(y[1],y[2])-mu2)%*%solve(sig2+E)%*%(c(y[1],y[2])-mu2))
  f = part1 - part2
  return(f)
}


## Obtain boundary for meVVV:
get.boundary.mevvv = function(y,simrun){
  param = simrun$res.mevvv$param
  tau1 = param$pro[1]
  tau2 = param$pro[2]
  sig1 = param$var$sigma[,,1]
  sig2 = param$var$sigma[,,2]
  mu1 = as.numeric(param$mean[,1])
  mu2 = as.numeric(param$mean[,2])
  
  part1 = tau1/sqrt(det(sig1))*exp((-1/2)*t(c(y[1],y[2])-mu1)%*%solve(sig1)%*%(c(y[1],y[2])-mu1))
  part2 = tau2/sqrt(det(sig2))*exp((-1/2)*t(c(y[1],y[2])-mu2)%*%solve(sig2)%*%(c(y[1],y[2])-mu2))
  f = part1 - part2
  return(f)
}


## Plot boundaries for MCME and meVVV
plot.boundary = function(simrun,seed){
  ## obtain sim results:
  res.mcme = simrun$res.mcme
  z.ini = simrun$z.ini
  rand.samples = simrun$rand.samples
  
  ## specify graphing parameters:
  col.blue = rgb(30,144,255,max=255)  
  xlow = min(rand.samples[,1])-0.5
  xup = max(rand.samples[,1])+0.5
  ylow = min(rand.samples[,2])-0.2
  yup = max(rand.samples[,2])+0.2

  ## evaluate boundary function over grid:
  x = seq(xlow, xup, len=250)
  y = seq(ylow, yup, len=250)
  z <- outer(x, y,
  FUN=function(x, y, ...){
  apply(cbind(x,y), 1, get.boundary, ...)
  }, simrun=simrun, k=0)
  
  kk = simrun$k
  zz <- outer(x, y,
  FUN=function(x, y, ...){
  apply(cbind(x,y), 1, get.boundary, ...)
  }, simrun=simrun, k=kk)
  
  ze <- outer(x, y,
  FUN=function(x, y, ...){
  apply(cbind(x,y), 1, get.boundary.mevvv, ...)
  }, simrun=simrun)
  

  ## add clustering result to contour plot
  z.mcme = res.mcme$z
  t = nrow(z.mcme)
  G = ncol(z.mcme)

  # Discretize membership matrix:
  z1 = (z.mcme[,1]>z.mcme[,2])

  # distinguish error-free and erroneous obs:
  index = simrun$index

  rs.new = cbind(rand.samples,z1,index,z.ini[,1])  

  #par(mfrow=c(1,2))
  #plot(rand.samples,pch=ifelse(rs.new[,4]==1,2,16),main=paste("True Grouping,seed=",seed,sep=""),
  #  col=ifelse(rs.new[,5]==1,col.blue,"red"),xlab="",ylab="",xlim=c(xlow,xup),ylim=c(ylow,yup),cex=1.6)
  
  #plot(rand.samples,pch=ifelse(rs.new[,4]==1,2,16),main=paste("MCLUST-ME,seed=",seed,sep=""),
  #  col=ifelse(rs.new[,3]==1,col.blue,"red"),xlab="",ylab="",xlim=c(xlow,xup),ylim=c(ylow,yup),cex=1.6)
  #plot(rand.samples,pch=ifelse(rs.new[,4]==1,2,16),main="MCLUST-ME",cex.main=2,
  #  col=ifelse(rs.new[,3]==1,col.blue,"red"),xlab="",ylab="",xlim=c(xlow,xup),ylim=c(ylow,yup),cex=1.7)
  plot(rand.samples,main="MCLUST-ME",cex.main=1,pch=ifelse(rs.new[,3]==1,22,24),
    bg=ifelse(rs.new[,4]==1,"black","white"), xlab="",ylab="",xlim=c(xlow,xup),ylim=c(ylow,yup),cex=0.7)
  
  # boundary for error-free obs:
  contour(x,y,z,level=0,add=TRUE,drawlabels=FALSE,lty="dashed",lwd=2)
  # boundary for erroneous obs:
  contour(x,y,zz,level=0,add=TRUE,drawlabels=FALSE,lty="dotted",lwd=2)
  # boundary for mevvv:
  contour(x,y,ze,level=0,add=TRUE,drawlabels=FALSE,lwd=2)
  
  
  z.mevvv = simrun$res.mevvv$z
  z2 = (z.mevvv[,1]>z.mevvv[,2])
  rs.new1 = cbind(rand.samples,z2,index)  
  #plot(rand.samples,pch=ifelse(rs.new1[,4]==1,2,16),main=paste("MCLUST,seed=",seed,sep=""),
  #  col=ifelse(rs.new1[,3]==1,col.blue,"red"),xlab="",ylab="",xlim=c(xlow,xup),ylim=c(ylow,yup),cex=1.6)
  plot(rand.samples,pch=ifelse(rs.new1[,3]==1,22,24),main="MCLUST",cex.main=1,
    bg=ifelse(rs.new1[,4]==1,"black","white"),xlab="",ylab="",xlim=c(xlow,xup),ylim=c(ylow,yup),cex=0.7)
  contour(x,y,ze,level=0,add=TRUE,drawlabels=FALSE,lwd=2)
  #par(mfrow=c(1,1))
  
} 


## Plot boundaries for MCME and meVVV with pointsize reflecting clustering uncertainty
plot.boundary.new = function(simrun,pointsize,pointsize2){
  ## obtain sim results:
  res.mcme = simrun$res.mcme
  z.ini = simrun$z.ini
  rand.samples = simrun$rand.samples
  
  ## specify graphing parameters:
  col.blue = rgb(30,144,255,max=255)  
  xlow = min(rand.samples[,1])-0.5
  xup = max(rand.samples[,1])+0.5
  ylow = min(rand.samples[,2])-0.2
  yup = max(rand.samples[,2])+0.2

  ## evaluate boundary function over grid:
  x = seq(xlow, xup, len=250)
  y = seq(ylow, yup, len=250)
  z <- outer(x, y,
  FUN=function(x, y, ...){
  apply(cbind(x,y), 1, get.boundary, ...)
  }, simrun=simrun, k=0)
  
  kk = simrun$k
  zz <- outer(x, y,
  FUN=function(x, y, ...){
  apply(cbind(x,y), 1, get.boundary, ...)
  }, simrun=simrun, k=kk)
  
  ze <- outer(x, y,
  FUN=function(x, y, ...){
  apply(cbind(x,y), 1, get.boundary.mevvv, ...)
  }, simrun=simrun)
  

  ## add clustering result to contour plot
  z.mcme = res.mcme$z
  t = nrow(z.mcme)
  G = ncol(z.mcme)

  # Discretize membership matrix:
  z1 = (z.mcme[,1]>z.mcme[,2])

  # distinguish error-free and erroneous obs:
  index = simrun$index

  rs.new = cbind(rand.samples,z1,index,z.ini[,1])  
 
  #par(mfrow=c(1,2))
  #plot(rand.samples,pch=ifelse(rs.new[,4]==1,2,16),main=paste("True Grouping,seed=",seed,sep=""),
  #  col=ifelse(rs.new[,5]==1,col.blue,"red"),xlab="",ylab="",xlim=c(xlow,xup),ylim=c(ylow,yup),cex=1.6)
  
  #plot(rand.samples,pch=ifelse(rs.new[,4]==1,2,16),main=paste("MCLUST-ME,seed=",seed,sep=""),
  #  col=ifelse(rs.new[,3]==1,col.blue,"red"),xlab="",ylab="",xlim=c(xlow,xup),ylim=c(ylow,yup),cex=1.6)
  #plot(rand.samples,pch=ifelse(rs.new[,4]==1,2,16),main="MCLUST-ME",cex.main=2,
  #  col=ifelse(rs.new[,3]==1,col.blue,"red"),xlab="",ylab="",xlim=c(xlow,xup),ylim=c(ylow,yup),cex=1.7)
  plot(rand.samples,main="MCLUST-ME",cex.main=1,pch=ifelse(rs.new[,3]==1,22,24),
    bg=ifelse(rs.new[,4]==1,"black","white"), xlab="",ylab="",xlim=c(xlow,xup),ylim=c(ylow,yup),cex=pointsize)
  
  # boundary for error-free obs:
  contour(x,y,z,level=0,add=TRUE,drawlabels=FALSE,lty="dashed",lwd=2)
  # boundary for erroneous obs:
  contour(x,y,zz,level=0,add=TRUE,drawlabels=FALSE,lty="dotted",lwd=2)
  # boundary for mevvv:
  contour(x,y,ze,level=0,add=TRUE,drawlabels=FALSE,lwd=2)
  
  
  z.mevvv = simrun$res.mevvv$z
  z2 = (z.mevvv[,1]>z.mevvv[,2])
  rs.new1 = cbind(rand.samples,z2,index)  
  #plot(rand.samples,pch=ifelse(rs.new1[,4]==1,2,16),main=paste("MCLUST,seed=",seed,sep=""),
  #  col=ifelse(rs.new1[,3]==1,col.blue,"red"),xlab="",ylab="",xlim=c(xlow,xup),ylim=c(ylow,yup),cex=1.6)
  plot(rand.samples,pch=ifelse(rs.new1[,3]==1,22,24),main="MCLUST",cex.main=1,
    bg=ifelse(rs.new1[,4]==1,"black","white"),xlab="",ylab="",xlim=c(xlow,xup),ylim=c(ylow,yup),cex=pointsize2)
  contour(x,y,ze,level=0,add=TRUE,drawlabels=FALSE,lwd=2)
  #par(mfrow=c(1,1))
  
} 

