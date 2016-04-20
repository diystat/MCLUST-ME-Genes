

  ## Plot stepwise clustering
  ## Current function only works on two target clusters.

plot.step = function(dat1, dat2, mu1, mu2, z.ini, errmat, itmax, filename){
  # Arguments: mevvv---clustering result object from stepwiseME()
  #            mcme---clustering result object from mcmeVVV()
  #            dat1,dat2---two groups of observations
  #            mu1,mu2---true means of each group
  #            z.ini---initial membership for all observations
  #            itmax---number of iterations
  #            filename---name of gif file produced
  
  col.blue = rgb(30,144,255,max=255)
  
  data = rbind(dat1,dat2)
  
  xlow = min(data[,1])-0.5
  xup = max(data[,1])+0.5
  ylow = min(data[,2])-0.2
  yup = max(data[,2])+0.2
  
  n1 = nrow(dat1)
  n2 = nrow(dat2)
  
  mevvv = stepwiseME(data,z.ini,itmax=itmax)
  mcme = mcmeVVV(data,z.ini,errmat,itmax=itmax)
  
  s = nrow(z.ini)
  G = ncol(z.ini)
  
  mb.mevvv = mevvv$member
  mb.mcme = mcme$member
  
  for(i in 1:nrow(mb.mevvv)){
    rowmax = max(mb.mevvv[i,])
    for(k in 1:G){
      mb.mevvv[i,k] = ifelse(mb.mevvv[i,k]==rowmax,1,0)
    }  
  }
  
  for(i in 1:nrow(mb.mcme)){
    rowmax = max(mb.mcme[i,])
    for(k in 1:G){
      mb.mcme[i,k] = ifelse(mb.mcme[i,k]==rowmax,1,0)
    }  
  }

  it = itmax
  mbary.mevvv = mbary.mcme = array(0,dim=c(s,G,(it+1)))
  for(i in 1:(it+1)){
    mbary.mevvv[,,i] = mb.mevvv[(i*s-s+1):(i*s),]
    mbary.mcme[,,i] = mb.mcme[(i*s-s+1):(i*s),]
  }


  # Estimated centers at each iteration
  cen.mevvv = mevvv$center
  cen.mcme = mcme$center
  p = ncol(data)
  cenary.mevvv = cenary.mcme = array(0,dim=c(p,G,(it+1)))
  for(i in 2:(it+1)){
    cenary.mevvv[,,i] = cen.mevvv[(i*p-p+1):(i*p),]
    cenary.mcme[,,i] = cen.mcme[(i*p-p+1):(i*p),]
  }
  
  cenary.mevvv[,,1] = cenary.mcme[,,1] = cbind(mu1,mu2)


  
  # Make an animation:
  library("animation")
  oopt <- ani.options(interval = 3) # set time between frames
  saveGIF({
  for (i in 1:(it+1)) {
        
    group1 = cbind(dat1,z.ini[1:n1,])
    group2 = cbind(dat2,z.ini[(n1+1):s,])
    
    group1.new = cbind(dat1,mbary.mevvv[1:n1,1,i],mbary.mcme[1:n1,1,i])
    group2.new = cbind(dat2,mbary.mevvv[(n1+1):s,1,i],mbary.mcme[(n1+1):s,1,i])
    
    #group1.mcme = cbind(dat1,mbary.mcme[1:n1,,i])
    #group2.mcme = cbind(dat2,mbary.mcme[(n1+1):s,,i])
    
    par(mfrow=c(1,3),cex=1.5)
    
    plot(dat1, xlim = c(xlow,xup), ylim = c(ylow,yup), pch=ifelse(group1[,3]==1,19,0),
      col=ifelse(group1[,3]==1,col.blue,"red") ,xlab="", ylab="",
      main="True Clusters")
    points(dat2, pch=ifelse(group2[,3]==1,19,0),
      col=ifelse(group2[,3]==1,col.blue,"red"))
    
    plot(dat1, xlim = c(xlow,xup), ylim = c(ylow,yup), pch=ifelse(group1.new[,3]==1,19,0),
      col=ifelse(group1.new[,3]==group1.new[,4],"gray",ifelse(group1.new[,3]==1,col.blue,"red")),
      xlab="", ylab="", main=paste("meVVV, iteration =",i-1))
    points(dat2, pch=ifelse(group2.new[,3]==1,19,0),
      col=ifelse(group2.new[,3]==group2.new[,4],"gray",ifelse(group2.new[,3]==1,col.blue,"red")))
    points(cenary.mevvv[1,,i],cenary.mevvv[2,,i],col=c("lawngreen","slateblue"),pch=4,cex=3,lwd=5)
    
    plot(dat1, xlim = c(xlow,xup), ylim = c(ylow,yup), pch=ifelse(group1.new[,4]==1,19,0),
      col=ifelse(group1.new[,3]==group1.new[,4],"gray",ifelse(group1.new[,4]==1,col.blue,"red")),
      xlab="", ylab="", main=paste("MCME, iteration =",i-1))
    points(dat2, pch=ifelse(group2.new[,4]==1,19,0),
      col=ifelse(group2.new[,3]==group2.new[,4],"gray",ifelse(group2.new[,4]==1,col.blue,"red")))
    points(cenary.mcme[1,,i],cenary.mcme[2,,i],col=c("lawngreen","slateblue"),pch=4,cex=3,lwd=5)
       
    par(mfrow=c(1,1),cex=1)
    
    ani.pause()
  }}, movie.name=paste(filename,".gif"), ani.width=1800, ani.height=600)
}
  