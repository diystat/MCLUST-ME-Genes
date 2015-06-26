

uncstep = function(data,my.result,filename){
  n = nrow(data)
  palette = c("magenta",heat.colors(10)[2:8],heat.colors(10)[10],"floralwhite")
    
  xlow = min(data[,1])
  xup = max(data[,1])
  ylow = min(data[,2])
  yup = max(data[,2])
  
  xl = seq(0,1,0.1)[1:10]
  yb = rep(0,10)
  xr = seq(0,1,0.1)[2:11]
  yt = rep(0.005,10)
  
  it = my.result$iteration
  
  library("animation")
  oopt <- ani.options(interval = 0.2) # set time between frames
  saveGIF({
  for (i in 1:it) {
    
    z = my.result$member[(n*i+1):(n*i+n),]
    unc = unc(z)
    obs.col = character()
    for(i in 1:n){
      for(k in 2:11){
        if(unc[i]>=(1.1-0.1*k) & unc[i]<(1.2-0.1*k)) obs.col[i] = palette[(k-1)]
      }
    }
    
    par(fig=c(0,1,0,0.8))
    plot(c(xlow,xup),c(ylow,yup),type="n",xlab="",ylab="")
    rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "black")
    points(data,col=obs.col,pch=16,cex=1.5)
  
    par(fig=c(0,1,0.55,1),new=T)
    plot(c(0,1),c(0,0.005),type="n",xlab="",ylab="",
      bty="n",axes=F)
    rect(xl,yb,xr,yt,col=rev(palette),border=NA)
    axis(3,at=seq(1,0,-0.1))

    
    ani.pause()
  }}, movie.name=paste(filename,".gif"), ani.width=1200, ani.height=600)
}  

  
  
  
  
  
  
  