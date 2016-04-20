
plot.boundary = function(param,seed){
  ## obtain sim results:
  set.seed(seed)
  sim.result = sim.run(param)
  res.mcme = sim.result$res.mcme
  z.ini = sim.result$z.ini
  rand.samples = sim.result$rand.samples
  
  ## specify graphing parameters:
  col.blue = rgb(30,144,255,max=255)  
  xlow = min(rand.samples[,1])-0.5
  xup = max(rand.samples[,1])+0.5
  ylow = min(rand.samples[,2])-0.2
  yup = max(rand.samples[,2])+0.2


  par(mfrow=c(1,2))
  ## draw two boundaries:
  x = seq(xlow, xup, len=250)
  y = seq(ylow, yup, len=250)
  z <- outer(x, y,
  FUN=function(x, y, ...){
  apply(cbind(x,y), 1, get.boundary, ...)
  }, sim.result=sim.result, k=0)
  # boundary for error-free obs:
  contour(x,y,z,level=0,drawlabels=FALSE)
  
  zz <- outer(x, y,
  FUN=function(x, y, ...){
  apply(cbind(x,y), 1, get.boundary, ...)
  }, sim.result=sim.result, k=9)
  # boundary for erroneous obs:
  contour(x,y,zz,level=0,add=TRUE,drawlabels=FALSE,lty="dashed")
  
  ze <- outer(x, y,
  FUN=function(x, y, ...){
  apply(cbind(x,y), 1, get.boundary.mevvv, ...)
  }, sim.result=sim.result)
  # boundary for erroneous obs:
  contour(x,y,ze,level=0,add=TRUE,drawlabels=FALSE,lty="dotted")


  ## add clustering result to contour plot
  z.mcme = res.mcme$z
  t = nrow(z.mcme)
  G = ncol(z.mcme)

  # Discretize membership matrix:
  z1 = (z.mcme[,1]>z.mcme[,2])

  # distinguish error-free and erroneous obs:
  index = sim.result$index

  rs.new = cbind(rand.samples,z1,index)  
  points(rand.samples,pch=ifelse(rs.new[,4]==1,2,16),
    col=ifelse(rs.new[,3]==1,col.blue,"red"))
  
  
  z.mevvv = sim.result$res.mevvv$z
  z2 = (z.mevvv[,1]>z.mevvv[,2])
  rs.new1 = cbind(rand.samples,z2,index)  
  points(rand.samples,pch=ifelse(rs.new[,4]==1,2,16),
    col=ifelse(rs.new[,3]==1,col.blue,"red"))
  contour(x,y,ze,level=0,add=TRUE,drawlabels=FALSE,lty="dotted")

  par(mfrow=c(1,1))
  
} 



