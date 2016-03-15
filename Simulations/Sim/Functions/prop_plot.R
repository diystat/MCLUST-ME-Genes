
prop.plot = function(tau,seed=4){
  N = 200
  mu1 = c(0,0)
  mu2 = c(8,0)
  sig1 = matrix(c(64,0,0,64),nrow=2)
  sig2 = matrix(c(16,0,0,16),nrow=2)
  k = 36
  p = 0.3

  data.par = sim.par(N,tau,mu1,mu2,sig1,sig2,k,p)
  set.seed(seed)
  simdata = sim.data(data.par)
  rand.samples = simdata$data
  z.ini = simdata$z.ini
  index = simdata$index

  col.blue = rgb(30,144,255,max=255)  
  xlow = min(rand.samples[,1])-0.5
  xup = max(rand.samples[,1])+0.5
  ylow = min(rand.samples[,2])-0.2
  yup = max(rand.samples[,2])+0.2
  
  rs.new = cbind(rand.samples,index,z.ini[,1])  
  plot(rand.samples,pch=ifelse(rs.new[,3]==1,2,16),main=paste("tau=",tau,sep=""),
    col=ifelse(rs.new[,4]==1,col.blue,"red"),xlab="",ylab="",xlim=c(xlow,xup),ylim=c(ylow,yup))
}



par(mfrow=c(2,3))
prop = c(0.1,0.3,0.5,0.7,0.9)
lapply(prop,prop.plot)
par(mfrow=c(1,1))



















