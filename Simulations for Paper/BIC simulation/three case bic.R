# Three well-separated clusters:
bic_wellsep = function(seed){
  N = 300
  tau1 = 0.3
  tau2 = 0.4
  mu1 = c(-16,0)
  mu2 = c(16,0)
  mu3 = c(0,24)
  sig1 = matrix(c(64,0,0,64),nrow=2)
  sig2 = matrix(c(16,0,0,16),nrow=2)
  sig3 = matrix(c(36,0,0,36),nrow=2)
  k = 6
  p = 0.5

  set.seed(seed)
  data.par = sim.par3(N,tau1,tau2,mu1,mu2,mu3,sig1,sig2,sig3,k,p)
  simdata = sim.data3(data.par)

  dat = simdata$data
  err = simdata$err
  z = simdata$z.ini
  group = numeric(N)
  for(i in 1:N){
    if(z[i,1]==1 & z[i,2]==0){
      group[i] = 1
    }else if(z[i,1]==0 & z[i,2]==1){
      group[i] = 2
    }else{
      group[i] = 3
    }
  }
  
  gp1 = dat[group==1,]
  gp2 = dat[group==2,]
  gp3 = dat[group==3,]
  
  xmin = min(dat[,1])
  xmax = max(dat[,1])
  ymin = min(dat[,2])
  ymax = max(dat[,2])
  
  plot(gp1,xlim=c(xmin,xmax),ylim=c(ymin,ymax),pch=0,
    main="Three well-separated clusters",xlab="",ylab="")
  points(gp2,pch=1)
  points(gp3,pch=2)
  
  out = test.3group(seed,N,tau1,tau2,mu1,mu2,mu3,sig1,sig2,sig3,k,p)
  return(out)
}
  
  
# Two clusters close, far away from third one:
bic_2close = function(seed){  
  N = 300
  tau1 = 0.3
  tau2 = 0.4
  mu1 = c(-6,0)
  mu2 = c(6,0)
  mu3 = c(0,36)
  sig1 = matrix(c(64,0,0,64),nrow=2)
  sig2 = matrix(c(16,0,0,16),nrow=2)
  sig3 = matrix(c(36,0,0,36),nrow=2)
  k = 25
  p = 0.5

  set.seed(seed)
  data.par = sim.par3(N,tau1,tau2,mu1,mu2,mu3,sig1,sig2,sig3,k,p)
  simdata = sim.data3(data.par)

  dat = simdata$data
  err = simdata$err
  z = simdata$z.ini
  group = numeric(N)
  for(i in 1:N){
    if(z[i,1]==1 & z[i,2]==0){
      group[i] = 1
    }else if(z[i,1]==0 & z[i,2]==1){
      group[i] = 2
    }else{
      group[i] = 3
    }
  }
  
  gp1 = dat[group==1,]
  gp2 = dat[group==2,]
  gp3 = dat[group==3,]
  
  xmin = min(dat[,1])
  xmax = max(dat[,1])
  ymin = min(dat[,2])
  ymax = max(dat[,2])
  
  plot(gp1,xlim=c(xmin,xmax),ylim=c(ymin,ymax),pch=0,
    main="Two clusters close, one far away",xlab="",ylab="")
  points(gp2,pch=1)
  points(gp3,pch=2)
  
  out = test.3group(seed,N,tau1,tau2,mu1,mu2,mu3,sig1,sig2,sig3,k,p)
  return(out)
}  
  
  
# All three clusters close to each other:  
bic_3close = function(seed){  
  N = 300
  tau1 = 0.3
  tau2 = 0.4
  mu1 = c(-8,0)
  mu2 = c(9,0)
  mu3 = c(0,20)
  sig1 = matrix(c(64,0,0,64),nrow=2)
  sig2 = matrix(c(16,0,0,16),nrow=2)
  sig3 = matrix(c(36,0,0,36),nrow=2)
  k = 36
  p = 0.5

  set.seed(seed)
  data.par = sim.par3(N,tau1,tau2,mu1,mu2,mu3,sig1,sig2,sig3,k,p)
  simdata = sim.data3(data.par)

  dat = simdata$data
  err = simdata$err
  z = simdata$z.ini
  group = numeric(N)
  for(i in 1:N){
    if(z[i,1]==1 & z[i,2]==0){
      group[i] = 1
    }else if(z[i,1]==0 & z[i,2]==1){
      group[i] = 2
    }else{
      group[i] = 3
    }
  }
  
  gp1 = dat[group==1,]
  gp2 = dat[group==2,]
  gp3 = dat[group==3,]
  
  xmin = min(dat[,1])
  xmax = max(dat[,1])
  ymin = min(dat[,2])
  ymax = max(dat[,2])
  
  plot(gp1,xlim=c(xmin,xmax),ylim=c(ymin,ymax),pch=0,
    main="All clusters close",xlab="",ylab="")
  points(gp2,pch=1)
  points(gp3,pch=2)
  
  out = test.3group(seed,N,tau1,tau2,mu1,mu2,mu3,sig1,sig2,sig3,k,p)
  return(out)
}  
  


plot.wellsep = function(seed,bicws){
  N = 300
  tau1 = 0.3
  tau2 = 0.4
  mu1 = c(-16,0)
  mu2 = c(16,0)
  mu3 = c(0,24)
  sig1 = matrix(c(64,0,0,64),nrow=2)
  sig2 = matrix(c(16,0,0,16),nrow=2)
  sig3 = matrix(c(36,0,0,36),nrow=2)
  k = 6
  p = 0.5

  set.seed(seed)
  data.par = sim.par3(N,tau1,tau2,mu1,mu2,mu3,sig1,sig2,sig3,k,p)
  simdata = sim.data3(data.par)

  dat = simdata$data
  err = simdata$err
  z = simdata$z.ini
  group = numeric(N)
  for(i in 1:N){
    if(z[i,1]==1 & z[i,2]==0){
      group[i] = 1
    }else if(z[i,1]==0 & z[i,2]==1){
      group[i] = 2
    }else{
      group[i] = 3
    }
  }
  
  gp1 = dat[group==1,]
  gp2 = dat[group==2,]
  gp3 = dat[group==3,]
  
  xmin = min(dat[,1])
  xmax = max(dat[,1])
  ymin = min(dat[,2])
  ymax = max(dat[,2])
  
  bic.mcme = bicws$out.mcme
  bic.mclust = bicws$out.mclust
  allbic = c(bic.mcme[,2],bic.mclust[1:7])
  ymin.bic = min(allbic)
  ymax.bic = max(allbic)

  #par(mfrow=c(1,2))
  plot(gp1,xlim=c(xmin,xmax),ylim=c(ymin,ymax),pch=0,
      main="Case 1: Three well-separated clusters",xlab="",ylab="")
  points(gp2,pch=1)
  points(gp3,pch=2)

  plot(bic.mcme,type="b",lwd=2,ylim=c(ymin.bic,ymax.bic),
    xlab="number of components",ylab="BIC",pch=0)
  lines(1:7,bic.mclust[1:7],type="b",lty="dashed",pch=2,lwd=2)
  legend("bottomright",legend=c("MCLUST-ME","MCLUST"),
    pch=c(0,2),lty=c("solid","dashed"),cex=0.5)
  #par(mfrow=c(1,1))
}



plot.2close = function(seed,bic2){
  N = 300
  tau1 = 0.3
  tau2 = 0.4
  mu1 = c(-6,0)
  mu2 = c(6,0)
  mu3 = c(0,36)
  sig1 = matrix(c(64,0,0,64),nrow=2)
  sig2 = matrix(c(16,0,0,16),nrow=2)
  sig3 = matrix(c(36,0,0,36),nrow=2)
  k = 25
  p = 0.5

  set.seed(seed)
  data.par = sim.par3(N,tau1,tau2,mu1,mu2,mu3,sig1,sig2,sig3,k,p)
  simdata = sim.data3(data.par)

  dat = simdata$data
  err = simdata$err
  z = simdata$z.ini
  group = numeric(N)
  for(i in 1:N){
    if(z[i,1]==1 & z[i,2]==0){
      group[i] = 1
    }else if(z[i,1]==0 & z[i,2]==1){
      group[i] = 2
    }else{
      group[i] = 3
    }
  }
  
  gp1 = dat[group==1,]
  gp2 = dat[group==2,]
  gp3 = dat[group==3,]
  
  xmin = min(dat[,1])
  xmax = max(dat[,1])
  ymin = min(dat[,2])
  ymax = max(dat[,2])
  
  bic.mcme = bic2$out.mcme
  bic.mclust = bic2$out.mclust
  allbic = c(bic.mcme[,2],bic.mclust[1:7])
  ymin.bic = min(allbic)
  ymax.bic = max(allbic)

  #par(mfrow=c(1,2))
  plot(gp1,xlim=c(xmin,xmax),ylim=c(ymin,ymax),pch=0,
      main="Case 2: Two clusters close, one far away",xlab="",ylab="")
  points(gp2,pch=1)
  points(gp3,pch=2)

  plot(bic.mcme,type="b",lwd=2,ylim=c(ymin.bic,ymax.bic),
    xlab="number of components",ylab="BIC",pch=0)
  lines(1:7,bic.mclust[1:7],type="b",lty="dashed",pch=2,lwd=2)
  legend("bottomright",legend=c("MCLUST-ME","MCLUST"),
    pch=c(0,2),lty=c("solid","dashed"),cex=0.5)
  #par(mfrow=c(1,1))
}


plot.3close = function(seed,bic3){
  N = 300
  tau1 = 0.3
  tau2 = 0.4
  mu1 = c(-8,0)
  mu2 = c(9,0)
  mu3 = c(0,20)
  sig1 = matrix(c(64,0,0,64),nrow=2)
  sig2 = matrix(c(16,0,0,16),nrow=2)
  sig3 = matrix(c(36,0,0,36),nrow=2)
  k = 25
  p = 0.5

  set.seed(seed)
  data.par = sim.par3(N,tau1,tau2,mu1,mu2,mu3,sig1,sig2,sig3,k,p)
  simdata = sim.data3(data.par)

  dat = simdata$data
  err = simdata$err
  z = simdata$z.ini
  group = numeric(N)
  for(i in 1:N){
    if(z[i,1]==1 & z[i,2]==0){
      group[i] = 1
    }else if(z[i,1]==0 & z[i,2]==1){
      group[i] = 2
    }else{
      group[i] = 3
    }
  }
  
  gp1 = dat[group==1,]
  gp2 = dat[group==2,]
  gp3 = dat[group==3,]
  
  xmin = min(dat[,1])
  xmax = max(dat[,1])
  ymin = min(dat[,2])
  ymax = max(dat[,2])
  
  bic.mcme = bic3$out.mcme
  bic.mclust = bic3$out.mclust
  allbic = c(bic.mcme[,2],bic.mclust[1:7])
  ymin.bic = min(allbic)
  ymax.bic = max(allbic)

  #par(mfrow=c(1,2))
  plot(gp1,xlim=c(xmin,xmax),ylim=c(ymin,ymax),pch=0,
      main="Case 3: Three clusters close",xlab="",ylab="")
  points(gp2,pch=1)
  points(gp3,pch=2)

  plot(bic.mcme,type="b",lwd=2,ylim=c(ymin.bic,ymax.bic),
    xlab="number of components",ylab="BIC",pch=0)
  lines(1:7,bic.mclust[1:7],type="b",lty="dashed",pch=2,lwd=2)
  legend("bottomright",legend=c("MCLUST-ME","MCLUST"),
    pch=c(0,2),lty=c("solid","dashed"),cex=0.5)
  #par(mfrow=c(1,1))
}





setwd("/Users/wzhang/Research Project/Simulations for Paper/BIC simulation")
source("bic simulation.R")

setwd("/Users/wzhang/Research Project/Simulations for Paper/BIC simulation/Results")
# Well-separated case
res0 = bic_wellsep(51)
res1 = bic_wellsep(68)
res2 = bic_wellsep(79)
res3 = bic_wellsep(86)
res4 = bic_wellsep(97)
res5 = bic_wellsep(101)
res6 = bic_wellsep(125)
res7 = bic_wellsep(137)
bicres = list(res0=res0,res1=res1,res2=res2,res3=res3,res4=res4,res5=res5,
  res6=res6,res7=res7)
save(bicres,file="bic_well_separated.RData")
  
  
# Other two cases:
bic2_51 = bic_2close(51)
save(bic2_51,file="bic2_51.RData")
bic3_51 = bic_3close(51)
save(bic3_51,file="bic3_51.RData")

bic2_68 = bic_2close(68)
save(bic2_68,file="bic2_68.RData")
bic3_68 = bic_3close(68)
save(bic3_68,file="bic3_68.RData")

bic2_79 = bic_2close(79)
save(bic2_79,file="bic2_79.RData")
bic3_79 = bic_3close(79)
save(bic3_79,file="bic3_79.RData")

bic2_86 = bic_2close(86)
save(bic2_86,file="bic2_86.RData")
bic3_86 = bic_3close(86)
save(bic3_86,file="bic3_86.RData")


### Plot the data and BIC side-by-side
setwd("/Users/wzhang/Research Project/Simulations/Sim/BIC/Results")
load("bic2_51.RData")
plot.2close(51,bic2_51)

load("bic2_68.RData")
plot.2close(68,bic2_68)

load("bic2_79.RData")
plot.2close(79,bic2_79)

load("bic2_86.RData")
plot.2close(86,bic2_86)


load("bic3_51.RData")
plot.3close(51,bic3_51)

load("bic3_68.RData")
plot.3close(68,bic3_68)

load("bic3_79.RData")
plot.3close(79,bic3_79)

load("bic3_86.RData")
plot.3close(86,bic3_86)


load("bic_well_separated.RData")
plot.wellsep(51,bicres$res0)
plot.wellsep(68,bicres$res1)
plot.wellsep(79,bicres$res2)
plot.wellsep(86,bicres$res3)
plot.wellsep(97,bicres$res4)
plot.wellsep(101,bicres$res5)
plot.wellsep(125,bicres$res6)
plot.wellsep(137,bicres$res7)


# Plot for paper:
layout(matrix(1:6, 3, 2, byrow = TRUE))
plot.wellsep(101,bicres$res5)
load("bic2_79.RData")
plot.2close(79,bic2_79)
load("bic3_51.RData")
plot.3close(51,bic3_51)




