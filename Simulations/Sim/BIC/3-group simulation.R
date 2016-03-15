

## Use BIC to select optimal model for 3-group data

# Set parameters for model:2 mixed clusters,1 well-separated from others
N = 300
tau1 = 0.3
tau2 = 0.4
mu1 = c(-3,0)
mu2 = c(3,0)
mu3 = c(0,24)
sig1 = matrix(c(64,0,0,64),nrow=2)
sig2 = matrix(c(16,0,0,16),nrow=2)
sig3 = matrix(c(36,0,0,36),nrow=2)
k = 6
p = 0.5

## seed = 51
res2.1 = test.3group(51,N,tau1,tau2,mu1,mu2,mu3,sig1,sig2,sig3,k,p)
#t1 = res$out.mcme
#t2 = res$out.mclust
#tp = cbind(1:7,t2[,1])
#rank.mcme = t1[order(t1[,2],decreasing = T),1]
#rank.mclust = tp[order(tp[,2],decreasing = T),1]
#rank.mcme == rank.mclust

## seed = 68
res2.2 = test.3group(68,N,tau1,tau2,mu1,mu2,mu3,sig1,sig2,sig3,k,p)

## seed = 79
res2.3 = test.3group(79,N,tau1,tau2,mu1,mu2,mu3,sig1,sig2,sig3,k,p)

## seed = 86
res2.4 = test.3group(86,N,tau1,tau2,mu1,mu2,mu3,sig1,sig2,sig3,k,p)

bicres = list(res1=res1,res2=res2,res3=res3,res4=res4)

setwd("/Users/wzhang/Research Project/Simulations/Sim/BIC/Results")
save(bicres,file="bicres_2.RData")




# Set parameters for model:2 mixed clusters,1 well-separated from others
N = 300
tau1 = 0.3
tau2 = 0.4
mu1 = c(-3,0)
mu2 = c(3,0)
mu3 = c(0,24)
sig1 = matrix(c(64,0,0,64),nrow=2)
sig2 = matrix(c(16,0,0,16),nrow=2)
sig3 = matrix(c(36,0,0,36),nrow=2)
k = 6
p = 0.5

## seed = 51
res2.1 = test.3group(51,N,tau1,tau2,mu1,mu2,mu3,sig1,sig2,sig3,k,p)
#t1 = res$out.mcme
#t2 = res$out.mclust
#tp = cbind(1:7,t2[,1])
#rank.mcme = t1[order(t1[,2],decreasing = T),1]
#rank.mclust = tp[order(tp[,2],decreasing = T),1]
#rank.mcme == rank.mclust

## seed = 68
res2.2 = test.3group(68,N,tau1,tau2,mu1,mu2,mu3,sig1,sig2,sig3,k,p)

## seed = 79
res2.3 = test.3group(79,N,tau1,tau2,mu1,mu2,mu3,sig1,sig2,sig3,k,p)

## seed = 86
res2.4 = test.3group(86,N,tau1,tau2,mu1,mu2,mu3,sig1,sig2,sig3,k,p)

bicres = list(res1=res1,res2=res2,res3=res3,res4=res4)

setwd("/Users/wzhang/Research Project/Simulations/Sim/BIC/Results")
save(bicres,file="bicres_2.RData")




























load("bicres.RData")


layout(matrix(1:8, 2, 4, byrow = TRUE))

plot(bicres$res0$out.mclust,legendArgs=list(plot=FALSE),lwd=2,col="magenta3")
lines(bicres$res0$out.mcme, type="b",col="cyan3",lty=2,lwd=2)
title(main="seed = 51")

plot(bicres$res1$out.mclust,legendArgs=list(plot=FALSE),lwd=2,col="magenta3")
lines(bicres$res1$out.mcme, type="b",col="cyan3",lty=2,lwd=2)
title(main="seed = 68")

plot(bicres$res2$out.mclust,legendArgs=list(plot=FALSE),lwd=2,col="magenta3")
lines(bicres$res2$out.mcme, type="b",col="cyan3",lty=2,lwd=2)
title(main="seed = 79")

plot(bicres$res3$out.mclust,legendArgs=list(plot=FALSE),lwd=2,col="magenta3")
lines(bicres$res3$out.mcme, type="b",col="cyan3",lty=2,lwd=2)
legend("bottomright",legend=c("MCLUST-ME","MCLUST"),col=c("cyan3","magenta3"),lty=c(2,1),lwd=c(2,2))
title(main="seed = 86")

