
## Vector of random seeds:
seedvec = c(91,92,93,94,95)





### Simulation 2 3d case:

## lambda = 1:
imgname50 = sapply(seedvec,paste,"-50.png",sep="")
res50 = list()
for(i in 1:5){
  res50[[i]] = sim2.3d(seedvec[i],n=50,lambda=1,imgname50[i])
}




imgname100 = sapply(seedvec,paste,"-100.png",sep="")
res100 = list()
for(i in 1:5){
  res100[[i]] = sim2.3d(seedvec[i],n=100,lambda=1,imgname100[i])
}



imgname200 = sapply(seedvec,paste,"-200.png",sep="")
res200 = list()
for(i in 1:5){
  res200[[i]] = sim2.3d(seedvec[i],n=200,lambda=1,imgname200[i])
}



imgname500 = sapply(seedvec,paste,"-500.png",sep="")
res500 = list()
for(i in 1:5){
  res500[[i]] = sim2.3d(seedvec[i],n=500,lambda=1,imgname500[i])
}



imgname800 = sapply(seedvec,paste,"-800.png",sep="")
res800 = list()
for(i in 1:5){
  res800[[i]] = sim2.3d(seedvec[i],n=800,lambda=1,imgname800[i])
}




out1 = list(res50=res50,res100=res100,res200=res200,
  res500=res500,res800=res800)
save(out1,file="sim2-3d-1")








## lambda = 0.1:
imgname50 = sapply(seedvec,paste,"-50.1.png",sep="")
res50.1 = list()
for(i in 1:5){
  res50.1[[i]] = sim2.3d(seedvec[i],n=50,lambda=0.1,imgname50[i])
}




imgname100 = sapply(seedvec,paste,"-100.1.png",sep="")
res100.1 = list()
for(i in 1:5){
  res100.1[[i]] = sim2.3d(seedvec[i],n=100,lambda=0.1,imgname100[i])
}



imgname200 = sapply(seedvec,paste,"-200.1.png",sep="")
res200.1 = list()
for(i in 1:5){
  res200.1[[i]] = sim2.3d(seedvec[i],n=200,lambda=0.1,imgname200[i])
}



imgname500 = sapply(seedvec,paste,"-500.1.png",sep="")
res500.1 = list()
for(i in 1:5){
  res500.1[[i]] = sim2.3d(seedvec[i],n=500,lambda=0.1,imgname500[i])
}



imgname800 = sapply(seedvec,paste,"-800.1.png",sep="")
res800.1 = list()
for(i in 1:5){
  res800.1[[i]] = sim2.3d(seedvec[i],n=800,lambda=0.1,imgname800[i])
}



out2 = list(res50=res50.1,res100=res100.1,res200=res200.1,
  res500=res500.1,res800=res800.1)
save(out2,file="sim2-3d-2")







    #-----------------------------------------------------------------#





### Simulation 2 with fixed centers:



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




out = list(res50=res50,res100=res100,res200=res200,res500=res500,res800=res800)
save(out,file="sim2res-1.fc")




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



out = list(res50=res50_2,res100=res100_2,res200=res200_2,res500=res500_2,res800=res800_2)
save(out,file="sim2res-2.fc")





