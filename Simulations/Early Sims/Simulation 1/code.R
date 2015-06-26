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
  z.ini = res500.ss[,i]$z.ini
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














png("50_ef.png",900,1500)
par(mfrow=c(5,3))
res50.ef = sapply(seedvec,sim1.ss,n=50,k=0)
par(mfrow=c(1,1))
dev.off()

png("50_ef_com.png",900,1500)
par(mfrow=c(5,3))
for(i in 1:5){
  z.ini = res50.ef[,i]$z.ini
  z1 = res50.ef[,i]$z1
  z2 = res50.ef[,i]$z2
  group1 = res50.ef[,i]$group1
  group2 = res50.ef[,i]$group2
  cp = sim1.complot(z.ini,z1,z2,group1,group2)
}
par(mfrow=c(1,1))
dev.off()



png("100_ef.png",900,1500)
par(mfrow=c(5,3))
res100.ef = sapply(seedvec,sim1.ss,n=100,k=0)
par(mfrow=c(1,1))
dev.off()

png("100_ef_com.png",900,1500)
par(mfrow=c(5,3))
for(i in 1:5){
  z.ini = res100.ef[,i]$z.ini
  z1 = res100.ef[,i]$z1
  z2 = res100.ef[,i]$z2
  group1 = res100.ef[,i]$group1
  group2 = res100.ef[,i]$group2
  cp = sim1.complot(z.ini,z1,z2,group1,group2)
}
par(mfrow=c(1,1))
dev.off()



png("200_ef.png",900,1500)
par(mfrow=c(5,3))
res200.ef = sapply(seedvec,sim1.ss,n=200,k=0)
par(mfrow=c(1,1))
dev.off()

png("200_ef_com.png",900,1500)
par(mfrow=c(5,3))
for(i in 1:5){
  z.ini = res200.ef[,i]$z.ini
  z1 = res200.ef[,i]$z1
  z2 = res200.ef[,i]$z2
  group1 = res200.ef[,i]$group1
  group2 = res200.ef[,i]$group2
  cp = sim1.complot(z.ini,z1,z2,group1,group2)
}
par(mfrow=c(1,1))
dev.off()



png("500_ef.png",900,1500)
par(mfrow=c(5,3))
res500.ef = sapply(seedvec,sim1.ss,n=500,k=0)
par(mfrow=c(1,1))
dev.off()

png("500_ef_com.png",900,1500)
par(mfrow=c(5,3))
for(i in 1:5){
  z.ini = res500.ef[,i]$z.ini
  z1 = res500.ef[,i]$z1
  z2 = res500.ef[,i]$z2
  group1 = res500.ef[,i]$group1
  group2 = res500.ef[,i]$group2
  cp = sim1.complot(z.ini,z1,z2,group1,group2)
}
par(mfrow=c(1,1))
dev.off()



png("800_ef.png",900,1500)
par(mfrow=c(5,3))
res800.ef = sapply(seedvec,sim1.ss,n=800,k=0)
par(mfrow=c(1,1))
dev.off()
  
png("800_ef_com.png",900,1500)
par(mfrow=c(5,3))
for(i in 1:5){
  z.ini = res800.ef[,i]$z.ini
  z1 = res800.ef[,i]$z1
  z2 = res800.ef[,i]$z2
  group1 = res800.ef[,i]$group1
  group2 = res800.ef[,i]$group2
  cp = sim1.complot(z.ini,z1,z2,group1,group2)
}
par(mfrow=c(1,1))
dev.off()


out = list(res50=res50.ef,res100=res100.ef,res200=res200.ef,res500=res500.ef,res800=res800.ef)
save(out,file="sim1_ef_res")



