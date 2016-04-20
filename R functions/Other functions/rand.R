
# Calculate rand indices
get.rand = function(simrun){
  
  z.ini = simrun$z.ini
  z.mcme = simrun$res.mcme$z
  z.mevvv = simrun$res.mevvv$z
  index = simrun$index
  
  N = nrow(z.mcme)
  G = ncol(z.mcme)
  
  # Discretize membership matrix:
  z1 = z2 = matrix(0,N,G)
  for(i in 1:N){
    rm.mcme = max(z.mcme[i,])
    rm.mevvv = max(z.mevvv[i,])
    for(k in 1:G){
      z1[i,k] = ifelse(z.mcme[i,k]==rm.mcme,1,0)
      z2[i,k] = ifelse(z.mevvv[i,k]==rm.mevvv,1,0)
    }  
  }
  
  z.ini1 = z.ini[,1] + 1
  z1 = z1[,1] + 1
  z2 = z2[,1] + 1
  
  ## Calculate Rand index:
  mcme = round(RRand(z.ini1,z1)$Rand,4)
  mcme.me = round(RRand(z.ini1[which(index==1)],z1[which(index==1)])$Rand,4)
  mcme.nn = round(RRand(z.ini1[which(index==0)],z1[which(index==0)])$Rand,4)
  mevvv = round(RRand(z.ini1,z2)$Rand,4)
  mevvv.me = round(RRand(z.ini1[which(index==1)],z2[which(index==1)])$Rand,4)
  mevvv.nn = round(RRand(z.ini1[which(index==0)],z2[which(index==0)])$Rand,4)
  mcme.f = round(fuzzyrand(t(z.ini),t(z.mcme)),4)
  mevvv.f = round(fuzzyrand(t(z.ini),t(z.mevvv)),4)
  
  out = c(mcme,mcme.me,mcme.nn,mevvv,mevvv.me,mevvv.nn,mcme.f,mevvv.f)
  return(out)
  
}