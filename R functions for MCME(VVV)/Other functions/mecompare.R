

## Function that compares outputs of two functions, meVVV and stepwiseME
mecompare = function(data,z){
  temp1 = meVVV(data,z)
  temp2 = stepwiseME(data,z)
  itmax = min(as.numeric(attributes(temp1)$info[1]),as.numeric(temp2$it)) # number of iterations meVVV used
  
  pro = mean = sigma = result = list()
  
  for(i in 1:itmax){ 
    res.step = stepwiseME(data,z,itmax=i)
    pro.step = res.step$par$pro
    mean.step = res.step$par$mean
    sigma.step = res.step$par$var$sigma
    
    res.mevvv = meVVV(data,z,control=emControl(itmax=i))
    pro.mevvv = res.mevvv$par$pro
    mean.mevvv = res.mevvv$par$mean
    sigma.mevvv = res.mevvv$par$var$sigma
    
    pro = list(pro.step=pro.step,pro.mevvv=pro.mevvv)
    mean = list(mean.step=mean.step,mean.mevvv=mean.mevvv)
    sigma = list(sigma.step=sigma.step,sigma.mevvv=sigma.mevvv)
    
    combined = list(iteration=i,pro=pro,mean=mean,sigma=sigma)
    
    result = append(result,combined)
  }
    return(result)
}

  