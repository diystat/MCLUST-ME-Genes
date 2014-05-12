

## Analysis of Di's RNA-Seq data

  print(load("../Di's data/rnaseq_data.RData"))
  # "coef" contains 500 observations of dimension 5. dim=c(500,5)
  # "err" contains 5*5 estimation error matrices for each observation. dim=c(5,5,500)
  
  ## obtain initial classification:
  library(mclust)
  
  hctree = hc("VVV", coef)
  temp = as.numeric(hclass(hctree, 8))
    
  ini.class = matrix(0, 500, 8)
    for(i in 1:500){
      g = temp[i]
      ini.class[i,g] = 1
    }
  
  
  ## our method:
  tmp = proc.time()
  my.class = ME.VVV.err(coef, ini.class, err)$z
  proc.time() - tmp
  
  
  ## mclust's method:
  tmp = proc.time()
  mc.class = meVVV(coef, ini.class)$z
  proc.time() - tmp
  
  
  
  
  
  
  






