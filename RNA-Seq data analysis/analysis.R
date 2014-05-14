

## Analysis of RNA-Seq data

  print(load("../Di's data/rnaseq_data.RData"))
  # "coef" contains 500 observations of dimension 5. dim=c(500,5)
  # "err" contains 5*5 estimation error matrices for each observation. dim=c(5,5,500)
  
  ## obtain initial classification:
  library(mclust)
  
  ## subset data:
  n = 10
  coef = coef[1:n,]
  err = err[,,1:n]
  
  G = 2 # number of clusters
  
  ## obtain an initial classification, using hierarchical clustering:
  hctree = hc("VVV", coef)
  temp = as.numeric(hclass(hctree, G))
    
  ini.class = matrix(0, n, G)
    for(i in 1:n){
      g = temp[i]
      ini.class[i,g] = 1
    }
  
  
  
  ## our method:
  tmp = proc.time()
  my.result = ME.VVV.err(coef, ini.class, err)
  proc.time() - tmp
  
  
  ## mclust's method:
  tmp = proc.time()
  mc.result = meVVV(coef, ini.class)
  proc.time() - tmp
  
