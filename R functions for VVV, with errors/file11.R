

## vector of sample size multipliers:
s = c(5, 10, 15, 20, 25, 30, 50, 100, 150, 200, 250)
#s = c(5, 6)
lapply(s, sim.goodini)

lapply(s, sim.badini)


foo = function(x){
  return(x^2)
}

lapply(s, foo)
