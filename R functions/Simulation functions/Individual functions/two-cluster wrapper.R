    
### Wrapper functions for two-cluster simulation

### Simulation 1 wrapper    
sim1 = function(p){
  N = 200
  tau = 0.5
  mu1 = c(0,0)
  mu2 = c(8,0)
  sig1 = matrix(c(64,0,0,64),nrow=2)
  sig2 = matrix(c(16,0,0,16),nrow=2)
  k = 36
  nseed = 100

  # Check if singularity occurs for any seed
  check.seed(N,tau,mu1,mu2,sig1,sig2,k,p)

  # Run the simulation
  out = sim.driver(N,tau,mu1,mu2,sig1,sig2,k,p,nseed)
  return(out)
}        
    
    

### Simulation 2 wrapper    
sim2 = function(p){
  N = 200
  tau = 0.5
  mu1 = c(0,0)
  mu2 = c(8,0)
  sig1 = matrix(c(64,0,0,64),nrow=2)
  sig2 = matrix(c(16,0,0,16),nrow=2)
  k = 9
  nseed = 100

  # Check if singularity occurs for any seed
  check.seed(N,tau,mu1,mu2,sig1,sig2,k,p)

  # Run the simulation
  out = sim.driver(N,tau,mu1,mu2,sig1,sig2,k,p,nseed)
  return(out)
}        
    