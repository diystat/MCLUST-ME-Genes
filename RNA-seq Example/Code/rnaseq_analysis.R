

#######################################################################
###------------------------ RNA-Seq Example ------------------------###
#######################################################################



# Load required packages
library(MASS)
library(gdata)
library(caret)
library(mclust)
library(phyclust)


###------------------- Importing Data ---------------------###

###################################################################
### Raw and processed data available in "RNA-seq Example/Data". ###
###################################################################

# Import raw data
load("rna_raw.RData")
obs = full$beta[top, 1:5] # data points
errary = v.beta[,,top] # measurement error covariances

# Export processed data
# save(obs,errary,file="rna_processed.RData")




###------------------- Clustering ---------------------###

# Cluster the data with MCLUST first:
res.mclust = Mclust(obs,modelNames="VVV") # Running time < 5min
summary(res.mclust)

# Define a wrapper function for MCLUST-ME:
mclust.run = function(G){
  # Generate initial membership matrix with hierarchical agglomeration
  N = nrow(obs)
  hcTree = hc(obs)
  cl = hclass(hcTree, G)
  z.ini = matrix(0,N,G)
  for(i in 1:N){
    for(j in 1:G){
      z.ini[i,j] = ifelse(cl[i]==j,1,0)
    }
  }
  # Group data with MCLUST-ME:
  res.mcme.full = mcmeVVV(obs,z.ini,errary)
  return(res.mcme.full)  
}


###########################################################################
### WARNING: Code from line 68 to line 75 are EXTREMELY time-consuming. ###
### The author strongly recommends using results in the folder          ###
### "RNA-seq Example/Results" for further analysis.                     ###
###########################################################################

# Run MCLUST-ME assuming 1~8 clusters:
run1 = mclust.run(1)
run2 = mclust.run(2)
run3 = mclust.run(3)
run4 = mclust.run(4)
run5 = mclust.run(5)
run6 = mclust.run(6)
run7 = mclust.run(7)
run8 = mclust.run(8)
####################################################################
### The BIC results above are stored in files named "run*.RData" ###
### under "RNA-seq Example/Results/BIC".                         ###
####################################################################


# Extract BIC values
b = numeric(8)
b[1] = run1$BIC
b[2] = run2$BIC
b[3] = run3$BIC
b[4] = run4$BIC
b[5] = run5$BIC
b[6] = run6$BIC
b[7] = run7$BIC
b[8] = run8$BIC

a = res.mclust$BIC[1:8]
ylow = min(c(a,b))
yup = max(c(a,b))

# Plot BIC values of both methods
plot(1:8,b,type="b",lwd=2,xlab="Number of components",ylab="BIC",pch=0,
  ylim=c(ylow,yup))
lines(1:8,a,type="b",lwd=2,pch=1,lty=2)
legend("bottomright",legend=c("MCLUST-ME","MCLUST"),
    pch=c(0,2),lty=c("solid","dashed"),lwd=2)

# BIC is greatest with 3 components, so keep this result:
res.mcme.full = run3
##########################################################
### The MLE and membership probabilities are stored in ###
### the RData file "cluster_res.RData", under          ###
### "RNA-seq Example/Results/MLE and membership".      ###
##########################################################





###---------------------- Analysis -----------------------###

### Obtain mean and covariance estimates
# MCLUST-ME
param = res.mcme.full$parameters
mu.mcme = param$muhat # centers
sigma.mcme = param$sigmahat # covariances

# MCLUST
param = res.mclust$parameters
mu.mclust = param$mean # centers
sigma.mclust = param$variance$sigma # covariances


### Predicted class labels:
# MCLUST-ME labels
zhat = res.mcme.full$z
predClass = numeric(1000)
for(i in 1:1000){
  temp = zhat[i,]
  predClass[i] = which(temp==max(temp))
}

# MCLUST labels
zhat.mclust = res.mclust$z
predClass.mclust = numeric(1000)
for(i in 1:1000){
  temp = zhat.mclust[i,]
  predClass.mclust[i] = which(temp==max(temp))
}






### Confusion matrix:
confusionMatrix(predClass,predClass.mclust,dnn=c("MCME","MCLUST"))

### Rand index:
phyclust::RRand(predClass,predClass.mclust)






    ##########################################################
    ###################### Figure 10 #########################
    ##########################################################

### Raw data scatterplot matrix:
pairs(obs,labels=c("10min","1h","3h","6h","12h"),pch=16,cex=0.6)






    ##########################################################
    ###################### Figure 11 #########################
    ##########################################################

### Partitioned data on 1h and 3h with confidence outlines
library(ellipse)
# Generate coordinates for confidence outlines(MCLUST-ME)
g1.coord = ellipse(sigma.mcme[2:3,2:3,1],centre=mu.mcme[2:3,1],level=0.9)
g2.coord = ellipse(sigma.mcme[2:3,2:3,2],centre=mu.mcme[2:3,2],level=0.9)
g3.coord = ellipse(sigma.mcme[2:3,2:3,3],centre=mu.mcme[2:3,3],level=0.9)
# Generate coordinates for confidence outlines(MCLUST)
g1.coord.mc = ellipse(sigma.mclust[2:3,2:3,1],centre=mu.mclust[2:3,1],level=0.9)
g2.coord.mc = ellipse(sigma.mclust[2:3,2:3,2],centre=mu.mclust[2:3,2],level=0.9)
g3.coord.mc = ellipse(sigma.mclust[2:3,2:3,3],centre=mu.mclust[2:3,3],level=0.9)

par(mfrow=c(2,2))
plot(1:8,b,main="BIC of MCLUST-ME",type="b",pch=0,lwd=2,xlab="Number of components",ylab="BIC")
lines(1:8,a,lty="dashed",type="b",pch=0,lwd=2)

plot(obs[,2],obs[,3],main="MCLUST-ME",pch=c(1,2,3)[predClass],cex=0.7,
  xlab="1h",ylab="3h",col="#00000070")
# Plot the confidence outlines
lines(g1.coord,type="l")
lines(g2.coord,type="l")
lines(g3.coord,type="l")

plot(obs[,2],obs[,3],main="MCLUST",pch=c(1,2,3)[predClass.mclust],cex=0.7,
  xlab="1h",ylab="3h",col="#00000070")
# Plot the confidence outlines
lines(g1.coord.mc,type="l")
lines(g2.coord.mc,type="l")
lines(g3.coord.mc,type="l")
lines(g1.coord,type="l",lty="dashed")
lines(g2.coord,type="l",lty="dashed")
lines(g3.coord,type="l",lty="dashed")
par(mfrow=c(1,1))






    ##########################################################
    ###################### Figure 12 #########################
    ##########################################################

### Line plot of grouped genes
g1 = obs[predClass==1,]
g2 = obs[predClass==2,]
g3 = obs[predClass==3,]
x = 1:5

ymin1 = min(g1)
ymax1 = max(g1)
ymin2 = min(g2)
ymax2 = max(g2)
ymin3 = min(g3)
ymax3 = max(g3)

mcme.g1 = which(predClass==1)
mcme.g2 = which(predClass==2)
mcme.g3 = which(predClass==3)

mclust.g1 = which(predClass.mclust==1)
mclust.g2 = which(predClass.mclust==2)
mclust.g3 = which(predClass.mclust==3)

g1m = obs[predClass.mclust==1,]
g2m = obs[predClass.mclust==2,]
g3m = obs[predClass.mclust==3,]

ymin1m = min(g1m)
ymax1m = max(g1m)
ymin2m = min(g2m)
ymax2m = max(g2m)
ymin3m = min(g3m)
ymax3m = max(g3m)

ymin = min(c(ymin1,ymin2,ymin3,ymin1m,ymin2m,ymin3m))
ymax = max(c(ymax1,ymax2,ymax3,ymax1m,ymax2m,ymax3m))


# Function for drawing line plots
bwlineplot = function(y,groupnum,method="MCLUST-ME"){
  x = 1:5
  d = rnorm(1000,0,0.01)
  eps1 = 0.01
  eps2 = 0.02
  title = paste(method," Group ",groupnum,sep="")

  # Individual lines with error bars
  plot(x,y[1,],xaxt="n",type="l",ylim=c(ymin,ymax),col="grey78",xlab="",
    ylab="Log2 fold change",main=title,lwd=0.5)
  sd1 = sqrt(diag(errary[,,1]))
  segments(x+d[1],y[1,]-sd1,x+d[1],y[1,]+sd1,lwd=0.5,col="grey78") # add 1 sd error bar
  segments(x+d[1]-eps1,y[1,]-sd1,x+d[1]+eps1,y[1,]-sd1,lwd=0.5,col="grey78")
  segments(x+d[1]-eps1,y[1,]+sd1,x+d[1]+eps1,y[1,]+sd1,lwd=0.5,col="grey78")
  axis(1,at=1:5,labels=c("10min","1h","3h","6h","12h"))
  for(i in 2:nrow(y)){
    lines(x,y[i,],col="grey78",lwd=0.5)
    sd1 = sqrt(diag(errary[,,i]))
    segments(x+d[i],y[i,]-sd1,x+d[i],y[i,]+sd1,lwd=0.5,col="grey78") # add 1 sd error bar
    segments(x+d[i]-eps1,y[i,]-sd1,x+d[i]+eps1,y[i,]-sd1,lwd=0.5,col="grey78")
    segments(x+d[i]-eps1,y[i,]+sd1,x+d[i]+eps1,y[i,]+sd1,lwd=0.5,col="grey78")
  }

  # group mean with error bar
  lines(x,mu.mcme[,groupnum],lwd=2) 
  sd1.mcme = sqrt(diag(sigma.mcme[,,groupnum]))
  mu1.mcme = mu.mcme[,groupnum]
  segments(x,mu1.mcme-sd1.mcme,x,mu1.mcme+sd1.mcme,lwd=2) # add 1 sd error bar
  segments(x-eps2,mu1.mcme-sd1.mcme,x+eps2,mu1.mcme-sd1.mcme,lwd=2)
  segments(x-eps2,mu1.mcme+sd1.mcme,x+eps2,mu1.mcme+sd1.mcme,lwd=2)
}

# Black and white line plots
par(mfrow=c(2,3))

bwlineplot(g1,1) # MCLUST-ME group 1
bwlineplot(g2,2) # MCLUST-ME group 2
bwlineplot(g3,3) # MCLUST-ME group 3
bwlineplot(g1m,1,"MCLUST") # MCLUST-ME group 1
bwlineplot(g2m,2,"MCLUST") # MCLUST-ME group 2
bwlineplot(g3m,3,"MCLUST") # MCLUST-ME group 3
  
par(mfrow=c(1,1))


