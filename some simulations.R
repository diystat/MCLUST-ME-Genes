sigma.y = matrix(c(2,1,1,4),nrow=2)
sigma.x = matrix(c(3,0.5,0.5,2),nrow=2)
sigma.e = matrix(c(5,1,1,5),nrow=2)
mu.y = c(2,3)
mu.x = c(7,9)
nsim = 50

library(MASS)
y = mvrnorm(nsim, mu.y, sigma.y)
x = mvrnorm(nsim, mu.x, sigma.y)

cov(rbind(y,x))

left.lim = min(min(x[,1],y[,1]))-1
right.lim = max(max(x[,1],y[,1]))+1
upper.lim = max(max(x[,2],y[,2]))+1
lower.lim = min(min(x[,2],y[,2]))-1

plot(y[,1],y[,2],xlim=c(left.lim,right.lim),ylim=c(lower.lim,upper.lim),
     col="blue")
points(x[,1],x[,2],col="red")

label = c(rep("y",nsim),rep("x",nsim))

data = data.frame(label,rbind(y,x))

mod = Mclust(data[,-1])
plot(mod, data[,-1], what="BIC")

coordProj(data[,-1], dimens=c(1,2), what="classification",
          classification=mod$classification,
          parameters=mod$parameters)

coordProj(data[,-1], dimens=c(1,2), what="errors",
          classification=mod$classification,
          parameters=mod$parameters,
          truth=data[,1])

