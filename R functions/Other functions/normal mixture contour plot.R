library(TeachingDemos)
# Bivariate Normal Density
# x: 2x1 vector, mu: 2x1 mean vector, Sigma: 2x2 covariance matrix
bivariate.normal <- function(x, mu1, Sigma1, mu2, Sigma2) {
    0.6*exp(-.5*t(x-mu1)%*%solve(Sigma1)%*%(x-mu1))/sqrt(2*pi*det(Sigma1))+
    0.4*exp(-.5*t(x-mu2)%*%solve(Sigma2)%*%(x-mu2))/sqrt(2*pi*det(Sigma2))
}
mu1 <- c(-1,0)
mu2 <- c(1,0)
Sigma1 <- Sigma2 <- matrix(c(1,0,0,1), nrow=2)
x <- y <- seq(-3, 3, len=250)
# Evaluate the bivariate normal density for each value of x and y
z <- outer(x, y,
FUN=function(x, y, ...){
apply(cbind(x,y), 1, bivariate.normal, ...)
}, mu1=mu1, Sigma1=Sigma1, mu2=mu2, Sigma2=Sigma2)
# Filled contour and surface plot of the bivariate normal density
contour(x,y,z, main="Bivariate Normal Density")
