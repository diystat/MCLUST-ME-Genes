



R = t(z.ini)
Q1 = t(q1)
Q2 = t(q2)

alpha = seq(0,10,0.1)

yy = sapply(alpha,wtdfuzzyrand,R=R,Q=Q2,err=errmat)
zz = sapply(alpha,wtdfuzzyrand,R=R,Q=Q1,err=errmat)

yup = max(c(yy,zz))+0.01
ylow = min(c(yy,zz))-0.01
plot(alpha,zz,type="l",ylim=c(ylow,yup),main="weighted fuzzy rand index",
  ylab="index")
lines(alpha,yy,lty="dashed")
legend("bottomright",legend=c("meVVV","MCME"),lty=c(1,2))



