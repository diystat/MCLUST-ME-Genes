
## Simulation 5 result
time = c(47.97, 227.30, 988.31, 3813.01, 5293.02)
size = c(1,2,5,10,15) * 10

d = data.frame(time=time,size=size)
d


ll = lm(time~size)
summary(ll)

Coefficients:
  Estimate Std. Error t value Pr(>|t|)   
(Intercept) -549.246    264.997  -2.073  0.12991   
size          39.745      3.145  12.638  0.00107 **
# so increasing sample size by one increases running time by about 40 seconds, when dimension = 2
  

ggplot(data=d,aes(y=time,x=size)) +
  geom_point() +
  geom_smooth(method="lm",se=FALSE) +
  annotate("text", x=40, y=4000, label="slope = 39.475") +
  xlab("sample size") +
  ylab("time in seconds")

ggsave("sample_size_mcr.png", width=8, height=5, units="in", dpi=400)

