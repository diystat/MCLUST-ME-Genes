
setwd("/Users/wzhang/Research Project/Simulations for Paper/BIC simulation/Results")
# Well-separated case
res0 = bic_wellsep(51)
res1 = bic_wellsep(68)
res2 = bic_wellsep(79)
res3 = bic_wellsep(86)
res4 = bic_wellsep(97)
res5 = bic_wellsep(101)
res6 = bic_wellsep(125)
res7 = bic_wellsep(137)
bicres = list(res0=res0,res1=res1,res2=res2,res3=res3,res4=res4,res5=res5,
  res6=res6,res7=res7)
save(bicres,file="bic_well_separated.RData")
  
  
# Other two cases:
bic2_51 = bic_2close(51)
save(bic2_51,file="bic2_51.RData")
bic3_51 = bic_3close(51)
save(bic3_51,file="bic3_51.RData")

bic2_68 = bic_2close(68)
save(bic2_68,file="bic2_68.RData")
bic3_68 = bic_3close(68)
save(bic3_68,file="bic3_68.RData")

bic2_79 = bic_2close(79)
save(bic2_79,file="bic2_79.RData")
bic3_79 = bic_3close(79)
save(bic3_79,file="bic3_79.RData")

bic2_86 = bic_2close(86)
save(bic2_86,file="bic2_86.RData")
bic3_86 = bic_3close(86)
save(bic3_86,file="bic3_86.RData")


### Plot the data and BIC side-by-side
setwd("/Users/wzhang/Research Project/Simulations/Sim/BIC/Results")
load("bic2_51.RData")
plot.2close(51,bic2_51)

load("bic2_68.RData")
plot.2close(68,bic2_68)

load("bic2_79.RData")
plot.2close(79,bic2_79)

load("bic2_86.RData")
plot.2close(86,bic2_86)


load("bic3_51.RData")
plot.3close(51,bic3_51)

load("bic3_68.RData")
plot.3close(68,bic3_68)

load("bic3_79.RData")
plot.3close(79,bic3_79)

load("bic3_86.RData")
plot.3close(86,bic3_86)


load("bic_well_separated.RData")
plot.wellsep(51,bicres$res0)
plot.wellsep(68,bicres$res1)
plot.wellsep(79,bicres$res2)
plot.wellsep(86,bicres$res3)
plot.wellsep(97,bicres$res4)
plot.wellsep(101,bicres$res5)
plot.wellsep(125,bicres$res6)
plot.wellsep(137,bicres$res7)


# Plot for paper:
layout(matrix(1:6, 3, 2, byrow = TRUE))
plot.wellsep(101,bicres$res5)
load("bic2_79.RData")
plot.2close(79,bic2_79)
load("bic3_51.RData")
plot.3close(51,bic3_51)




