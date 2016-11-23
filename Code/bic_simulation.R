

#######################################################################
###------------------------ BIC Simulation -------------------------###
#######################################################################


## Source all necessary functions
source("core_functions.R")
source("simulation_functions.R")

## Load required packages
library(MASS)
library(gdata)
library(caret)
library(mclust)
library(phyclust)


## Three cases:
# Case 1: 3 clusters well separated
bic1_101 = bic_wellsep(101)

# Case 2: 2 clusters close
bic2_79 = bic_2close(79)

# Case 3: all 3 clusters close
bic3_51 = bic_3close(51)

##################################################################
### All simulations results are available in "Results" folder. ###
##################################################################




    ##########################################################
    ###################### Figure 9 ##########################
    ##########################################################

## Plot the optimal clusters and BIC values for each case
#png("bic_paper.png", width = 12, height = 14,units="in",res=300)
tiff("Fig9.tiff",width=1580,height=1580,compression="lzw",res=240,pointsize=7.5)
par(cex.main=2.3,cex.lab=2.3,cex.axis=2)
layout(matrix(1:6, 3, 2, byrow = TRUE))

load("Results/bic1_multiseeds.RData")
bic1_101 = bicres$res5
plot.wellsep(101,bic1_101)

load("Results/bic2_79.RData")
plot.2close(79,bic2_79)

load("Results/bic3_51.RData")
plot.3close(51,bic3_51)

par(cex.lab=1,cex.axis=1,cex.main=1)
dev.off()
