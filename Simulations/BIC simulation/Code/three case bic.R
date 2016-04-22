

#######################################################################
###------------------------ BIC Simulation -------------------------###
#######################################################################


## Three cases:
# Case 1: 3 clusters well separated
bic1_101 = bic_wellsep(101)

# Case 2: 2 clusters close
bic2_79 = bic_2close(79)

# Case 3: all 3 clusters close
bic3_51 = bic_3close(51)

#########################################################
### All simulations results are available in:         ###
### "Simulations/BIC Simulation/Sim 2/Results".       ###
#########################################################




    ##########################################################
    ###################### Figure 9 ##########################
    ##########################################################

## Plot the optimal clusters and BIC values for each case
layout(matrix(1:6, 3, 2, byrow = TRUE))
plot.wellsep(101,bic1_101)
plot.2close(79,bic2_79)
plot.3close(51,bic3_51)
par(mfrow=c(1,1))