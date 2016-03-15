# another option
makePairs <- function(data) 
{
  temp = c("10min","1h","3h","6h","12h")
  grid <- expand.grid(x = 1:ncol(data), y = 1:ncol(data))
  grid <- subset(grid, x != y)
  all <- do.call("rbind", lapply(1:nrow(grid), function(i) {
    xcol <- grid[i, "x"]
    ycol <- grid[i, "y"]
    data.frame(xvar = temp[ycol], yvar = temp[xcol], 
               x = data[, xcol], y = data[, ycol], data)
  }))
  
  all$xvar <- factor(all$xvar, levels = c("10min","1h","3h","6h","12h"))
  all$yvar <- factor(all$yvar, levels = c("10min","1h","3h","6h","12h"))
  densities <- do.call("rbind", lapply(1:ncol(data), function(i) {
    data.frame(xvar = temp[i], yvar = temp[i], x = data[, i])
  }))
  list(all=all, densities=densities)
}
 
# expand iris data frame for pairs plot
pc = as.character(predClass)
obs.comb = data.frame(label=pc,tenmin=obs[,1],oneh=obs[,2],
  threeh=obs[,3],sixh=obs[,4],twh=obs[,5])
gg1 = makePairs(obs.comb[,-1])
 
# new data frame mega iris
mega_obs = data.frame(gg1$all, Groups=rep(obs.comb$label, length=nrow(gg1$all)))
 
# pairs plot
png("realdata_mcme2.png",1500,1200)
ggplot(mega_obs,aes_string(x="x",y="y")) + 
  facet_grid(xvar ~ yvar, scales = "free") + 
  geom_point(aes(colour=Groups), na.rm = TRUE, alpha=0.8)
dev.off()

# produces high-resolution plot for LATEX
ggsave("realdata_mcme.png", width=7, height=4, units="in", dpi=300) 




