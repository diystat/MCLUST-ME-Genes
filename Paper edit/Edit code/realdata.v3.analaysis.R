#######################################################################
###------------------------ RNA-Seq Example ------------------------###
#######################################################################

## Linux working directory
## setwd("/home/yanming/ongoing/Model-based-Clustering-Research");

## Windows working directory
setwd("D:/ongoing/Model-based-Clustering-Research");

library(mclust);
library(phyclust)
library(psych)

## library(MASS)
## Load required packages
## library(gdata)
## library(caret)
## library(openintro)

## Load raw data
load("Data/rna_raw.RData")

## Load sampled data and the mclust results
m = 1000;
seed = 999;

compare.with.old = FALSE;

if (compare.with.old) {

  file.mclust = sprintf("Results/res_real2d_v2_m%d_seed%d_mclust.Rdata", m, seed);
  list(load(file.mclust));
  res.mclust.old = res.mclust;

  file.results = sprintf("Results/res_real2d_v3_m%d_seed%d.Rdata", m, seed);
  list(load(file.results));
  res3.old = res3;

}

date = "2020-01-15";
file.results = sprintf("Results/res_real2d_v3_m%d_seed%d.%s.Rdata", m, seed, date);
print(load(file.results));

if (compare.with.old) {
  identical(res.mclust, res.mclust.old);
  identical(res3, res3.old);
}

plot(res.mclust, "classification");


## Extract clustering uncertainty
unc.mcme = res0$uncertainty
unc.mclust = 1-apply(res.mclust$z, 1, max)

plot(unc.mclust, unc.mcme)

abline(0,1,lty="dashed")

## Summmarize error variation
## Use trace()
st = numeric(1000)
for(i in 1:1000){
  st[i] = tr(errary[,,i])
}

## Use det()
sd = numeric(1000)
for(i in 1:1000){
  sd[i] = det(errary[,,i])
}

s = st;

plot(log(s), res.mclust$z[,1]);

## Thu Jan 16 16:45:10 PST 2020
## Table and Figures for the paper revision

## 1. Contingency table for the paper revision (res4 vs res5 results were used)
table(res4$classification, res5$classification);

table(res4$classification, res0$classification);
table(res5$classification, res0$classification);


## 2. Comparing two sets of simulated error covariance estimates
## 2020-02-03
## 2. Comparing two sets of simulated standard errors of the estimated log fold changes
## used in the paper revision.
file.eps = "Graphs/real_data_simulated_se_beta1.eps";
postscript(file.eps, paper="special", width=6, height=6);

plot(sqrt(err4[1,1,]), sqrt(err5[1,1,]),
     main = "Simulated standard errors for the estimated \n log fold changes at 1h",
     xlim = c(0, 1),
     ylim = c(0, 1),
     xlab = "Set 1",
     ylab = "Set 2");

dev.off();


file.eps = "Graphs/real_data_simulated_se_beta2.eps";
postscript(file.eps, paper="special", width=6, height=6);

plot(sqrt(err4[2,2,]), sqrt(err5[2,2,]),
     main = "Simulated standard errors for the estimated \n log fold changes at 3h",
     xlim = c(0, 1),
     ylim = c(0, 1),
     xlab = "Set 1",
     ylab = "Set 2");

dev.off();


## 3. Compare two MCLUST-ME runs 
file.eps = "Graphs/real_data_change_in_z_mclust_vs_mclust-me.eps"; 
postscript(file.eps, paper="special", width=6, height=6);

plot(res.mclust$z[,1], res0$z[,1],
     main = "Estimated membership probability to the non-DE cluster",
     xlab = "MCLUST",
     ylab = "MCLUST-ME");

dev.off();

file.eps = "Graphs/real_data_change_in_z_mclust-me-resimulated.v.beta.eps"; 
postscript(file.eps, paper="special", width=6, height=6);
plot(res4$z[,1], res5$z[,1],
     main = "Estimated membership probability to the non-DE cluster",
     xlab = "MCLUST-ME with resimulated covariances (seed 1)",
     ylab = "MCLUST-ME with resimulated covariances (seed 2)");
dev.off();

## 4. A table summarizes the data and MCLUST/MCLUST-ME membership probabilities
plot(res.mclust, "classification");
table(res.mclust$classification, res0$classification);

## Summarize log fold changes, standard errors at timepoints 1h and 3h, and clustering results.
se1 = errary[1,1,];
se2 = errary[2,2,];
s = cbind(obs, se1, se2, res.mclust$z[,1], res0$z[,1]);
rownames(s) = rownames(nb.data$counts)[ss];
s[1:5,];

## Identify observations that are classified differently by the two methods
s12 = (1:m)[res.mclust$classification == 1 & res0$classification==2];
s21 = (1:m)[res.mclust$classification == 2 & res0$classification==1];

s[1:10,];
s[s12,];
s[s21,];


## For the paper, we will take 5 rows from each set:
summary.table = round(rbind(s[1:5,], s[s12[1:5],], s[s21[1:5],]), 3);

library(xtable);

xtable(summary.table);
sink(file = "Results/intro.summary.table.txt");
xtable(summary.table);

cat (sprintf("%s & %.3f (%.3f) & %.3f (%.3f) & %.3f & %.3f \\\\ \n",
        rownames(summary.table),
        summary.table[,1], summary.table[,3],
        summary.table[,2], summary.table[,4],
        summary.table[,5], summary.table[,6]));

sink();

## Look closer
library(NBPSeq);
nb.data[ss,][s12,];
nb.data[ss,][s21,];


## Additional results not used in the paper
## More detailed look at how the membership probability to cluster 1 has changed.
plot.change.in.z = function(res.a, res.b) {
  ## Use arrows to show estimated membership probability to cluster 1.
  unc.a = apply(res.a$z, 1, min);
  unc.b = apply(res.b$z, 1, min);
  par(mar=c(5.1, 4.1, 4.1, 4.1), xpd=TRUE)
  plot(s, res.a$z[,1], xlab="magnitude of error variance", ylab="changes to membership probability to cluster 1", 
       main="Changes in membership probability to cluster 1: MCLUST --> MCLUST-ME", type="n")
  ## points(s, res.b$z[,1]);
  for(i in 1:m){
    arrows(s[i], res.b$z[i,1], s[i], res.a$z[i,1], lwd=1.3, length=0.1,
           col=ifelse(unc.a[i]>unc.b[i],"red","dodgerblue"))
  }
  ## legend("topright", legend=c("decrease", "increase"), col=c("red","dodgerblue"),lty=1, inset=c(-0.13,0), cex=0.7)

  segments(-3, 0.5, 103, 0.5, lty=2);
  invisible();
}

plot.change.in.z(res.mclust, res0);
plot.change.in.z(res4, res0);
plot.change.in.z(res4, res5);
