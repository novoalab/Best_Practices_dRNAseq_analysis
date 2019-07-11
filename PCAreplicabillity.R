#PCA replicability
basecaller <- "GU_3.0.3"
mapper <- "graphmap"
dataname1="UNM_C"
dataname2="hm5C"
plot_PCA_replicates <- function(dataname){
rep1 <- read.delim(paste(basecaller, "_", mapper, "/replicate1/", dataname, "_per_site.var.csv.slided.onekmer.oneline.5mer.filtered.csv", sep=""), sep=",")
rep2 <- read.delim(paste(basecaller, "_", mapper, "/replicate2/", dataname, "_per_site.var.csv.slided.onekmer.oneline.5mer.filtered.csv", sep=""), sep=",")

rep1$modification <- "rep1"
rep2$modification <- "rep2"
data <- rbind(rep1, rep2)

data$modification <- as.factor(data$modification)
pca = prcomp(data[,c(5,6,7,8,9,10,11,12,13,14,20,21,22,23,24)])
eigs <- pca$sdev^2
p1 <- eigs[1] / sum(eigs)
p2 <- eigs[2] / sum(eigs)
p <- plot(pca$x, col=data$modification, ylab=paste("PC2 (",round(p2,2)*100,"%)"), xlab=paste("PC1 (",round(p1,2)*100,"%)"), pch=19, xaxt='n', yaxt='n', main=paste(dataname, "(rep1 vs rep2)"))
return(p)
}

plot_PCA_non_replicates <- function(dataname1, dataname2){
  rep1 <- read.delim(paste(basecaller, "_", mapper, "/replicate1/", dataname1, "_per_site.var.csv.slided.onekmer.oneline.5mer.filtered.csv", sep=""), sep=",")
  rep2 <- read.delim(paste(basecaller, "_", mapper, "/replicate1/", dataname2, "_per_site.var.csv.slided.onekmer.oneline.5mer.filtered.csv", sep=""), sep=",")
  rep1$modification <- "rep1"
  rep2$modification <- "rep2"
  data <- rbind(rep1, rep2)
  data$modification <- as.factor(data$modification)
  pca = prcomp(data[,c(5,6,7,8,9,10,11,12,13,14,20,21,22,23,24)])
  eigs <- pca$sdev^2
  p1 <- eigs[1] / sum(eigs)
  p2 <- eigs[2] / sum(eigs)
  p <- plot(pca$x, col=data$modification, ylab=paste("PC2 (",round(p2,2)*100,"%)"), xlab=paste("PC1 (",round(p1,2)*100,"%)"), pch=19, xaxt='n', yaxt='n', main=paste(dataname1, "(rep1) vs", dataname2, "(rep1)"))
  return(p)
}

pdf(paste(dataname1, "_", dataname2, "_", basecaller, "_", mapper,"_PCA_replicability.pdf", sep=""), width=9, height=3)
par(mfrow = c(1,3))
plot_PCA_replicates(dataname1)
plot_PCA_replicates(dataname2)
plot_PCA_non_replicates(dataname1,dataname2)
dev.off()

