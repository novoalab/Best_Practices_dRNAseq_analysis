library(optparse)

parser <- OptionParser()
parser <- add_option(parser, opt_str=c("-u", "--unmodified"), type="character",
                     dest='unmodified',
                     help="Unmodified .csv file")
parser <- add_option(parser, opt_str=c("-m", "--modified"), type="character",
                     dest='modified',
                     help="Modified .csv file")
parser <- add_option(parser, opt_str=c("-n", "--names"), type="character",
                     dest='names',
                     help="Names")
parser <- add_option(parser, opt_str=c("-c", "--columns"), type="character",
                     dest='columns',
                     help="columns")
parser <- add_option(parser, opt_str=c("-s", "--scaling"), type="boolean",
                     action="store_true",
                     dest='scaling',
                     help="scaling")

options=parse_args(parser)
unmodified=options$unmodified
modified=options$modified
scaling=options$scaling
names=options$names
columns=options$columns

# unmodified <- "GU_3.0.3_graphmap/replicate1/UNM_A_per_site.var.csv.slided.onekmer.oneline.5mer.filtered.csv"
# modified <- "GU_3.0.3_graphmap/replicate1/m6A_per_site.var.csv.slided.onekmer.oneline.5mer.filtered.csv"
# names <- "unm,m6A"
# columns <- "5,6,7,8,9,10,11,12,13,14,20,21,22,23,24"
# scaling=TRUE

n <- strsplit(names, ",")[[1]]
c <- as.numeric(strsplit(columns, ",")[[1]])

unmodified <- read.delim("GU_3.0.3_graphmap/replicate1/UNM_A_per_site.var.csv.slided.onekmer.oneline.5mer.filtered.csv", sep = ",")
modified <- read.delim("GU_3.0.3_graphmap/replicate1/m6A_per_site.var.csv.slided.onekmer.oneline.5mer.filtered.csv", sep = ",")

unmodified$modification <- "UNM"
modified$modification <- "MOD"

data <- rbind(unmodified, modified)
data$modification <- as.factor(data$modification)

pca = prcomp(data[,c], scale=scaling)
pdf("scree_plot.pdf")
plot(1:length(pca$sdev), pca$sdev, type="l", xlab="Index", 
     ylab="Eigenvalue", main="Scree Plot")
dev.off() 

pdf("pca_individuals.pdf")
plot(pca$x, col=data$modification, pch=19, xlab=paste("PC1 (", 100*round(pca$sdev[1]/sum(pca$sdev),2), "%)", sep=""),
     ylab=paste("PC2 (", 100*round(pca$sdev[2]/sum(pca$sdev),2), "%)", sep=""),) 
dev.off()

proj_ind <- as.matrix(data[,c(5,6,7,8,9,10,11,12,13,14,20,21,22,23,24)])%*%pca$rotation
proj_var <- cor(as.matrix(data[,c(5,6,7,8,9,10,11,12,13,14,20,21,22,23,24)]), proj_ind)

pdf("pca_variables.pdf")
ang = seq(from=0, to=2*pi, by = pi/1000)
Uang <-  cbind(cos(ang), sin(ang))
plot(Uang[,1], Uang[,2], type="l", 
     xlab=paste("PC1 (", 100*round(pca$sdev[1]/sum(pca$sdev),2), "%)", sep=""),
     ylab=paste("PC2 (", 100*round(pca$sdev[2]/sum(pca$sdev),2), "%)", sep=""),
     main="Variables Factor Map (PCA)", col="gray40")
abline(0, 0, lty=2, col="gray40")
abline(v=0, lty=2, col="gray40")
text(proj_var, labels=rownames(proj_var),cex= 0.7, pos=3)
for(i in 1:nrow(proj_var)){
  arrows(0,0,proj_var[i,1], proj_var[i,2], length = 0.1)
}
dev.off()



