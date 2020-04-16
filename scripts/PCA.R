
# Loading packages
library(optparse)
library(FactoMineR)
library(factoextra)

#Parsing
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


options=parse_args(parser)
unmodified=options$unmodified
modified=options$modified
names=options$names
columns=options$columns


#Pre-processing

n <- strsplit(names, ",")[[1]]
c <- as.numeric(strsplit(columns, ",")[[1]])

unmodified <- read.delim(unmodified, sep = ",")
modified <- read.delim(modified, sep = ",")

unmodified$modification <- "UNM"
modified$modification <- "MOD"

data <- rbind(unmodified, modified)
data$modification <- as.factor(data$modification)

#PCA
pca_f <- PCA(data[,c])

pdf(paste("PCA/", names, "_scree_plot.pdf", sep=""))
plot(pca_f$eig[,1], type="l", ylab="Eigenvalue", main="Scree Plot")
dev.off() 



pdf(paste("PCA/", names, "_pca_individuals.pdf", sep=""))
plot(pca_f$ind$coord, pch=19, col=data$modification, xlab=paste("PC1 (", round(pca_f$eig[1,2],2), "%)", sep=""),
     ylab=paste("PC2 (", round(pca_f$eig[2,2],2), "%)", sep=""))
dev.off()


pdf(paste("PCA/", names, "_pca_variables.pdf", sep=""))
fviz_pca_var(pca_f, xlab=paste("PC1 (", round(pca_f$eig[1,2],2), "%)", sep=""),
             ylab=paste("PC2 (", round(pca_f$eig[2,2],2), "%)", sep=""))
dev.off()

