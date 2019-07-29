#Loading packages
required.packages <- c("gridExtra","grDevices","KernSmooth")
new.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

suppressMessages(library(ggplot2))
library(optparse)
library(gridExtra)
library(grDevices)
options(warn=-1)


#KernSmooth is a dependency for denscols

parser <- OptionParser()
parser <- add_option(parser, opt_str=c("-i", "--input"), type="character",
                     dest='input',
                     help="Input Directory")
parser <- add_option(parser, opt_str=c("-n", "--names"), type="character",
                     dest='names',
                     help="Names")
parser <- add_option(parser, opt_str=c("-k", "--threshold"), type="double",
                     default=0.1,
                     dest='threshold',
                     help="threshold that will be used for removing rows with lower coverage than mean_coverage * this_threshold [0.1]")


options=parse_args(parser)
dir=options$input
names=options$names
threshold=options$threshold
names="1,2"
dir="~/Documentos/Uni/Internship/Project/Files/scripts/output"
threshold=0.1

n <- strsplit(names, ",")


filesmis <- list.files(path=dir, pattern = ".*mismatches")
filesstats <- list.files(path=dir, pattern = ".*STATS")


lmis <- list()
lstats <- list()
for (i in 1:length(filesmis)){
  assign(filesmis[i],read.delim(paste(dir, filesmis[i], sep="/")))
  assign(filesstats[i],read.delim(paste(dir, filesstats[i], sep="/")))
  lmis[[i]] <- get(filesmis[i])
  lstats[[i]] <- get(filesstats[i])
  lstats[[i]]$dataset <- n[[1]][i]
}

datasets <- list()
# We merge .mismatches and .STATS files 
for (i in 1:length(n[[1]])){
  assign(n[[1]][i], merge(lmis[[i]], 
                          lstats[[i]], 
                          by = c("chr", "pos", "ref_nuc"), 
                          sort = F))
  datasets[[i]] <- get(n[[1]][i])}


# We remove  positions with low coverage
remove_low_coverage <- function(i){
  cov <- mean(i$coverage) # We compute the coverage mean
  threshold <- cov*threshold # We create a threshold for removing values lower than 10% of the computed mean
  i <- i[i$coverage>threshold,]
}

datasets <- lapply(datasets, remove_low_coverage)

# We create a function for computing the total number of mismatches
mismatches_function <- function(k){ 
  mismatches <- c()  
  for (i in c(1:nrow(k))){   
    base <- k[i,3]
    a <- sum(k[i,4:7])-k[i,toString(base)]
    mismatches <- c(mismatches, a)
  }
  k <- cbind(k, mismatches)
}

datasets <- lapply(datasets, mismatches_function)

datasets <- lapply(datasets, transform, mismatches_freq = log(mismatches/coverage))
datasets <- lapply(datasets, transform, del_freq = log(del/coverage))


# We merge all the data
data <- merge(datasets[[1]], datasets[[2]], by=c("chr","pos", "ref_nuc"))




#quality
a <-cor(data$meanQual.x, data$meanQual.y, method = "pearson")
a <- paste("r=", round(a,3), sep="")
g1 <- ggplot(data, aes(meanQual.x, meanQual.y)) + geom_point(col=densCols(data$meanQual.x, data$meanQual.y, colramp = colorRampPalette(blues9[-(1:3)]))) + theme_bw() +
  theme(legend.position="none", plot.title = element_text(hjust = 0.5)) + 
  geom_text(aes(x=-Inf, y=+Inf, label=a), color="black", size=5, hjust = -0.5, vjust = 2) +
  xlab(paste("base quality", n[[1]][1])) + ylab(paste("base quality", n[[1]][2]))
m <- min(layer_scales(g1)$x$range$range[[1]], layer_scales(g1)$y$range$range[[1]])
M <- max(layer_scales(g1)$x$range$range[[2]], layer_scales(g1)$y$range$range[[2]])
g1 <- g1 + expand_limits(x=c(m,M), y=c(m,M))

#mismatch freq
datamis <- data[data$mismatches_freq.x!=-Inf & data$mismatches_freq.y!=-Inf,]
b <-cor(datamis$mismatches_freq.x, datamis$mismatches_freq.y, method = "pearson")
b <- paste("r=", round(b,3), sep="")
g2 <- ggplot(datamis, aes(mismatches_freq.x, mismatches_freq.y)) + geom_point(col=densCols(datamis$mismatches_freq.x, datamis$mismatches_freq.y, colramp = colorRampPalette(blues9[-(1:3)]))) + theme_bw() +
  theme(legend.position="none", plot.title = element_text(hjust = 0.5)) +
  geom_text(aes(x=-Inf, y=+Inf, label=b), color="black", size=5, hjust = -0.5, vjust = 2) +
  xlab(paste("log mismatch freq", n[[1]][1])) + ylab(paste("log mismatch freq", n[[1]][2]))
m <- min(layer_scales(g2)$x$range$range[[1]], layer_scales(g2)$y$range$range[[1]])
M <- max(layer_scales(g2)$x$range$range[[2]], layer_scales(g2)$y$range$range[[2]])
g2 <- g2 + expand_limits(x=c(m,M), y=c(m,M))

#deletion freq
datadel <- data[data$del_freq.x!=-Inf & data$del_freq.y!=-Inf,]
c <-cor(datadel$del_freq.x, datadel$del_freq.y, method = "pearson")
c <- paste("r=", round(c,3), sep="")
g3 <- ggplot(datadel, aes(del_freq.x, del_freq.y)) + geom_point(col=densCols(datadel$del_freq.x, datadel$del_freq.y, colramp = colorRampPalette(blues9[-(1:3)]))) + theme_bw() +
  theme(legend.position="none", plot.title = element_text(hjust = 0.5)) +
  geom_text(aes(x=-Inf, y=+Inf, label=c), color="black", size=5, hjust = -0.5, vjust = 2) +
  xlab(paste("log deletion freq", n[[1]][1])) + ylab(paste("log deletion freq", n[[1]][2]))
m <- min(layer_scales(g3)$x$range$range[[1]], layer_scales(g3)$y$range$range[[1]])
M <- max(layer_scales(g3)$x$range$range[[2]], layer_scales(g3)$y$range$range[[2]])
g3 <- g3 + expand_limits(x=c(m,M), y=c(m,M))




g <- grid.arrange(g1, g2, g3, nrow=1, ncol=3)

ggsave(paste("replicability", n[[1]][1], n[[1]][2], ".pdf", sep="_"), g, path = dir, device = "pdf", width = 10, height = 9)

