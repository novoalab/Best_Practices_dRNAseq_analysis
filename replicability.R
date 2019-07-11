library(ggplot2)
library(gridExtra)
library(grDevices)
#install.packages("KernSmooth") dependenci for denscols

name1 <- "hm5C"
name2 <- "C"
basecaller="AL_2.1.7" 
mapper="graphmap"
long_basecallername="albacore2.1.7"

dir_rep1 <- paste(mapper, name1, sep="/", end="")
dir_rep2 <- paste(mapper, "/", name1, "_rep2", sep="", end="/")
dir_rep3 <- paste(mapper, "/UNM_", name2, sep="", end="/")
dir_rep4 <- paste(mapper, "/UNM_", name2, "_rep2", sep="", end="/")


rep1_S <- list.files(path=dir_rep1, pattern = "albacore2.1.7.*STATS")
rep1_S <- read.delim(paste(dir_rep1, rep1_S, sep="/"))
rep1_M <- list.files(path=dir_rep1, pattern = "albacore2.1.7.*mismatches")
rep1_M <- read.delim(paste(dir_rep1, rep1_M, sep="/"))

rep2_S <- list.files(path=dir_rep2, pattern = "albacore_2.1.7.*STATS")
rep2_S <- read.delim(paste(dir_rep2, rep2_S, sep="/"))
rep2_M <- list.files(path=dir_rep2, pattern = "albacore_2.1.7.*mismatches")
rep2_M <- read.delim(paste(dir_rep2, rep2_M, sep="/"))

rep3_S <- list.files(path=dir_rep3, pattern = "albacore2.1.7.*STATS")
rep3_S <- read.delim(paste(dir_rep3, rep3_S, sep="/"))
rep3_M <- list.files(path=dir_rep3, pattern = "albacore2.1.7.*mismatches")
rep3_M <- read.delim(paste(dir_rep3, rep3_M, sep="/"))

rep4_S <- list.files(path=dir_rep4, pattern = "albacore_2.1.7.*STATS")
rep4_S <- read.delim(paste(dir_rep4, rep4_S, sep="/"))
rep4_M <- list.files(path=dir_rep4, pattern = "albacore_2.1.7.*mismatches")
rep4_M <- read.delim(paste(dir_rep4, rep4_M, sep="/"))


# We merge .mismatches and .STATS files 
names <- c("rep1", "rep2", "rep3", "rep4")
for (i in 1:length(names)){
  assign(names[i], merge(get(paste(names[i],"_S",sep="")), 
                         get(paste(names[i],"_M",sep="")), 
                         by = c("chr", "pos", "ref_nuc"), 
                         sort = F))}

# We remove  positions with low coverage

l <- list(rep1, rep2, rep3, rep4)
remove_low_coverage <- function(i){
  cov <- mean(i$coverage) # We compute the coverage mean
  threshold <- cov*0.1 # We create a threshold for removing values lower than 10% of the computed mean
  i <- i[i$coverage>threshold,]
}

l <- lapply(l, remove_low_coverage)

# We create a function for computing the total number of mismatches
mismatches_function <- function(k){ 
  mismatches <- c()  
  for (i in c(1:nrow(k))){   
    base <- k[i,3]
    a <- sum(k[i,10:13])-k[i,toString(base)]
    mismatches <- c(mismatches, a)
  }
  k <- cbind(k, mismatches)
}

l <- lapply(l, mismatches_function)

rep1 <- l[[1]]
rep2 <- l[[2]]
rep3 <- l[[3]]
rep4 <- l[[4]]
# We compute the log mismatch frequency
rep1$mismatches_freq <- log(rep1$mismatches/rep1$coverage)
rep2$mismatches_freq <- log(rep2$mismatches/rep2$coverage)
rep3$mismatches_freq <- log(rep3$mismatches/rep3$coverage)
rep4$mismatches_freq <- log(rep4$mismatches/rep4$coverage)
rep1$del_freq <- log(rep1$del/rep1$coverage)
rep2$del_freq <- log(rep2$del/rep2$coverage)
rep3$del_freq <- log(rep3$del/rep3$coverage)
rep4$del_freq <- log(rep4$del/rep4$coverage)


rep1$modification <- "rep1"
rep2$modification <- "rep2"
rep3$modification <- "rep3"
rep4$modification <- "rep4"

data <- merge(rep1,rep2, by=c("chr","pos", "ref_nuc"))
data <- data[data$ref_nuc==name2,]
data2 <- merge(rep3,rep4, by=c("chr","pos", "ref_nuc"))
data2 <- data2[data2$ref_nuc==name2,]


#########m6A
#quality
a <-cor(data$meanQual.x, data$meanQual.y, method = "pearson")
a <- paste("r=", round(a,3), sep="")
g1 <- ggplot(data, aes(meanQual.x, meanQual.y)) + geom_point(col=densCols(data$meanQual.x, data$meanQual.y, colramp = colorRampPalette(blues9[-(1:3)]))) + theme_bw() +
  theme(legend.position="none", plot.title = element_text(hjust = 0.5)) + 
  geom_text(aes(x=-Inf, y=+Inf, label=a), color="black", size=5, hjust = -0.5, vjust = 2) +
  ggtitle(name1) + xlab("base quality (rep1)") + ylab("base quality (rep2)")
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
  ggtitle(name1) + xlab("log mismatch freq (rep1)") + ylab("log mismatch freq (rep2)")
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
  ggtitle(name1) + xlab("log deletion freq (rep1)") + ylab("log deletion freq (rep2)")
m <- min(layer_scales(g3)$x$range$range[[1]], layer_scales(g3)$y$range$range[[1]])
M <- max(layer_scales(g3)$x$range$range[[2]], layer_scales(g3)$y$range$range[[2]])
g3 <- g3 + expand_limits(x=c(m,M), y=c(m,M))


############UNM
#quality
d <-cor(data2$meanQual.x, data2$meanQual.y, method = "pearson")
d <- paste("r=", round(d,3), sep="")
g4 <- ggplot(data2, aes(meanQual.x, meanQual.y)) + geom_point(col=densCols(data2$meanQual.x, data2$meanQual.y, colramp = colorRampPalette(blues9[-(1:3)]))) + theme_bw() +
  theme(legend.position="none", plot.title = element_text(hjust = 0.5)) + 
  geom_text(aes(x=-Inf, y=+Inf, label=d), color="black", size=5, hjust = -0.5, vjust = 2) +
  ggtitle(name2) + xlab("base quality (rep1)") + ylab("base quality (rep2)")
m <- min(layer_scales(g4)$x$range$range[[1]], layer_scales(g4)$y$range$range[[1]])
M <- max(layer_scales(g4)$x$range$range[[2]], layer_scales(g4)$y$range$range[[2]])
g4 <- g4 + expand_limits(x=c(m,M), y=c(m,M))

#mismatch freq
datamis <- data2[data2$mismatches_freq.x!=-Inf & data2$mismatches_freq.y!=-Inf,]
e <-cor(datamis$mismatches_freq.x, datamis$mismatches_freq.y, method = "pearson")
e <- paste("r=", round(e,3), sep="")
g5 <- ggplot(datamis, aes(mismatches_freq.x, mismatches_freq.y)) + geom_point(col=densCols(datamis$mismatches_freq.x, datamis$mismatches_freq.y, colramp = colorRampPalette(blues9[-(1:3)]))) + theme_bw() +
  theme(legend.position="none", plot.title = element_text(hjust = 0.5)) +
  geom_text(aes(x=-Inf, y=+Inf, label=e), color="black", size=5, hjust = -0.5, vjust = 2) +
  ggtitle(name2) + xlab("log mismatch freq (rep1)") + ylab("log mismatch freq (rep2)")
m <- min(layer_scales(g5)$x$range$range[[1]], layer_scales(g5)$y$range$range[[1]])
M <- max(layer_scales(g5)$x$range$range[[2]], layer_scales(g5)$y$range$range[[2]])
g5 <- g5 + expand_limits(x=c(m,M), y=c(m,M))

#deletion freq
datadel <- data2[data2$del_freq.x!=-Inf & data2$del_freq.y!=-Inf,]
f <-cor(datadel$del_freq.x, datadel$del_freq.y, method = "pearson")
f <- paste("r=", round(f,3), sep="")
g6 <- ggplot(datadel, aes(del_freq.x, del_freq.y)) + geom_point(col=densCols(datadel$del_freq.x, datadel$del_freq.y, colramp = colorRampPalette(blues9[-(1:3)]))) + theme_bw() +
  theme(legend.position="none", plot.title = element_text(hjust = 0.5)) +
  geom_text(aes(x=-Inf, y=+Inf, label=f), color="black", size=5, hjust = -0.5, vjust = 2) +
  ggtitle(name2) + xlab("log deletion freq (rep1)") + ylab("log deletion freq (rep2)")
m <- min(layer_scales(g6)$x$range$range[[1]], layer_scales(g6)$y$range$range[[1]])
M <- max(layer_scales(g6)$x$range$range[[2]], layer_scales(g6)$y$range$range[[2]])
g6 <- g6 + expand_limits(x=c(m,M), y=c(m,M))

######## UNM vs m6A

data3 <- merge(rep1,rep3, by=c("chr","pos", "ref_nuc"))
data3 <- data3[data3$ref_nuc==name2,]
data4 <- merge(rep2,rep4, by=c("chr","pos", "ref_nuc"))
data4 <- data4[data4$ref_nuc==name2,]

#quality
g <-cor(data3$meanQual.x, data3$meanQual.y, method = "pearson")
g <- paste("r=", round(g,3), sep="")
g7 <- ggplot(data3, aes(meanQual.x, meanQual.y)) + geom_point(col=densCols(data3$meanQual.x, data3$meanQual.y, colramp = colorRampPalette(blues9[-(1:3)]))) +
  theme_bw() +
  theme(legend.position="none", plot.title = element_text(hjust = 0.5)) + 
  geom_text(aes(x=-Inf, y=+Inf, label=g), color="black", size=5, hjust = -0.5, vjust = 2) +
  ggtitle(paste(name2, "vs", name1)) + xlab(paste("base quality (", name1, " rep1)", sep="")) + ylab(paste("base quality (", name2, " rep1)", sep=""))
m <- min(layer_scales(g7)$x$range$range[[1]], layer_scales(g7)$y$range$range[[1]])
M <- max(layer_scales(g7)$x$range$range[[2]], layer_scales(g7)$y$range$range[[2]])
g7 <- g7 + expand_limits(x=c(m,M), y=c(m,M))


#mismatch freq
datamis <- data3[data3$mismatches_freq.x!=-Inf & data3$mismatches_freq.y!=-Inf,]
h <-cor(datamis$mismatches_freq.x, datamis$mismatches_freq.y, method = "pearson")
h <- paste("r=", round(h,3), sep="")
g8 <- ggplot(datamis, aes(mismatches_freq.x, mismatches_freq.y)) + geom_point(col=densCols(datamis$mismatches_freq.x, datamis$mismatches_freq.y, colramp = colorRampPalette(blues9[-(1:3)]))) + theme_bw() +
  theme(legend.position="none", plot.title = element_text(hjust = 0.5)) +
  geom_text(aes(x=-Inf, y=+Inf, label=h), color="black", size=5, hjust = -0.5, vjust = 2) +
  ggtitle(paste(name2, "vs", name1)) + xlab(paste("log mismatch freq (", name1, " rep1)", sep="")) + ylab(paste("log mismatch freq (", name2, " rep1)", sep=""))
m <- min(layer_scales(g8)$x$range$range[[1]], layer_scales(g8)$y$range$range[[1]])
M <- max(layer_scales(g8)$x$range$range[[2]], layer_scales(g8)$y$range$range[[2]])
g8 <- g8 + expand_limits(x=c(m,M), y=c(m,M))

#deletion freq
datadel <- data3[data3$del_freq.x!=-Inf & data3$del_freq.y!=-Inf,]
i <-cor(datadel$del_freq.x, datadel$del_freq.y, method = "pearson")
i <- paste("r=", round(i,3), sep="")
g9 <- ggplot(datadel, aes(del_freq.x, del_freq.y)) + geom_point(col=densCols(datadel$del_freq.x, datadel$del_freq.y, colramp = colorRampPalette(blues9[-(1:3)]))) + theme_bw() +  theme(legend.position="none", plot.title = element_text(hjust = 0.5)) +
  geom_text(aes(x=-Inf, y=+Inf, label=i), color="black", size=5, hjust = -0.5, vjust = 2) +
  ggtitle(paste(name2, "vs", name1)) + xlab(paste("log deletion freq (", name1, " rep1)", sep="")) + ylab(paste("log deletion freq (", name2, " rep1)", sep=""))
m <- min(layer_scales(g9)$x$range$range[[1]], layer_scales(g9)$y$range$range[[1]])
M <- max(layer_scales(g9)$x$range$range[[2]], layer_scales(g9)$y$range$range[[2]])
g9 <- g9 + expand_limits(x=c(m,M), y=c(m,M))

#quality
j <-cor(data4$meanQual.x, data4$meanQual.y, method = "pearson")
j <- paste("r=", round(j,3), sep="")
g10 <- ggplot(data4, aes(meanQual.x, meanQual.y)) + geom_point(col=densCols(data4$meanQual.x, data4$meanQual.y, colramp = colorRampPalette(blues9[-(1:3)]))) + theme_bw() +
  theme(legend.position="none", plot.title = element_text(hjust = 0.5)) + 
  geom_text(aes(x=-Inf, y=+Inf, label=j), color="black", size=5, hjust = -0.5, vjust = 2) +
  ggtitle(paste(name2, "vs", name1)) + xlab(paste("base quality (", name1, " rep2)", sep="")) + ylab(paste("base quality (", name2, " rep2)", sep=""))
m <- min(layer_scales(g10)$x$range$range[[1]], layer_scales(g10)$y$range$range[[1]])
M <- max(layer_scales(g10)$x$range$range[[2]], layer_scales(g10)$y$range$range[[2]])

g10 <- g10 + expand_limits(x=c(m,M), y=c(m,M))

#mismatch freq
datamis <- data4[data4$mismatches_freq.x!=-Inf & data4$mismatches_freq.y!=-Inf,]
k <-cor(datamis$mismatches_freq.x, datamis$mismatches_freq.y, method = "pearson")
k <- paste("r=", round(k,3), sep="")
g11 <- ggplot(datamis, aes(mismatches_freq.x, mismatches_freq.y)) + geom_point(col=densCols(datamis$mismatches_freq.x, datamis$mismatches_freq.y, colramp = colorRampPalette(blues9[-(1:3)]))) + theme_bw() +
  theme(legend.position="none", plot.title = element_text(hjust = 0.5)) +
  geom_text(aes(x=-Inf, y=+Inf, label=k), color="black", size=5, hjust = -0.5, vjust = 2) +
  ggtitle(paste(name2, "vs", name1)) + xlab(paste("log mismatch freq (", name1, " rep2)", sep="")) + ylab(paste("log mismatch freq (", name2, " rep2)", sep=""))
m <- min(layer_scales(g11)$x$range$range[[1]], layer_scales(g11)$y$range$range[[1]])
M <- max(layer_scales(g11)$x$range$range[[2]], layer_scales(g11)$y$range$range[[2]])
g11 <- g11 + expand_limits(x=c(m,M), y=c(m,M))

#deletion freq
datadel <- data4[data4$del_freq.x!=-Inf & data4$del_freq.y!=-Inf,]
l <-cor(datadel$del_freq.x, datadel$del_freq.y, method = "pearson")
l <- paste("r=", round(l,3), sep="")
g12 <- ggplot(datadel, aes(del_freq.x, del_freq.y)) + geom_point(col=densCols(datadel$del_freq.x, datadel$del_freq.y, colramp = colorRampPalette(blues9[-(1:3)]))) + theme_bw() +
  theme(legend.position="none", plot.title = element_text(hjust = 0.5)) +
  geom_text(aes(x=-Inf, y=+Inf, label=l), color="black", size=5, hjust = -0.5, vjust = 2) +
  ggtitle(paste(name2, "vs", name1)) + xlab(paste("log deletion freq (", name1, " rep2)", sep="")) + ylab(paste("log deletion freq (", name2, " rep2)", sep=""))
m <- min(layer_scales(g12)$x$range$range[[1]], layer_scales(g12)$y$range$range[[1]])
M <- max(layer_scales(g12)$x$range$range[[2]], layer_scales(g12)$y$range$range[[2]])
g12 <- g12 + expand_limits(x=c(m,M), y=c(m,M))

g <- grid.arrange(g4, g1, g7, g10, g5, g2, g8, g11, g6, g3, g9, g12, nrow=3, ncol=4)

ggsave(paste("replicability/replicability", name2, name1, basecaller, mapper, ".pdf", sep="_"), g, device = "pdf", width = 10, height = 9)

