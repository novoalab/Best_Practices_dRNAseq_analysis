library(optparse)

## Parsing

parser <- OptionParser()
parser <- add_option(parser, opt_str=c("-b", "--basecaller"), type="character",
                     dest='basecaller',
                     help="basecaller")
parser <- add_option(parser, opt_str=c("-s", "--basecaller_short"), type="character",
                     dest='basecaller_short',
                     help="basecaller_short")
options=parse_args(parser)
basecaller=options$basecaller
basecaller_short=options$basecaller_short

#basecaller_short=GU_3.0.3
#basecaller=guppy3.0.3
#mapper=graphmap

UNM = read.delim(paste(basecaller_short, "_graphmap/UNM_C_per_site.var.csv.slided.onekmer.oneline.5mer.filtered.csv"), sep=",")
m5C = read.delim(paste(basecaller_short, "_graphmap/m5C_per_site.var.csv.slided.onekmer.oneline.5mer.filtered.csv"), sep=",")
hm5C =  read.delim(paste(basecaller_short, "_graphmap/5hmC_per_site.var.csv.slided.onekmer.oneline.5mer.filtered.csv"), sep=",")

UNM$modification <- "UNM"
m5C$modification <- "m5C"
hm5C$modification <- "hm5C"

data <- rbind(UNM, m5C, hm5C)

data$modification <- as.factor(data$modification)


#UNM vs m5C:

data1 <- rbind(UNM, m5C)
data1$modification <- as.factor(data1$modification)


pdf(paste("PCAs/C_", basecaller, "_graphmap/UNM_m5C_multiples.pdf"))
pca = prcomp(data1[,c(5,6,7,8,9,10,11,12,13,14,20,21,22,23,24)])
plot(pca$x, col=data1$modification) 
dev.off()



## UNM vs 5hmC

data2 <- rbind(UNM, hm5C)
data2$modification <- as.factor(data2$modification)

pca = prcomp(data2[,c(5,6,7,8,9,10,11,12,13,14,20,21,22,23,24)])
pdf(paste("PCAs/C_", basecaller, "_graphmap/UNM_5hmC_multiples.pdf"))
plot(pca$x, col=data2$modification) 
dev.off()
 


## UNM vs m5C vs 5hmC

data <- rbind(UNM, m5C, hm5C)
data$modification <- as.factor(data$modification)
data$modification <- factor(data$modification, levels=c("UNM","m5C","hm5C"))

pca = prcomp(data[,c(7,12,22)])  #q3, mis3, del3
pdf(paste("PCAs/C_", basecaller, "_graphmap/UNM_m5C_5hmC_singles.pdf"))
plot(pca$x, col=data$modification, pch=19)  # 1 is black, 2 is red, 3 is green. 1 is m5C, 0 is UNM
dev.off()
# q1,q2,q3,q4,q5, mis1,mis2,mis3,mis4,mis5, del1,del2,del3,del4,del5
pca = prcomp(data[,c(5,6,7,8,9,10,11,12,13,14,20,21,22,23,24)])
pdf(paste("PCAs/C_", basecaller,"_graphmap/UNM_m5C_5hmC_multiples_1.pdf"))
plot(pca$x, col=data$modification, pch = 19) #PC 1 vs PC 2   now 3 is gree is UNM, we should change this, 2 is red m5c, 1 is black 5hmC
dev.off()
pdf(paste("PCAs/C_", basecaller, "_graphmap/UNM_m5C_5hmC_multiples_2.pdf"))
plot(pca$x[,c(3,4)], col=data$modification, pch=19) #nothing separates m5C and 5hmC.
dev.off()
# the same but scaling
pca = prcomp(data[,c(5,6,7,8,9,10,11,12,13,14,20,21,22,23,24)], scale=TRUE)
plot(pca$x[,c(3,4)], col=data$modification) 



## m5C vs 5hmC

data3 <- rbind(m5C, hm5C)
data3$modification <- as.factor(data3$modification)

# q1,q2,q3,q4,q5, mis1,mis2,mis3,mis4,mis5, del1,del2,del3,del4,del5
pca = prcomp(data3[,c(5,6,7,8,9,10,11,12,13,14,20,21,22,23,24)])
pdf(paste("PCAs/C_", basecaller, "_graphmap/m5C_5hmC_multiples.pdf"))
plot(pca$x, col=data3$modification) #PC 1 vs PC 2   now 3 is gree is UNM, we should change this, 2 is red m5c, 1 is black 5hmC
dev.off()
plot(pca$x[,c(3,4)], col=data3$modification) #nothing separates m5C and 5hmC.






#m6A

UNM = read.delim(paste(basecaller_short, "_graphmap/UNM_A_per_site.var.csv.slided.onekmer.oneline.5mer.filtered.csv"), sep=",")
m6A = read.delim(paste(basecaller_short, "_graphmap/m6A_per_site.var.csv.slided.onekmer.oneline.5mer.filtered.csv"), sep=",")

UNM$modification <- "UNM"
m6A$modification <- "m5C"

data <- rbind(UNM, m6A)

data$modification <- as.factor(data$modification)


pca = prcomp(data[,c(7,12,22)])  #q3, mis3, del3
pdf(paste("PCAs/A_", basecaller, "_graphmap/UNM_m6A_singles.pdf"))
plot(pca$x, col=data$modification, pch=19)  # 1 is black, 2 is red, 3 is green. 1 is m5C, 0 is UNM
dev.off()

pca = prcomp(data[,c(7,12,22)], scale = TRUE)  #q3, mis3, del3
plot(pca$x, col=data$modification, pch=19)  # 1 is black, 2 is red, 3 is green. 1 is m5C, 0 is UNM

# q1,q2,q3,q4,q5, mis1,mis2,mis3,mis4,mis5, del1,del2,del3,del4,del5
pdf(paste("PCAs/A_", basecaller, "_graphmap/UNM_m6A_multiples.pdf"))
pca = prcomp(data[,c(5,6,7,8,9,10,11,12,13,14,20,21,22,23,24)])
plot(pca$x, col=data$modification, pch=19) 
dev.off()
# the same but scaling
pca = prcomp(data[,c(5,6,7,8,9,10,11,12,13,14,20,21,22,23,24)], scale=TRUE)
plot(pca$x, col=data$modification, pch=19) 

#with insertions, no difference
pca = prcomp(data[,c(5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24)]) #with insertions
pdf(paste("PCAs/A_", basecaller, "_graphmap/UNM_m6A_multiples_insertions.pdf"))
plot(pca$x, col=data$modification, pch=19)
dev.off()


#pU

UNM = read.delim(paste(basecaller_short, "_graphmap/UNM_T_per_site.var.csv.slided.onekmer.oneline.5mer.filtered.csv"), sep=",")
pU = read.delim(paste(basecaller_short, "_graphmap/pseudoU_per_site.var.csv.slided.onekmer.oneline.5mer.filtered.csv"), sep=",")
UNM$modification <- "UNM"
pU$modification <- "pU"
data <- rbind(UNM, pU)
data$modification <- as.factor(data$modification)

pdf(paste("PCAs/T_", basecaller, "_graphmap/UNM_pU_multiples.pdf"))
pca = prcomp(data[,c(5,6,7,8,9,10,11,12,13,14,20,21,22,23,24)])
plot(pca$x, col=data$modification, pch=19) 
dev.off()





