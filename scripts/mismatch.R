#####################################################################################
#This script takes as input a directory containing files finishing with .STATS and 
# and with .mismatches (output of previous step) and it compares the mismatch frequency 
# per base-caller/dataset and nucleotide, the mismatch pattern, the bases more modified 
# than a given frequency and it outputs ternary diagrams per base-caller/dataset.
#The mismatch pattern and the ternary diagrams are computed using the k-mers 
#that contain -m in its central position.  
#The parameters needed are an input directory -i, a string with the names of each dataset
#separated by commas -n, an optional title name -t that will be added in the plots,
#a base that is being modificied to focus on it -m, an option -e that studies the input
#as epinano output if it is present (taken k-mers with only -m base in its central position),
#a threshold k for removing positions with lower coverage than this percentage of the mean
#coverage (default=0.1), and a zoom -z for plotting the mismatch pattern (default=0.8)
######################################################################################


#implement it in my .sh correctly, add to github, finish next R scripts that are going to be in the .sh too

## Loading libraries

required.packages <- c("ggplot2","ggExtra","optparse","reshape2","ggtern","lattice","latticeExtra")
new.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

suppressMessages(library(ggplot2))
library(ggpubr)
library(optparse)
library(reshape2)
suppressMessages(library(ggtern))
library(lattice)
suppressMessages(library(latticeExtra))
options(warn=-1)

## Parsing

parser <- OptionParser()
parser <- add_option(parser, opt_str=c("-i", "--input"), type="character",
                     dest='input',
                     help="Input Directory")
parser <- add_option(parser, opt_str=c("-n", "--names"), type="character",
                     dest='names',
                     help="Names")
parser <- add_option(parser, opt_str=c("-t", "--title"), type="character",
                     default="",
                     dest='title',
                     help="Title")
parser <- add_option(parser, opt_str=c("-m", "--modification"), type="character",
                     default="",
                     dest='mod',
                     help="mod")
parser <- add_option(parser, opt_str=c("-e", "--epinano"), type="character",
                     action="store_true",
                     dest='epinano',
                     help="if -e is present, everything is computed as epinano output")
parser <- add_option(parser, opt_str=c("-k", "--threshold"), type="double",
                     default=0.1,
                     dest='threshold',
                     help="threshold that will be used for removing rows with lower coverage than mean_coverage * this_threshold [0.1]")
parser <- add_option(parser, opt_str=c("-z", "--zoom"), type="double",
                     default=0.8,
                     dest='zoom',
                     help="lower limit in y axis for mismatch_pattern_100.pdf [0.8]")

options=parse_args(parser)
dir=options$input
names=options$names
title=options$title
threshold=options$threshold
epinano=options$epinano
mod=options$mod
zoom=options$zoom
#dir="~/Documentos/Uni/Internship/Project/Files/scripts/output"
#names="minimap2,graphmap"
#epinano=TRUE
#title=""
#mod="A"
#threshold=0.1
#zoom=0.8
n <- strsplit(names, ",")

# In the input directory there should be a file ending in .mismatches and another one ending in .STATS per dataset/base-caller.
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

datasets <- lapply(datasets, transform, mismatches_freq = mismatches/coverage)


# We merge all the data
data <- c()
for (i in 1:length(n[[1]])){
  assign(n[[1]][i], datasets[[i]])
  data <- rbind(data, datasets[[i]])
}


## Plotting

#The first plot compares the global mismatch frequency per dataset/base-caller
#colors=c("cyan3","cyan4", "red", "red3")

g <- ggplot(data, aes(dataset, log(mismatches_freq), color=dataset, fill=dataset)) + 
  geom_violin(alpha = 0.3) +
  geom_boxplot(alpha = 0.4, width=0.2, show.legend = FALSE) +
 # scale_color_manual(values = colors) + 
#  scale_fill_manual(values = colors) +
  ylab("Mismatch Frequency (log)") + theme_bw() + ggtitle(paste("Mismatch Frequency per dataset:", title)) + 
  theme(plot.title = element_text(size = 10, face = "bold"))

ggsave("mismatch_freq_per_basecaller.pdf", g, device = "pdf", path = dir, width = 4, height = 5)

#The next set of plots evaluate the number of bases with mismatch frequency higher than a given threshold. First with the total number and then with the proportion
plotting_threshold <- function(threshold){
  ggplot(data[data$mismatches_freq>threshold,], aes(dataset, fill=ref_nuc)) + 
    geom_bar()  + labs(fill = "Reference Nucleotide") + theme_bw() + ggtitle(paste("Threshold:", threshold*100, "%")) +
    geom_text(stat="count", aes(label=..count..), position = "stack", vjust=1.5, size=3.5) + 
    #scale_x_discrete(labels= c("AL\n2.1.7","AL\n2.3.4","GU\n2.3.1","GU\n3.0.3")) +
    theme(plot.title = element_text(hjust = 0.5)) + xlab("")
}

plotting_threshold_proportion <- function(threshold){
  ggplot(data[data$mismatches_freq>threshold,], aes(dataset, fill=ref_nuc)) + 
    geom_bar(position = "fill")  + labs(fill = "Reference Nucleotide") + theme_bw() + ggtitle(paste("Threshold:", threshold*100, "%")) +
    #scale_x_discrete(labels= c("AL\n2.1.7","AL\n2.3.4","GU\n2.3.1","GU\n3.0.3")) +
    theme(plot.title = element_text(hjust = 0.5)) + xlab("")
}


g1 <- plotting_threshold(0.05)
g2 <- plotting_threshold(0.1)
g3 <- plotting_threshold(0.2)
g4 <- plotting_threshold(0.5)

l <- list(g1, g2, g3, g4)
figure <- ggarrange(plotlist=l, ncol=4, nrow=1, common.legend = TRUE, legend = "bottom") 
f <- annotate_figure(figure, top = text_grob(paste("Number of Bases with Mismatch Frequency Higher than a Threhold:", title), face = "bold"))

ggsave("bases_above_threshold.pdf", f, device = "pdf", path = dir, width = 9, height = 7)

g1 <- plotting_threshold_proportion(0.05)
g2 <- plotting_threshold_proportion(0.1)
g3 <- plotting_threshold_proportion(0.2)
g4 <- plotting_threshold_proportion(0.5)

l <- list(g1, g2, g3, g4)
figure <- ggarrange(plotlist=l, ncol=4, nrow=1, common.legend = TRUE, legend = "bottom") 
f <- annotate_figure(figure, top = text_grob(paste("Proportion of Bases with Mismatch Frequency Higher than a Threhold:", title), face = "bold"))

ggsave("bases_above_threshold_proportion.pdf", f, device = "pdf", path = dir, width = 9, height = 7)


# These plots now compute the mismatch frequencies per nucleotide

plotting_mismatch_frequencies <- function(i, n){
  ggplot(data = i, mapping = aes(ref_nuc, log(mismatches_freq), color=ref_nuc, fill=ref_nuc)) +
    geom_violin(alpha = 0.3) +
    geom_boxplot(alpha = 0.5, show.legend = FALSE, width = 0.2) +
    ggtitle(n) + ylab("Mismatch Frequency (log)") + xlab("Nucleotide") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
}

l <- list()
for (i in 1:length(n[[1]])){
  l[[n[[1]][i]]] <- plotting_mismatch_frequencies(get(n[[1]][i]), n[[1]][i])
}

figure <- ggarrange(plotlist=l, ncol=length(n[[1]]), nrow=1, common.legend = TRUE, legend = "bottom") 
f <- annotate_figure(figure, top = text_grob(paste("Mismatch Frequency per Nucleotide:", title), face = "bold"))

ggsave("mismatch_per_nucleotide.pdf", f, device = "pdf", path = dir, width = 9, height = 5)

#levelplots

A <- data[data$ref_nuc=="A",]
C <- data[data$ref_nuc=="C",]
G <- data[data$ref_nuc=="G",]
T <- data[data$ref_nuc=="T",]

m <- matrix(0,4, length(n[[1]]))
for (i in 1:length(n[[1]])){
  a <- log(A[A$dataset==n[[1]][i],]$mismatches_freq)
  a <- a[a!=-Inf]
  m[1,i]<- median(a)
  c <- log(C[C$dataset==n[[1]][i],]$mismatches_freq)
  c <- c[c!=-Inf]
  m[2,i]<- median(c)
  g <- log(G[G$dataset==n[[1]][i],]$mismatches_freq)
  g <- g[g!=-Inf]
  m[3,i]<- median(g)
  t <- log(T[T$dataset==n[[1]][i],]$mismatches_freq)
  t <- t[t!=-Inf]
  m[4,i]<- median(t)
}

colnames(m) <- n[[1]]
rownames(m) <- c("A","C","G","T")

pdf(paste(dir, "levelplot.pdf", sep="/"))
levelplot(m, xlab="", ylab="", margin=FALSE, 
          col.regions = heat.colors(100)[1:length(heat.colors(90))], 
          scales = list(tck = c(0,0)), main="Levelplot of Mismatch Frequencies")
dev.off()

##########
# Now we compute the mismatch signature per dataset/base-caller using 5-mers centered on each nucleotide that doesn't contain mod in positions -2,-1,1,2. This is thought to be used in EpiNano outputs

col <- list("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF")
names(col) <- c("A", "C", "G", "T")

if (epinano==TRUE){
A <- c()
for (i in 3:(nrow(data)-2)){
  if(data[i,]$ref_nuc=="A" && data[i-1,]$ref_nuc!=mod && data[i-2,]$ref_nuc!=mod && data[i+1,]$ref_nuc!=mod && data[i+2,]$ref_nuc!=mod){
    A <- rbind(A, data[i,])}}
C <- c()
for (i in 3:(nrow(data)-2)){
  if(data[i,]$ref_nuc=="C" && data[i-1,]$ref_nuc!=mod && data[i-2,]$ref_nuc!=mod && data[i+1,]$ref_nuc!=mod && data[i+2,]$ref_nuc!=mod){
    C <- rbind(C, data[i,])}}
G <- c()
for (i in 3:(nrow(data)-2)){
  if(data[i,]$ref_nuc=="G" && data[i-1,]$ref_nuc!=mod && data[i-2,]$ref_nuc!=mod && data[i+1,]$ref_nuc!=mod && data[i+2,]$ref_nuc!=mod){
    G <- rbind(G, data[i,])}}
T <- c()
for (i in 3:(nrow(data)-2)){
  if(data[i,]$ref_nuc=="T" && data[i-1,]$ref_nuc!=mod && data[i-2,]$ref_nuc!=mod && data[i+1,]$ref_nuc!=mod && data[i+2,]$ref_nuc!=mod){
    T <- rbind(T, data[i,])}}
} else {
  A <- data[,ref_nuc=="A"]
  C <- data[,ref_nuc=="C"]
  G <- data[,ref_nuc=="G"]
  T <- data[,ref_nuc=="T"]
}

mismatch_pattern <- function(input, number){
  
  dat <- A[A$dataset==input,]
  dat <- suppressMessages({melt(dat[,c("ref_nuc", "A", "T", "C", "G")])})
  dat <- tapply(dat$value, list(dat$ref_nuc,dat$variable), sum)[1,] 
  datA <- as.data.frame(cbind(c("A","A","A","A"), c("A","T","C","G"), dat))
  colnames(datA) <- c("ref_nuc", "variable","value")
  datA$variable <- factor(datA$variable, levels = as.character(datA[order(as.numeric(as.character(datA$value))),]$variable))
  datA$value <- as.numeric(as.character(datA$value))
  
  dat <- C[C$dataset==input,]
  dat <- suppressMessages({melt(dat[,c("ref_nuc", "A", "T", "C", "G")])})
  dat <- tapply(dat$value, list(dat$ref_nuc,dat$variable), sum)[2,]
  datC <- as.data.frame(cbind(c("C","C","C","C"), c("A","T","C","G"), dat))
  colnames(datC) <- c("ref_nuc", "variable","value")
  datC$variable <- factor(datC$variable, levels = as.character(datC[order(as.numeric(as.character(datC$value))),]$variable))
  datC$value <- as.numeric(as.character(datC$value))
  
  dat <- G[G$dataset==input,]
  dat <- suppressMessages({melt(dat[,c("ref_nuc", "A", "T", "C", "G")])})
  dat <- tapply(dat$value, list(dat$ref_nuc,dat$variable), sum)[3,]
  datG <- as.data.frame(cbind(c("G","G","G","G"), c("A","T","C","G"), dat))
  colnames(datG) <- c("ref_nuc", "variable","value")
  datG$variable <- factor(datG$variable, levels = as.character(datG[order(as.numeric(as.character(datG$value))),]$variable))
  datG$value <- as.numeric(as.character(datG$value))
  
  dat <- T[T$dataset==input,]
  dat <- suppressMessages({melt(dat[,c("ref_nuc", "A", "T", "C", "G")])})
  dat <- tapply(dat$value, list(dat$ref_nuc,dat$variable), sum)[4,]
  datT <- as.data.frame(cbind(c("T","T","T","T"), c("A","T","C","G"), dat))
  colnames(datT) <- c("ref_nuc", "variable","value")
  datT$variable <- factor(datT$variable, levels = as.character(datT[order(as.numeric(as.character(datT$value))),]$variable))
  datT$value <- as.numeric(as.character(datT$value))
  
  first <- as.character(datA[order(datA$value),]$variable)[1]
  second <- as.character(datA[order(datA$value),]$variable)[2]
  third <- as.character(datA[order(datA$value),]$variable)[3]
  fourth <- as.character(datA[order(datA$value),]$variable)[4]
  
  ggplot() + 
    geom_bar(data=datA, mapping=aes(ref_nuc, value, fill=variable), stat="identity", position="fill") +
    geom_bar(data=datC, mapping=aes(ref_nuc, value, fill=variable), stat="identity", position="fill") +
    geom_bar(data=datG, mapping=aes(ref_nuc, value, fill=variable), stat="identity", position="fill") +
    geom_bar(data=datT, mapping=aes(ref_nuc, value, fill=variable), stat="identity", position="fill") +
    scale_fill_manual("Base-called Nucleotide", 
                      breaks = c(first, second, third, fourth),
                      values = c(col[[first]], col[[second]], col[[third]], col[[fourth]]),
                      guide = guide_legend(reverse = TRUE)) +
    xlab("Reference Nucleotide") + ylab("count") +
    coord_cartesian(ylim=c(number,1)) +
    theme_bw() +
    ggtitle(input) +
    theme(plot.title = element_text(size = 12, hjust = 0.5), axis.text=element_text(size=6), axis.title=element_text(size=9))
}

l <- list()
for (i in 1:length(n[[1]])){
  l[[n[[1]][i]]] <- mismatch_pattern(n[[1]][i], zoom)
}


figure <- ggarrange(plotlist=l, ncol=length(n[[1]]), nrow=1, common.legend = TRUE, legend = "bottom") 
f <- annotate_figure(figure, top = text_grob(paste("Mismatch Frequency per Nucleotide:", title), face = "bold"))

ggsave("mismatch_pattern_zoom.pdf", f, device = "pdf", path = dir, width = 8, height = 4)


l <- list()
for (i in 1:length(n[[1]])){
  l[[n[[1]][i]]] <- mismatch_pattern(n[[1]][i], 0)
}


figure <- ggarrange(plotlist=l, ncol=length(n[[1]]), nrow=1, common.legend = TRUE, legend = "bottom") 
f <- annotate_figure(figure, top = text_grob(paste("Mismatch Frequency per Nucleotide:", title), face = "bold"))

ggsave("mismatch_pattern_100.pdf", f, device = "pdf", path = dir, width = 8, height = 4)

# Finally we compute the ternary diagrams for plotting the mismatch signature per nucleotide

  
ternaries <- function(input){
  datA <- A[A$dataset==input,]
g1 <- ggtern(datA, aes(C,G,T)) +
  geom_point(size = 0.05, aes(color = log(coverage))) +
  ggtitle("Ternary Diagram for A") + theme_bw(base_size = 6) +
  scale_color_continuous("Coverage (log)") +
  theme(plot.title = element_text(size = 8, hjust = 0.5)) 
legend <- get_legend(g1)
g1 <- g1 + theme(legend.position="none")

datC <- C[C$dataset==input,]
g2 <- ggtern(datC, aes(A,G,T)) +
  geom_point(size = 0.05, aes(color = log(coverage))) +
  ggtitle("Ternary Diagram for C") + theme_bw(base_size = 6) +
  scale_color_continuous("Coverage (log)") +
  theme(legend.position="none", plot.title = element_text(size = 8, hjust = 0.5))

datG <- G[G$dataset==input,]
g3 <- ggtern(datG, aes(A,C,T)) +
  geom_point(size = 0.05, aes(color = log(coverage))) +
  ggtitle("Ternary Diagram for G") + theme_bw(base_size = 6) +
  scale_color_continuous("Coverage (log)") +
  theme(legend.position="none", plot.title = element_text(size = 8, hjust = 0.5))

datT <- T[T$dataset==input,]
g4 <- ggtern(datT, aes(A,C,G)) +
  geom_point(size = 0.05, aes(color = log(coverage))) +
  ggtitle("Ternary Diagram for T") + theme_bw(base_size = 6) +
  scale_color_continuous("Coverage (log)") +
  theme(legend.position="none", plot.title = element_text(size = 8, hjust = 0.5))

g<- grid.arrange(g1,g2,g3,g4, nrow=1, ncol=5, legend, top=paste("Ternary Diagrams in", input))
ggsave(paste("ternary_diagram_", input, ".pdf", sep=""), g, device = "pdf", path = dir, width = 7, height = 3)
}

for (i in n[[1]]){
  ternaries(i)
}
