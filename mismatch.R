## Loading libraries
library(ggplot2)
library(ggpubr)
library(optparse)
library(reshape2)
library(ggtern)



## Parsing

parser <- OptionParser()
parser <- add_option(parser, opt_str=c("-i", "--input"), type="character",
                     dest='input',
                     help="Input Directory")
parser <- add_option(parser, opt_str=c("-n", "--name"), type="character",
                     dest='name',
                     help="Name")
options=parse_args(parser)
dir=options$input
name=options$name

#dir="RNAAB091808/"
#name="pU"

# In the input directory there should be a file ending in .mismatches and another one ending in .STATS per base-caller.
AL_2.1.7_S <- list.files(path=dir, pattern = "albacore2.1.7.sorted.bam.STATS")
AL_2.1.7_M <- list.files(path=dir, pattern = "albacore2.1.7.sorted.bam.mismatches")

AL_2.3.4_S <- list.files(path=dir, pattern = "albacore2.3.4.sorted.bam.STATS")
AL_2.3.4_M <- list.files(path=dir, pattern = "albacore2.3.4.sorted.bam.mismatches")

GU_2.3.1_S <- list.files(path=dir, pattern = "guppy2.3.1.sorted.bam.STATS")
GU_2.3.1_M <- list.files(path=dir, pattern = "guppy2.3.1.sorted.bam.mismatches")

GU_3.0.3_S <- list.files(path=dir, pattern = "guppy3.0.3.sorted.bam.STATS")
GU_3.0.3_M <- list.files(path=dir, pattern = "guppy3.0.3.sorted.bam.mismatches")

names <- c("AL_2.1.7_S", "AL_2.1.7_M", "AL_2.3.4_S", "AL_2.3.4_M", "GU_2.3.1_S", "GU_2.3.1_M", "GU_3.0.3_S", "GU_3.0.3_M")
l <- list(AL_2.1.7_S, AL_2.1.7_M, AL_2.3.4_S, AL_2.3.4_M, GU_2.3.1_S, GU_2.3.1_M, GU_3.0.3_S, GU_3.0.3_M)

for (i in 1:length(l)){
  assign(names[i],read.delim(paste(dir, l[[i]], sep="/")))}

# We merge .mismatches and .STATS files 
names <- c("AL_2.1.7", "AL_2.3.4", "GU_2.3.1", "GU_3.0.3")
for (i in 1:length(names)){
  assign(names[i], merge(get(paste(names[i],"_S",sep="")), 
                         get(paste(names[i],"_M",sep="")), 
                         by = c("chr", "pos", "ref_nuc"), 
                         sort = F))}


# We remove  positions with low coverage

l <- list(AL_2.1.7, AL_2.3.4, GU_2.3.1, GU_3.0.3)
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

AL_2.1.7 <- l[[1]]
AL_2.3.4 <- l[[2]]
GU_2.3.1 <- l[[3]]
GU_3.0.3 <- l[[4]]

# We compute the mismatch frequency
AL_2.1.7$mismatches_freq <- AL_2.1.7$mismatches/AL_2.1.7$coverage
AL_2.3.4$mismatches_freq <- AL_2.3.4$mismatches/AL_2.3.4$coverage
GU_2.3.1$mismatches_freq <- GU_2.3.1$mismatches/GU_2.3.1$coverage
GU_3.0.3$mismatches_freq <- GU_3.0.3$mismatches/GU_3.0.3$coverage

# We add a variable to identify the base-caller
AL_2.1.7$Basecaller <- "AL 2.1.7"
AL_2.3.4$Basecaller <- "AL 2.3.4"
GU_2.3.1$Basecaller <- "GU 2.3.1"
GU_3.0.3$Basecaller <- "GU 3.0.3"

# We merge the data from all the base-callers
data <- rbind(AL_2.1.7, AL_2.3.4, GU_2.3.1, GU_3.0.3)


## Plotting

#The first plot compares the global mismatch frequency per base-caller
colors=c("cyan3","cyan4", "red", "red3")

g <- ggplot(data, aes(Basecaller, log(mismatches_freq), color=Basecaller, fill=Basecaller)) + 
  geom_violin(alpha = 0.3) +
  geom_boxplot(alpha = 0.4, width=0.2, show.legend = FALSE) +
  scale_color_manual(values = colors) + 
  scale_fill_manual(values = colors) +
  ylab("Mismatch Frequency (log)") + theme_bw() + ggtitle(paste("Mismatch Frequency per Base-Caller:", name)) + 
  theme(plot.title = element_text(size = 10, face = "bold"))

ggsave("mismatch_freq_per_basecaller.pdf", g, device = "pdf", path = dir, width = 4, height = 5)

#The next set of plots evaluate the number of bases with mismatch frequency higher than a given threshold. First with the total number and then with the proportion
plotting_threshold <- function(threshold){
  ggplot(data[data$mismatches_freq>threshold,], aes(Basecaller, fill=ref_nuc)) + 
    geom_bar()  + labs(fill = "Reference Nucleotide") + theme_bw() + ggtitle(paste("Threshold:", threshold*100, "%")) +
    geom_text(stat="count", aes(label=..count..), position = "stack", vjust=1.5, size=3.5) + 
    scale_x_discrete(labels= c("AL\n2.1.7","AL\n2.3.4","GU\n2.3.1","GU\n3.0.3")) +
    theme(plot.title = element_text(hjust = 0.5))
}

plotting_threshold_proportion <- function(threshold){
  ggplot(data[data$mismatches_freq>threshold,], aes(Basecaller, fill=ref_nuc)) + 
    geom_bar(position = "fill")  + labs(fill = "Reference Nucleotide") + theme_bw() + ggtitle(paste("Threshold:", threshold*100, "%")) +
    scale_x_discrete(labels= c("AL\n2.1.7","AL\n2.3.4","GU\n2.3.1","GU\n3.0.3")) +
    theme(plot.title = element_text(hjust = 0.5))
}


g1 <- plotting_threshold(0.05)
g2 <- plotting_threshold(0.1)
g3 <- plotting_threshold(0.2)
g4 <- plotting_threshold(0.5)

l <- list(g1, g2, g3, g4)
figure <- ggarrange(plotlist=l, ncol=4, nrow=1, common.legend = TRUE, legend = "bottom") 
f <- annotate_figure(figure, top = text_grob(paste("Number of Bases with Mismatch Frequency Higher than a Threhold in", name), face = "bold"))

ggsave("bases_above_threshold.pdf", f, device = "pdf", path = dir, width = 9, height = 7)

g1 <- plotting_threshold_proportion(0.05)
g2 <- plotting_threshold_proportion(0.1)
g3 <- plotting_threshold_proportion(0.2)
g4 <- plotting_threshold_proportion(0.5)

l <- list(g1, g2, g3, g4)
figure <- ggarrange(plotlist=l, ncol=4, nrow=1, common.legend = TRUE, legend = "bottom") 
f <- annotate_figure(figure, top = text_grob(paste("Proportion of Bases with Mismatch Frequency Higher than a Threhold in", name), face = "bold"))

ggsave("bases_above_threshold_proportion.pdf", f, device = "pdf", path = dir, width = 9, height = 7)


# These plots now compute the mismatch frequencies per nucleotide

plotting_mismatch_frequencies <- function(i){
  ggplot(data = i, mapping = aes(ref_nuc, log(mismatches_freq), color=ref_nuc, fill=ref_nuc)) +
    geom_violin(alpha = 0.3) +
    geom_boxplot(alpha = 0.5, show.legend = FALSE, width = 0.2) +
    ggtitle(deparse(substitute(i))) + ylab("Mismatch Frequency (log)") + xlab("Nucleotide") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
}

g1 <- plotting_mismatch_frequencies(AL_2.1.7)
g2 <- plotting_mismatch_frequencies(AL_2.3.4)
g3 <- plotting_mismatch_frequencies(GU_2.3.1)
g4 <- plotting_mismatch_frequencies(GU_3.0.3)

l <- list(g1,g2,g3,g4)
figure <- ggarrange(plotlist=l, ncol=4, nrow=1, common.legend = TRUE, legend = "bottom") 
f <- annotate_figure(figure, top = text_grob(paste("Mismatch Frequency per Nucleotide in", name), face = "bold"))

ggsave("mismatch_per_nucleotide.pdf", f, device = "pdf", path = dir, width = 9, height = 5)


##########
# Now we compute the mismatch signature per base-caller

col <- list("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF")
names(col) <- c("A", "C", "G", "T")

mismatch_pattern <- function(input){
dat <- melt(input[,c("ref_nuc", "A", "T", "C", "G")])
dat <- tapply(dat$value, list(dat$ref_nuc,dat$variable), sum) #rows are ref_nuc, columns are read nucleotides
dat <- melt(dat)
colnames(dat) <- c("ref_nuc", "Variable", "value")

A <- dat[c(1,5,9,13),]
C <- dat[c(2,6,10,14),]
G <- dat[c(3,7,11,15),]
T <- dat[c(4,8,12,16),]

A$Variable <- factor(A$Variable, levels = as.character(A[order(A$value),]$Variable))
C$Variable <- factor(C$Variable, levels = as.character(C[order(C$value),]$Variable))
G$Variable <- factor(G$Variable, levels = as.character(G[order(G$value),]$Variable))
T$Variable <- factor(T$Variable, levels = as.character(T[order(T$value),]$Variable))

first <- as.character(A[order(A$value),]$Variable)[1]
second <- as.character(A[order(A$value),]$Variable)[2]
third <- as.character(A[order(A$value),]$Variable)[3]
fourth <- as.character(A[order(A$value),]$Variable)[4]

ggplot() + 
  geom_bar(data=A, mapping=aes(ref_nuc, value, fill=Variable), stat="identity", position="fill") +
  geom_bar(data=C, mapping=aes(ref_nuc, value, fill=Variable), stat="identity", position="fill") +
  geom_bar(data=G, mapping=aes(ref_nuc, value, fill=Variable), stat="identity", position="fill") +
  geom_bar(data=T, mapping=aes(ref_nuc, value, fill=Variable), stat="identity", position="fill") +
  coord_cartesian(ylim=c(0.9,1)) +
  scale_fill_manual("Base-called Nucleotide", 
                    breaks = c(first, second, third, fourth),
                    values = c(col[[first]], col[[second]], col[[third]], col[[fourth]]),
                    guide = guide_legend(reverse = TRUE)) +
  xlab("Reference Nucleotide") + ylab("count") +
  theme_bw() +
  ggtitle(levels(factor(input$Basecaller))) +
  theme(plot.title = element_text(size = 12, hjust = 0.5), axis.text=element_text(size=6), axis.title=element_text(size=9))
}

g1 <- mismatch_pattern(AL_2.1.7)
g2 <- mismatch_pattern(AL_2.3.4)
g3 <- mismatch_pattern(GU_2.3.1)
g4 <- mismatch_pattern(GU_3.0.3)

l <- list(g1,g2,g3,g4)
figure <- ggarrange(plotlist=l, ncol=4, nrow=1, common.legend = TRUE, legend = "bottom") 
f <- annotate_figure(figure, top = text_grob(paste("Mismatch Frequency per Nucleotide in", name), face = "bold"))

ggsave("mismatch_pattern_zoom.pdf", f, device = "pdf", path = dir, width = 8, height = 4)


mismatch_pattern_100 <- function(input){
  dat <- melt(input[,c("ref_nuc", "A", "T", "C", "G")])
  dat <- tapply(dat$value, list(dat$ref_nuc,dat$variable), sum) #rows are ref_nuc, columns are read nucleotides
  dat <- melt(dat)
  colnames(dat) <- c("ref_nuc", "Variable", "value")
  
  A <- dat[c(1,5,9,13),]
  C <- dat[c(2,6,10,14),]
  G <- dat[c(3,7,11,15),]
  T <- dat[c(4,8,12,16),]
  
  A$Variable <- factor(A$Variable, levels = as.character(A[order(A$value),]$Variable))
  C$Variable <- factor(C$Variable, levels = as.character(C[order(C$value),]$Variable))
  G$Variable <- factor(G$Variable, levels = as.character(G[order(G$value),]$Variable))
  T$Variable <- factor(T$Variable, levels = as.character(T[order(T$value),]$Variable))
  
  first <- as.character(A[order(A$value),]$Variable)[1]
  second <- as.character(A[order(A$value),]$Variable)[2]
  third <- as.character(A[order(A$value),]$Variable)[3]
  fourth <- as.character(A[order(A$value),]$Variable)[4]
  
  ggplot() + 
    geom_bar(data=A, mapping=aes(ref_nuc, value, fill=Variable), stat="identity", position="fill") +
    geom_bar(data=C, mapping=aes(ref_nuc, value, fill=Variable), stat="identity", position="fill") +
    geom_bar(data=G, mapping=aes(ref_nuc, value, fill=Variable), stat="identity", position="fill") +
    geom_bar(data=T, mapping=aes(ref_nuc, value, fill=Variable), stat="identity", position="fill") +
    scale_fill_manual("Base-called Nucleotide", 
                      breaks = c(first, second, third, fourth),
                      values = c(col[[first]], col[[second]], col[[third]], col[[fourth]]),
                      guide = guide_legend(reverse = TRUE)) +
    xlab("Reference Nucleotide") + ylab("count") +
    theme_bw() +
    ggtitle(levels(factor(input$Basecaller))) +
    theme(plot.title = element_text(size = 12, hjust = 0.5), axis.text=element_text(size=6), axis.title=element_text(size=9))
}

g1 <- mismatch_pattern_100(AL_2.1.7)
g2 <- mismatch_pattern_100(AL_2.3.4)
g3 <- mismatch_pattern_100(GU_2.3.1)
g4 <- mismatch_pattern_100(GU_3.0.3)

l <- list(g1,g2,g3,g4)
figure <- ggarrange(plotlist=l, ncol=4, nrow=1, common.legend = TRUE, legend = "bottom") 
f <- annotate_figure(figure, top = text_grob(paste("Mismatch Frequency per Nucleotide in", name), face = "bold"))

ggsave("mismatch_pattern_100.pdf", f, device = "pdf", path = dir, width = 8, height = 4)

# Finally we compute the ternary diagrams for plotting the mismatch signature per nucleotide

ternaries <- function(data){
A <- data[data$ref_nuc=="A",]
C <- data[data$ref_nuc=="C",]
G <- data[data$ref_nuc=="G",]
T <- data[data$ref_nuc=="T",]

g1 <- ggtern(A, aes(C,G,T)) +
  geom_point(size = 0.05, aes(color = log(coverage))) +
  ggtitle("Ternary Diagram for A") + theme_bw(base_size = 6) +
  scale_color_continuous("Coverage (log)") +
  theme(plot.title = element_text(size = 8, hjust = 0.5)) 
legend <- get_legend(g1)
g1 <- g1 + theme(legend.position="none")

g2 <- ggtern(C, aes(A,G,T)) +
  geom_point(size = 0.05, aes(color = log(coverage))) +
  ggtitle("Ternary Diagram for C") + theme_bw(base_size = 6) +
  scale_color_continuous("Coverage (log)") +
  theme(legend.position="none", plot.title = element_text(size = 8, hjust = 0.5))
g3 <- ggtern(G, aes(A,C,T)) +
  geom_point(size = 0.05, aes(color = log(coverage))) +
  ggtitle("Ternary Diagram for G") + theme_bw(base_size = 6) +
  scale_color_continuous("Coverage (log)") +
  theme(legend.position="none", plot.title = element_text(size = 8, hjust = 0.5))
g4 <- ggtern(T, aes(A,C,G)) +
  geom_point(size = 0.05, aes(color = log(coverage))) +
  ggtitle("Ternary Diagram for T") + theme_bw(base_size = 6) +
  scale_color_continuous("Coverage (log)") +
  theme(legend.position="none", plot.title = element_text(size = 8, hjust = 0.5))

g<- grid.arrange(g1,g2,g3,g4, nrow=1, ncol=5, legend, top=paste("Ternary Diagrams in", levels(factor(data$Basecaller))))
ggsave(paste("ternary_diagram", levels(factor(data$Basecaller)), ".pdf", sep="_"), g, device = "pdf", path = dir, width = 7, height = 3)
}

ternaries(AL_2.1.7)
ternaries(AL_2.3.4)
ternaries(GU_2.3.1)
ternaries(GU_3.0.3)
  
