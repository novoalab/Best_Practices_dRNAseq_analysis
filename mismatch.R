#Loading libraries
library(ggplot2)
library(ggpubr)
library(optparse)


#Parsing

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
print(dir)
dir="RNA080420/"
name="5hmC"

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


names <- c("AL_2.1.7", "AL_2.3.4", "GU_2.3.1", "GU_3.0.3")
for (i in 1:length(names)){
  assign(names[i], merge(get(paste(names[i],"_S",sep="")), 
                         get(paste(names[i],"_M",sep="")), 
                         by = c("chr", "pos", "ref_nuc"), 
                         sort = F))}

#I don't have NAs by now, so I skip the step in which I was removing NAs
# I am now removing positions with low coverage by now, because mismatch frequency may be very inflated if the coverage is 1 and the mismatch frequency in that position is 100 %

l <- list(AL_2.1.7, AL_2.3.4, GU_2.3.1, GU_3.0.3)
remove_low_coverage <- function(i){
  cov <- mean(i$coverage) # We compute the coverage mean
  threshold <- cov*0.1 #We create a threshold for removing values lower than 10% of the computed mean
  i <- i[i$coverage>threshold,]
}

l <- lapply(l, remove_low_coverage)


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

AL_2.1.7$mismatches_freq <- AL_2.1.7$mismatches/AL_2.1.7$coverage
AL_2.3.4$mismatches_freq <- AL_2.3.4$mismatches/AL_2.3.4$coverage
GU_2.3.1$mismatches_freq <- GU_2.3.1$mismatches/GU_2.3.1$coverage
GU_3.0.3$mismatches_freq <- GU_3.0.3$mismatches/GU_3.0.3$coverage

AL_2.1.7$Basecaller <- "AL 2.1.7"
AL_2.3.4$Basecaller <- "AL 2.3.4"
GU_2.3.1$Basecaller <- "GU 2.3.1"
GU_3.0.3$Basecaller <- "GU 3.0.3"

data <- rbind(AL_2.1.7, AL_2.3.4, GU_2.3.1, GU_3.0.3)

#Plotting

colors=c("cyan3","cyan4", "red", "red3")

g <- ggplot(data, aes(Basecaller, log(mismatches_freq), color=Basecaller)) + 
  geom_boxplot(coef=3, fill=colors, alpha=0.2, show.legend = FALSE)+
  geom_point(alpha=0.2, position=position_jitter(width=0.2), show.legend = FALSE) +
  scale_color_manual(values = colors) + 
  ylab("Mismatch Frequency (log)") + theme_bw() + ggtitle(paste("Mismatch Frequency per Base-Caller:", name))

ggsave("mismatch_freq_per_basecaller.png", g, device = "png", path = dir, width = 5, height = 7)

# 

plotting_threshold <- function(threshold){
  ggplot(data[data$mismatches_freq>threshold,], aes(Basecaller, fill=ref_nuc)) + 
    geom_bar()  + labs(fill = "Reference Nucleotide") + theme_bw() + ggtitle(paste("Threshold:", threshold*100, "%")) +
    geom_text(stat="count", aes(label=..count..), position = "stack", vjust=1.5, size=3.5) + 
    scale_x_discrete(labels= c("AL\n2.1.7","AL\n2.3.4","GU\n2.3.1","GU\n3.0.3"))
}

g1 <- plotting_threshold(0.05)
g2 <- plotting_threshold(0.1)
g3 <- plotting_threshold(0.2)
g4 <- plotting_threshold(0.5)

l <- list(g1, g2, g3, g4)
figure <- ggarrange(plotlist=l, ncol=4, nrow=1, common.legend = TRUE, legend = "bottom") 
f <- annotate_figure(figure, top = text_grob(paste("Number of Bases with Mismatch Frequency Higher than a Threhold in", name), face = "bold"))

ggsave("bases_above_threshold.png", f, device = "png", path = dir, width = 8, height = 7)

#

plotting_mismatch_frequencies <- function(i){
ggplot(data = i, mapping = aes(ref_nuc, log(mismatches_freq))) +
  geom_boxplot(coef=3, alpha=0.4, show.legend = FALSE, aes(fill=factor(ref_nuc))) + geom_point(aes(colour=factor(ref_nuc)), show.legend = FALSE, alpha =0.6, position = position_jitter(width=0.15)) +
  ggtitle(deparse(substitute(i))) + ylab("Mismatch Frequency (log)") + xlab("Nucleotide") +
  theme_bw() 
}

g1 <- plotting_mismatch_frequencies(AL_2.1.7)
g2 <- plotting_mismatch_frequencies(AL_2.3.4)
g3 <- plotting_mismatch_frequencies(GU_2.3.1)
g4 <- plotting_mismatch_frequencies(GU_3.0.3)

l <- list(g1,g2,g3,g4)
figure <- ggarrange(plotlist=l, ncol=4, nrow=1) 
f <- annotate_figure(figure, top = text_grob(paste("Mismatch Frequency per Nucleotide in", name), face = "bold"))

ggsave("mismatch_per_nucleotide.png", f, device = "png", path = dir, width = 9, height = 7)


#i=GU_3.0.3
#mean(i$mismatches_freq)
#mean(i[i$ref_nuc=="C",]$mismatches_freq)
#mean(i[i$ref_nuc=="T",]$mismatches_freq)
#It makes sense that C has not the highest mismatch frequency, because as I have seen in IGV, the whole k-mer is affected


