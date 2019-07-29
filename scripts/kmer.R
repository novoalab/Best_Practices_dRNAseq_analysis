## Loading libraries

suppressMessages(library(ggplot2))
library(ggpubr)
library(optparse)
options(warn=-1)

## Parsing

parser <- OptionParser()
parser <- add_option(parser, opt_str=c("-u", "--unm"), type="character",
                     dest='unm',
                     help="unm file")
parser <- add_option(parser, opt_str=c("-m", "--mod"), type="character",
                     dest='mod',
                     help="mod file")
parser <- add_option(parser, opt_str=c("-o", "--output"), type="character",
                     dest='output',
                     help="Output Directory")
parser <- add_option(parser, opt_str=c("-n", "--names"), type="character",
                     default = "UNM,MOD",
                     dest='names',
                     help="Names")
options=parse_args(parser)
unm=options$unm
mod=options$mod
output=options$output
names=options$names

n <- strsplit(names, ",")

## Creating global variables

cols <- rainbow(5)

## Creating functions

# Quality Plots
quality <- function(u, name){
  g <- ggplot(u) + geom_boxplot(aes("q1", q1, color=cols[1], fill=cols[1], alpha=0.6)) +
    geom_boxplot(aes("q2",q2, color=cols[2], fill=cols[2], alpha= 0.6)) +
    geom_boxplot(aes("q3",q3, color=cols[3], fill=cols[3], alpha=0.6)) +
    geom_boxplot(aes("q4",q4, color=cols[4], fill=cols[4], alpha=0.6)) +
    geom_boxplot(aes("q5",q5, color=cols[5], fill=cols[5], alpha=0.6)) +
    ylab("Quality") +
    xlab("") + 
    ggtitle(paste(name, "Quality")) +
    guides(fill=FALSE, color=FALSE, alpha=FALSE) +
    theme_bw()
  
  ggsave(paste(output, name, "_qual.pdf", sep=""), g, device = "pdf", width = 3.5, height = 2.7)
}

#Mimsatches Plot
mismatches <- function(u, name){
  g <- ggplot(u) + geom_boxplot(aes("mis1", log(mis1), color=cols[1], fill=cols[1], alpha=0.6)) +
    geom_boxplot(aes("mis2", log(mis2), color=cols[2], fill=cols[2], alpha= 0.6)) +
    geom_boxplot(aes("mis3", log(mis3), color=cols[3], fill=cols[3], alpha=0.6)) +
    geom_boxplot(aes("mis4", log(mis4), color=cols[4], fill=cols[4], alpha=0.6)) +
    geom_boxplot(aes("mis5", log(mis5), color=cols[5], fill=cols[5], alpha=0.6)) +
    ylab("Mismatch frequency (log)") +
    xlab("") + 
    ggtitle(paste(name, "Mismatches")) +
    guides(fill=FALSE, color=FALSE, alpha=FALSE) +
    theme_bw()
  
  ggsave(paste(output, name, "_mis.pdf", sep=""), g, device = "pdf", width = 3.5, height = 2.7)
}

#Final Plot
final_plot <- function(dat){
  g1 <- ggplot(data=dat) + geom_boxplot(aes(condition, q1, fill=condition, alpha = 0.6, color=condition)) + theme_bw() + xlab("") + theme(legend.position="none")
  g2 <- ggplot(data=dat) + geom_boxplot(aes(condition, q2, fill=condition, alpha = 0.6, color=condition)) + theme_bw() + xlab("") + theme(legend.position="none")
  g3 <- ggplot(data=dat) + geom_boxplot(aes(condition, q3, fill=condition, alpha = 0.6, color=condition)) + theme_bw() + xlab("") + theme(legend.position="none")
  g4 <- ggplot(data=dat) + geom_boxplot(aes(condition, q4, fill=condition, alpha = 0.6, color=condition)) + theme_bw() + xlab("") + theme(legend.position="none")
  g5 <- ggplot(data=dat) + geom_boxplot(aes(condition, q5, fill=condition, alpha = 0.6, color=condition)) + theme_bw() + xlab("") + theme(legend.position="none")
  
  g6 <- ggplot(data=dat) + geom_boxplot(aes(condition, log(mis1), fill=condition, alpha = 0.6, color=condition)) + theme_bw() + xlab("") + theme(legend.position="none")
  g7 <- ggplot(data=dat) + geom_boxplot(aes(condition, log(mis2), fill=condition, alpha = 0.6, color=condition)) + theme_bw() + xlab("") + theme(legend.position="none")
  g8 <- ggplot(data=dat) + geom_boxplot(aes(condition, log(mis3), fill=condition, alpha = 0.6, color=condition)) + theme_bw() + xlab("") + theme(legend.position="none")
  g9 <- ggplot(data=dat) + geom_boxplot(aes(condition, log(mis4), fill=condition, alpha = 0.6, color=condition)) + theme_bw() + xlab("") + theme(legend.position="none")
  g10 <- ggplot(data=dat) + geom_boxplot(aes(condition, log(mis5), fill=condition, alpha = 0.6, color=condition)) + theme_bw() + xlab("") + theme(legend.position="none")
  
  g11 <- ggplot(data=dat) + geom_boxplot(aes(condition, log(del1), fill=condition, alpha = 0.6, color=condition)) + theme_bw() + xlab("") + theme(legend.position="none")
  g12 <- ggplot(data=dat) + geom_boxplot(aes(condition, log(del2), fill=condition, alpha = 0.6, color=condition)) + theme_bw() + xlab("") + theme(legend.position="none")
  g13 <- ggplot(data=dat) + geom_boxplot(aes(condition, log(del3), fill=condition, alpha = 0.6, color=condition)) + theme_bw() + xlab("") + theme(legend.position="none")
  g14 <- ggplot(data=dat) + geom_boxplot(aes(condition, log(del4), fill=condition, alpha = 0.6, color=condition)) + theme_bw() + xlab("") + theme(legend.position="none")
  g15 <- ggplot(data=dat) + geom_boxplot(aes(condition, log(del5), fill=condition, alpha = 0.6, color=condition)) + theme_bw() + xlab("") + theme(legend.position="none")
  
  #g16 <- ggplot(data) + geom_boxplot(aes(condition, log(ins1), fill=condition, alpha = 0.6, color=condition)) + theme_bw() + xlab("") + theme(legend.position="none")
  #g17 <- ggplot(data) + geom_boxplot(aes(condition, log(ins2), fill=condition, alpha = 0.6, color=condition)) + theme_bw() + xlab("") + theme(legend.position="none")
  #g18 <- ggplot(data) + geom_boxplot(aes(condition, log(ins3), fill=condition, alpha = 0.6, color=condition)) + theme_bw() + xlab("") + theme(legend.position="none")
  #g19 <- ggplot(data) + geom_boxplot(aes(condition, log(ins4), fill=condition, alpha = 0.6, color=condition)) + theme_bw() + xlab("") + theme(legend.position="none")
  #g20 <- ggplot(data) + geom_boxplot(aes(condition, log(ins5), fill=condition, alpha = 0.6, color=condition)) + theme_bw() + xlab("") + theme(legend.position="none")
  
  l <- list(g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12,g13,g14,g15)
  figure <- ggarrange(plotlist=l, ncol=5, nrow=3)
  return(figure)
}



## m6A

UNM <- read.delim(unm, sep=",")
MOD <- read.delim(mod, sep=",")
quality(UNM, n[[1]][1])
quality(MOD, n[[1]][2])
mismatches(UNM, n[[1]][1])
mismatches(MOD, n[[1]][2])

UNM$condition <- n[[1]][1]
MOD$condition <- n[[1]][2]

data <- rbind(UNM, MOD)
data$condition <- factor(data$condition, levels = c(n[[1]][1],n[[1]][2]))

figure <- final_plot(data)

f <- annotate_figure(figure, top = text_grob(paste(n[[1]][1], "vs", n[[1]][2]), face = "bold"))
ggsave(paste(output, n[[1]][1], "_", n[[1]][2], ".pdf", sep=""), f, device = "pdf", width = 8, height = 8)

# statistical test: Mann-Whitney U test:
print("P values of Mann-Whitney U comparing UNM vs MOD")
m <- matrix(1:15, nrow = 3, ncol = 5)
colnames(m) <- c("pos1", "pos2", "pos3", "pos4", "pos5")
rownames(m) <- c("qualities", "mismatches", "deletions")
m[1,1] <- wilcox.test(q1 ~ condition, data=data)$p.value
m[1,2] <- wilcox.test(q2 ~ condition, data=data)$p.value
m[1,3] <- wilcox.test(q3 ~ condition, data=data)$p.value
m[1,4] <- wilcox.test(q4 ~ condition, data=data)$p.value
m[1,5] <- wilcox.test(q5 ~ condition, data=data)$p.value

m[2,1] <- wilcox.test(mis1 ~ condition, data=data)$p.value
m[2,2] <- wilcox.test(mis2 ~ condition, data=data)$p.value
m[2,3] <- wilcox.test(mis3 ~ condition, data=data)$p.value
m[2,4] <- wilcox.test(mis4 ~ condition, data=data)$p.value 
m[2,5] <- wilcox.test(mis5 ~ condition, data=data)$p.value

m[3,1] <- wilcox.test(del1 ~ condition, data=data)$p.value
m[3,2] <- wilcox.test(del2 ~ condition, data=data)$p.value
m[3,3] <- wilcox.test(del3 ~ condition, data=data)$p.value
m[3,4] <- wilcox.test(del4 ~ condition, data=data)$p.value
m[3,5] <- wilcox.test(del5 ~ condition, data=data)$p.value

print(m)
