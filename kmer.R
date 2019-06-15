## Loading libraries

library(ggplot2)
library(ggpubr)
library(optparse)

## Parsing

parser <- OptionParser()
parser <- add_option(parser, opt_str=c("-i", "--input"), type="character",
                     dest='input',
                     help="Input Directory")
options=parse_args(parser)
dir=options$input
#dir="GU_2.3.1_graphmap/"

## Creating global variables

cols <- rainbow(5)

## Creating functions

# Quality Plots
quality <- function(u){
  g <- ggplot(u) + geom_boxplot(aes("q1", q1, color=cols[1], fill=cols[1], alpha=0.6)) +
    geom_boxplot(aes("q2",q2, color=cols[2], fill=cols[2], alpha= 0.6)) +
    geom_boxplot(aes("q3",q3, color=cols[3], fill=cols[3], alpha=0.6)) +
    geom_boxplot(aes("q4",q4, color=cols[4], fill=cols[4], alpha=0.6)) +
    geom_boxplot(aes("q5",q5, color=cols[5], fill=cols[5], alpha=0.6)) +
    ylab("Quality") +
    xlab("") + 
    ggtitle(paste(deparse(substitute(u)), "Quality")) +
    guides(fill=FALSE, color=FALSE, alpha=FALSE) +
    theme_bw()
  
  ggsave(paste(dir, deparse(substitute(u)), "_qual.pdf", sep=""), g, device = "pdf", width = 3.5, height = 2.7)
}

#Mimsatches Plot
mismatches <- function(u, mod){
  g <- ggplot(u) + geom_boxplot(aes("mis1", log(mis1), color=cols[1], fill=cols[1], alpha=0.6)) +
    geom_boxplot(aes("mis2", log(mis2), color=cols[2], fill=cols[2], alpha= 0.6)) +
    geom_boxplot(aes("mis3", log(mis3), color=cols[3], fill=cols[3], alpha=0.6)) +
    geom_boxplot(aes("mis4", log(mis4), color=cols[4], fill=cols[4], alpha=0.6)) +
    geom_boxplot(aes("mis5", log(mis5), color=cols[5], fill=cols[5], alpha=0.6)) +
    ylab("Mismatch frequency (log)") +
    xlab("") + 
    ggtitle(paste(deparse(substitute(u)), "Mismatches")) +
    guides(fill=FALSE, color=FALSE, alpha=FALSE) +
    theme_bw()
  
  ggsave(paste(dir, deparse(substitute(u)), "_mis.pdf", sep=""), g, device = "pdf", width = 3.5, height = 2.7)
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

UNM_A <- read.delim(paste(dir, "UNM_A_per_site.var.csv.slided.onekmer.oneline.5mer.filtered.csv", sep=""), sep=",")
m6A <- read.delim(paste(dir, "m6A_per_site.var.csv.slided.onekmer.oneline.5mer.filtered.csv", sep=""), sep=",")

quality(UNM_A)
quality(m6A)
mismatches(UNM_A)
mismatches(m6A)

UNM_A$condition <- "UNM"
m6A$condition <- "m6A"

data <- rbind(UNM_A, m6A)
data$condition <- factor(data$condition, levels = c("UNM","m6A"))

figure <- final_plot(data)

f <- annotate_figure(figure, top = text_grob("UNM vs m6A", face = "bold"))
ggsave(paste(dir, "UNM_m6A.pdf", sep=""), f, device = "pdf", width = 8, height = 8)

# statistical test: Mann-Whitney U test:
wilcox.test(q1 ~ condition, data=data) 
wilcox.test(q2 ~ condition, data=data) 
wilcox.test(q3 ~ condition, data=data) 
wilcox.test(q4 ~ condition, data=data) 
wilcox.test(q5 ~ condition, data=data) 

wilcox.test(mis1 ~ condition, data=data) 
wilcox.test(mis2 ~ condition, data=data) 
wilcox.test(mis3 ~ condition, data=data) 
wilcox.test(mis4 ~ condition, data=data) 
wilcox.test(mis5 ~ condition, data=data) 

wilcox.test(del1 ~ condition, data=data) 
wilcox.test(del2 ~ condition, data=data) 
wilcox.test(del3 ~ condition, data=data) 
wilcox.test(del4 ~ condition, data=data) 
wilcox.test(del5 ~ condition, data=data) 

## m5C and 5hmC

UNM_C <- read.delim(paste(dir, "UNM_C_per_site.var.csv.slided.onekmer.oneline.5mer.filtered.csv", sep=""), sep=",")
m5C <- read.delim(paste(dir, "m5C_per_site.var.csv.slided.onekmer.oneline.5mer.filtered.csv", sep=""), sep=",")
hm5C <- read.delim(paste(dir, "5hmC_per_site.var.csv.slided.onekmer.oneline.5mer.filtered.csv", sep=""), sep=",")

quality(UNM_C)
quality(m5C)
quality(hm5C)

mismatches(UNM_C)
mismatches(m5C)
mismatches(hm5C)

UNM_C$condition <- "UNM"
m5C$condition <- "m5C"
hm5C$condition <- "5hmC"

data <- rbind(UNM_C, m5C, hm5C)
data$condition <- factor(data$condition, levels = c("UNM","m5C","5hmC"))

figure <- final_plot(data)

f <- annotate_figure(figure, top = text_grob("UNM vs m5C vs 5hmC", face = "bold"))
ggsave(paste(dir, "UNM_m5C_5hmC.pdf", sep=""), f, device = "pdf", width = 8.3, height = 8)


## pU

UNM_T <- read.delim(paste(dir, "UNM_T_per_site.var.csv.slided.onekmer.oneline.5mer.filtered.csv", sep=""), sep=",")
pU <- read.delim(paste(dir, "pU_per_site.var.csv.slided.onekmer.oneline.5mer.filtered.csv", sep=""), sep=",")

quality(UNM_T)
quality(pU)
mismatches(UNM_T)
mismatches(pU)

UNM_T$condition <- "UNM"
pU$condition <- "pU"

data <- rbind(UNM_T, pU)
data$condition <- factor(data$condition, levels = c("UNM","pU"))

figure <- final_plot(data)

f <- annotate_figure(figure, top = text_grob("UNM vs pU", face = "bold"))
ggsave(paste(dir, "UNM_pU.pdf", sep=""), f, device = "pdf", width = 8, height = 8)



