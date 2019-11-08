######################################################################################################################################
#This script reads n files with a name starting by "OUTPUT_", containing for each base-called read its ID, its length and its mean quality
#It outputs a boxplot comparing the mean read qualities between approaches and another one comparing the read lengths per base-caller.
#It also compares how the reads base-called by AL 2.3.4 but not by AL 2.1.7 look regarding qualities and lengths. 
#And this same comparison for the reads base-called by GU 3.0.3 and not by GU 2.3.1
######################################################################################################################################


#Loading libraries


required.packages <- c("ggplot2","ggExtra","optparse")
new.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

suppressMessages(library(ggplot2))
library(ggExtra)
library(optparse)
source("scripts/marginal_plot_edited_colors.R")

#Parsing

parser <- OptionParser()
parser <- add_option(parser, opt_str=c("-i", "--input"), type="character",
                     dest='input',
                     help="Input Directory")
parser <- add_option(parser, opt_str=c("-t", "--title"), type="character",
                     default = "",
                     dest='title',
                     help="Title")
parser <- add_option(parser, opt_str=c("-n", "--names"), type="character",
                     default = "1,2,3,4",
                     dest='names',
                     help="Names of the datasets")
options=parse_args(parser)
dir=options$input 
title=options$title
names=options$names

print(names)
#Loading data:
n <- strsplit(names, ",")
files <- list.files(path=dir, pattern = "OUTPUT_.*")

l <- list()
table <- c()
for (i in 1:length(files)){
  assign(files[i],read.delim(paste(dir, files[i], sep="/")))
  l[[i]] <- get(files[i])
  l[[i]]$dataset <- n[[1]][i]
  table <- rbind(table, l[[i]])
}


#Plotting read lengths

title=paste("Read Length", title)

g <- ggplot(table, aes(dataset, log(read_length))) + 
  geom_boxplot(outlier.size = 0.7, colour=c(colorRampPalette(c("red", "blue"))(length(files))), 
               fill=c(colorRampPalette(c("red", "blue"))(length(files))), alpha=0.4) + 
  ylab("read length (log)") + ggtitle(title) + xlab("") + theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, size = 9), 
        axis.text=element_text(size=4), 
        axis.title=element_text(size=6))
ggsave("read_length.png", g, device = "png", path = dir, width = 1.8, height = 2.5)

#Plotting qualities

title=paste("Mean Quality", title)
g <- ggplot(table, aes(dataset, log(mean_quality))) + 
  geom_boxplot(outlier.size = 0.7, colour=c(colorRampPalette(c("red", "blue"))(length(files))), 
               fill=c(colorRampPalette(c("red", "blue"))(length(files))), alpha=0.4) + 
  ylab("mean quality (log)") + ggtitle(title) + xlab("") + theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, size = 9), 
        axis.text=element_text(size=4), 
        axis.title=element_text(size=6))
ggsave("quality.png", g, device = "png", path = dir, width = 1.8, height = 2.5)



# Plotting comparisons between versions

# The following plots will do pairwise comparisons, comparing the reads that are present in one fastq file but not in the other and the common reads

for (i in 1:length(files)){
  for (j in 1:length(files)){
    if (i!=j){
      l[[i]]$Condition <- l[[i]]$read_id %in% l[[j]]$read_id
      if(nrow(l[[i]][l[[i]]$Condition==FALSE,])!=0){
        l[[i]][l[[i]]$Condition==FALSE,]$Condition <- "NEW"}
      l[[i]][l[[i]]$Condition==TRUE,]$Condition <- "COMMON"
      png(paste(dir, "/",n[[1]][i], "_against_", n[[1]][j], ".png", sep=""))
      marginal_plot(x = read_length, y = mean_quality, group = Condition, data = l[[i]], alpha = 0.5, bw = "nrd", lm_formula = NULL, 
                    xlab = "Read Length (log)", ylab = "Mean Quality (log)", pch = 15, cex = 0.5, log="xy")
      dev.off()
}}}

