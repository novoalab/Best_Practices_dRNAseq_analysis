######################################################################################################################################
#This script reads 4 files, one per base-caller, containing for each base-called read its ID, its length and its mean quality
#It outputs a boxplot comparing the mean read qualities between approaches and another one compare the read lengths per base-caller.
#It also compares how the reads base-called by AL 2.3.4 but not by AL 2.1.7 look regarding qualities and lengths. 
#And this same comparison for the reads base-called by GU 3.0.3 and not by GU 2.3.1
######################################################################################################################################


#Loading libraries
library(ggplot2)
library(ggExtra)
library(optparse)
source("marginal_plot_edited_colors.R")


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

#Loading data:

AL_O <- read.delim(paste(dir, "/ALBACORE_O_OUTPUT", sep = ""))
AL_N <- read.delim(paste(dir, "/ALBACORE_N_OUTPUT", sep = ""))
GU_O <- read.delim(paste(dir, "/GUPPY_O_OUTPUT", sep = ""))
GU_N <- read.delim(paste(dir, "/GUPPY_N_OUTPUT", sep = ""))

AL_O$basecaller <- "Albacore 2.1.7"
AL_N$basecaller <- "Albacore 2.3.4"
GU_O$basecaller <- "Guppy 2.3.1"
GU_N$basecaller <- "Guppy 3.0.3"

table <- rbind(AL_O, AL_N, GU_O, GU_N)

#Plotting read lengths

title=paste("Read Length", name)

g <- ggplot(table, aes(basecaller, log(read_length))) + 
  geom_boxplot(outlier.size = 0.7, colour=c("cyan3","cyan4", "red", "red3"), 
               fill=c("cyan3","cyan4", "red", "red3"), alpha=0.4) + 
  ylab("read length (log)") + ggtitle(title) + xlab("") + theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, size = 9), 
        axis.text=element_text(size=3.5), 
        axis.title=element_text(size=6))
ggsave("read_length.png", g, device = "png", path = dir, width = 1.8, height = 2.5)

#Plotting qualities

title=paste("Mean Quality", name)
g <- ggplot(table, aes(basecaller, log(mean_quality))) + 
  geom_boxplot(outlier.size = 0.7, colour=c("cyan3","cyan4", "red", "red3"), 
               fill=c("cyan3","cyan4", "red", "red3"), alpha=0.4) + 
  ylab("Mean Quality (log)") + ggtitle(title) + xlab("") + theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, size = 9), 
        axis.text=element_text(size=3.5), 
        axis.title=element_text(size=6))
ggsave("quality.png", g, device = "png", path = dir, width = 1.8, height = 2.5)

# Plotting comparisons between versions

#The following plots will compare the reads base-called by both versions of each base-caller against the reads only 
#base-called by the new version.

# Albacore

#we check whether a read in the new version of Albacore was also base-called by the old version or it is new.
AL_N$Condition <- AL_N$read_id %in% AL_O$read_id
if(nrow(AL_N[AL_N$Condition==FALSE,])!=0){
  AL_N[AL_N$Condition==FALSE,]$Condition <- "NEW"}
AL_N[AL_N$Condition==TRUE,]$Condition <- "COMMON"

png(paste(dir, "/Albacore_2.3.4.png", sep=""))
marginal_plot(x = read_length, y = mean_quality, group = Condition, data = AL_N, alpha = 0.5, bw = "nrd", lm_formula = NULL, 
              xlab = "Read Length (log)", ylab = "Mean Quality (log)", pch = 15, cex = 0.5, log="xy")
dev.off()


# Guppy

GU_N$Condition <- GU_N$read_id %in% GU_O$read_id
if(nrow(GU_N[GU_N$Condition==FALSE,])!=0){GU_N[GU_N$Condition==FALSE,]$Condition <- "NEW"}
GU_N[GU_N$Condition==TRUE,]$Condition <- "COMMON"

png(paste(dir, "/Guppy_3.0.3.png", sep=""))
marginal_plot(x = read_length, y = mean_quality, group = Condition, data = GU_N, alpha = 0.5, bw = "nrd", lm_formula = NULL, 
              xlab = "Read Length (log)", ylab = "Mean Quality (log)", pch = 15, cex = 0.5, log="xy")
dev.off()


#Short reads in GU 3.0.3
  
GU_N_short <- GU_N[GU_N$read_length<10,]
