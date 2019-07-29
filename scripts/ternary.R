## Loading libraries
library(ggplot2)
library(ggpubr)
library(optparse)
library(reshape2)
library(ggtern)
library(lattice)
library(latticeExtra)

## Parsing

parser <- OptionParser()
parser <- add_option(parser, opt_str=c("-m", "--mapper"), type="character",
                     dest='mapper',
                     help="mapper")
parser <- add_option(parser, opt_str=c("-b", "--basecaller"), type="character",
                     dest='basecaller',
                     help="basecaller")
options=parse_args(parser)
mapper=options$mapper
basecaller=options$basecaller


#mapper="graphmap"
#basecaller="guppy2.3.1" 

UNM_S <- read.delim(paste(mapper, "UNM_C/RNA081120181_", basecaller, ".sorted.bam.STATS", sep="/"))
UNM_M <- read.delim(paste(mapper, "UNM_C/RNA081120181_", basecaller, ".sorted.bam.mismatches", sep="/"))
m5C_S <- read.delim(paste(mapper, "m5C/RNA010220191_", basecaller, ".sorted.bam.STATS", sep="/"))
m5C_M <- read.delim(paste(mapper, "m5C/RNA010220191_", basecaller, ".sorted.bam.mismatches", sep="/"))
hm5C_S <- read.delim(paste(mapper, "5hmC/RNA080420_", basecaller, ".sorted.bam.STATS", sep="/"))
hm5C_M <- read.delim(paste(mapper, "5hmC/RNA080420_", basecaller, ".sorted.bam.mismatches", sep="/"))

names <- c("UNM", "m5C", "hm5C")
for (i in 1:length(names)){
  assign(names[i], merge(get(paste(names[i],"_S",sep="")), 
                         get(paste(names[i],"_M",sep="")), 
                         by = c("chr", "pos", "ref_nuc"), 
                         sort = F))}


# We remove  positions with low coverage
l <- list(UNM, m5C, hm5C)
remove_low_coverage <- function(i){
  cov <- mean(i$coverage) # We compute the coverage mean
  threshold <- cov*0.1 # We create a threshold for removing values lower than 10% of the computed mean
  i <- i[i$coverage>threshold,]}
l <- lapply(l, remove_low_coverage)

# We create a function for computing the total number of mismatches
mismatches_function <- function(k){ 
  mismatches <- c()  
  for (i in c(1:nrow(k))){   
    base <- k[i,3]
    a <- sum(k[i,10:13])-k[i,toString(base)]
    mismatches <- c(mismatches, a)}
  k <- cbind(k, mismatches)}
l <- lapply(l, mismatches_function)

UNM <- l[[1]]
m5C <- l[[2]]
hm5C <- l[[3]]

# We compute the mismatch frequency
UNM$mismatches_freq <- UNM$mismatches/UNM$coverage
m5C$mismatches_freq <- m5C$mismatches/m5C$coverage
hm5C$mismatches_freq <- hm5C$mismatches/hm5C$coverage

# We add a variable to identify the modification
UNM$modification <- "UNM"
m5C$modification <- "m5C"
hm5C$modification <- "hm5C"



C_UNM <- c()
for (i in 3:(nrow(UNM)-2)){
  if(UNM[i,]$ref_nuc=="C" && UNM[i-1,]$ref_nuc!="C" && UNM[i-2,]$ref_nuc!="C" && UNM[i+1,]$ref_nuc!="C" && UNM[i+2,]$ref_nuc!="C"){
    C_UNM <- rbind(C_UNM, UNM[i,])}}

C_m5C <- c()
for (i in 3:(nrow(m5C)-2)){
  if(m5C[i,]$ref_nuc=="C" && m5C[i-1,]$ref_nuc!="C" && m5C[i-2,]$ref_nuc!="C" && m5C[i+1,]$ref_nuc!="C" && m5C[i+2,]$ref_nuc!="C"){
    C_m5C <- rbind(C_m5C, m5C[i,])}}

C_5hmC <- c()
for (i in 3:(nrow(hm5C)-2)){
  if(hm5C[i,]$ref_nuc=="C" && hm5C[i-1,]$ref_nuc!="C" && hm5C[i-2,]$ref_nuc!="C" && hm5C[i+1,]$ref_nuc!="C" && hm5C[i+2,]$ref_nuc!="C"){
    C_5hmC <- rbind(C_5hmC, hm5C[i,])}}


g1 <- ggtern(C_UNM, aes(A,G,T)) +
  geom_point(size = 0.05, aes(color = log(coverage))) +
  ggtitle("Ternary Diagram for unmodified C") + theme_bw(base_size = 6) +
  scale_color_continuous("Coverage (log)", low = 'cyan2', high = 'blue4') +
  theme(plot.title = element_text(size = 8, hjust = 0.5))
legend <- get_legend(g1)
g1 <- g1 + theme(legend.position="none")
g2 <- ggtern(C_m5C, aes(A,G,T)) +
  geom_point(size = 0.05, aes(color = log(coverage))) +
  ggtitle("Ternary Diagram for m5C") + theme_bw(base_size = 6) +
  scale_color_continuous("Coverage (log)", low = 'cyan2', high = 'blue4') +
  theme(legend.position="none", plot.title = element_text(size = 8, hjust = 0.5))

g3 <- ggtern(C_5hmC, aes(A,G,T)) +
  geom_point(size = 0.05, aes(color = log(coverage))) +
  ggtitle("Ternary Diagram for 5hmC") + theme_bw(base_size = 6) +
  scale_color_continuous("Coverage (log)", low = 'cyan2', high = 'blue4') +
  theme(legend.position="none", plot.title = element_text(size = 8, hjust = 0.5))

g<- grid.arrange(g1,g2,g3, nrow=1, ncol=4, legend, top=paste("Ternary Diagrams for Cs"))

ggsave(paste("ternary_diagrams_C", mapper, ".pdf", sep="_"), g, device = "pdf", width = 7, height = 3)











##### m6A


UNM_S <- read.delim(paste(mapper, "UNM_A/RNA081120181_", basecaller, ".sorted.bam.STATS", sep="/"))
UNM_M <- read.delim(paste(mapper, "UNM_A/RNA081120181_", basecaller, ".sorted.bam.mismatches", sep="/"))
m6A_S <- read.delim(paste(mapper, "m6A/RNA081120182_", basecaller, ".sorted.bam.STATS", sep="/"))
m6A_M <- read.delim(paste(mapper, "m6A/RNA081120182_", basecaller, ".sorted.bam.mismatches", sep="/"))

names <- c("UNM", "m6A")
for (i in 1:length(names)){
  assign(names[i], merge(get(paste(names[i],"_S",sep="")), 
                         get(paste(names[i],"_M",sep="")), 
                         by = c("chr", "pos", "ref_nuc"), 
                         sort = F))}

l <- list(UNM, m6A)
l <- lapply(l, remove_low_coverage)
l <- lapply(l, mismatches_function)

UNM <- l[[1]]
m6A <- l[[2]]



mod="A"
C_UNM <- c()
for (i in 3:(nrow(UNM)-2)){
  if(UNM[i,]$ref_nuc==mod && UNM[i-1,]$ref_nuc!=mod && UNM[i-2,]$ref_nuc!=mod && UNM[i+1,]$ref_nuc!=mod && UNM[i+2,]$ref_nuc!=mod){
    C_UNM <- rbind(C_UNM, UNM[i,])}}

C_m6A <- c()
for (i in 3:(nrow(m6A)-2)){
  if(m6A[i,]$ref_nuc==mod && m6A[i-1,]$ref_nuc!=mod && m6A[i-2,]$ref_nuc!=mod && m6A[i+1,]$ref_nuc!=mod && m6A[i+2,]$ref_nuc!=mod){
    C_m6A <- rbind(C_m6A, m6A[i,])}}



g1 <- ggtern(C_UNM, aes(C,G,T)) +
  geom_point(size = 0.05, aes(color = log(coverage))) +
  ggtitle(paste("Ternary Diagram for unmodified", mod)) + theme_bw(base_size = 6) +
  scale_color_continuous("Coverage (log)", low = 'cyan2', high = 'blue4') +
  theme(plot.title = element_text(size = 8, hjust = 0.5))
legend <- get_legend(g1)
g1 <- g1 + theme(legend.position="none")
g2 <- ggtern(C_m6A, aes(C,G,T)) +
  geom_point(size = 0.05, aes(color = log(coverage))) +
  ggtitle("Ternary Diagram for m6A") + theme_bw(base_size = 6) +
  scale_color_continuous("Coverage (log)", low = 'cyan2', high = 'blue4') +
  theme(legend.position="none", plot.title = element_text(size = 8, hjust = 0.5))

g<- grid.arrange(g1,g2, nrow=1, ncol=3, legend, top=paste("Ternary Diagrams for As"))

ggsave(paste("ternary_diagrams_A", mapper, ".pdf", sep="_"), g, device = "pdf", width = 7, height = 3)









##### pU


UNM_S <- read.delim(paste(mapper, "UNM_T/RNA081120181_", basecaller, ".sorted.bam.STATS", sep="/"))
UNM_M <- read.delim(paste(mapper, "UNM_T/RNA081120181_", basecaller, ".sorted.bam.mismatches", sep="/"))
pU_S <- read.delim(paste(mapper, "pU/RNAAB091808_", basecaller, ".sorted.bam.STATS", sep="/"))
pU_M <- read.delim(paste(mapper, "pU/RNAAB091808_", basecaller, ".sorted.bam.mismatches", sep="/"))

names <- c("UNM", "pU")
for (i in 1:length(names)){
  assign(names[i], merge(get(paste(names[i],"_S",sep="")), 
                         get(paste(names[i],"_M",sep="")), 
                         by = c("chr", "pos", "ref_nuc"), 
                         sort = F))}

l <- list(UNM, pU)
l <- lapply(l, remove_low_coverage)
l <- lapply(l, mismatches_function)

UNM <- l[[1]]
pU <- l[[2]]


mod="T"
C_UNM <- c()
for (i in 3:(nrow(UNM)-2)){
  if(UNM[i,]$ref_nuc==mod && UNM[i-1,]$ref_nuc!=mod && UNM[i-2,]$ref_nuc!=mod && UNM[i+1,]$ref_nuc!=mod && UNM[i+2,]$ref_nuc!=mod){
    C_UNM <- rbind(C_UNM, UNM[i,])}}

C_pU <- c()
for (i in 3:(nrow(pU)-2)){
  if(pU[i,]$ref_nuc==mod && pU[i-1,]$ref_nuc!=mod && pU[i-2,]$ref_nuc!=mod && pU[i+1,]$ref_nuc!=mod && pU[i+2,]$ref_nuc!=mod){
    C_pU <- rbind(C_pU, pU[i,])}}

g1 <- ggtern(C_UNM, aes(A,C,G)) +
  geom_point(size = 0.05, aes(color = log(coverage))) +
  ggtitle(paste("Ternary Diagram for unmodified U")) + theme_bw(base_size = 6) +
  scale_color_continuous("Coverage (log)", low = 'cyan2', high = 'blue4') +
  theme(plot.title = element_text(size = 8, hjust = 0.5))
legend <- get_legend(g1)
g1 <- g1 + theme(legend.position="none")

g2 <- ggtern(C_pU, aes(A,C,G)) +
  geom_point(size = 0.05, aes(color = log(coverage))) +
  ggtitle("Ternary Diagram for pU") + theme_bw(base_size = 6) +
  scale_color_continuous("Coverage (log)", low = 'cyan2', high = 'blue4') +
  theme(legend.position="none", plot.title = element_text(size = 8, hjust = 0.5))

g<- grid.arrange(g1,g2, nrow=1, ncol=3, legend, top=paste("Ternary Diagrams for Us"))

ggsave(paste("ternary_diagrams_U", mapper, ".pdf", sep="_"), g, device = "pdf", width = 7, height = 3)

