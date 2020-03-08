#####################################################################################
#This script reads from different directories .STATS and .mismatches files, one of each
#per base-caller. In total, 8 files should be in each given directory
######################################################################################


## Loading libraries
library(lattice)
library(optparse)

## Parsing

parser <- OptionParser()
parser <- add_option(parser, opt_str=c("-i", "--input"), type="character",
                      dest='input',
                      help="Input Directories")
parser <- add_option(parser, opt_str=c("-n", "--names"), type="character",
                    dest='names',
                    help="Base-caller names")
parser <- add_option(parser, opt_str=c("-d", "--datasets"), type="character",
                     dest='datasets',
                     help="Data set names")


options=parse_args(parser)
dir=options$input
names=options$names
datasets=options$datasets

#dir="minimap2/m6A,minimap2/m5C"
dirs <- strsplit(dir, ",")[[1]]
datasets <- strsplit(datasets, ",")[[1]]


# In the input directories there should be one file ending in .mismatches and one in .STATS per base-caller.
preprocessing <- function(dir){
  AL_2.1.7_S <- list.files(path=dir, pattern = "*albacore2.1.7.*STATS")
  AL_2.1.7_M <- list.files(path=dir, pattern = "albacore2.1.7.*mismatches")

  AL_2.3.4_S <- list.files(path=dir, pattern = "albacore2.3.4.*STATS")
  AL_2.3.4_M <- list.files(path=dir, pattern = "albacore2.3.4.*mismatches")
  
  GU_2.3.1_S <- list.files(path=dir, pattern = "guppy2.3.1.*STATS")
  GU_2.3.1_M <- list.files(path=dir, pattern = "guppy2.3.1.*mismatches")
  
  GU_3.0.3_S <- list.files(path=dir, pattern = "guppy3.0.3.*STATS")
  GU_3.0.3_M <- list.files(path=dir, pattern = "guppy3.0.3.*mismatches")
  
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
  
  AL_2.1.7_plain <- AL_2.1.7
  AL_2.3.4_plain <- AL_2.3.4
  GU_2.3.1_plain <- GU_2.3.1
  GU_3.0.3_plain <- GU_3.0.3
  
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
  return(data)
}


#matrices

constructing_matrix <- function(data){
  A <- data[data$ref_nuc=="A",]
  C <- data[data$ref_nuc=="C",]
  G <- data[data$ref_nuc=="G",]
  T <- data[data$ref_nuc=="T",]
  a <- log(A[A$Basecaller=="AL 2.1.7",]$mismatches_freq)
  a <- a[a!=-Inf]
  b <- log(C[C$Basecaller=="AL 2.1.7",]$mismatches_freq)
  b <- b[b!=-Inf]
  c <- log(G[G$Basecaller=="AL 2.1.7",]$mismatches_freq)
  c <- c[c!=-Inf]
  d <- log(T[T$Basecaller=="AL 2.1.7",]$mismatches_freq)
  d <- d[d!=-Inf]
  e <- log(A[A$Basecaller=="AL 2.3.4",]$mismatches_freq)
  e <- e[e!=-Inf]
  f <- log(C[C$Basecaller=="AL 2.3.4",]$mismatches_freq)
  f <- f[f!=-Inf]
  g <- log(G[G$Basecaller=="AL 2.3.4",]$mismatches_freq)
  g <- g[g!=-Inf]
  h <- log(T[T$Basecaller=="AL 2.3.4",]$mismatches_freq)
  h <- h[h!=-Inf]
  i <- log(A[A$Basecaller=="GU 2.3.1",]$mismatches_freq)
  i <- i[i!=-Inf]
  j <- log(C[C$Basecaller=="GU 2.3.1",]$mismatches_freq)
  j <- j[j!=-Inf]
  k <- log(G[G$Basecaller=="GU 2.3.1",]$mismatches_freq)
  k <- k[k!=-Inf]
  l <- log(T[T$Basecaller=="GU 2.3.1",]$mismatches_freq)
  l <- l[l!=-Inf]
  m <- log(A[A$Basecaller=="GU 3.0.3",]$mismatches_freq)
  m <- m[m!=-Inf]
  n <- log(C[C$Basecaller=="GU 3.0.3",]$mismatches_freq)
  n <- n[n!=-Inf]
  o <- log(G[G$Basecaller=="GU 3.0.3",]$mismatches_freq)
  a <- a[a!=-Inf]
  p <- log(T[T$Basecaller=="GU 3.0.3",]$mismatches_freq)
  o <- o[o!=-Inf]
  test <- matrix(c(median(a),median(b),
                   median(c),median(d),
                   median(e),median(f),
                   median(g),median(h),
                   median(i),median(j),
                   median(k),median(l),
                   median(m),median(n),
                   median(o),median(p)), 4, 4)
  colnames(test) <- c("AL 2.1.7", "AL 2.3.4", "GU 2.3.1", "GU 3.0.3")
  rownames(test) <- c("A","C","G","T")
  return(test)
}

matrix <- matrix(1:4, nrow=1)
for (i in dirs){
  data <- preprocessing(i)
  small_matrix <- constructing_matrix(data)
  matrix <- rbind(matrix, small_matrix)
}
matrix <- matrix[-1,]

g <-levelplot(matrix, xlab="", ylab="", margin=FALSE, 
          col.regions = heat.colors(100)[length(heat.colors(90)):1], 
          scales = list(tck = c(0,0)), panel=function(...){
            panel.levelplot(...) 
            panel.abline(v=4.5)
            panel.abline(v=8.5)
            panel.abline(v=12.5)
            panel.abline(v=16.5)
            panel.abline(v=20.5)}, main=datasets)
pdf(paste("levelplot.pdf"))
print(g)
dev.off()
