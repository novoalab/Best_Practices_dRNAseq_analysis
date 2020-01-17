#This script is placed in a folder with a file called sites.bed and two subfolders called nanopolish and tailfindr

### Libraries
# We install the libraries that are not already installed
packages <- c("ggplot2", "grDevices", "dplyr", "lattice", "ggpubr", "limma", "seqsetvis") 
new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

# We load them
lapply(packages, library, character.only = TRUE)


### Parsing

clean_dataframe_t <- function(d, n, t){ #We create a function for cleaning the tailfindr inputs
  d <- na.omit(d)
  d <- d[,c(1,2,6)]
  colnames(d) <- c("readname", "gene_name", "polya_length")
  d$groups <- "tailfindr"
  d$rep <- n
  d$type <- t
  d <- distinct(d)
  return(d)}

for(a in c("ko1", "ko2", "ko3", "wt1", "wt2", "wt3")){
  assign(paste(a, "_tailfindr", sep=""), read.delim(paste("tailfindr/", toupper(a), "_GU_3.0.3_tailfindr_gene_names.tsv", sep="")))
  assign(paste(a, "_tailfindr", sep=""), clean_dataframe_t(get(paste(a, "_tailfindr", sep="")), paste(a, "_tailfindr", sep=""), substr(a,1,2)))
  assign(paste("n", a, "_tailfindr", sep=""), nrow(get(paste(a, "_tailfindr", sep=""))))
}

clean_dataframe_n <- function(d, n, t){ # We create a function for cleaning the nanopolish inputs
  d[,1] <- sub("_Basecall_.*", "", d[,1])
  d <- d[,c(1,9,10)]
  d <- na.omit(d)
  d$groups <- "nanopolish"
  d$rep <- n
  d$type <- t
  d <- distinct(d)
  return(d)}

for(a in c("ko1", "ko2", "ko3", "wt1", "wt2", "wt3")){
  assign(paste(a, "_nanopolish", sep=""), read.delim(paste("nanopolish/", toupper(a), "_GU_3.0.3_nanopolish_gene_names.tsv", sep="")))
  assign(paste(a, "_nanopolish", sep=""), clean_dataframe_n(get(paste(a, "_nanopolish", sep="")), paste(a, "_nanopolish", sep=""), substr(a,1,2)))
  assign(paste("n", a, "_nanopolish", sep=""), nrow(get(paste(a, "_nanopolish", sep=""))))
}


#Plots

dir.create("plots")
for (i in c("general_stats", "Nanopolish", "Tailfindr", "Nanopolsih_vs_Tailfindr", "final_plots", "m6A_vs_not_m6A")){
  dir.create(paste("plots", i, sep="/"))}

nanopolish <- rbind(ko1_nanopolish, ko2_nanopolish, ko3_nanopolish, wt1_nanopolish, wt2_nanopolish, wt3_nanopolish)
g <- ggplot(nanopolish, aes(rep, log(polya_length), fill=type, col=type)) + geom_boxplot(alpha=0.4, show.legend = FALSE) + theme_bw()+ xlab("") +
  scale_fill_manual(values=c("cyan", "green")) + scale_color_manual(values=c("cyan4", "green4")) + ggtitle("Nanopolish") +
  scale_y_continuous(breaks=c(-2.5,0,2.5,5.0,7.5), labels=c(round(exp(-2.5),0), round(exp(0),0), round(exp(2.5),0), round(exp(5),0), round(exp(7.55),0))) +
  ylab("polya length (log)") + scale_x_discrete(labels=c("ko1_nanopolish" = "KO1", "ko2_nanopolish" = "KO2", "ko3_nanopolish" = "KO3", "wt1_nanopolish" = "WT1", "wt2_nanopolish" = "WT2", "wt3_nanopolish" = "WT3"))
ggsave("plots/general_stats/nanopolish_boxplots.pdf", g, device = "pdf", width = 5, height = 4)

g <- ggplot(nanopolish, aes(type, log(polya_length), fill=type, col=type)) + geom_boxplot(alpha=0.4, show.legend = FALSE) + theme_bw()+ xlab("") +
  scale_fill_manual(values=c("cyan", "green")) + scale_color_manual(values=c("cyan4", "green4")) + ggtitle("Nanopolish") +
  scale_y_continuous(breaks=c(-2.5,0,2.5,5.0,7.5), labels=c(round(exp(-2.5),0), round(exp(0),0), round(exp(2.5),0), round(exp(5),0), round(exp(7.55),0))) +
  ylab("polya length (log)") + scale_x_discrete(labels=c("ko1_nanopolish" = "KO1", "ko2_nanopolish" = "KO2", "ko3_nanopolish" = "KO3", "wt1_nanopolish" = "WT1", "wt2_nanopolish" = "WT2", "wt3_nanopolish" = "WT3"))
ggsave("plots/general_stats/nanopolish_wt_vs_ko_boxplots.pdf", g, device = "pdf", width = 5, height = 4)

tailfindr <- rbind(ko1_tailfindr, ko2_tailfindr, ko3_tailfindr, wt1_tailfindr, wt2_tailfindr, wt3_tailfindr)
g <- ggplot(tailfindr, aes(rep, log(polya_length), fill=type, col=type)) + geom_boxplot(alpha=0.4, show.legend = FALSE) + theme_bw()+ xlab("") +
  scale_fill_manual(values=c("cyan", "green")) + scale_color_manual(values=c("cyan4", "green4")) + ggtitle("Tailfindr") +
  scale_y_continuous(breaks=c(0,2,4,6,8), labels=c(round(exp(0),0), round(exp(2),0), round(exp(4),0), round(exp(6),0), round(exp(8),0))) +
  ylab("polya length (log)") + scale_x_discrete(labels=c("ko1_tailfindr" = "KO1", "ko2_tailfindr" = "KO2", "ko3_tailfindr" = "KO3", "wt1_tailfindr" = "WT1", "wt2_tailfindr" = "WT2", "wt3_tailfindr" = "WT3"))
ggsave("plots/general_stats/tailfindr_boxplots.pdf", g, device = "pdf", width = 5, height = 4)

g <- ggplot(tailfindr, aes(type, log(polya_length), fill=type, col=type)) + geom_boxplot(alpha=0.4, show.legend = FALSE) + theme_bw()+ xlab("") +
  scale_fill_manual(values=c("cyan", "green")) + scale_color_manual(values=c("cyan4", "green4")) + ggtitle("Tailfindr") +
  scale_y_continuous(breaks=c(0,2,4,6,8), labels=c(round(exp(0),0), round(exp(2),0), round(exp(4),0), round(exp(6),0), round(exp(8),0))) +
  ylab("polya length (log)") + scale_x_discrete(labels=c("ko1_tailfindr" = "KO1", "ko2_tailfindr" = "KO2", "ko3_tailfindr" = "KO3", "wt1_tailfindr" = "WT1", "wt2_tailfindr" = "WT2", "wt3_tailfindr" = "WT3"))
ggsave("plots/general_stats/tailfindr_wt_vs_ko_boxplots.pdf", g, device = "pdf", width = 5, height = 4)
  

genes <- read.table("sites.bed")
genes <- genes$V6

nanopolish_filtered <- nanopolish[nanopolish$gene_name %in% genes,]
tailfindr_filtered <- tailfindr[tailfindr$gene_name %in% genes,]

g <- ggplot(nanopolish_filtered, aes(rep, log(polya_length), fill=type, col=type)) + geom_boxplot(alpha=0.4,show.legend = FALSE) + theme_bw()+ xlab("") +
  scale_fill_manual(values=c("cyan", "green")) + scale_color_manual(values=c("cyan4", "green4")) + ggtitle("Nanopolish 1308") +
  scale_y_continuous(breaks=c(-2.5,0,2.5,5.0,7.5), labels=c(round(exp(-2.5),0), round(exp(0),0), round(exp(2.5),0), round(exp(5),0), round(exp(7.55),0))) +
  ylab("polya length (log)") + scale_x_discrete(labels=c("ko1_nanopolish" = "KO1", "ko2_nanopolish" = "KO2", "ko3_nanopolish" = "KO3", "wt1_nanopolish" = "WT1", "wt2_nanopolish" = "WT2", "wt3_nanopolish" = "WT3"))
ggsave("plots/general_stats/nanopolish_boxplots_1308.pdf", g, device = "pdf", width = 5, height = 4)


g <- ggplot(tailfindr_filtered, aes(rep, log(polya_length), fill=type, col=type)) + geom_boxplot(alpha=0.4, show.legend = FALSE) + theme_bw()+ xlab("") +
  scale_fill_manual(values=c("cyan", "green")) + scale_color_manual(values=c("cyan4", "green4")) + ggtitle("Tailfindr 1308") +
  scale_y_continuous(breaks=c(2,4,6), labels=c(round(exp(2),0), round(exp(4),0), round(exp(6),0))) +
  ylab("polya length (log)") + scale_x_discrete(labels=c("ko1_tailfindr" = "KO1", "ko2_tailfindr" = "KO2", "ko3_tailfindr" = "KO3", "wt1_tailfindr" = "WT1", "wt2_tailfindr" = "WT2", "wt3_tailfindr" = "WT3"))
ggsave("plots/general_stats/tailfindr_boxplots_1308.pdf", g, device = "pdf", width = 5, height = 4)



#replicability:

data <- rbind(ko1_nanopolish, ko1_tailfindr, ko2_nanopolish, ko2_tailfindr, ko3_nanopolish, ko3_tailfindr, wt1_nanopolish, wt1_tailfindr, wt2_nanopolish, wt2_tailfindr, wt3_nanopolish, wt3_tailfindr)
data <- rbind(nanopolish, tailfindr)
#barplot comparing number of reads that have had an assigned tail length
g <- ggplot(data, aes(rep, log(polya_length), color=groups, fill=groups)) + 
  geom_boxplot(alpha=0.3) + theme_bw() + xlab("") +
  ylab("polya length (log)") + scale_x_discrete(labels=c("ko1_nanopolish" = "KO1", "ko1_tailfindr" = "KO1", "ko2_nanopolish" = "KO2", "ko2_tailfindr" = "KO2", "ko3_nanopolish" = "KO3", "ko3_tailfindr" = "KO3", "wt1_nanopolish" = "WT1", "wt1_tailfindr" = "WT1", "wt2_nanopolish" = "WT2", "wt2_tailfindr" = "WT2", "wt3_nanopolish" = "WT3", "wt3_tailfindr" = "WT3"))
g <- g + scale_y_continuous(breaks=c(as.numeric(ggplot_build(g)$layout$panel_params[[1]]$y.labels)), labels=c(round(exp(as.numeric(ggplot_build(g)$layout$panel_params[[1]]$y.labels)),0)))
ggsave("plots/Nanopolish_vs_Tailfindr/nanopolish_vs_tailfindr_boxplots.pdf", g, device = "pdf", width = 6, height = 4)


merged_ko2 <- merge(ko2_tailfindr, ko2_nanopolish, by = "readname")
nmerged2 <- nrow(merged_ko2)

#e <-cor(merged_ko2$polya_length.x, merged_ko2$polya_length.y, method = "pearson")
#e <- paste("r=", round(e,3), sep="")

# g <- ggplot(merged_ko2, aes(log(polya_length.x), log(polya_length.y))) + theme_bw() +
#    geom_point(col=densCols(merged_ko2$polya_length.x, merged_ko2$polya_length.y, colramp = colorRampPalette(blues9[-(1:3)]))) +
#    geom_text(aes(x=-Inf, y=+Inf, label=e), color="black", size=5, hjust = -0.5, vjust = 2) +
#    xlab("tailfindr") + ylab("nanopolish") + ggtitle("kO2 Tail Lengths") 
# m <- min(layer_scales(g)$x$range$range[[1]], layer_scales(g)$y$range$range[[1]])
# M <- max(layer_scales(g)$x$range$range[[2]], layer_scales(g)$y$range$range[[2]])
# g <- g + expand_limits(x=c(m,M), y=c(m,M))
# g <- g + scale_y_continuous(breaks=c(as.numeric(ggplot_build(g)$layout$panel_params[[1]]$y.labels)), labels=c(round(exp(as.numeric(ggplot_build(g)$layout$panel_params[[1]]$y.labels)),0))) +
#   scale_x_continuous(breaks=c(as.numeric(ggplot_build(g)$layout$panel_params[[1]]$x.labels)), labels=c(round(exp(as.numeric(ggplot_build(g)$layout$panel_params[[1]]$x.labels)),0)))
#ggsave("nanopolish_vs_tailfindr_KO2.pdf", g, device = "pdf", width = 4, height = 4)


merged_ko3 <- merge(ko3_tailfindr, ko3_nanopolish, by = "readname")
nmerged3 <- nrow(merged_ko3)

merged_ko1 <- merge(ko1_tailfindr, ko1_nanopolish, by = "readname")
nmerged1 <- nrow(merged_ko1)
merged_wt1 <- merge(wt1_tailfindr, wt1_nanopolish, by = "readname")
nmerged4 <- nrow(merged_wt1)
merged_wt2 <- merge(wt2_tailfindr, wt2_nanopolish, by = "readname")
nmerged5 <- nrow(merged_wt2)
merged_wt3 <- merge(wt3_tailfindr, wt3_nanopolish, by = "readname")
nmerged6 <- nrow(merged_wt3)

numbers <- rbind(nko1_tailfindr, nko1_nanopolish, nmerged1, nko2_tailfindr, nko2_nanopolish, nmerged2, nko3_tailfindr, nko3_nanopolish, nmerged3, nwt1_tailfindr, nwt1_nanopolish, nmerged4, nwt2_tailfindr, nwt2_nanopolish, nmerged5, nwt3_tailfindr, nwt3_nanopolish, nmerged6)
names <- rbind("ko1_t", "ko1_n", "ko1_zc", "ko2_t", "ko2_n", "ko2_zc", "ko3_t", "ko3_n", "ko3_zc", "wt1_t", "wt1_n", "wt1_zc", "wt2_t", "wt2_n", "wt2_zc", "wt3_t", "wt3_n", "wt3_zc")
groups <- cbind(rep(c("tailfindr","nanopolish","common"), times = 6))
numbersdata <- data.frame(numbers, names, groups)
numbersdata$groups <- as.factor(numbersdata$groups)
g <- ggplot(numbersdata, aes(names, numbers, fill=groups)) + geom_col() + theme_bw() + xlab("") + ylab("Number of estimated tails") +
  scale_x_discrete(labels=c("ko1_t" = "KO1", "ko1_n" = " ", "ko1_zc" = " ", "ko2_t" = "KO2", "ko2_n" = " ", "ko2_zc" = " ", "ko3_t" = "KO3", "ko3_n" = " ", "ko3_zc" = " ", "wt1_t" = "WT1", "wt1_n" = " ", "wt1_zc" = " ", "wt2_t" = "WT2", "wt2_n" = " ", "wt2_zc" = " ", "wt3_t" = "WT3", "wt3_n" = " ", "wt3_zc" = " ")) +
  theme(axis.text.x = element_text(size = 12))
ggsave("plots/general_stats/estimated_tails.pdf", g, device = "pdf", width = 5, height = 4)


#Replicability in Nanopolish

ko2_medians <- aggregate(polya_length~gene_name, ko2_nanopolish, median)
ko3_medians <- aggregate(polya_length~gene_name, ko3_nanopolish, median)
merged <- merge(ko2_medians, ko3_medians, by = "gene_name")

e <-cor(merged$polya_length.x, merged$polya_length.y, method = "pearson")
e <- paste("r=", round(e,3), sep="")
g <- ggplot(merged, aes(polya_length.x, polya_length.y)) + geom_point(col=densCols(merged$polya_length.x, merged$polya_length.y, colramp = colorRampPalette(blues9[-(1:3)]))) +
  theme_bw() + ggtitle("ko2 vs ko3 in Nanopolish") + xlab("polyA tail lengths in KO2") + ylab("polyA tail lengths in KO3") +
  geom_text(aes(x=-Inf, y=+Inf, label=e), color="black", size=5, hjust = -0.5, vjust = 2) 
ggsave("plots/general_stats/replicability_nanopolish_KO2_KO3.pdf", g, device = "pdf", width = 5, height = 4)

  
filtering <- function(input, n){
  test <- input %>% 
    group_by(gene_name) %>% 
    filter(n() >= n)
  return(test)
}

#min 30


ko2_filtered <- filtering(ko2_nanopolish, 30)
ko3_filtered <- filtering(ko3_nanopolish, 30)

ko2_medians <- aggregate(polya_length~gene_name, ko2_filtered, median)
ko3_medians <- aggregate(polya_length~gene_name, ko3_filtered, median)

merged <- merge(ko2_medians, ko3_medians, by = "gene_name")

e <-cor(merged$polya_length.x, merged$polya_length.y, method = "pearson")
e <- paste("r=", round(e,3), sep="")
g <- ggplot(merged, aes(polya_length.x, polya_length.y)) + geom_point(col=densCols(merged$polya_length.x, merged$polya_length.y, colramp = colorRampPalette(blues9[-(1:3)]))) +
  theme_bw() + ggtitle("KO2 vs KO3 min30 reads/gene") + xlab("polyA tail lengths in KO2") + ylab("polyA tail lengths in KO3") +
  geom_text(aes(x=-Inf, y=+Inf, label=e), color="black", size=5, hjust = -0.5, vjust = 2)  
ggsave("plots/general_stats/replicability_nanopolish_KO2_KO3_min30.pdf", g, device = "pdf", width = 5, height = 4)


# min 50

ko2_filtered <- filtering(ko2_nanopolish, 50)
nko2_filtered <- nrow(ko2_filtered)
ko3_filtered <- filtering(ko3_nanopolish, 50)
nko3_filtered <- nrow(ko3_filtered)

ko2_medians <- aggregate(polya_length~gene_name, ko2_filtered, median)
ko3_medians <- aggregate(polya_length~gene_name, ko3_filtered, median)

merged <- merge(ko2_medians, ko3_medians, by = "gene_name")

e <-cor(merged$polya_length.x, merged$polya_length.y, method = "pearson")
e <- paste("r=", round(e,3), sep="")
g <- ggplot(merged, aes(polya_length.x, polya_length.y)) + geom_point(col=densCols(merged$polya_length.x, merged$polya_length.y, colramp = colorRampPalette(blues9[-(1:3)]))) +
  theme_bw() + ggtitle("KO2 vs KO3 min50 reads/gene") + xlab("polyA tail lengths in KO2") + ylab("polyA tail lengths in KO3") +
  geom_text(aes(x=-Inf, y=+Inf, label=e), color="black", size=5, hjust = -0.5, vjust = 2) 
ggsave("plots/general_stats/replicability_nanopolish_KO2_KO3_min50.pdf", g, device = "pdf", width = 5, height = 4)


ko1_filtered <- filtering(ko1_nanopolish, 50)
nko1_filtered <- nrow(ko1_filtered)
wt1_filtered <- filtering(wt1_nanopolish, 50)
nwt1_filtered <- nrow(wt1_filtered)
wt2_filtered <- filtering(wt2_nanopolish, 50)
nwt2_filtered <- nrow(wt2_filtered)
wt3_filtered <- filtering(wt3_nanopolish, 50)
nwt3_filtered <- nrow(wt3_filtered)

ko1_filtered_t <- filtering(ko1_tailfindr, 50)
nko1_filtered_t <- nrow(ko1_filtered_t)
ko2_filtered_t <- filtering(ko2_tailfindr, 50)
nko2_filtered_t <- nrow(ko2_filtered_t)
ko3_filtered_t <- filtering(ko3_tailfindr, 50)
nko3_filtered_t <- nrow(ko3_filtered_t)
wt1_filtered_t <- filtering(wt1_tailfindr, 50)
nwt1_filtered_t <- nrow(wt1_filtered_t)
wt2_filtered_t <- filtering(wt2_tailfindr, 50)
nwt2_filtered_t <- nrow(wt2_filtered_t)
wt3_filtered_t <- filtering(wt3_tailfindr, 50)
nwt3_filtered_t <- nrow(wt3_filtered_t)

######################
nanopolish_filtered <- rbind(ko1_filtered, ko2_filtered, ko3_filtered, wt1_filtered, wt2_filtered, wt3_filtered)
g <- ggplot(nanopolish_filtered, aes(rep, log(polya_length), fill=type, col=type)) + geom_boxplot(alpha=0.4, show.legend = FALSE) + theme_bw()+ xlab("") +
  scale_fill_manual(values=c("cyan", "green")) + scale_color_manual(values=c("cyan4", "green4")) + ggtitle("Nanopolish min 50 reads/gene") + ylab("polyA length (log)") +
  scale_x_discrete(labels=c("ko1_nanopolish" = "KO1", "ko2_nanopolish" = "KO2", "ko3_nanopolish" = "KO3", "wt1_nanopolish" = "WT1", "wt2_nanopolish" = "WT2", "wt3_nanopolish" = "WT3"))
g <- g + scale_y_continuous(breaks=c(as.numeric(ggplot_build(g)$layout$panel_params[[1]]$y.labels)), labels=c(round(exp(as.numeric(ggplot_build(g)$layout$panel_params[[1]]$y.labels)),0)))
ggsave("plots/general_stats/nanopolish_boxplots_min50.pdf", g, device = "pdf", width = 5, height = 4)


tailfindr_filtered <- rbind(ko1_filtered_t, ko2_filtered_t, ko3_filtered_t, wt1_filtered_t, wt2_filtered_t, wt3_filtered_t)
g <- ggplot(tailfindr_filtered, aes(rep, log(polya_length), fill=type, col=type)) + geom_boxplot(alpha=0.4, show.legend = FALSE) + theme_bw()+ xlab("") +
  scale_fill_manual(values=c("cyan", "green")) + scale_color_manual(values=c("cyan4", "green4")) + ggtitle("Tailfindr min 50 reads/gene") +
  scale_x_discrete(labels=c("ko1_tailfindr" = "KO1", "ko2_tailfindr" = "KO2", "ko3_tailfindr" = "KO3", "wt1_tailfindr" = "WT1", "wt2_tailfindr" = "WT2", "wt3_tailfindr" = "WT3"))
g <- g + scale_y_continuous(breaks=c(as.numeric(ggplot_build(g)$layout$panel_params[[1]]$y.labels)), labels=c(round(exp(as.numeric(ggplot_build(g)$layout$panel_params[[1]]$y.labels)),0)))
ggsave("plots/general_stats/tailfindr_boxplots_min50.pdf", g, device = "pdf", width = 5, height = 4)

#now the table but with genes that have more than 50 reads. 
nko1_tailfindr_filtered <- nrow(ko1_filtered_t)
merged_ko1 <- merge(ko1_filtered_t, ko1_filtered, by = "readname")
nmerged1 <- nrow(merged_ko1)

nko2_tailfindr_filtered <- nrow(ko2_filtered_t)
merged_ko2 <- merge(ko2_filtered_t, ko2_filtered, by = "readname")
nmerged2 <- nrow(merged_ko2)
nko3_tailfindr_filtered <- nrow(ko3_filtered_t)
merged_ko3 <- merge(ko3_filtered_t, ko3_filtered, by = "readname")
nmerged3 <- nrow(merged_ko3)
nwt1_tailfindr_filtered <- nrow(wt1_filtered_t)
merged_wt1 <- merge(wt1_filtered_t, wt1_filtered, by = "readname")
nmerged4 <- nrow(merged_wt1)
nwt2_tailfindr_filtered <- nrow(wt2_filtered_t)
merged_wt2 <- merge(wt2_filtered_t, wt2_filtered, by = "readname")
nmerged5 <- nrow(merged_wt2)
nwt3_tailfindr_filtered <- nrow(wt3_filtered_t)
merged_wt3 <- merge(wt3_filtered_t, wt3_filtered, by = "readname")
nmerged6 <- nrow(merged_wt3)

numbers <- rbind(nko1_tailfindr_filtered, nko1_filtered, nmerged1, nko2_tailfindr_filtered, nko2_filtered, nmerged2, nko3_tailfindr_filtered, nko3_filtered, nmerged3, nwt1_tailfindr_filtered, nwt1_filtered, nmerged4, nwt2_tailfindr_filtered, nwt2_filtered, nmerged5, nwt3_tailfindr_filtered, nwt3_filtered, nmerged6)
names <- rbind("ko1_t", "ko1_n", "ko1_zc", "ko2_t", "ko2_n", "ko2_zc", "ko3_t", "ko3_n", "ko3_zc", "wt1_t", "wt1_n", "wt1_zc", "wt2_t", "wt2_n", "wt2_zc", "wt3_t", "wt3_n", "wt3_zc")
groups <- cbind(rep(c("tailfindr","nanopolish","common"), times = 6))
numbersdata <- data.frame(numbers, names, groups)
numbersdata$groups <- as.factor(numbersdata$groups)
g <- ggplot(numbersdata, aes(names, numbers, fill=groups)) + geom_col() + theme_bw() + xlab("") + ylab("Number of estimated tails") + ggtitle("min50") +
  scale_x_discrete(labels=c("ko1_t" = "KO1", "ko1_n" = " ", "ko1_zc" = " ", "ko2_t" = "KO2", "ko2_n" = " ", "ko2_zc" = " ", "ko3_t" = "KO3", "ko3_n" = " ", "ko3_zc" = " ", "wt1_t" = "WT1", "wt1_n" = " ", "wt1_zc" = " ", "wt2_t" = "WT2", "wt2_n" = " ", "wt2_zc" = " ", "wt3_t" = "WT3", "wt3_n" = " ", "wt3_zc" = " ")) +
  theme(axis.text.x = element_text(size = 12))
ggsave("plots/general_stats/estimated_tails_min50.pdf", g, device = "pdf", width = 5, height = 4)


nanopolish_filtered_x2 <- nanopolish_filtered[nanopolish_filtered$gene_name %in% genes,]
tailfindr_filtered_x2 <- tailfindr_filtered[tailfindr_filtered$gene_name %in% genes,]

g <- ggplot(nanopolish_filtered_x2, aes(rep, log(polya_length), fill=type, col=type)) + geom_boxplot(alpha=0.4, show.legend = FALSE) + theme_bw()+ xlab("") +
  scale_fill_manual(values=c("cyan", "green")) + scale_color_manual(values=c("cyan4", "green4")) + ggtitle("Nanopolish min 50 from 1308 genes") +
  ylab("polya length (log)") + scale_x_discrete(labels=c("ko1_nanopolish" = "KO1", "ko2_nanopolish" = "KO2", "ko3_nanopolish" = "KO3", "wt1_nanopolish" = "WT1", "wt2_nanopolish" = "WT2", "wt3_nanopolish" = "WT3"))
g <- g + scale_y_continuous(breaks=c(as.numeric(ggplot_build(g)$layout$panel_params[[1]]$y.labels)), labels=c(round(exp(as.numeric(ggplot_build(g)$layout$panel_params[[1]]$y.labels)),0)))
ggsave("plots/general_stats/nanopolish_boxplots_min50_1308.pdf", g, device = "pdf", width = 5, height = 4)

    
g <- ggplot(tailfindr_filtered_x2, aes(rep, log(polya_length), fill=type, col=type)) + geom_boxplot(alpha=0.4, show.legend = FALSE) + theme_bw()+ xlab("") +
  scale_fill_manual(values=c("cyan", "green")) + scale_color_manual(values=c("cyan4", "green4")) + ggtitle("Tailfindr min 50 from 1308 genes") +
  ylab("polya length (log)") + scale_x_discrete(labels=c("ko1_tailfindr" = "KO1", "ko2_tailfindr" = "KO2", "ko3_tailfindr" = "KO3", "wt1_tailfindr" = "WT1", "wt2_tailfindr" = "WT2", "wt3_tailfindr" = "WT3"))
g <- g + scale_y_continuous(breaks=c(as.numeric(ggplot_build(g)$layout$panel_params[[1]]$y.labels)), labels=c(round(exp(as.numeric(ggplot_build(g)$layout$panel_params[[1]]$y.labels)),0)))
ggsave("plots/general_stats/tailfindr_boxplots_min50_1308.pdf", g, device = "pdf", width = 5, height = 4)


### MEDIANS 

ko1_medians <- aggregate(polya_length~gene_name, ko1_filtered, median) #ko1_filtered has reads with genes that appear at leasy 50 times
ko1_medians$rep <- "KO1"
ko1_medians$type <- "ko"
ko2_medians$rep <- "KO2"
ko2_medians$type <- "ko"
ko3_medians$rep <- "KO3"
ko3_medians$type <- "ko"
wt1_medians <- aggregate(polya_length~gene_name, wt1_filtered, median)
wt1_medians$rep <- "WT1"
wt1_medians$type <- "wt"
wt2_medians <- aggregate(polya_length~gene_name, wt2_filtered, median)
wt2_medians$rep <- "WT2"
wt2_medians$type <- "wt"
wt3_medians <- aggregate(polya_length~gene_name, wt3_filtered, median)
wt3_medians$rep <- "WT3"
wt3_medians$type <- "wt"
#########################
nanopolish_filtered <- rbind(ko1_medians, ko2_medians, ko3_medians, wt1_medians, wt2_medians, wt3_medians)

ko1_medians <- aggregate(polya_length~gene_name, ko1_filtered_t, median)
ko1_medians$rep <- "KO1"
ko1_medians$type <- "ko"
ko2_medians <- aggregate(polya_length~gene_name, ko2_filtered_t, median)
ko2_medians$rep <- "KO2"
ko2_medians$type <- "ko"
ko3_medians <- aggregate(polya_length~gene_name, ko3_filtered_t, median)
ko3_medians$rep <- "KO3"
ko3_medians$type <- "ko"
wt1_medians <- aggregate(polya_length~gene_name, wt1_filtered_t, median)
wt1_medians$rep <- "WT1"
wt1_medians$type <- "wt"
wt2_medians <- aggregate(polya_length~gene_name, wt2_filtered_t, median)
wt2_medians$rep <- "WT2"
wt2_medians$type <- "wt"
wt3_medians <- aggregate(polya_length~gene_name, wt3_filtered_t, median)
wt3_medians$rep <- "WT3"
wt3_medians$type <- "wt"
tailfindr_filtered <- rbind(ko1_medians, ko2_medians, ko3_medians, wt1_medians, wt2_medians, wt3_medians)

nanopolish_filtered_x2 <- nanopolish_filtered[nanopolish_filtered$gene_name %in% genes,]
tailfindr_filtered_x2 <- tailfindr_filtered[tailfindr_filtered$gene_name %in% genes,]

g <- ggplot(nanopolish_filtered_x2, aes(rep, log(polya_length), fill=type, col=type)) + geom_boxplot(alpha=0.4, show.legend = FALSE) + theme_bw()+ xlab("") +
  scale_fill_manual(values=c("cyan", "green")) + scale_color_manual(values=c("cyan4", "green4")) + ggtitle("Nanopolish min 50 from 1308 genes") +
  ylab("median of polya length (log) per gene") + scale_x_discrete(labels=c("ko1_nanopolish" = "KO1", "ko2_nanopolish" = "KO2", "ko3_nanopolish" = "KO3", "wt1_nanopolish" = "WT1", "wt2_nanopolish" = "WT2", "wt3_nanopolish" = "WT3"))
g <- g + scale_y_continuous(breaks=c(as.numeric(ggplot_build(g)$layout$panel_params[[1]]$y.labels)), labels=c(round(exp(as.numeric(ggplot_build(g)$layout$panel_params[[1]]$y.labels)),0)))
ggsave("plots/general_stats/nanopolish_boxplots_median_min50_1308.pdf", g, device = "pdf", width = 5, height = 4)


g <- ggplot(tailfindr_filtered_x2, aes(rep, log(polya_length), fill=type, col=type)) + geom_boxplot(alpha=0.4, show.legend = FALSE) + theme_bw()+ xlab("") +
  scale_fill_manual(values=c("cyan", "green")) + scale_color_manual(values=c("cyan4", "green4")) + ggtitle("Tailfindr min 50 from 1308 genes") +
  ylab("median of polya length (log) per gene") #+ scale_x_discrete(labels=c("ko1_tailfindr" = "KO1", "ko2_tailfindr" = "KO2", "ko3_tailfindr" = "KO3", "wt1_tailfindr" = "WT1", "wt2_tailfindr" = "WT2", "wt3_tailfindr" = "WT3"))
g <- g + scale_y_continuous(breaks=c(as.numeric(ggplot_build(g)$layout$panel_params[[1]]$y.labels)), labels=c(round(exp(as.numeric(ggplot_build(g)$layout$panel_params[[1]]$y.labels)),0)))
ggsave("plots/general_stats/tailfindr_boxplots_median_min50_1308.pdf", g, device = "pdf", width = 5, height = 4)


#Now comparing x-y wt vs ko median of genes.

nanopolish_filtered_x2_ko <- nanopolish_filtered_x2[nanopolish_filtered_x2$type=="ko",]
nanopolish_filtered_x2_ko <- aggregate(polya_length~gene_name, nanopolish_filtered_x2_ko, mean)
nanopolish_filtered_x2_wt <- nanopolish_filtered_x2[nanopolish_filtered_x2$type=="wt",]
nanopolish_filtered_x2_wt <- aggregate(polya_length~gene_name, nanopolish_filtered_x2_wt, mean)

nanopolish <- merge(nanopolish_filtered_x2_ko, nanopolish_filtered_x2_wt, by="gene_name")
e <-cor(nanopolish$polya_length.x, nanopolish$polya_length.y, method = "pearson")
e <- paste("r=", round(e,3), sep="")
g <- ggplot(nanopolish) + geom_point(aes(polya_length.x,polya_length.y), col=densCols(nanopolish$polya_length.x, nanopolish$polya_length.y, colramp = colorRampPalette(blues9[-(1:3)]))) + theme_bw() +
  xlab("KO polyA lengths (mean=38.9)") + ylab("WT polyA lengths (mean=45.4)") + ggtitle(paste("Nanopolish: genes (from 1308) with more than 50 reads (n=",nrow(nanopolish),")", sep="")) +
  geom_text(aes(x=-Inf, y=+Inf, label=e), color="black", size=5, hjust = -0.5, vjust = 2) +
  theme(plot.title = element_text(size=10)) + geom_abline(slope = 1)
m <- min(layer_scales(g)$x$range$range[[1]], layer_scales(g)$y$range$range[[1]])
M <- max(layer_scales(g)$x$range$range[[2]], layer_scales(g)$y$range$range[[2]])
g <- g + expand_limits(x=c(m,M), y=c(m,M))
ggsave("plots/Nanopolish/nanopolish_KO_vs_WT_min50_1308.pdf", g, device = "pdf", width = 5, height = 4)

mean(nanopolish$polya_length.x)
mean(nanopolish$polya_length.y)


tailfindr_filtered_x2_ko <- tailfindr_filtered_x2[tailfindr_filtered_x2$type=="ko",]
tailfindr_filtered_x2_ko <- aggregate(polya_length~gene_name, tailfindr_filtered_x2_ko, mean)
tailfindr_filtered_x2_wt <- tailfindr_filtered_x2[tailfindr_filtered_x2$type=="wt",]
tailfindr_filtered_x2_wt <- aggregate(polya_length~gene_name, tailfindr_filtered_x2_wt, mean)

tailfindr <- merge(tailfindr_filtered_x2_ko, tailfindr_filtered_x2_wt, by="gene_name")
f <-cor(tailfindr$polya_length.x, tailfindr$polya_length.y, method = "pearson")
f <- paste("r=", round(f,3), sep="")
g <- ggplot(tailfindr) + geom_point(aes(polya_length.x,polya_length.y), col=densCols(tailfindr$polya_length.x, tailfindr$polya_length.y, colramp = colorRampPalette(blues9[-(1:3)]))) + 
  xlab("KO polyA_lengths (mean=37.1)") + ylab("WT polyA_lengths (mean=38.9)") + theme_bw() + ggtitle(paste("Tailfindr: genes (from 1308) with more than 50 reads (n=",nrow(tailfindr),")", sep="")) +
  geom_text(aes(x=-Inf, y=+Inf, label=f), color="black", size=5, hjust = -0.5, vjust = 2) +
  theme(plot.title = element_text(size=10)) + geom_abline(slope = 1)
m <- min(layer_scales(g)$x$range$range[[1]], layer_scales(g)$y$range$range[[1]])
M <- max(layer_scales(g)$x$range$range[[2]], layer_scales(g)$y$range$range[[2]])
g <- g + expand_limits(x=c(m,M), y=c(m,M))
ggsave("plots/Tailfindr/tailfindr_KO_vs_WT_min50_1308.pdf", g, device = "pdf", width = 5, height = 4)

mean(tailfindr$polya_length.x)
mean(tailfindr$polya_length.y)

#Same but for ALL genes

nanopolish_filtered_x2_ko <- nanopolish_filtered[nanopolish_filtered$type=="ko",]
nanopolish_filtered_x2_ko <- aggregate(polya_length~gene_name, nanopolish_filtered_x2_ko, mean)
nanopolish_filtered_x2_wt <- nanopolish_filtered[nanopolish_filtered$type=="wt",]
nanopolish_filtered_x2_wt <- aggregate(polya_length~gene_name, nanopolish_filtered_x2_wt, mean)

nanopolish <- merge(nanopolish_filtered_x2_ko, nanopolish_filtered_x2_wt, by="gene_name")
t <- t.test(nanopolish$polya_length.x, nanopolish$polya_length.y, paired = TRUE, alternative = "two.sided")
p_value <- t$p.value  #Pired t-test
a<-strsplit(format(p_value, scientific=T),"e")[[1]]
p_value <- paste("p.value=", round(p_value, -as.integer(a[2]) +1), sep="")
e <-cor(nanopolish$polya_length.x, nanopolish$polya_length.y, method = "pearson")
e <- paste("r=", round(e,3), sep="")
g <- ggplot(nanopolish) + geom_point(aes(polya_length.x,polya_length.y), col=densCols(nanopolish$polya_length.x, nanopolish$polya_length.y, colramp = colorRampPalette(blues9[-(1:3)]))) + theme_bw() + 
  xlab("KO polyA lengths (mean=38)") + ylab("WT polyA lengths (mean=44)") + ggtitle(paste("Nanopolish: genes with more than 50 reads (n=",nrow(nanopolish),")", sep="")) +
  geom_text(aes(x=-Inf, y=+Inf, label=e), color="black", size=5, hjust = -0.5, vjust = 2) +
  geom_text(aes(x=+Inf, y=-Inf, label=p_value), color="black", size=5, hjust = +1.2, vjust = -2) +
  theme(plot.title = element_text(size=10)) + geom_abline(slope = 1)
m <- min(layer_scales(g)$x$range$range[[1]], layer_scales(g)$y$range$range[[1]])
M <- max(layer_scales(g)$x$range$range[[2]], layer_scales(g)$y$range$range[[2]])
g <- g + expand_limits(x=c(m,M), y=c(m,M))
ggsave("plots/Nanopolish/nanopolish_KO_vs_WT_min50.pdf", g, device = "pdf", width = 5, height = 4)

tailfindr_filtered_x2_ko <- tailfindr_filtered[tailfindr_filtered$type=="ko",]
tailfindr_filtered_x2_ko <- aggregate(polya_length~gene_name, tailfindr_filtered_x2_ko, mean)
tailfindr_filtered_x2_wt <- tailfindr_filtered[tailfindr_filtered$type=="wt",]
tailfindr_filtered_x2_wt <- aggregate(polya_length~gene_name, tailfindr_filtered_x2_wt, mean)

tailfindr <- merge(tailfindr_filtered_x2_ko, tailfindr_filtered_x2_wt, by="gene_name")
t <- t.test(tailfindr$polya_length.x, tailfindr$polya_length.y, paired = TRUE, alternative = "two.sided")
p_value <- t$p.value  #Pired t-test
a<-strsplit(format(p_value, scientific=T),"e")[[1]]
p_value <- paste("p.value=", round(p_value, -as.integer(a[2]) +1), sep="")
e <-cor(tailfindr$polya_length.x, tailfindr$polya_length.y, method = "pearson")
e <- paste("r=", round(e,3), sep="")
g <- ggplot(tailfindr) + geom_point(aes(polya_length.x,polya_length.y), col=densCols(tailfindr$polya_length.x, tailfindr$polya_length.y, colramp = colorRampPalette(blues9[-(1:3)]))) + theme_bw() + 
  xlab("KO polyA lengths (mean=38)") + ylab("WT polyA lengths (mean=44)") + ggtitle(paste("Tailfindr: genes with more than 50 reads (n=",nrow(tailfindr),")", sep="")) +
  geom_text(aes(x=-Inf, y=+Inf, label=e), color="black", size=5, hjust = -0.5, vjust = 2) +
  geom_text(aes(x=+Inf, y=-Inf, label=p_value), color="black", size=5, hjust = +1.2, vjust = -2) +
  theme(plot.title = element_text(size=10)) + geom_abline(slope = 1)
m <- min(layer_scales(g)$x$range$range[[1]], layer_scales(g)$y$range$range[[1]])
M <- max(layer_scales(g)$x$range$range[[2]], layer_scales(g)$y$range$range[[2]])
g <- g + expand_limits(x=c(m,M), y=c(m,M))
ggsave("plots/Tailfindr/tailfindr_KO_vs_WT_min50.pdf", g, device = "pdf", width = 5, height = 4)

#p_value table
m <- data.frame("KO1" = 1:6, "KO2" = 1:6, "KO3" = 1:6, "WT1" = 1:6, "WT2" = 1:6, "WT3" = 1:6)
rownames(m) <- c("KO1", "KO2", "KO3", "WT1","WT2","WT3")
for (i in c("KO1", "KO2", "KO3", "WT1","WT2","WT3")){
  datai <- nanopolish_filtered[nanopolish_filtered$rep==i,]
  for (j in c("KO1", "KO2", "KO3", "WT1","WT2","WT3")){
    dataj <- nanopolish_filtered[nanopolish_filtered$rep==j,]
    merged <- merge(datai, dataj, by="gene_name")
    t <- t.test(merged$polya_length.x, merged$polya_length.y, paired = TRUE, alternative = "two.sided")
    m[j,i] <- t$p.value
  }}

pdf("plots/Nanopolish/nanopolish_p_values.pdf", width = 4, height = 6) 
levelplot(log(as.matrix(m)), colorkey=list(space="bottom", tck=c(0,0)) , main="Nanopolish p.values", xlab="log(p.value)", ylab="", margin=FALSE, 
          col.regions = heat.colors(100)[1:length(heat.colors(70))],at=seq(-630,0,by=9))
dev.off()


m <- data.frame("KO1" = 1:6, "KO2" = 1:6, "KO3" = 1:6, "WT1" = 1:6, "WT2" = 1:6, "WT3" = 1:6)
rownames(m) <- c("KO1", "KO2", "KO3", "WT1","WT2","WT3")
for (i in c("KO1", "KO2", "KO3", "WT1","WT2","WT3")){
  datai <- tailfindr_filtered[tailfindr_filtered$rep==i,]
  for (j in c("KO1", "KO2", "KO3", "WT1","WT2","WT3")){
    dataj <- tailfindr_filtered[tailfindr_filtered$rep==j,]
    merged <- merge(datai, dataj, by="gene_name")
    t <- t.test(merged$polya_length.x, merged$polya_length.y, paired = TRUE, alternative = "two.sided")
    m[j,i] <- t$p.value
  }}

pdf("plots/Tailfindr/tailfindr_p_values.pdf", width = 4, height = 6) 
levelplot(log(as.matrix(m)), colorkey=list(space="bottom", tck=c(0,0)) , main="Tailfindr p.values", xlab="log(p.value)", ylab="", margin=FALSE, 
          col.regions = heat.colors(100)[1:length(heat.colors(70))],at=seq(-630,0,by=9))
dev.off()


#KO1 nanopolish vs tailfindr

tailfindr_filtered_x2_ko <- tailfindr_filtered[tailfindr_filtered$type=="ko",]
nanopolish_filtered_x2_ko <- nanopolish_filtered[nanopolish_filtered$type=="ko",]
tailfindr_filtered_x2_wt <- tailfindr_filtered[tailfindr_filtered$type=="wt",]
nanopolish_filtered_x2_wt <- nanopolish_filtered[nanopolish_filtered$type=="wt",]

nanopolish_vs_tailfindr<- function(data1, data2, string){
  test1 <- data1[data1$rep==string,]
  test2 <- data2[data2$rep==string,]
  merged <- merge(test1, test2, by=c("gene_name","rep","type"))
  g <-cor(merged$polya_length.x, merged$polya_length.y, method = "pearson")
  g <- paste("r=", round(g,3), sep="")
  plot <- ggplot(merged) + geom_point(aes(polya_length.x,polya_length.y), col=densCols(merged$polya_length.x, merged$polya_length.y, colramp = colorRampPalette(blues9[-(1:3)]))) + ylim(20,65) + xlim(20,65) +
    xlab(paste("Nanopolish", string,"polyA lengths")) + ylab(paste("Tailfindr", string,"polyA lengths")) + theme_bw() + ggtitle(paste(string, ": genes with more than 50 reads (n=",nrow(merged),")", sep="")) +
    geom_text(aes(x=-Inf, y=+Inf, label=g), color="black", size=5, hjust = -0.5, vjust = 2) +
    theme(plot.title = element_text(size=10)) + geom_abline(slope = 1)
  ggsave(paste("plots/Nanopolish_vs_Tailfindr/", string, "_nanopolish_vs_tailfindr_min50.pdf",sep=""), plot, device = "pdf", width = 5, height = 4)
}
nanopolish_vs_tailfindr(nanopolish_filtered_x2_ko, tailfindr_filtered_x2_ko, "KO1")
nanopolish_vs_tailfindr(nanopolish_filtered_x2_ko, tailfindr_filtered_x2_ko, "KO2")
nanopolish_vs_tailfindr(nanopolish_filtered_x2_ko, tailfindr_filtered_x2_ko, "KO3")
nanopolish_vs_tailfindr(nanopolish_filtered_x2_wt, tailfindr_filtered_x2_wt, "WT1")
nanopolish_vs_tailfindr(nanopolish_filtered_x2_wt, tailfindr_filtered_x2_wt, "WT2")
nanopolish_vs_tailfindr(nanopolish_filtered_x2_wt, tailfindr_filtered_x2_wt, "WT3")

nanopolish_vs_tailfindr_samedata<- function(data, string1, string2, string3){
  test1 <- data[data$rep==string1,]
  test2 <- data[data$rep==string2,]
  merged <- merge(test1, test2, by="gene_name")
  t <- t.test(merged$polya_length.x, merged$polya_length.y, paired = TRUE, alternative = "two.sided")
  p_value <- t$p.value  #Paired t-test
  a<-strsplit(format(p_value, scientific=T),"e")[[1]]
  p_value <- paste("p.value=", round(p_value, -as.integer(a[2]) +1), sep="")
  g <-cor(merged$polya_length.x, merged$polya_length.y, method = "pearson")
  g <- paste("r=", round(g,3), sep="")
  plot <- ggplot(merged) + geom_point(aes(polya_length.x,polya_length.y), col=densCols(merged$polya_length.x, merged$polya_length.y, colramp = colorRampPalette(blues9[-(1:3)]))) + ylim(15,65) + xlim(15,65) +
    xlab(paste(string3, string1,"polyA lengths")) + ylab(paste(string3, string2,"polyA lengths")) + theme_bw() + ggtitle(paste(string3, ", ", string1, " vs ", string2, ": genes with more than 50 reads (n=",nrow(merged),")", sep="")) +
    geom_text(aes(x=-Inf, y=+Inf, label=g), color="black", size=5, hjust = -0.5, vjust = 2) +
    geom_text(aes(x=+Inf, y=-Inf, label=p_value), color="black", size=5, hjust = +1.2, vjust = -2) +
    theme(plot.title = element_text(size=10)) + geom_abline(slope=1)
  ggsave(paste("plots/", string3, "/", string3, "_", string1, "_vs_", string2, "_min50.pdf",sep=""), plot, device = "pdf", width = 5, height = 4)
}

strings <- levels(as.factor(nanopolish_filtered$rep))
for (i in 1:5){
  for (j in (i+1):6){
    if (i==j){next}
    nanopolish_vs_tailfindr_samedata(nanopolish_filtered, strings[i], strings[j], "Nanopolish")
    nanopolish_vs_tailfindr_samedata(tailfindr_filtered, strings[i], strings[j], "Tailfindr")
  }}


### Now, what we are going to be comparing are: KO1 vs WT1, KO2 vs WT2, KO3 vs WT3 and KO vs WT.

colored_plots <- function(rep, string1, string2, string3){
  dist <- sqrt((rep$polya_length.x - rep$polya_length.y)^2 / 2)
  threshold <- 2.5*sd(dist)
  genes_u <- dist>threshold
  genes <- rep[genes_u,]$gene_name
  upper <- rep[rep$gene_name %in% genes,]$polya_length.x < rep[rep$gene_name %in% genes,]$polya_length.y
  genes_upper <- genes[upper]
  genes_lower <- genes[!upper]
  t <- t.test(rep$polya_length.x, rep$polya_length.y, paired = TRUE, alternative = "two.sided")
  p_value <- t$p.value  #Paired t-test
  a<-strsplit(format(p_value, scientific=T),"e")[[1]]
  p_value <- paste("p.value=", round(p_value, -as.integer(a[2]) +1), sep="")
  g <-cor(rep$polya_length.x, rep$polya_length.y, method = "pearson")
  g <- paste("r=", round(g,3), sep="")
  plot <- ggplot(rep, aes(polya_length.x,polya_length.y)) + geom_point(col=densCols(rep$polya_length.x, rep$polya_length.y, colramp = colorRampPalette(blues9[-(1:3)])), 
                                                                       size=1) + ylim(15,65) + xlim(15,65) +
    xlab(paste(string1, "polyA lengths")) + ylab(paste(string2, "polyA lengths")) + theme_bw() + ggtitle(paste(string3, " (nº genes = ", nrow(rep), ")", sep="")) +
    geom_text(aes(x=-Inf, y=+Inf, label=g), size=4, hjust = -0.2, vjust = 1.8) +
    geom_text(aes(x=+Inf, y=-Inf, label=p_value), size=4, hjust = +1.1, vjust = -1) +
    theme(plot.title = element_text(size=10)) + geom_abline(slope=1) + 
    geom_point(data=rep[rep$gene_name %in% genes,], aes(polya_length.x,polya_length.y), col="red", size=1) +
    theme(plot.title = element_text(hjust = 0.5, size=12, face = "bold"))
  l <- list(plot, genes, genes_upper, genes_lower)
  return(l)
}

rep1 <- merge(nanopolish_filtered[nanopolish_filtered$rep=="KO1",], nanopolish_filtered[nanopolish_filtered$rep=="WT1",], by="gene_name")
out1 <- colored_plots(rep1, "KO1", "WT1", "Rep1")
rep2 <- merge(nanopolish_filtered[nanopolish_filtered$rep=="KO2",], nanopolish_filtered[nanopolish_filtered$rep=="WT2",], by="gene_name")
out2 <- colored_plots(rep2, "KO2", "WT2", "Rep2")
rep3 <- merge(nanopolish_filtered[nanopolish_filtered$rep=="KO3",], nanopolish_filtered[nanopolish_filtered$rep=="WT3",], by="gene_name")
out3 <- colored_plots(rep3, "KO3", "WT3", "Rep3")

n_upper <- Reduce(intersect, list(out1[[3]],out2[[3]],out3[[3]])) 
n_lower <- Reduce(intersect, list(out1[[4]],out2[[4]],out3[[4]])) 

rep4 <- merge(tailfindr_filtered[tailfindr_filtered$rep=="KO1",], tailfindr_filtered[tailfindr_filtered$rep=="WT1",], by="gene_name")
out4 <- colored_plots(rep4, "KO1", "WT1", "Rep1")
rep5 <- merge(tailfindr_filtered[tailfindr_filtered$rep=="KO2",], tailfindr_filtered[tailfindr_filtered$rep=="WT2",], by="gene_name")
out5 <- colored_plots(rep5, "KO2", "WT2", "Rep2")
rep6 <- merge(tailfindr_filtered[tailfindr_filtered$rep=="KO3",], tailfindr_filtered[tailfindr_filtered$rep=="WT3",], by="gene_name")
out6 <- colored_plots(rep6, "KO3", "WT3", "Rep3")

t_upper <- Reduce(intersect, list(out4[[3]],out5[[3]],out6[[3]])) 
t_lower <- Reduce(intersect, list(out4[[4]],out5[[4]],out6[[4]])) 

u12 <- intersect(out1[[3]],out2[[3]])
u13 <- intersect(out1[[3]],out3[[3]])
u23 <- intersect(out2[[3]],out3[[3]])
total_upper <- c(u12,u13,u23)
total_upper <- unique(total_upper) #upper differentially expressed genes in at least two of the replicates
write.table(total_upper, "plots/final_plots/nanopolish_differentially_expressed_genes_at_least_2_replicates.txt", append = FALSE, quote = FALSE,
            row.names = FALSE, col.names = FALSE)

u12 <- intersect(out4[[3]],out5[[3]])
u13 <- intersect(out4[[3]],out6[[3]])
u23 <- intersect(out5[[3]],out6[[3]])
total_upper <- c(u12,u13,u23)
total_upper <- unique(total_upper) #upper differentially expressed genes in at least two of the replicates
write.table(total_upper, "plots/final_plots/tailfindr_differentially_expressed_genes_at_least_2_replicates.txt", append = FALSE, quote = FALSE,
            row.names = FALSE, col.names = FALSE)

venn_nanopolish_upper <- ssvFeatureVenn(list(out1[[3]], out2[[3]], out3[[3]]), group_names = c("Rep1", "Rep2", "Rep3"),
                                        line_width = 2) + ggtitle("Overlap of upper genes") +
  theme(legend.position="right", legend.text = element_text(size = 15), plot.title = element_text(hjust = 0.5, size=12, face = "bold"))
venn_nanopolish_lower <- ssvFeatureVenn(list(out1[[4]], out2[[4]], out3[[4]]), group_names = c("Rep1", "Rep2", "Rep3"),
                                        line_width = 2) + ggtitle("Overlap of lower genes") +
  theme(legend.position="right", legend.text = element_text(size = 15), plot.title = element_text(hjust = 0.5, size=12, face = "bold"))

venn_tailfindr_upper <- ssvFeatureVenn(list(out4[[3]], out5[[3]], out6[[3]]), group_names = c("Rep1", "Rep2", "Rep3"),
                                       line_width = 2) + ggtitle("Overlap of upper genes") +
  theme(legend.position="right", legend.text = element_text(size = 15), plot.title = element_text(hjust = 0.5, size=12, face = "bold"))
venn_tailfindr_lower <- ssvFeatureVenn(list(out4[[4]], out5[[4]], out6[[4]]), group_names = c("Rep1", "Rep2", "Rep3"),
                                       line_width = 2) + ggtitle("Overlap of lower genes") +
  theme(legend.position="right", legend.text = element_text(size = 15), plot.title = element_text(hjust = 0.5, size=12, face = "bold"))

g <- ggarrange(out1[[1]], out2[[1]], out3[[1]], venn_nanopolish_upper, venn_nanopolish_lower, out4[[1]], out5[[1]], out6[[1]], venn_tailfindr_upper, venn_tailfindr_lower, nrow=2, ncol=5, common.legend = TRUE, legend = "right")
g <- annotate_figure(g, left = text_grob("          Tailfindr                                      Nanopolish", face = "bold", size = 14, rot=90),)
ggsave("plots/final_plots/xy_differential_genes_min50.pdf", g, device = "pdf", width = 15, height = 6)

#Same for n=30
ko1_filtered <- filtering(ko1_nanopolish, 30)
ko2_filtered <- filtering(ko2_nanopolish, 30)
ko3_filtered <- filtering(ko3_nanopolish, 30)
wt1_filtered <- filtering(wt1_nanopolish, 30)
wt2_filtered <- filtering(wt2_nanopolish, 30)
wt3_filtered <- filtering(wt3_nanopolish, 30)

ko1_filtered_t <- filtering(ko1_tailfindr, 30)
ko2_filtered_t <- filtering(ko2_tailfindr, 30)
ko3_filtered_t <- filtering(ko3_tailfindr, 30)
wt1_filtered_t <- filtering(wt1_tailfindr, 30)
wt2_filtered_t <- filtering(wt2_tailfindr, 30)
wt3_filtered_t <- filtering(wt3_tailfindr, 30)

ko1_medians <- aggregate(polya_length~gene_name, ko1_filtered, median) #ko1_filtered has reads with genes that appear at least 50 times
ko1_medians$rep <- "KO1"
ko2_medians <- aggregate(polya_length~gene_name, ko2_filtered, median)
ko2_medians$rep <- "KO2"
ko3_medians <- aggregate(polya_length~gene_name, ko3_filtered, median)
ko3_medians$rep <- "KO3"
wt1_medians <- aggregate(polya_length~gene_name, wt1_filtered, median)
wt1_medians$rep <- "WT1"
wt2_medians <- aggregate(polya_length~gene_name, wt2_filtered, median)
wt2_medians$rep <- "WT2"
wt3_medians <- aggregate(polya_length~gene_name, wt3_filtered, median)
wt3_medians$rep <- "WT3"
nanopolish_filtered_30 <- rbind(ko1_medians, ko2_medians, ko3_medians, wt1_medians, wt2_medians, wt3_medians)

ko1_medians <- aggregate(polya_length~gene_name, ko1_filtered_t, median)
ko1_medians$rep <- "KO1"
ko2_medians <- aggregate(polya_length~gene_name, ko2_filtered_t, median)
ko2_medians$rep <- "KO2"
ko3_medians <- aggregate(polya_length~gene_name, ko3_filtered_t, median)
ko3_medians$rep <- "KO3"
wt1_medians <- aggregate(polya_length~gene_name, wt1_filtered_t, median)
wt1_medians$rep <- "WT1"
wt2_medians <- aggregate(polya_length~gene_name, wt2_filtered_t, median)
wt2_medians$rep <- "WT2"
wt3_medians <- aggregate(polya_length~gene_name, wt3_filtered_t, median)
wt3_medians$rep <- "WT3"
tailfindr_filtered_30 <- rbind(ko1_medians, ko2_medians, ko3_medians, wt1_medians, wt2_medians, wt3_medians)

rep1 <- merge(nanopolish_filtered_30[nanopolish_filtered_30$rep=="KO1",], nanopolish_filtered_30[nanopolish_filtered_30$rep=="WT1",], by="gene_name")
out1 <- colored_plots(rep1, "KO1", "WT1", "Rep1")
rep2 <- merge(nanopolish_filtered_30[nanopolish_filtered_30$rep=="KO2",], nanopolish_filtered_30[nanopolish_filtered_30$rep=="WT2",], by="gene_name")
out2 <- colored_plots(rep2, "KO2", "WT2", "Rep2")
rep3 <- merge(nanopolish_filtered_30[nanopolish_filtered_30$rep=="KO3",], nanopolish_filtered_30[nanopolish_filtered_30$rep=="WT3",], by="gene_name")
out3 <- colored_plots(rep3, "KO3", "WT3", "Rep3")

rep4 <- merge(tailfindr_filtered_30[tailfindr_filtered_30$rep=="KO1",], tailfindr_filtered_30[tailfindr_filtered_30$rep=="WT1",], by="gene_name")
out4 <- colored_plots(rep4, "KO1", "WT1", "Rep1")
rep5 <- merge(tailfindr_filtered_30[tailfindr_filtered_30$rep=="KO2",], tailfindr_filtered_30[tailfindr_filtered_30$rep=="WT2",], by="gene_name")
out5 <- colored_plots(rep5, "KO2", "WT2", "Rep2")
rep6 <- merge(tailfindr_filtered_30[tailfindr_filtered_30$rep=="KO3",], tailfindr_filtered_30[tailfindr_filtered_30$rep=="WT3",], by="gene_name")
out6 <- colored_plots(rep6, "KO3", "WT3", "Rep3")

u12 <- intersect(out1[[3]],out2[[3]])
u13 <- intersect(out1[[3]],out3[[3]])
u23 <- intersect(out2[[3]],out3[[3]])
total_upper <- c(u12,u13,u23)
total_upper <- unique(total_upper) #upper differentially expressed genes in at least two of the replicates
write.table(total_upper, "plots/final_plots/nanopolish_differentially_expressed_genes_at_least_2_replicates_min30.txt", append = FALSE, quote = FALSE,
            row.names = FALSE, col.names = FALSE)

u12 <- intersect(out4[[3]],out5[[3]])
u13 <- intersect(out4[[3]],out6[[3]])
u23 <- intersect(out5[[3]],out6[[3]])
total_upper <- c(u12,u13,u23)
total_upper <- unique(total_upper) #upper differentially expressed genes in at least two of the replicates
write.table(total_upper, "plots/final_plots/tailfindr_differentially_expressed_genes_at_least_2_replicates_min30.txt", append = FALSE, quote = FALSE,
            row.names = FALSE, col.names = FALSE)

venn_nanopolish_upper <- ssvFeatureVenn(list(out1[[3]], out2[[3]], out3[[3]]), group_names = c("Rep1", "Rep2", "Rep3"),
                                        line_width = 2) + ggtitle("Overlap of upper genes") +
  theme(legend.position="right", legend.text = element_text(size = 15), plot.title = element_text(hjust = 0.5, size=12, face = "bold"))
venn_nanopolish_lower <- ssvFeatureVenn(list(out1[[4]], out2[[4]], out3[[4]]), group_names = c("Rep1", "Rep2", "Rep3"),
                                        line_width = 2) + ggtitle("Overlap of lower genes") +
  theme(legend.position="right", legend.text = element_text(size = 15), plot.title = element_text(hjust = 0.5, size=12, face = "bold"))

venn_tailfindr_upper <- ssvFeatureVenn(list(out4[[3]], out5[[3]], out6[[3]]), group_names = c("Rep1", "Rep2", "Rep3"),
                                       line_width = 2) + ggtitle("Overlap of upper genes") +
  theme(legend.position="right", legend.text = element_text(size = 15), plot.title = element_text(hjust = 0.5, size=12, face = "bold"))
venn_tailfindr_lower <- ssvFeatureVenn(list(out4[[4]], out5[[4]], out6[[4]]), group_names = c("Rep1", "Rep2", "Rep3"),
                                       line_width = 2) + ggtitle("Overlap of lower genes") +
  theme(legend.position="right", legend.text = element_text(size = 15), plot.title = element_text(hjust = 0.5, size=12, face = "bold"))

g <- ggarrange(out1[[1]], out2[[1]], out3[[1]], venn_nanopolish_upper, venn_nanopolish_lower, out4[[1]], out5[[1]], out6[[1]], venn_tailfindr_upper, venn_tailfindr_lower, nrow=2, ncol=5, common.legend = TRUE, legend = "right")
g <- annotate_figure(g, left = text_grob("          Tailfindr                                      Nanopolish", face = "bold", size = 14, rot=90),)
ggsave("plots/final_plots/xy_differential_genes_min30.pdf", g, device = "pdf", width = 15, height = 6)

#Same with KO vs WT in general

merged_plot <- function(merged, string){
dist <- sqrt((merged$polya_length.x - merged$polya_length.y)^2 / 2)
threshold <- 2.5*sd(dist)
genes_u <- dist>threshold
diff <- merged[genes_u,]
upper <- diff$polya_length.x < diff$polya_length.y
genes_upper <- diff[upper,]
genes_lower <- diff[!upper,]
t <- t.test(merged$polya_length.x, merged$polya_length.y, paired = TRUE, alternative = "two.sided")
p_value <- t$p.value  #Paired t-test
a<-strsplit(format(p_value, scientific=T),"e")[[1]]
p_value <- paste("p.value=", round(p_value, -as.integer(a[2]) +1), sep="")
e <-cor(merged$polya_length.x, merged$polya_length.y, method = "pearson")
e <- paste("r=", round(e,3), sep="")
g <- ggplot(merged) + geom_point(aes(polya_length.x,polya_length.y), col=densCols(merged$polya_length.x, merged$polya_length.y, colramp = colorRampPalette(blues9[-(1:3)]))) + theme_bw() +
  xlab("KO polyA lengths") + ylab("WT polyA lengths") + ggtitle(paste(string, ": rep1 + rep2 + rep3 (n=",nrow(merged),")", sep="")) +
  geom_text(aes(x=-Inf, y=+Inf, label=e), color="black", size=5, hjust = -0.5, vjust = 2) +
  geom_text(aes(x=+Inf, y=-Inf, label=p_value), color="black", size=5, hjust = +1.2, vjust = -2) +
  theme(plot.title = element_text(size=10)) + geom_abline(slope = 1) +
  geom_point(data=genes_upper, aes(polya_length.x,polya_length.y), col="red") +
  geom_point(data=genes_lower, aes(polya_length.x,polya_length.y), col="red")
m <- min(layer_scales(g)$x$range$range[[1]], layer_scales(g)$y$range$range[[1]])
M <- max(layer_scales(g)$x$range$range[[2]], layer_scales(g)$y$range$range[[2]])
g <- g + expand_limits(x=c(m,M), y=c(m,M))
ggsave(paste("plots/final_plots/", string, "_KO_vs_WT.pdf", sep=""), g, device = "pdf", width = 6, height = 6)
return(genes_upper)
}

merged <- rbind(rep1, rep2, rep3)
genes_upper_nanopolish <- merged_plot(merged, "Nanopolish")
merged <- rbind(rep4, rep5, rep6)
genes_upper_tailfindr <- merged_plot(merged, "Tailfindr")

# Check GO enrichment of the two files and the intersection

write.table(genes_upper_tailfindr$gene_name, "plots/final_plots/tailfindr_differentially_expressed_genes_global_min50.txt", append = FALSE, quote = FALSE,
            row.names = FALSE, col.names = FALSE)
write.table(genes_upper_nanopolish$gene_name, "plots/final_plots/nanopolish_differentially_expressed_genes_global_min50.txt", append = FALSE, quote = FALSE,
            row.names = FALSE, col.names = FALSE)

common <- intersect(genes_upper_tailfindr$gene_name, genes_upper_nanopolish$gene_name)

write.table(common, "plots/final_plots/common_differentially_expressed_genes_global_min50.txt", append = FALSE, quote = FALSE,
            row.names = FALSE, col.names = FALSE)

# 5) do genes that are thought to contain m6A (the 1308 genes) in general have 
# shorter/longer pA lengths than those that are thought not to have it? And if 
# this is the case, is it visible in WT but not in KO?

tailfindr_filtered_genes <- tailfindr_filtered[tailfindr_filtered$gene_name %in% genes,]
tailfindr_filtered_not_genes <- tailfindr_filtered[!tailfindr_filtered$gene_name %in% genes,]

tailfindr_filtered_genes$group <- "m6A"
tailfindr_filtered_not_genes$group <- "not_m6A"
data <- rbind(tailfindr_filtered_genes, tailfindr_filtered_not_genes)
g <- ggplot(data) + geom_boxplot(aes(group, polya_length, fill=group, alpha=0.4, color=group), show.legend = FALSE) + ggtitle("Tailfindr min50") + theme_bw() + xlab("")
ggsave("plots/m6A_vs_not_m6A/Tailfindr.pdf", g, device = "pdf", width = 4, height = 4)

tailfindr_filtered_genes_ko <- tailfindr_filtered_genes[tailfindr_filtered_genes$type=="ko",]
tailfindr_filtered_genes_ko$group <- "m6A_KO"
tailfindr_filtered_genes_wt <- tailfindr_filtered_genes[tailfindr_filtered_genes$type=="wt",]
tailfindr_filtered_genes_wt$group <- "m6A_WT"

tailfindr_filtered_not_genes_ko <- tailfindr_filtered_not_genes[tailfindr_filtered_not_genes$type=="ko",]

tailfindr_filtered_not_genes_ko$group <- "not_m6A_KO"
tailfindr_filtered_not_genes_wt <- tailfindr_filtered_not_genes[tailfindr_filtered_not_genes$type=="wt",]
tailfindr_filtered_not_genes_wt$group <- "not_m6A_WT"

data <- rbind(tailfindr_filtered_genes_ko, tailfindr_filtered_genes_wt, tailfindr_filtered_not_genes_ko, tailfindr_filtered_not_genes_wt)
g <- ggplot(data) + geom_boxplot(aes(group, polya_length, fill=group, alpha=0.4, color=group), show.legend = FALSE) + ggtitle("Tailfindr min50") + theme_bw() + xlab("")
ggsave("plots/m6A_vs_not_m6A/Tailfindr_KO_WT.pdf", g, device = "pdf", width = 5, height = 4)

nanopolish_filtered_genes <- nanopolish_filtered[nanopolish_filtered$gene_name %in% genes,]
nanopolish_filtered_not_genes <- nanopolish_filtered[!nanopolish_filtered$gene_name %in% genes,]

nanopolish_filtered_genes$group <- "m6A"
nanopolish_filtered_not_genes$group <- "not_m6A"
data <- rbind(nanopolish_filtered_genes, nanopolish_filtered_not_genes)
g <- ggplot(data) + geom_boxplot(aes(group, polya_length, fill=group, alpha=0.4, color=group), show.legend = FALSE) + ggtitle("Nanopolish min50") + xlab("") + theme_bw()
ggsave("plots/m6A_vs_not_m6A/Nanopolish.pdf", g, device = "pdf", width = 4, height = 4)

nanopolish_filtered_genes_ko <- nanopolish_filtered_genes[nanopolish_filtered_genes$type=="ko",]
nanopolish_filtered_genes_ko$group <- "m6A_KO"
nanopolish_filtered_genes_wt <- nanopolish_filtered_genes[nanopolish_filtered_genes$type=="wt",]
nanopolish_filtered_genes_wt$group <- "m6A_WT"

nanopolish_filtered_not_genes_ko <- nanopolish_filtered_not_genes[nanopolish_filtered_not_genes$type=="ko",]
nanopolish_filtered_not_genes_ko$group <- "not_m6A_KO"
nanopolish_filtered_not_genes_wt <- nanopolish_filtered_not_genes[nanopolish_filtered_not_genes$type=="wt",]
nanopolish_filtered_not_genes_wt$group <- "not_m6A_WT"

data <- rbind(nanopolish_filtered_genes_ko, nanopolish_filtered_genes_wt, nanopolish_filtered_not_genes_ko, nanopolish_filtered_not_genes_wt)
g <- ggplot(data) + geom_boxplot(aes(group, polya_length, fill=group, alpha=0.4, color=group), show.legend = FALSE) + ggtitle("Nanopolish min50") + xlab("") + theme_bw()
ggsave("plots/m6A_vs_not_m6A/Nanopolish_KO_WT.pdf", g, device = "pdf", width = 5, height = 4)

# 8) for invididual READS (not just genes), can you also make the comparison 
# (like the plots above) between tailfindr and nanopolish?

tailfindr <- rbind(ko1_tailfindr, ko2_tailfindr, ko3_tailfindr, wt1_tailfindr, wt2_tailfindr, wt3_tailfindr)
nanopolish <- rbind(ko1_nanopolish, ko2_nanopolish, ko3_nanopolish, wt1_nanopolish, wt2_nanopolish, wt3_nanopolish)

head(tailfindr)
head(nanopolish)
merged <- merge(nanopolish, tailfindr, by="readname")

reads <- function(data, string){
  merged_temp <- data[data$rep.x==paste(string, "_nanopolish", sep=""),]
  f <-cor(merged_temp$polya_length.x, merged_temp$polya_length.y, method = "pearson")
  f <- paste("r=", round(f,3), sep="")
  g <- ggplot(merged_temp) + geom_point(aes(log(polya_length.x),log(polya_length.y)), col=densCols(merged_temp$polya_length.x, merged_temp$polya_length.y, colramp = colorRampPalette(blues9[-(1:3)]))) + 
    xlab("Nanopolish polya lengths (log)") + ylab("Tailfindr polya lengths (log)") + theme_bw() + ggtitle(paste("Nanopolish vs Tailfindr: ", toupper(string), " (n=",nrow(merged_temp),")", sep="")) +
    geom_text(aes(x=-Inf, y=+Inf, label=f), color="black", size=5, hjust = -0.5, vjust = 2) +
    theme(plot.title = element_text(size=10))
  g <- g + scale_y_continuous(breaks=c(as.numeric(ggplot_build(g)$layout$panel_params[[1]]$y.labels)), labels=c(round(exp(as.numeric(ggplot_build(g)$layout$panel_params[[1]]$y.labels)),0))) +
    scale_x_continuous(breaks=c(as.numeric(ggplot_build(g)$layout$panel_params[[1]]$x.labels)), labels=c(round(exp(as.numeric(ggplot_build(g)$layout$panel_params[[1]]$x.labels)),0)))
  
  ggsave(paste("plots/reads/", toupper(string), "_tailfindr_vs_nanopolish_all_reads.pdf", sep=""), g, device = "pdf", width = 5, height = 4)
  
}

reads(merged, "ko1")
reads(merged, "ko2")
reads(merged, "ko3")
reads(merged, "wt1")
reads(merged, "wt2")
reads(merged, "wt3")

# 9)Can we make a plot of the predicted lengths (boxplots) for each of these 
# subgroups of genes? I wonder if nanopolish common for KO1 is significantly 
# different than nanopolish (not in tailfindr) for KO1 for example. And same 
# for tailfindr, I wonder if the length of tailfindr specific cases is different 
# than than of those that are tailfindr-common5) do genes that are thought to 
# contain m6A (the 1308 genes) in general have shorter/longer pA lengths than 
# those that are thought not to have it? And if this is the case, is it visible 
# in WT but not in KO?

#for all reads.

ko1_merged <- merge(ko1_nanopolish, ko1_tailfindr, by=c("readname", "gene_name"))
ko1_common_nanopolish <- ko1_merged[,c("readname","gene_name","polya_length.x","rep.x", "type.x")]
colnames(ko1_common_nanopolish) <- c("readname","gene_name","polya_length","rep","type")
ko1_common_nanopolish$groups <- "common_nanopolish"

ko1_common_tailfindr <- ko1_merged[,c("readname","gene_name","polya_length.y","rep.y", "type.y")]
colnames(ko1_common_tailfindr) <- c("readname","gene_name","polya_length","rep", "type")
ko1_common_tailfindr$groups <- "common_taifindr"

ko1_only_nanopolish <- setdiff(ko1_nanopolish$readname, ko1_merged$readname)
ko1_only_nanopolish <- ko1_nanopolish[ko1_nanopolish$readname %in% ko1_only_nanopolish,]
ko1_only_nanopolish$groups <- "only_in_nanopolish"
ko1_only_tailfindr <- setdiff(ko1_tailfindr$readname, ko1_merged$readname)
ko1_only_tailfindr <- ko1_tailfindr[ko1_tailfindr$readname %in% ko1_only_tailfindr,]
ko1_only_tailfindr$groups <- "only_in_tailfindr"

ko1_merged <- rbind(ko1_only_nanopolish, ko1_only_tailfindr, ko1_common_nanopolish, ko1_common_tailfindr)

g <- ggplot(ko1_merged) + geom_boxplot(aes(groups, log(polya_length), fill=groups, alpha=0.4, color=groups), show.legend = FALSE) 
g <- g + scale_y_continuous(breaks=c(as.numeric(ggplot_build(g)$layout$panel_params[[1]]$y.labels)), labels=c(round(exp(as.numeric(ggplot_build(g)$layout$panel_params[[1]]$y.labels)),0))) +
  theme_bw() + xlab("") + ylab("polyA length (log)")
ggsave("plots/general_stats/common_unique_boxplots.pdf", g, device = "pdf", width = 6, height = 4)





