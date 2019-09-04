#Loading libraries
library(tailfindr)
library(optparse)

#Parsing
parser <- OptionParser()
parser <- add_option(parser, opt_str=c("-i", "--input"), type="character",
                     dest='input',
                     help="Input Directory")
parser <- add_option(parser, opt_str=c("-o", "--output"), type="character",
                     dest='output',
                     help="Output Dir")
parser <- add_option(parser, opt_str=c("-n", "--name"), type="character",
                     dest='name',
                     help="Name")
options=parse_args(parser)
dir=options$input
output=options$output
name=options$name

find_tails(fast5_dir= dir, save_dir = output, csv_filename = paste(name, ".csv", sep=""), num_cores = 8)
