###################################################
#
# script written by Samarth T. Setty and Anke Busch
#
###################################################

library(R.utils)
library(Biostrings)
library(ShortRead)

# get input arguments
args=commandArgs(asValues = TRUE)

# check the args are properly assigned
print(args)

# assign input arguments
input.folder  <- args$inputPath 
output.folder <- args$outputPath
fastqfile     <- args$fastqfile
libraryname   <- args$lib
read.number   <- as.character(args$read)

# check if input arguments are correct
if ((read.number != "1") && (read.number != "2")) {
  cat("ERROR! The read number should be either 1 or 2, specify by --read=1 or --read=2\n")
  quit(save="no")
}

# read in fastq file
subset.fastq <- readFastq(dirPath = input.folder, pattern = fastqfile)

# split ids to get the barcodes
id.split <- strsplit(as.character(id(subset.fastq)), paste0(read.number,':N:0:1:'))

# select second entry of each list element
barcodes <- sapply(id.split, "[", 2) 

# split sequences based on their barcode into separate elements of a list,
# where each list element holds all sequences with a certain barcode
split.fastq <- split(subset.fastq, barcodes)

# write reads to separate fastq files (per barcode)
lapply(names(split.fastq), 
       function(barcodename){ 
         writeFastq(split.fastq[[barcodename]], file = paste0(output.folder,"/",libraryname,'_',read.number,'_', barcodename,'.fastq.gz'))})

# write summary stats table (count per barcode)
bc.count <-sapply(split.fastq, function(i){ return(length(i)) })
write.table(bc.count, file=paste0(input.folder,"/barcode.counts.read",read.number,".tsv"), sep="\t", quote=FALSE, col.names=FALSE)



