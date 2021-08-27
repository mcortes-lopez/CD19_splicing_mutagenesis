###################################################
#
# script adapted by Mariela Cortés-López
# from the original written by Samarth T. Setty and Anke Busch
#
###################################################


library(R.utils)
library(Biostrings)
library(ShortRead)

## get input arguments
args <- commandArgs(asValues = TRUE)

## check the args are properly assigned
print(args)

## assign input arguments
input.folder <- args$inputfolder
output.folder <- args$outputfolder
fastqfile <- args$fastqfile

libraryname <- args$lib


## read in fastq file
print(paste0(input.folder, "/", fastqfile))
subset.fastq <- readFastq(dirPath = paste0(input.folder, "/", fastqfile))


## split ids to get the barcodes
id.split <- strsplit(as.character(id(subset.fastq)), ":")


barcodes <- sapply(id.split, "[", 2) # returns a vector


## split sequences based on their barcode into separate elements of a list,
## where each list element holds all sequences with a certain barcode
split.fastq <- split(subset.fastq, barcodes)
head(split.fastq)

# Correction of the names to
change_ids <- function(x) {
  ShortReadQ(
    sread(split.fastq[[x]]), quality(split.fastq[[x]]),
    BStringSet(paste0(
      gsub(paste0("ccs:", x), "", id(split.fastq[[x]])),
      "0_", width(split.fastq[[x]])
    ))
  )
}


split.fastq2 <- lapply(names(split.fastq), change_ids)
names(split.fastq2) <- names(split.fastq)
rm(split.fastq)
lapply(names(split.fastq2), function(barcodename) {
  writeFastq(split.fastq2[[barcodename]],
    file = paste(output.folder, libraryname, "_", barcodename, ".fastq.gz", sep = "")
  )
})



## write summary stats table (count per barcode)
bc.count <- sapply(split.fastq2, function(i) {
  return(length(i))
})
write.table(bc.count, file = paste0(input.folder, "/barcode.counts.read.tsv"), sep = "\t", quote = FALSE, col.names = FALSE)
