###################################################
#
# script written by Samarth T. Setty and Anke Busch
#
###################################################

####### load libraries
library(R.utils)
library(Biostrings)
library(ShortRead)
options(stringsAsFactors=FALSE)

# read arguments
args <- commandArgs(asValues = TRUE)

# check if args are properly assigned
print(args)


###### get input arguments
input.folder   <- args$inputPath
file.read1     <- args$Read1
file.read2     <- args$Read2
output.folder  <- args$outputPath
barcode.length <- as.integer(args$barcodelength)
minlen.read1   <- as.integer(args$minlengthread1)
minlen.read2   <- as.integer(args$minlengthread2)


# anchor sequences around barcode
lseq    <- "TGCAGAATTC"
rseq    <- "GGATCC"

# random N at the beginning of read1 and read2
read1.N <- 10
read2.N <- 10

# number of extra positions to be trimmed off of read2 after read2.N (to account for a few more random bases than the intended read2.N)
read2.bonus.trimmed <- 2


#######  input
read1.fastq <- readFastq(dirPath=paste0(input.folder,"/",file.read1))
read2.fastq <- readFastq(dirPath=paste0(input.folder,"/",file.read2))

# Stats:
total.readcount <- length(read1.fastq)
results.table   <- data.frame(category = "Total reads", number = total.readcount, percent = 100, stringsAsFactors=FALSE)

# filter reads based on length
length_filter_read1 <- which( width(read1.fastq) >= minlen.read1 )
length_filter_read2 <- which( width(read2.fastq) >= minlen.read2 )

length_filter_total <- intersect(length_filter_read1, length_filter_read2)

# Stats:
results.table[nrow(results.table)+1, ] <- c("Removed due to trimmed length of Read1", 
                                            number = total.readcount - length(length_filter_read1), 
                                            percent = sprintf("%.1f", (total.readcount - length(length_filter_read1))/total.readcount*100))
results.table[nrow(results.table)+1, ] <- c("Removed due to trimmed length of Read2", 
                                            number = total.readcount - length(length_filter_read2), 
                                            percent = sprintf("%.1f", (total.readcount - length(length_filter_read2))/total.readcount*100))
results.table[nrow(results.table)+1, ] <- c("Total removed due to trimmed length", 
                                            number = total.readcount - length(length_filter_total),
                                            percent = sprintf("%.1f", (total.readcount - length(length_filter_total))/total.readcount*100))

pos <- length_filter_total

## keep reads fulfilling length requirements
read1.fastq <- read1.fastq[pos]
read2.fastq <- read2.fastq[pos]



#############################################
# READ1 stats and filtering                 #
#############################################

##########################
########  subset the reads from 26 to 51+bc nts (for barcode search in this part)
# (subset region: flank + lseq + bc + rseq + flank)

# subset limits
flanklength <- 5
aimed.barcode.start.pos <- 41 # position in read1 where the barcode is expected to start

subset.start.pos <- aimed.barcode.start.pos - nchar(lseq) - flanklength
subset.end.pos   <- aimed.barcode.start.pos + barcode.length + nchar(rseq) + flanklength - 1

# subset
bc.subset.fastq <- subseq(sread(read1.fastq), start=subset.start.pos, end=subset.end.pos)


##########################
########  matchLRPattern function (search for barcode between anchor regions)

# number of mismatches allowed in anchor sequences
max.lseq.mismatch <- 1
max.rseq.mismatch <- 1

# can mismatch be an indel?
indel.lseq <- TRUE
indel.rseq <- TRUE

# find pattern between anchors, report start and end pos. of string (incl. anchors)
match.pattern = sapply(bc.subset.fastq,
                       function(a){matchLRPatterns(Lpattern = lseq,
                                                   Rpattern = rseq,
                                                   max.gaplength = barcode.length + 5,
                                                   subject = a,
                                                   max.Lmismatch = max.lseq.mismatch,
                                                   max.Rmismatch = max.rseq.mismatch,
                                                   with.Lindels = indel.lseq, 
                                                   with.Rindels = indel.rseq)})

# delete no longer need bc.subset.fastq object
rm(bc.subset.fastq)


##########################
###### Trash items with no or multiple views

# find positions with >1 view
pos <- which(elementNROWS(match.pattern)>1)

# select view to keep: 
# width needs to be at least the length of the bc plus the two anchors minus the one
# allowed mismatch in each anchor in case the mismatch can be an indel
# (if more than one view fulfills requirements, it will be removed later)
temp.list <- sapply(match.pattern[pos], 
                    function(a) a[which(width(a) >= barcode.length + width(lseq) - ifelse(indel.lseq,max.lseq.mismatch,0) + width(rseq) - ifelse(indel.rseq,max.rseq.mismatch,0))])
match.pattern[pos] <- temp.list

# Stats: 
# no.views
no.views <- sum(elementNROWS(match.pattern) == 0)
results.table[nrow(results.table)+1, ] <- c("No barcode found", number = no.views, percent = sprintf("%.1f", no.views/total.readcount*100))

## write reads without view, i.e. without found anchor regions, to fastq.gz files
no.view.pos <- which(elementNROWS(match.pattern)==0)

file.read1.nobc <- sub(".fastq.gz$",".noBC.fastq.gz",file.read1)
file.read2.nobc <- sub(".fastq.gz$",".noBC.fastq.gz",file.read2)
writeFastq(read1.fastq[no.view.pos],paste0(output.folder,"/",file.read1.nobc))
writeFastq(read2.fastq[no.view.pos],paste0(output.folder,"/",file.read2.nobc))

# remove no longer needed large object
rm(temp.list)

## remove all views for entries, which still have >1 view, even after width filter
pos <- which(elementNROWS(match.pattern)>1)
if (length(pos) > 0) { match.pattern[pos] <- Views('') }

# Stats: reads with still more views even after width filter
more.views <- length(pos)
results.table[nrow(results.table)+1, ] = c("Multiple barcodes found", number = more.views, percent = sprintf("%.1f", more.views/total.readcount*100))


##########################
######## trim anchors off, barcodes should be left

barcodes <- sapply(match.pattern, function(a){ trimLRPatterns(Lpattern = lseq,
                                                              Rpattern = rseq,
                                                              subject = as(a,"DNAStringSet"),
                                                              max.Lmismatch = max.lseq.mismatch,
                                                              max.Rmismatch = max.rseq.mismatch,
                                                              with.Lindels = indel.lseq, 
                                                              with.Rindels = indel.rseq)})

# Stats: barcode lengths
barcode.width <- sapply(barcodes, width)

# replace those entries w/o length w/ length 0
sel <- which(sapply(barcode.width, length) == 0)
barcode.width[sel] <- 0
barcode.width <- unlist(barcode.width)

# write all barcodes to a file
str.barc <- sapply(barcodes, as.character)
writeLines(unlist(str.barc), con = paste(output.folder,'barcodes.all.txt', sep='/'))

# remove barcodes longer/shorter than barcode.length +/- 1
pos.to.remove <- which( barcode.width < (barcode.length-1) | barcode.width > (barcode.length+1) )

# Stats: unwanted barcodes filter
length.pos.to.remove <- length(pos.to.remove)
removed.for.barcode.length <- length.pos.to.remove - (no.views + more.views)

results.table[nrow(results.table)+1, ] = c(paste("Barcode length not in ",(barcode.length-1),":",(barcode.length+1),sep=""), number = removed.for.barcode.length, percent = sprintf("%.1f", removed.for.barcode.length/total.readcount*100))
results.table[nrow(results.table)+1, ] = c("Total removed for barcode reasons", number = length.pos.to.remove, percent = sprintf("%.1f", length.pos.to.remove/total.readcount*100))

newbarcodes <- barcodes[-pos.to.remove]
str.newbarcodes <- sapply(newbarcodes, as.character)

# write to txt file
writeLines(unlist(str.newbarcodes), con = paste0(output.folder,'/barcodes.',(barcode.length-1),"to",(barcode.length+1),'.txt'))

# remove views of invalid barcodes
match.pattern[pos.to.remove] <- Views('')


##########################
######## get start and end positions of barcodes

# starts
barcode.starts1 <- sapply(match.pattern, function(a){start(a)})

# replace those w/o entries w/ NA
sel <- which(sapply(barcode.starts1, length) == 0)
barcode.starts1[sel] <- NA
barcode.starts <- unlist(barcode.starts1)

# ends
barcode.ends1 <- sapply(match.pattern, function(a){end(a)})

# replace those w/o entries w/ NA
sel <- which(sapply(barcode.ends1, length) == 0)
barcode.ends1[sel] <- NA
barcode.ends <- unlist(barcode.ends1)


# Stats: barcode shifts 
shifts.starts <- sum(barcode.starts!= aimed.barcode.start.pos-width(lseq)-(subset.start.pos-1), na.rm=TRUE)


# find start positions of the trimmed reads (after removing random, primer, anchors and barcode areas)
# (trimmed reads should start after the second anchor (rseq))
trimmed_read_starts <- barcode.ends + subset.start.pos  # start of trimmed read = subset.start.pos - 1 + barcode.ends + 1


##########################
######## filter and trim reads 

# remove not needed fastq entries
read1.fastq.kept    <- read1.fastq[-pos.to.remove]
read2.fastq.kept    <- read2.fastq[-pos.to.remove]
trimmed_read_starts <- trimmed_read_starts[-pos.to.remove]

# trim start of reads
read1.subset.fastq.trimmed <- narrow(read1.fastq.kept, start = trimmed_read_starts)
read2.subset.fastq.trimmed <- narrow(read2.fastq.kept, start = read2.N + read2.bonus.trimmed + 1)

total.kept <- length(read1.subset.fastq.trimmed)
results.table[nrow(results.table)+1, ] = c("Total kept", number = total.kept, percent = sprintf("%.1f", total.kept/total.readcount*100))

results.table[nrow(results.table)+1, ] = c("Shifted barcodes", number = shifts.starts, percent = sprintf("%.1f", shifts.starts/total.kept*100))

# concatenate ids with barcodes
ids.read1 <- mapply(function(myid,mybarcode) paste(myid,mybarcode,sep=":"), id(read1.subset.fastq.trimmed),newbarcodes)
ids.read2 <- mapply(function(myid,mybarcode) paste(myid,mybarcode,sep=":"), id(read2.subset.fastq.trimmed),newbarcodes)
ids.read1 <- as(ids.read1, "BStringSet")
ids.read2 <- as(ids.read2, "BStringSet")

# extract quality strings
qual.read1 <- as(quality(read1.subset.fastq.trimmed), "PhredQuality")
qual.read2 <- as(quality(read2.subset.fastq.trimmed), "PhredQuality")

# combine reads, qualities, ids to ShortReadQ object
Seq.read1 <- ShortReadQ(sread(read1.subset.fastq.trimmed), qual.read1, ids.read1)
Seq.read2 <- ShortReadQ(sread(read2.subset.fastq.trimmed), qual.read2, ids.read2)

# write trimmed reads to fastq files
file.read1.output <- sub(".fastq.gz$",".trimmedBC.fastq.gz",file.read1)
file.read2.output <- sub(".fastq.gz$",".trimmedBC.fastq.gz",file.read2)
writeFastq(Seq.read1, file = paste(output.folder,file.read1.output, sep='/'))
writeFastq(Seq.read2, file = paste(output.folder,file.read2.output, sep='/'))

# write stats table
write.table(results.table, file=paste(output.folder, "stats_table.txt", sep="/"), quote=FALSE, row.names=FALSE, sep="\t")








