###################################################
#
# script adapted by Mariela Cortés-López based on
# first version written by Samarth T. Setty and Anke Busch
#
###################################################

#######  load the library and declare the flanking RE sites as strings for matching
library(R.utils)
library(Biostrings)
library(ShortRead)
options(stringsAsFactors = FALSE)

# read arguments
args <- commandArgs(asValues = TRUE)

# check if args are properly assigned
print(args)


###### get input arguments
input.folder <- args$inputPath
file.read <- args$Read
rds.folder <- args$rdsPath
output.folder <- args$outputPath
barcode.length <- as.integer(args$barcodelength)
minlen.read <- as.integer(args$minlengthread)
maxlen.read <- as.integer(args$maxlengthread)
lseq <- args$LAnchor
rseq <- args$RAnchor


###### global variables (WHERE IS seqtype USED?)
seqtype <- "DNASeq"

# anchor sequences around barcode


# Read the fastq file
read.fastq <- readFastq(dirPath = paste0(input.folder, "/", file.read))


# Stats1:
total.readcount <- length(read.fastq)
results.table <- data.frame(
  category = "Total reads",
  number = total.readcount,
  percent = 100, stringsAsFactors = FALSE
)


## filter reads based on length
length_filter_read <- which(width(read.fastq) >= minlen.read & width(read.fastq) <= maxlen.read)

results.table[nrow(results.table) + 1, ] <- c("Removed due to trimmed length of Read",
  number = total.readcount - length(length_filter_read),
  percent = sprintf("%.1f", (total.readcount - length(length_filter_read)) / total.readcount * 100)
)


pos <- length_filter_read
pos

saveRDS(object = pos, file = paste0(rds.folder, "pos.removed.due.to.read.length.rds"))

## remove reads not fulfilling length requirements
read.fastq <- read.fastq[pos]


# Stats2:
thrown.out <- total.readcount - length(read.fastq)
results.table[nrow(results.table) + 1, ] <- c("Total removed due to trimmed length",
  number = thrown.out,
  percent = sprintf("%.1f", thrown.out / total.readcount * 100)
)

## check 1
if (length(pos) != length(read.fastq)) {
  print("pos:")(length(pos)("read1:")(length(read.fastq)))
}
head(read.fastq)




#############################################
# READ1 stats and filtering                 #
#############################################

########  subset the reads from 26 to 51+bc nts (for barcode search in this part)
## (subset region: flank - lseq - bc - rseq - flank)

## subset limits
window.to.search <- 50
## subset


bc.subset.fastq <- subseq(sread(read.fastq), -window.to.search)
bc.subset.fastq

########  matchLRPattern function (search for barcode between anchor regions)

## number of mismatches allowed in anchor sequences
max.lseq.mismatch <- 1
max.rseq.mismatch <- 1

## can mismatch be an indel?
indel.lseq <- TRUE
indel.rseq <- TRUE

## find pattern between anchors, report start and end pos. of string (incl. anchors)
# Only allowing 3 or more nucleotides of match
match.pattern <- sapply(
  bc.subset.fastq,
  function(a) {
    matchLRPatterns(
      Lpattern = lseq,
      Rpattern = rseq,
      max.gaplength = barcode.length + 3,
      subject = a,
      max.Lmismatch = max.lseq.mismatch,
      max.Rmismatch = max.rseq.mismatch,
      with.Lindels = indel.lseq,
      with.Rindels = indel.rseq
    )
  }
)

## check the list
head(match.pattern)

###### Trash items with no or multiple views

# Complementary strand search for reads without a match in the sense strand-------------------

no.view.pos <- which(elementNROWS(match.pattern) == 0)
no.views <- length(no.view.pos)
no.views


read.fastq[no.view.pos] <- reverseComplement(read.fastq[no.view.pos])

bc.revcom.subset.fastq <- subseq(sread(read.fastq[no.view.pos]), -window.to.search)

match.pattern.revcomp <- sapply(
  bc.revcom.subset.fastq,
  function(a) {
    matchLRPatterns(
      Lpattern = lseq,
      Rpattern = rseq,
      max.gaplength = barcode.length + 3,
      subject = a,
      max.Lmismatch = max.lseq.mismatch,
      max.Rmismatch = max.rseq.mismatch,
      with.Lindels = indel.lseq,
      with.Rindels = indel.rseq
    )
  }
)


match.pattern[no.view.pos] <- match.pattern.revcomp

saveRDS(object = match.pattern, file = paste(rds.folder, "match.pattern.rds", sep = "/"))

pos <- which(elementNROWS(match.pattern) > 1)
more.views <- length(pos)
more.views

temp.list <- sapply(
  match.pattern[pos],
  function(a) a[which(width(a) >= barcode.length + width(lseq) - ifelse(indel.lseq, max.lseq.mismatch, 0) + width(rseq) - ifelse(indel.rseq, max.rseq.mismatch, 0))]
)

new.match.pattern <- match.pattern
new.match.pattern[pos] <- temp.list


# Stats 3: reads with more views (before correction)
more.views <- length(pos)

# find the view to keep (what if several fulfil requirements? ---> check and remove those later)
# seq.length (width(a) needs to be at least the length of the bc plus the two anchor minus the one
# allowed mismatch in each anchor in case the mismatch can be an indel)
# temp.list = sapply(match.pattern[pos], function(a)a[which(width(a)>=barcode.length+10)])

no.views <- sum(elementNROWS(new.match.pattern) == 0)
no.views

results.table[nrow(results.table) + 1, ] <- c("No barcode found",
  number = no.views,
  percent = sprintf("%.1f", no.views / total.readcount * 100)
)


results.table
## write reads without view aka without found anchor regions to fastq.gz files
no.view.pos <- which(elementNROWS(new.match.pattern) == 0)
no.view.pos

file.read.nobc <- sub(".fastq$", ".noBC.fastq", file.read)
writeFastq(read.fastq[no.view.pos], paste0(output.folder, "/", file.read.nobc))


# check 3
if (length(new.match.pattern) != length(match.pattern)) {
  print("new.match.pattern:")(length(new.match.pattern)("match.pattern:")(length(match.pattern)))
}
head(new.match.pattern)

### REMEMBER TO REMOVE LARGE UNNEEDED OBJECTS (e.g. match.pattern, temp.list)
rm(match.pattern, match.pattern.revcomp, temp.list)

## remove all views for entries, which still have >1 view, even after width filter
pos <- which(elementNROWS(new.match.pattern) > 1)
if (length(pos) != 0) {
  new.match.pattern[pos] <- Views("")
}
# Stats 5: reads with still more views even after width filter!!!
more.plus.views <- length(pos)
results.table[nrow(results.table) + 1, ] <- c("Multiple barcodes found", number = more.plus.views, percent = sprintf("%.1f", more.plus.views / total.readcount * 100))

results.table
saveRDS(object = new.match.pattern, file = paste(rds.folder, "new.match.pattern.rds", sep = "/"))


##########################
######## trim anchors off and barcodes should be left

barcodes <- sapply(new.match.pattern, function(a) {
  trimLRPatterns(
    Lpattern = lseq,
    Rpattern = rseq,
    subject = as(a, "DNAStringSet"),
    max.Lmismatch = max.lseq.mismatch,
    max.Rmismatch = max.rseq.mismatch,
    with.Lindels = indel.lseq,
    with.Rindels = indel.rseq
  )
})

# check 4:
head(barcodes)

# Stats 6: barcode lengths
barcode.width <- sapply(barcodes, width)

### ADDED CODE TO CATCH integer(0)
sel <- which(sapply(barcode.width, length) == 0)
barcode.width[sel] <- 0
barcode.width <- unlist(barcode.width)

saveRDS(object = barcodes, file = paste(rds.folder, "barcodes.rds", sep = "/"))

## all barcodes list into one file
str.barc <- sapply(barcodes, as.character)
writeLines(unlist(str.barc), con = paste(output.folder, "barcodes.all.txt", sep = "/"))

### remove barcodes longer/shorter than barcode.length+-1
pos.to.remove <- which(barcode.width < (barcode.length - 1) | barcode.width > (barcode.length + 1))

# Stats 7: unwanted barcodes filter
length.pos.to.remove <- length(pos.to.remove)
removed.for.barcode.length <- length.pos.to.remove - (no.views + more.plus.views)

results.table[nrow(results.table) + 1, ] <- c(paste("Barcode length not in ", (barcode.length - 1), ":", (barcode.length + 1), sep = ""), number = removed.for.barcode.length, percent = sprintf("%.1f", removed.for.barcode.length / total.readcount * 100))
results.table[nrow(results.table) + 1, ] <- c("Total removed for barcode reasons", number = length.pos.to.remove, percent = sprintf("%.1f", length.pos.to.remove / total.readcount * 100))

newbarcodes <- barcodes[-pos.to.remove]
# save to file
saveRDS(object = newbarcodes, file = paste0(rds.folder, "/barcodes.", (barcode.length - 1), "to", (barcode.length + 1), ".rds"))

str.newbarcodes <- sapply(newbarcodes, as.character)
# write to txt file
writeLines(unlist(str.newbarcodes), con = paste0(output.folder, "/barcodes.", (barcode.length - 1), "to", (barcode.length + 1), ".txt"))

# so that starts and ends, shifts which are not in barcodelengths 9:11 can be counted as no views
new.match.pattern[pos.to.remove] <- Views("")


# for read2 subsetting
saveRDS(object = pos.to.remove, file = paste(rds.folder, "pos.to.remove.rds", sep = "/"))



######## save start end positions of barcodes for later MSA

# starts
barcode.starts1 <- sapply(new.match.pattern, function(a) {
  start(a)
})

## added later check again?
### CAPTURE integer(0)
sel <- which(sapply(barcode.starts1, length) == 0 | is.na(sapply(barcode.starts1, length)))
barcode.starts1[sel] <- NA

barcode.starts <- unlist(barcode.starts1)

# check 5:
if (length(barcode.starts) != length(barcode.starts1)) {
  print(paste("barcode.starts:", length(barcode.starts), "barcodes.starts1:", length(barcode.starts1)))
} #### THIS WILL GIVE FALSE
head(barcode.starts)

# Stats 8: shifted barcodes - is this necessary?
barcode.starts_shifts <- table(barcode.starts)
barcode.starts_shifts

new.barcode.starts <- barcode.starts
#### new.barcode.starts[pos.to.remove]=NA ### THIS WILL NOT WORK SINCE THE UNLIST ABOVE REMOVED EMPTY ELEMENTS
saveRDS(object = new.barcode.starts, file = paste(rds.folder, "barcode.starts.rds", sep = "/"))

# end
barcode.ends1 <- sapply(new.match.pattern, function(a) {
  end(a)
})

### CAPTURE integer(0)
sel <- which(sapply(barcode.ends1, length) == 0 | is.na(sapply(barcode.starts1, length)))
barcode.ends1[sel] <- NA

barcode.ends <- unlist(barcode.ends1)

# check 6:
if (length(barcode.ends) != length(barcode.ends1)) {
  print(paste("barcode.ends:", length(barcode.ends), "barcode.ends1:", length(barcode.ends1)))
}
head(barcode.ends)

## new.barcode.ends=barcode.ends
saveRDS(object = barcode.ends, file = paste(rds.folder, "barcode.ends.rds", sep = "/"))



###### Trimming --------------------------------------------------------------


trimmed_read_ends <- barcode.starts - 1 # start of trimmed read = subset.start.pos - 1 + barcode.ends + 1

# saveRDS(object = trimmed_read_starts,file =paste(rds.folder,'trimmed_read_starts.rds', sep='/'))


# check 8:
if (length(trimmed_read_ends) != length(new.match.pattern)) {
  print(paste("trimmed_read_starts:", length(trimmed_read_ends), "new.match.pattern:", length(new.match.pattern)))
}
head(trimmed_read_ends)

## remove unneeded fastq entries
read.fastq.kept <- read.fastq[-pos.to.remove]

trimmed_read_ends <- trimmed_read_ends[-pos.to.remove]
read.subset.fastq.trimmed <- narrow(read.fastq.kept, end = -trimmed_read_ends - 1)

total.kept <- length(read.subset.fastq.trimmed)

total.kept
results.table[nrow(results.table) + 1, ] <- c("Total kept", number = total.kept, percent = sprintf("%.1f", total.kept / total.readcount * 100))

## concatenate ids with barcodes
ids.read <- mapply(
  function(myid, mybarcode) paste(myid, mybarcode, sep = ":"),
  ShortRead::id(read.subset.fastq.trimmed), newbarcodes
)
ids.read <- as(ids.read, "BStringSet")

qual.read <- as(quality(read.subset.fastq.trimmed), "PhredQuality")


Seq.read <- ShortReadQ(sread(read.subset.fastq.trimmed), qual.read, ids.read)



# check 9: (should always return TRUE, i.e. not equal and print)
if (length(Seq.read) != length(new.match.pattern)) {
  cat("Seq.read1:", length(Seq.read), "new.match.pattern:", length(new.match.pattern))
}

# check 10:
if (length(Seq.read) != length(read.subset.fastq.trimmed)) {
  cat("Seq.read1:", length(Seq.read), "read1.subset.fastq.trimmed:", length(read1.subset.fastq.trimmed))
}

rds.file.read <- sub(".fastq.gz$|.fastq$", ".trimmedBC.rds", file.read)

saveRDS(object = Seq.read, file = paste(rds.folder, rds.file.read, sep = "/"))


file.read.output <- sub(".fastq.gz$|.fastq$", ".trimmedBC.fastq", file.read)

writeFastq(Seq.read, file = paste(output.folder, file.read.output, sep = "/"), compress = T)


write.table(results.table, file = paste(output.folder, "stats_table.txt", sep = "/"), quote = FALSE, row.names = FALSE, sep = "\t")
