library(ShortRead)
fq<-readFastq("~/Seafile/CD19/NOVARTIS_paper_BAMs/fastq/10S_RNA_1.fastq.gz")


library(Biostrings)

Biostrings::chartr("N", "", sread(fq))
Ns<-vmatchPattern("N", sread(fq))


sread(fq[1][-Ns[[1]]])


min_ranges<-lapply(Ns, reduce)


table(sapply(min_ranges, length) )
