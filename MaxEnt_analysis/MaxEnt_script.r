
minigene_seq<-seqinr::read.fasta("~/Seafile/CD19/CD19_WT.minigene.fa")

minigene_seq<-toupper(paste0(minigene_seq$chr, collapse = "") )


library(zoo)

five_ss9mers<-sapply(strsplit(minigene_seq, ""), rollapplyr, 9, paste, collapse = "")[,1]
three_ss23mers<-sapply(strsplit(minigene_seq, ""), rollapplyr, 23, paste, collapse = "")[,1]


library(reticulate)
reticulate:::conda_list()
use_condaenv(condaenv = "anaconda3", required = T)
py_module_available("maxentpy")


## From https://stackoverflow.com/questions/57226087/making-a-list-of-all-mutations-of-a-sequence-dna
## makeSub - create an arbitrary number of substitutions in each 1:n positions
##   with the IUPAC bases specified by 'symbol'
##   return a character vector with all possible variants
##
makeSub <- function(seq, n, symbol = "N")
{
  # IUPAC codes for ambiguous bases
  iupac <- c(N = "ACGT", A = "A", C = "C", G = "G", T = "T", M = "AC", R = "AG",
             W = "AT", S = "CG", Y = "CT", K = "GT", V = "ACG", H = "ACT",
             D = "AGT", B = "CGT")
  
  # accept only a single value for 'seq'
  cseq <- as.character(seq)
  cseq <- unlist(strsplit(cseq[1], ""))
  nseq <- length(cseq)
  
  # simple argument checks
  if (!is(n, "numeric")) stop("'n' must be an integer")
  symbol <- toupper(symbol)
  if (nchar(symbol) != 1 | !symbol %in% names(iupac))
    stop("'symbol' must be a single valid IUPAC symbol")
  if (n == 0) return(paste(cseq, collapse = ""))
  if (n > nseq) stop("too many substitutions for ", nseq, " bases")
  
  # which bases are to be used for the substitution?
  ACGT <- strsplit(iupac[symbol], "")[[1]]
  
  # create all possible combinations of positions to be changed in 'index'
  index <- lapply(seq_len(n), function(j) combn(nseq, j, simplify = FALSE))
  index <- unlist(index, recursive = FALSE)
  
  # for each numeric vector in index, create as many variants as
  # alternative bases are needed, collect in 'ans'
  ans <- lapply(index, function(idx) {
    bases <- lapply(cseq[idx], function(v) setdiff(ACGT, v))
    bases <- bases[sapply(bases, length) > 0] # defensive 
    bases <- expand.grid(bases, stringsAsFactors = FALSE)
    bases <- as.matrix(bases)
    nvars <- nrow(bases)
    
    vars <- do.call(rbind, rep(list(cseq), nvars))
    vars[ ,idx] <- bases
    if (!is.null(vars))
      return(split(vars, seq_len(nvars)))
  })
  ans <- unlist(ans, recursive = FALSE)
  ans <- sapply(ans, paste, collapse = "")
  ans <- unique(ans) # remove duplicates
  return(ans)
}



mutant_variants_9mers<-lapply(five_ss9mers, function(x) makeSub(x,1))
mutant_variants_23mers<-lapply(three_ss23mers, function(x) makeSub(x,1))

maxent<-reticulate::import(module = "maxentpy")


# Calculation scores functions
max_ent5<-function(seq){
  maxent$maxent_fast$score5(seq)
}

max_ent3<-function(seq){
  maxent$maxent_fast$score3(seq)
}


maxent5_scores<-lapply(five_ss9mers, max_ent5)
maxent3_scores<-lapply(three_ss23mers, max_ent3)


mut_variants_5ss <- lapply(mutant_variants_9mers, function(x){
  sapply(x, max_ent5)})

mut_variants_3ss<- lapply(mutant_variants_23mers, function(x){
  sapply(x, max_ent3)})





names(maxent3_scores)<-three_ss23mers
names(mut_variants_3ss)<-three_ss23mers
names(maxent5_scores)<-five_ss9mers
names(mut_variants_5ss)<-five_ss9mers

datalist = list(ss3=maxent3_scores, ss5=maxent5_scores, mut_3ss=mut_variants_3ss, mut5ss=mut_variants_5ss)
saveRDS(datalist, "~/Seafile/CD19/MaxEnt_values_with_single_mutant_scores.rds")
rm(list=ls())


