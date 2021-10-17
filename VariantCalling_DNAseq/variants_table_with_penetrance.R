# Create barcodes table 

library(data.table)
library(dplyr)
vcf_var_count<-fread("/fsimb/groups/imb-koeniggr/mariela/CD19/scripts/merged_vcf_variant_count.txt", header = F,col.names = c("barcode", "variant_number"))

non_empty_vcfs<-(vcf_var_count)%>%filter(variant_number>0)

library("vcfR", lib.loc = "~/R_libs/")
library("memuse", lib.loc = "~/R_libs/")
vcf_dir<-"/fsimb/groups/imb-koeniggr/mariela/CD19/scripts/merged_vcfs/"

barcode_counts1<-fread("/fsimb/groups/imb-koeniggr/mariela/CD19/scripts/dna_demultiplex_101018/results/barcode.counts.read.tsv", header =F)
barcode_counts2<-fread("/fsimb/groups/imb-koeniggr/mariela/CD19/scripts/dna_demultiplex_311018/results/barcode.counts.read.tsv", header =F)

barcode_counts<-rbind(barcode_counts1, barcode_counts2)
colnames(barcode_counts)<-c("barcode", "reads")
barcode_counts<-barcode_counts[,list(reads=sum(reads)), by=barcode]
barcode_counts_f<-filter(barcode_counts,reads>3)

valid_bc<-merge(barcode_counts_f, non_empty_vcfs, by="barcode")

vcf_info_list<-function(vcf,bc){
  fvcf<-read.vcfR(vcf,verbose = F)
  sfdata<-vcfR2tidy(fvcf, single_frame = T)$dat%>%
    mutate(BARCODE=bc)%>%
    dplyr::select(c("BARCODE","POS","REF", "ALT", "AF", "QUAL","gt_AD", "gt_DP","gt_GT_alleles"))
  return(sfdata)
}



vcf_table<-lapply(seq(1:nrow(valid_bc)), function(i){
  print(i)
  vcf<-paste0(vcf_dir, valid_bc[i,"barcode"], ".vcf.gz")
  vcf_info_list(vcf,valid_bc[i,"barcode"])
})%>%data.table::rbindlist()


#saveRDS(vcf_table, file="variants_table_prefilters_Nov29.rds")
saveRDS(vcf_table, file="variants_table_prefilters_Dec10.rds")


