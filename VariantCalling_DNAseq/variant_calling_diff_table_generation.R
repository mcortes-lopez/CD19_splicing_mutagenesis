# Create barcodes table 

library(data.table)
library(dplyr)

#complicated.vcf<-read.vcf("/fsimb/groups/imb-koeniggr/mariela/CD19/scripts/merged_vcfs_2/ACCTTCTAATCCACT.vcf.gz")

variant_count<-fread("/fsimb/groups/imb-koeniggr/mariela/CD19/scripts/variant_count.diffpar.comparison.tab")
variant_count<-variant_count[,-3]
colnames(variant_count)<-c("barcode","variants.default","variants.kmer.ploidy")

variant_count$diff<-variant_count[,"variants.kmer.ploidy"]-variant_count[,"variants.default"]


diff.bc<-(variant_count)%>% filter(variants.kmer.ploidy>0 & variants.default>0)
diff.bc<-diff.bc%>% filter(barcode!="ACCTTCTAATCCACT")


library("vcfR", lib.loc = "~/R_libs/")
library("memuse", lib.loc = "~/R_libs/")
library("memuse", lib.loc = "~/R_libs/")

vcf1_dir<-"/fsimb/groups/imb-koeniggr/mariela/CD19/scripts/merged_vcfs/"
vcf2_dir<-"/fsimb/groups/imb-koeniggr/mariela/CD19/scripts/merged_vcfs_2/"



vcf_info_list<-function(vcf,bc){
  fvcf<-read.vcfR(vcf,verbose = T,cols = c(1:10), nrows = 30)
  sfdata<-vcfR2tidy(fvcf, single_frame = T )$dat%>%
    mutate(BARCODE=bc)%>%
    dplyr::select(c("BARCODE","POS","REF", "ALT", "AF", "QUAL","gt_AD", "gt_DP","gt_GT_alleles"))
  return(sfdata)
}





vcf_table1<-lapply(seq(1:nrow(diff.bc)), function(i){
  print(i)
  vcf1<-paste0(vcf1_dir, diff.bc[i,"barcode"], ".vcf.gz")
  vcf_info_list(vcf1,diff.bc[i,"barcode"])
  
})%>%data.table::rbindlist()




# 4 READS BARCODE FILTER

barcode_counts1<-fread("/fsimb/groups/imb-koeniggr/mariela/CD19/scripts/dna_demultiplex_101018/results/barcode.counts.read.tsv", header =F)
barcode_counts2<-fread("/fsimb/groups/imb-koeniggr/mariela/CD19/scripts/dna_demultiplex_311018/results/barcode.counts.read.tsv", header =F)

dnaseq.ccs.count<-rbind(barcode_counts1, barcode_counts2)
colnames(dnaseq.ccs.count)<-c("barcode", "reads")
dnaseq.ccs.count<-dnaseq.ccs.count[,list(reads=sum(reads)), by=barcode]

total_bc<-nrow(dnaseq.ccs.count)
bc_4m<-filter(dnaseq.ccs.count, reads>3) %>% nrow()

bc_4m


vcf2.4reads_cutoff<-(diff.bc)%>% filter(barcode%in%filter(dnaseq.ccs.count, reads>3)$barcode)

vcf2.4reads_cutoff<-vcf2.4reads_cutoff%>% filter(barcode!="ACCTTCTAATCCACT")%>% unique()

vcf_table2<-lapply(seq(1,nrow(vcf2.4reads_cutoff)), function(i){
  print(i)
  vcf2<-paste0(vcf2_dir, vcf2.4reads_cutoff[i,"barcode"], ".vcf.gz")
  vcf_info_list(vcf2,vcf2.4reads_cutoff[i,"barcode"])
})%>%data.table::rbindlist()


library(VariantAnnotation)

complicated.vcf<-readVcf("/fsimb/groups/imb-koeniggr/mariela/CD19/scripts/merged_vcfs_2/ACCTTCTAATCCACT.vcf.gz")







variants_second_method<-rbind(vcf_table2%>%filter(!is.na(AF)),data.frame(BARCODE="ACCTTCTAATCCACT",
              POS=ranges(complicated.vcf@rowRanges)%>% start(),
              REF=fixed(complicated.vcf)$REF,
              ALT=as.character(fixed(complicated.vcf)$ALT%>% unlist()),
              AF=info(complicated.vcf)$AF%>% unlist() %>% as.vector(), 
              QUAL=complicated.vcf@fixed$QUAL,
              gt_AD=unstrsplit(CharacterList(geno(complicated.vcf)$AD), ","),
              gt_DP=geno(complicated.vcf)$DP[,1]%>% as.vector(),
              gt_GT_alleles=unlist(genotypeToSnpMatrix(complicated.vcf)[[2]]["allele.2"]) %>% CharacterList() %>% unstrsplit()%>% as.vector(),
              stringsAsFactors = F))


library(splitstackshape)
var_called_filter<-variants_second_method %>% 
  mutate(gt_AD=gsub("^\\d+\\,","",gt_AD)) %>%
  as.data.table()%>%
  cSplit(., splitCols = c("ALT", "AF", "gt_AD"), sep = ",", direction = "long",makeEqual = F) %>%
  group_by(BARCODE, POS) %>% 
  mutate(REC_AF=as.numeric(gt_AD)/ max(as.numeric(gt_DP)))%>%
  filter(REC_AF == max(REC_AF))%>%
  ungroup()




saveRDS(var_called_filter,file = "/fsimb/groups/imb-koeniggr/mariela/CD19/scripts/2018_09/rds/variants_2nd_method.rds")

bc_to_search<-seq(1:nrow(diff.bc))
bc_to_search.list<-split(bc_to_search, ceiling(seq_along(bc_to_search)/1000))
specific_mut_tab<-lapply(bc_to_search.list, function(index){
        lapply(index,function(i){
          print(diff.bc[i,"barcode"])
          vcf1<-paste0(vcf1_dir, diff.bc[i,"barcode"], ".vcf.gz")
          vcf2<-paste0(vcf2_dir, diff.bc[i,"barcode"], ".vcf.gz")
          var_set1<-vcf_info_list(vcf1,diff.bc[i,"barcode"])
          var_set2<-vcf_info_list(vcf2,diff.bc[i,"barcode"])
          rbind(anti_join(var_set1, var_set2, by=c("BARCODE", "POS", "REF")) %>% mutate(specific="default"),
               anti_join(var_set2, var_set1, by=c("BARCODE", "POS", "REF")) %>% mutate(specific="kmer10ploidy1"))
          })  %>%data.table::rbindlist()})
  
  



specific_mut_tab.final<-specific_mut_tab%>%data.table::rbindlist()

View(specific_mut_tab.final)


saveRDS(specific_mut_tab.final,file = "/fsimb/groups/imb-koeniggr/mariela/CD19/scripts/2018_09/rds/differential_variants.rds")



library(ggplot2)
colnames(specific_mut_tab.final)
ggplot(specific_mut_tab.final, aes(POS, fill=specific))+ 
  geom_histogram(binwidth =10)
