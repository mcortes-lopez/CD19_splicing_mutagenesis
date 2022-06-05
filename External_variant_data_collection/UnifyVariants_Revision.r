library(data.table)
library(dplyr)
library(tidyr)

hg38_coord <- seq(28931871,28933144) 
hg38_coord_st <- paste0("16:", hg38_coord) 

library("BSgenome.Hsapiens.UCSC.hg38")
seq_cd19 <- strsplit(as.character(Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg38, 
                                         "chr16",28931871, 28933144)), "")[[1]]

names(seq_cd19) <- hg38_coord 
# Generate a single variant table

# ENSEMBL -------------------------------------------------------------------------------------------
ensembl <- fread("data/SNPs/Annotated_variants/ensembl-export.csv")
colnames(ensembl) <- make.unique(colnames(ensembl))

ensembl <- ensembl %>% 
  dplyr::select(`Location`,`Variant ID`, `Conseq. Type`, vf_allele, Alleles) %>% 
  filter( Alleles != vf_allele & Location %in% hg38_coord_st )  %>% 
  unique() %>% 
  separate_rows(Alleles) %>% 
  separate_rows(vf_allele) %>% 
  unique()%>% 
  mutate(POS = gsub("16:", "", Location)) %>% 
  filter(Alleles != vf_allele) %>% 
  #dplyr::rename(REF = vf_allele, 
   #      ALT = Alleles) %>% 
  dplyr::select(-Location) %>%
  mutate(REF = seq_cd19[as.character(POS)]) %>% 
  filter(!(Alleles != REF  & nchar(Alleles) == 1)) %>%
  dplyr::select(-REF) %>%
  dplyr::rename(REF = Alleles, 
        ALT = vf_allele) 

# COSMIC ---------------------------------------------------------------------------------------------

cosmic_coding <- fread("data/SNPs/Annotated_variants/V95_38_MUTANT_CD19_COSMIC.csv")
cosmic_noncoding <- fread("data/SNPs/Annotated_variants/V95_COSMIC_NCV_CD19_reg.csv", sep = "\t")

cosmic_coding <-cosmic_coding %>%
  dplyr::select(MUTATION_GENOME_POSITION, GENOMIC_MUTATION_ID, HGVSG)

cosmic_noncoding <-cosmic_noncoding %>%
  dplyr::select(genome_position, GENOMIC_MUTATION_ID, HGVSG) %>% 
  dplyr::rename(MUTATION_GENOME_POSITION = genome_position)

cosmic <- rbind(cosmic_coding, cosmic_noncoding) %>% 
  unique()

cosmic_vcf<- vcfR::read.vcfR("data/2022-03_variants/vcf_files/COSMIC_vcf.vcf")
cosmic_vcf <- cosmic_vcf@fix %>% 
  as.data.frame() %>% 
  dplyr::select(POS,ID,REF,ALT)

cosmic_vep <- fread("data/2022-03_variants/COSMIC_VEP.txt") 

cosmic <- left_join(cosmic, cosmic_vep, by = c("HGVSG" = "Uploaded_variation")) %>%
  dplyr::select(MUTATION_GENOME_POSITION, GENOMIC_MUTATION_ID, HGVSG, Consequence,Allele , CLIN_SIG, Existing_variation) %>% 
  unique()

cosmic <- left_join(cosmic, cosmic_vcf, by = c("HGVSG" = "ID"))

cosmic <- cosmic %>%
  dplyr::select(GENOMIC_MUTATION_ID, Consequence, REF, ALT, POS, Existing_variation) %>% 
  dplyr::rename(`Variant ID` = GENOMIC_MUTATION_ID,
                `Conseq. Type` = Consequence)
  

# GNOMAD -------------------------------------------------------------------------------------------------------

gnomad <- fread("data/SNPs/Annotated_variants/gnomAD_v3.1.1_ENSG00000177455_2021_09_09_12_33_20.csv")

#  `Variant ID` `Conseq. Type`      REF   ALT   POS 

gnomad <- gnomad  %>% 
  dplyr::select(rsIDs, `VEP Annotation`, Reference, Alternate, Position, `ClinVar Clinical Significance`) %>% 
  dplyr::rename(`Variant ID` = rsIDs, 
                `Conseq. Type` = `VEP Annotation`,
                REF = Reference, 
                ALT = Alternate, 
                POS = Position, 
                `Clinical Significance` = `ClinVar Clinical Significance`)%>% 
  unique() 


# ClinVar ------------------------------------------------------------------------------------------------------

clinvar <- fread("data/2022-03_variants/CLINVAR_VEP.txt")
clinvar <- clinvar %>% 
  dplyr::select(Location, Uploaded_variation,  Consequence,Allele , CLIN_SIG, Existing_variation) %>% 
  mutate(Uploaded_variation=as.character(Uploaded_variation)) %>% 
  unique()

clinvar_vcf <- vcfR::read.vcfR("data/2022-03_variants/vcf_files/clinvar_vcf.vcf")
clinvar_vcf <- clinvar_vcf@fix %>% 
  as.data.frame() %>% 
  dplyr::select(POS,ID,REF,ALT) %>%
  mutate(ID=as.character(ID))

clinvar <- left_join(clinvar, clinvar_vcf, by = c("Uploaded_variation" = "ID"))

clinvar <- clinvar %>%
  dplyr::select(Uploaded_variation, Consequence, REF, ALT, POS, Existing_variation, CLIN_SIG) %>% 
  dplyr::rename(`Variant ID` = Uploaded_variation,
                `Conseq. Type` = Consequence, 
                `Clinical Significance` = CLIN_SIG)

#  `Variant ID` `Conseq. Type`      REF   ALT   POS 


# ORLANDO ------------------------------------------------------------------------------------

orlando <- fread("data/2022-03_variants/ORLANDO_VEP.txt")
orlando <- orlando %>% 
  dplyr::select(Location, Uploaded_variation,  Consequence,Allele , CLIN_SIG, Existing_variation) %>% 
  mutate(Uploaded_variation=as.character(Uploaded_variation)) %>%
  unique()

orlando_vcf <- vcfR::read.vcfR("data/2022-03_variants/vcf_files/ORLANDO_VCF.vcf")
orlando_vcf <- orlando_vcf@fix %>% 
  as.data.frame() %>% 
  dplyr::select(POS,ID,REF,ALT) %>%
  mutate(ID=as.character(ID)) 

orlando <- left_join(orlando, orlando_vcf, by = c("Uploaded_variation" = "ID"))

orlando <- orlando %>%
  dplyr::select(Uploaded_variation, Consequence, REF, ALT, POS, Existing_variation, CLIN_SIG) %>% 
  dplyr::rename(`Variant ID` = Uploaded_variation,
                `Conseq. Type` = Consequence, 
                `Clinical Significance` = CLIN_SIG)



# TARGET --------------------------------------------------------------------------------------

library(readxl)
target <- read_xls("../Supplementary_Tables/data/CD19_variants_04FEB22.xls")



target <- target %>%
  dplyr::select(ID, variant_type, REF, ALT, POS, variant_type, variant_cat) %>% 
  dplyr::rename(`Variant ID` = ID,
                `Conseq. Type` = variant_type, 
                `Clinical Significance` = variant_cat) %>%
  unique()


merged_variants <- rbindlist(list("Ensembl" = ensembl, 
                                  "ClinVar" = clinvar, 
                                  "gnomAD" = gnomad, 
                                  "COSMIC" = cosmic, 
                                  "TARGET" = target, 
                                  "Orlando et al." = orlando),
                             idcol = "set", fill = T ) %>% 
  filter(POS != ""  & POS %in% hg38_coord)


merged_variants <- merged_variants %>% 
  group_by(POS, REF, ALT) %>% 
  summarize( set = stringr::str_c(unique(set), collapse = ","), 
             `Variant ID`= stringr::str_c(unique(`Variant ID`), collapse = ","), 
             `Conseq. Type` = stringr::str_c(unique(`Conseq. Type`), collapse = ",")) %>% 
  ungroup()

writexl::write_xlsx(merged_variants, "output/Resource_CD19_mutations_Exons1-3.xlsx")


### Numbers ---------------------------------------------------------------------
# Some variants were manually cleaned 
library(readxl)

merged_variants <- read_xlsx("data/Resource_CD19_mutations_Exons1-3.xlsx")

# Indels
merged_variants %>% 
  filter(nchar(as.character(REF))>1 | nchar(as.character(ALT)) >1) %>% 
  nrow()

# SNPs
merged_variants %>% 
  filter(nchar(as.character(REF)) == 1 & nchar(as.character(ALT)) == 1) %>% 
  nrow()

# Indels


indels_gen <- merged_variants %>% 
  filter(nchar(as.character(REF))>1 | nchar(as.character(ALT)) >1) %>% 
  pull(set) %>% 
  strsplit(., ",")%>% unlist() %>% 
  table() %>% 
  as.data.frame()

indels_plus5 <- merged_variants %>% 
  filter(nchar(as.character(REF))>5 | nchar(as.character(ALT)) >5) %>% 
  pull(set) %>% 
  strsplit(., ",")%>% unlist() %>% 
  table() %>%
  as.data.frame()

snvs <- merged_variants %>% 
  filter(nchar(as.character(REF)) == 1 & nchar(as.character(ALT)) ==  1) %>% 
  pull(set) %>% 
  strsplit(., ",")%>% unlist() %>% 
  table() %>% 
  as.data.frame()

all_vars <- merged_variants %>% 
  pull(set) %>% 
  strsplit(., ",")%>% unlist() %>% 
  table() %>% 
  as.data.frame( ) 

colnames(all_vars) <- c("set", "Total")
colnames(snvs) <- c("set", "SNVs")
colnames(indels_gen) <- c("set", "Indel (>1)")
colnames(indels_plus5) <- c("set", "Indel (>5)")

summary_tab <- Reduce( function(...) merge(...,all = T, by="set"), list( all_vars, snvs, indels_gen, indels_plus5))

sum_tab_2save<- summary_tab %>% 
  mutate(SNVs_p= round(SNVs/Total*100,1),
         `Indel (>1) p` = round(`Indel (>1)`/Total *100, 1), 
         `Indel (>5) p` = round(`Indel (>5)`/Total *100,1)) %>% 
  mutate( SNVs = paste0(SNVs, " (", SNVs_p, "%)") , 
          `Indel (>1)` = paste0(`Indel (>1)`, " (", `Indel (>1) p`, "%)") ,
          `Indel (>5)` = paste0(`Indel (>5)`, " (", `Indel (>5) p`, "%)")) %>% 
  dplyr::select(-contains("p")) 

fwrite(sum_tab_2save, "output/Resource_summary_tab.txt", sep = "\t", quote = F)

sum_tab_2save
dim(merged_variants)
