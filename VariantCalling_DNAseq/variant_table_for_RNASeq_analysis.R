# Generate variants table usindg the second method------------------------------------------------------

# Read file
var_m2<-readRDS("/fsimb/groups/imb-koeniggr/mariela/CD19/scripts/2018_09/rds/variants_3rd_method.rds")


# Remove low frequency variants
library(dplyr)
library(tidyr)


# Calculate bad count ratio 
var_m2_rec_af<-var_m2%>%mutate(LOW_AF= (AF < 0.8))%>%
  group_by(BARCODE, LOW_AF)%>%
  summarise(n = n()) %>%
  spread(LOW_AF,n) %>% dplyr::select(BARCODE,`TRUE`) %>% rename("n_low_p_var"=`TRUE`) %>%
  mutate_if(is.numeric , replace_na, replace = 0)%>%
  left_join(var_m2, ., by="BARCODE") %>%
  rename("n_low_p_var.rec"=n_low_p_var.x) %>%
  mutate(bad_count_ratio=n_low_p_var.rec/n_variants) %>%
  dplyr::select(-n_low_p_var.y)


no_wt_bc<-var_m2_rec_af$BARCODE%>%unique()%>%length()

total_called_var<-var_m2_rec_af%>%nrow() 

# Select all the low penetrance barcodes
low_p_bc<-var_m2_rec_af%>%
  arrange(desc(n_low_p_var.rec)) %>% 
  filter(bad_count_ratio >= 0.25  |  n_low_p_var.rec>1)

# Discart low penetrance barcodes

var_m2_rec_af %>% filter(!BARCODE%in%low_p_bc$BARCODE)


var_m1<-readRDS("/fsimb/groups/imb-koeniggr/mariela/CD19/scripts/2018_09/rds/variants_1st_method.rds")



differential_variants<-rbind(anti_join(var_m1, var_m2, by=c("BARCODE","POS","REF","ALT")) %>% mutate(specific="default"),
                             anti_join(var_m2, var_m1, by=c("BARCODE","POS","REF","ALT"))%>% mutate(specific="kmer"))




# Adding the missing variants from the first method 

missing_first<-differential_variants%>% filter(specific=="default") %>% dplyr::select(-specific)


# select the columns that we want

library(Biostrings)

final_mut_table<-rbind(var_m2_rec_af %>% filter(!BARCODE%in%low_p_bc$BARCODE) %>% 
  select(BARCODE,POS,REF, ALT, AF),
missing_first%>% 
  select(BARCODE,POS,REF, ALT, AF)) %>% 
  arrange(BARCODE)

final_mut_table$RNA_BARCODE<-as.character(reverseComplement(DNAStringSet(final_mut_table$BARCODE)))

final_mut_table<-final_mut_table %>% rename(AF="PENETRANCE")


write.table(final_mut_table,"mutations_individual_variants_110319.tab", row.names = F,quote = F, sep = "\t")

