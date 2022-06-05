library(data.table)
library(dplyr)
combined_score <- readRDS("data/combined_scores_tidy.rds")
hg38_coord <- seq(28931871,28933144) 

prevalence_variants <- combined_score %>%
  filter(combined_score > 0.25) %>%
  dplyr::select(POS, ALT, REF) %>%
  unique() %>% 
  mutate(hg38_POS = hg38_coord[POS], 
          mut_minigene_code = paste0(ALT,POS,REF))

prevalence_variants

isoform_df <- readRDS("../Mutation_effects/Softmax_modelling/distribution_cutoffs.rds")

sig_mut <- isoform_df %>% 
  mutate(effective_mutation = ifelse(reps_cutoff_empirical == 0, "no", "yes")) %>% 
  filter(effective_mutation == "yes") %>% 
  select(POS, REF, ALT) %>%
  unique() %>%
  mutate(hg38_POS = hg38_coord[POS], 
         mut_minigene_code = paste0(ALT,POS,REF))

sig_mut


collapsed_var <- fread("data/SNPs/Annotated_variants/Collapsed_variants_10092021.tab")


prevalence_variants <- left_join(prevalence_variants, collapsed_var, by = c("hg38_POS" = "POS", "REF" = "REF", "ALT" = "ALT"))
sig_mut <- left_join(sig_mut, collapsed_var, by = c("hg38_POS" = "POS", "REF" = "REF", "ALT" = "ALT"))


prevalence_variants_annotated <- prevalence_variants %>% 
  filter(set != "") %>%  
  select(hg38_POS,REF, ALT, mut_minigene_code, set, Variation_ID) %>% 
  mutate(hg38_POS = paste0("chr16:", hg38_POS))

sig_mut_annotated <- sig_mut %>% 
  filter(set != "") %>%  
  select(hg38_POS, REF, ALT, mut_minigene_code, set, Variation_ID) %>% 
  mutate(hg38_POS = paste0("chr16:", hg38_POS))


annotated_variants <- rbindlist(list("high_prevalence_score" = prevalence_variants_annotated, 
          "effective_mutation" = sig_mut_annotated), idcol = "dataset")

library(writexl)
write_xlsx(annotated_variants, "output/Table_Sxx_01_AnnotatedVariants.xlsx")

