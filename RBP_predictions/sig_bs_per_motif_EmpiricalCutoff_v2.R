# Quantification of significant sites in binding motifs
library(data.table)

filtered_mut <- readRDS("../Mutation_effects/Softmax_modelling/distribution_cutoffs.rds")

library(dplyr)
library(tidyr)
library(GenomicRanges)

# Get model data from the minigene --------------------------------



# Filter only single point mutations
filtered_mut <- filtered_mut %>%
  filter(nchar(ALT) == 1 & nchar(REF) == 1)

# Get the significant ones (cutoff empirical > 0)
sig_pos_scores <- filtered_mut %>%
  filter(reps_cutoff_empirical>0) %>%
  select(isoform, ALT, POS) %>%
  unique()

isofact<-levels(filtered_mut$isoform)

sig_pos_scores$isoform <- factor(sig_pos_scores$isoform,
                                 levels = isofact,
                                 labels = isofact,
                                 ordered = T
)


sig_pos.gr <- GRanges(
  seqnames = "minigene",
  ranges = IRanges(
    start = sig_pos_scores$POS,
    end = sig_pos_scores$POS
  ),
  strand = "+",
  ALT = sig_pos_scores$ALT,
  isoforms = sig_pos_scores$isoform
)


# chr16:28943192-28944465:+ (minigene)
cd19_nts_hg19 <- seq(28943192, 28944465)

### AtTRACT binding sites -------------------------------------------------------------------------------------------------------------------------

attract_motifs <- fread("data/CD19_merged_bs_wSeq.bed")
colnames(attract_motifs) <- c("chr", "hg38_start", "hg38_end", "RBP_motif", "merged_motifs", "strand")



attract_motifs <- attract_motifs %>%
  separate(RBP_motif, c("RBP", "motif")) %>%
  mutate(
    minigene_st = match(hg38_start, cd19_nts_hg19),
    minigene_end = match(hg38_end, cd19_nts_hg19)
  )




attract.gr <- GRanges(
  seqnames = "minigene",
  ranges = IRanges(start = attract_motifs$minigene_st, end = attract_motifs$minigene_end),
  strand = "+",
  RBP = attract_motifs$RBP
)

ovp_bs <- findOverlapPairs(attract.gr, sig_pos.gr)

ovp_df <- cbind(
  (first(ovp_bs) %>%
     as.data.frame()),
  second(ovp_bs) %>%
    as.data.frame()
)

colnames(ovp_df) <- paste0(colnames(ovp_df), c(rep(c("_attract", ""), each = 6), ""))

ovp_df <- ovp_df %>%
  select(start_attract, end_attract, width_attract, RBP_attract, start, ALT, isoforms)

not_ov_df <- as.data.frame(attract.gr[!(attract.gr %in% first(ovp_bs))] )
colnames(not_ov_df) <- paste0(colnames(not_ov_df), "_attract")
not_ov_df <- not_ov_df %>%
  select(start_attract, end_attract, width_attract, RBP_attract)

ovp_df <- rbindlist(list(ovp_df, not_ov_df), fill = T)
ovp_df <- ovp_df %>% 
  mutate(isoforms = forcats::fct_explicit_na(isoforms, na_level = "No significant"))

# Only considering the positions, not the mutation (A site can be significant in up to 3 different mutations)

attract_sig_per_siteiso <- ovp_df %>%
  select(-ALT) %>%
  unique() %>%
  dplyr::count(start_attract, end_attract, width_attract, RBP_attract, isoforms, name = "sig_sites") %>%
  mutate(sig_sites = ifelse(isoforms == "No significant", 0, sig_sites))


attract_sig_per_site <- ovp_df %>%
  mutate(sig_bs = ifelse(isoforms != "No significant", "Sig", "NoSig")) %>% 
  select(-ALT, -isoforms, -start) %>%
  unique()


# Counts how many are significant per isoform
attract_sig_per_RBPiso <- ovp_df %>%
  select(-ALT) %>%
  unique() %>%
  dplyr::count(RBP_attract, isoforms, name = "sig_bs") 


attract_total_bs <- ovp_df %>% 
  select(start_attract, end_attract, width_attract, RBP_attract) %>% 
  unique() %>% 
  count(RBP_attract, name = "total_bs")


attract_sig_per_RBPiso <- left_join(attract_sig_per_RBPiso, attract_total_bs, by = "RBP_attract")

# Significant and not significant per RBP 
attract_sig_per_RBP <- ovp_df %>%
  mutate(sig_bs = ifelse(isoforms != "No significant", "Total_Any", "Total_NoSig")) %>% 
  select(-ALT, -isoforms, -start) %>%
  unique() %>%
  dplyr::count(RBP_attract, sig_bs, name = "total_sig")  %>% 
  pivot_wider(names_from = sig_bs, values_from = total_sig, values_fill = 0)


colnames(attract_sig_per_RBP) <- gsub("_attract", "", colnames(attract_sig_per_RBP))
colnames(attract_sig_per_RBPiso) <- gsub("_attract", "", colnames(attract_sig_per_RBPiso))
colnames(attract_sig_per_site) <- gsub("_attract", "", colnames(attract_sig_per_site))
colnames(attract_sig_per_siteiso) <- gsub("_attract", "", colnames(attract_sig_per_siteiso))

fwrite(attract_sig_per_site, "output/tab/motif_sig_pos/Empirical_Cutoff/AtTRACT_sig_pos_per_site.tab")
fwrite(attract_sig_per_siteiso, "output/tab/motif_sig_pos/Empirical_Cutoff/AtTRACT_sig_pos_per_site_iso.tab")
fwrite(attract_sig_per_RBPiso, "output/tab/motif_sig_pos/Empirical_Cutoff/AtTRACT_sig_pos_per_RBP.tab")
fwrite(attract_sig_per_RBP, "output/tab/motif_sig_pos/Empirical_Cutoff/AtTRACT_sig_pos_per_RBPiso.tab")


## oRNAment binding sites ---------------------------------------------------------------------------------------------


ornament_motifs <- fread("data/motifs_exon1-3.cleaned.tsv")
ornament_motifs <- ornament_motifs %>%
  mutate(RBP = gsub("\\s\\(.*\\)", "", motif))

ornament.gr <- GRanges(
  seqnames = "minigene",
  ranges = IRanges(start = ornament_motifs$minigene_pos_st, end = ornament_motifs$minigene_pos_end),
  strand = "+",
  RBP = ornament_motifs$RBP
)

ovp_bs <- findOverlapPairs(ornament.gr, sig_pos.gr)

ovp_df <- cbind(
  (first(ovp_bs) %>%
     as.data.frame()),
  second(ovp_bs) %>%
    as.data.frame()
)

colnames(ovp_df) <- paste0(colnames(ovp_df), c(rep(c("_ornament", ""), each = 6), ""))


ovp_df <- ovp_df %>%
  select(start_ornament, end_ornament, width_ornament, RBP_ornament, start, ALT, isoforms)

not_ov_df <- as.data.frame(ornament.gr[!(ornament.gr %in% first(ovp_bs))] )
colnames(not_ov_df) <- paste0(colnames(not_ov_df), "_ornament")
not_ov_df <- not_ov_df %>%
  select(start_ornament, end_ornament, width_ornament, RBP_ornament)

ovp_df <- rbindlist(list(ovp_df, not_ov_df), fill = T)
ovp_df <- ovp_df %>% 
  mutate(isoforms = forcats::fct_explicit_na(isoforms, na_level = "No significant"))

# Only considering the positions, not the mutation (A site can be significant in up to 3 different mutations)

ornament_sig_per_siteiso <- ovp_df %>%
  select(-ALT) %>%
  unique() %>%
  dplyr::count(start_ornament, end_ornament, width_ornament, RBP_ornament, isoforms, name = "sig_sites") %>%
  mutate(sig_sites = ifelse(isoforms == "No significant", 0, sig_sites))


ornament_sig_per_site <- ovp_df %>%
  mutate(sig_bs = ifelse(isoforms != "No significant", "Sig", "NoSig")) %>% 
  select(-ALT, -isoforms, -start) %>%
  unique()


# Counts how many are significant per isoform
ornament_sig_per_RBPiso <- ovp_df %>%
  select(-ALT) %>%
  unique() %>%
  dplyr::count(RBP_ornament, isoforms, name = "sig_bs") 


ornament_total_bs <- ovp_df %>% 
  select(start_ornament, end_ornament, width_ornament, RBP_ornament) %>% 
  unique() %>% 
  count(RBP_ornament, name = "total_bs")


ornament_sig_per_RBPiso <- left_join(ornament_sig_per_RBPiso, ornament_total_bs, by = "RBP_ornament")

# Significant and not significant per RBP 
ornament_sig_per_RBP <- ovp_df %>%
  mutate(sig_bs = ifelse(isoforms != "No significant", "Total_Any", "Total_NoSig")) %>% 
  select(-ALT, -isoforms, -start) %>%
  unique() %>%
  dplyr::count(RBP_ornament, sig_bs, name = "total_sig")  %>% 
  pivot_wider(names_from = sig_bs, values_from = total_sig, values_fill = 0)


colnames(ornament_sig_per_RBP) <- gsub("_ornament", "", colnames(ornament_sig_per_RBP))
colnames(ornament_sig_per_RBPiso) <- gsub("_ornament", "", colnames(ornament_sig_per_RBPiso))
colnames(ornament_sig_per_site) <- gsub("_ornament", "", colnames(ornament_sig_per_site))
colnames(ornament_sig_per_siteiso) <- gsub("_ornament", "", colnames(ornament_sig_per_siteiso))

fwrite(ornament_sig_per_site, "output/tab/motif_sig_pos/Empirical_Cutoff/oRNAment_sig_pos_per_site.tab")
fwrite(ornament_sig_per_siteiso, "output/tab/motif_sig_pos/Empirical_Cutoff/oRNAment_sig_pos_per_site_iso.tab")
fwrite(ornament_sig_per_RBPiso, "output/tab/motif_sig_pos/Empirical_Cutoff/oRNAment_sig_pos_per_RBP.tab")
fwrite(ornament_sig_per_RBP, "output/tab/motif_sig_pos/Empirical_Cutoff/oRNAment_sig_pos_per_RBPiso.tab")



####  Significant sites in DeepRiPe predictions---------------------------------------------------------############

deepripe <- fread("data/rbp_variant_scores_DeepRiPe.csv") %>%
  mutate(POS = gsub("[A-Z]([0-9]+)[A-Z]", "\\1", variant)) %>%
  separate(ref_alt, c("REF", "ALT"))

deepripe$POS <- as.numeric(deepripe$POS)


deepripe_sig_mut <- left_join(deepripe,
                         sig_pos_scores %>%
                           mutate(sig_mut = 1),
                         by = c("POS", "ALT")
) 

deepripe_bs <- deepripe_bs <- left_join(deepripe,
                                        sig_pos_scores %>%
                                          mutate(sig_mut = 1),
                                        by = c("POS", "ALT")
                                        ) %>% 
  group_by(RBP, set) %>%
  summarise(
    total_pos = n(),
    total_sig_pos = sum(sig_mut, na.rm = T)
  ) %>%
  ungroup()


# Vector with unique significant variants in the dataset
deepripe_sig <- left_join(deepripe,
                          sig_pos_scores %>%
                            mutate(sig_mut = 1),
                          by = c("POS", "ALT")
) %>%
  filter(sig_mut == 1) %>%
  select(isoform, variant) %>%
  unique()

deepripe_sig.list <- split(deepripe_sig$variant, deepripe_sig$isoform)



# GRanges with all the positions
deepripe.grl <- deepripe %>%
  mutate(strand = ifelse(score > 0, "+", "-")) %>%
  select(POS, variant, set, strand, RBP) %>%
  unique() %>%
  mutate(
    seqnames = "minigene",
    start = POS,
    end = POS,
    rbp_set = paste0(RBP, "_", set)
  ) %>%
  makeGRangesListFromDataFrame(.,
                               keep.extra.columns = T,
                               split.field = "rbp_set",
                               names.field = "rbp_set"
  )

deepripe.grl

# Granges reduced (to get clusters)
deepripe_red.grl <- reduce(deepripe.grl, drop.empty.ranges = F, with.revmap = T)


# Getting the information of the records that where collapsed


deepripe_red.grl_var <- mapply(FUN = function(set, red) {
  variants <- sapply(red$revmap, function(ind) {
    paste(set[ind]$variant, collapse = ",")
  })
  
  red$variants <- variants
  red
}, red = deepripe_red.grl, set = deepripe.grl, USE.NAMES = T)


deepripe_red.grl_var <- GenomicRangesList(deepripe_red.grl_var) %>%
  unlist()




# Filtering motifs for any cluster with at least 4 consecutive positions
deepripe_motifs <- deepripe_red.grl_var[width(deepripe_red.grl_var) >= 4]

# Length of variants in a motif
deepripe_motifs$n_variants <- sapply(
  deepripe_motifs$variants,
  function(x) {
    stringr::str_dplyr::count(x, ",")
  }
)

deepripe_motifs$RBP_set <- names(deepripe_motifs)


# Convert the motif GR to a data frame

deepripe_motifs.df <- as_tibble(deepripe_motifs) %>%
  select(-seqnames, -revmap) %>%
  separate(RBP_set, into = c("RBP", "set"), sep = "_", extra = "merge")


# Quantification of total motifs and total significant motifs


sig_var_iso <- lapply(deepripe_sig.list, function(iso) {
  
  # Getting the significant variants that overlap the detected with DeepRiPe
  sig_list <- lapply(
    strsplit(deepripe_motifs.df$variants, ","),
    function(v) {
      int_var <- intersect(iso, v)
      sig_var <- paste(int_var, collapse = ",")
      n_sig <- length(int_var)
      n_var <- length(v)
      return(list(sig_var, n_sig, n_var))
    }
  )
  
  # Adding them as new columns
  deepripe_motifs.df$n_sig_variants <- sapply(sig_list, `[[`, 2)
  deepripe_motifs.df$sig_variants <- sapply(sig_list, `[[`, 1)
  deepripe_motifs.df$n_variants <- sapply(sig_list, `[[`, 3)
  
  deepripe_motifs.df
})

# Convert the list to dataframe
deepripe_motifs_iso.df <- rbindlist(sig_var_iso, use.names = T, idcol = "isoform")


# dplyr::count significant motifs per binding site and isoform

deepripe_motifs_tab <- left_join(deepripe_motifs_iso.df %>%
                                   dplyr::count(RBP, isoform, set, name = "n_motifs"),
                                 deepripe_motifs_iso.df %>%
                                   filter(n_sig_variants > 0) %>%
                                   dplyr::count(RBP, isoform, set, name = "n_motifs_with_sig_var"),
                                 by = c("RBP", "set", "isoform")
)
# Adding variant ratio
deepripe_motifs_iso.df <- deepripe_motifs_iso.df %>%
  mutate(
    total_pos_var = width * 3,
    raw_var_ratio = n_sig_variants / n_variants,
    norm_var_ratio = n_sig_variants / total_pos_var
  )


var_ratio <- deepripe_motifs_iso.df %>%
  group_by(RBP, set) %>%
  summarise(
    max_raw_var_ratio = max(raw_var_ratio),
    min_raw_var_ratio = min(raw_var_ratio),
    max_norm_var_ratio = max(norm_var_ratio),
    min_norm_var_ratio = min(norm_var_ratio)
  ) %>%
  ungroup()


# Merge with the individual variants
# deepripe_bs<-left_join(
# deepripe_bs,
# left_join(var_ratio, deepripe_motifs_tab, by=c("RBP", "set", "isoform")),
# by=c("RBP", "set", "isoform"))

fwrite(deepripe_sig_mut, "output/tab/motif_sig_pos/Empirical_Cutoff/DeepRiPe_predictions_sig_mut.tab")
fwrite(deepripe_motifs_iso.df, "output/tab/motif_sig_pos/Empirical_Cutoff/DeepRiPe_predictions_sig_var_and_isoforms.tab")
fwrite(deepripe_bs, "output/tab/motif_sig_pos/Empirical_Cutoff/DeepRiPe_predictions_sig_var.tab")
fwrite(deepripe_motifs.df, "output/tab/motif_sig_pos/Empirical_Cutoff/DeepRiPe_predictions_sig_var_perbs.tab")
