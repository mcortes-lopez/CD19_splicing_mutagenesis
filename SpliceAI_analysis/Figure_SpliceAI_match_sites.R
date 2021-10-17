# Paper figure 4B: Splice AI sites 
library(data.table)
library(dplyr)

spliceai_vars <- fread("output/SpliceAI_unfiltered_variants.tab")
barcode_ss <- fread("output/max_ss_usage_per_barcode.tab")


single_vars <- spliceai_vars %>% 
  filter(set == "minigene_mask") %>% 
  filter(gain_score > 0 ) 




experimental_ss <- barcode_ss %>% 
  select(splice_sites, ss_type) %>% 
  unique() %>% 
  mutate(found_in_minigene = 1)


single_vars <- left_join(single_vars, experimental_ss, by = c("POS_gain_minigene" = "splice_sites", "ss_type" = "ss_type"))
single_vars$found_in_minigene <- ifelse(is.na(single_vars$found_in_minigene), F, T)

single_vars_unique <- single_vars %>% 
  select(POS_gain_minigene, gain_score, found_in_minigene, ss_type) %>% 
  group_by(POS_gain_minigene) %>% 
  arrange(desc(gain_score)) %>% 
  slice(1) %>% 
  ungroup()

cd19_exons<-data.frame(start=c(69,476,1041),
                       end=c(218,742,1244),
                       start_hg38=c(28931939,28932346,28932911),
                       end_hg38=c(28932088,28932612,28933114))


library(ggplot2)
library(ggthemes)
library(extrafont)
font_import(prompt = F)


pdf("out_plots/Fig_SpliceAI_splice_site_overlap.pdf", useDingbats = F, width = 10, height = 6, family = "ArialMT")
ggplot()+
  geom_rect(data=cd19_exons, aes(xmin=start, xmax=end, ymin=0, ymax=1),  alpha=0.3)+
  geom_bar(data=single_vars_unique, aes(x=POS_gain_minigene, y=gain_score, fill=found_in_minigene), stat="identity")+
  xlab("Splice Sites")+
  ylab("Max. SpliceAI score")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))+
  scale_x_continuous(limits = c(1, 1244)) +
  theme_classic()+
  scale_y_continuous(expand = c(0,0))+
  geom_hline(yintercept = 0.20, linetype = 2, colour = "grey50" )+
  labs(fill = "Overlapping \nsplice site")+ 
  scale_fill_manual(values = colorblind_pal()(8)[c(7,4)])+ 
  facet_wrap(~ss_type,ncol = 1, strip.position = "right")+
  theme(strip.background =  element_blank())
  
dev.off()

