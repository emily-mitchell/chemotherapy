
SV_sum <- read.table("/Users/mp29/Documents/Main PhD/Blood/Metadata/SV_summary.txt", header=T,check.names =F, sep="\t",quote = "", row.names=1)

SV_sum <- tibble::rownames_to_column(SV_sum, "PDid")

CN_sum <- read.table("/Users/mp29/Documents/Main PhD/Blood/Metadata/CN_summary.txt", header=T,check.names =F, sep="\t",quote = "", row.names=1)

CN_sum <- tibble::rownames_to_column(CN_sum, "PDid")


### Structural Variants Plots

PD_order_chemo <- c("PD50308 (27F)", "PD50307 (40F)", "PD37580 (43F)", "PD47703 (48F)","PD47537 (61F)", "PD44579 (63F)")
p_SV_chemo <- SV_sum[SV_sum$Status == "Chemotherapy",] %>% 
  pivot_longer(cols = matches("^F"), names_to = "group", values_to = "proportion") %>%
  mutate(PDid = factor(PDid, level= PD_order_chemo)) %>%
  ggplot(aes(x= PDid, y = proportion, fill = group)) + ylim (0,0.08) +
  geom_bar(position= "stack",stat = "identity", width = 0.9, color = "black") +
  #geom_text(aes(label = SV_number, hjust = ifelse(SV_number < 0, 1.5, -1),vjust = 0.5), size = 3) + 
  ylab("Indipendently acquired SV per cell") + 
  scale_fill_manual(values = c("#7a4988","#c3b1e1","#bf40bf","#66023c"), labels = c("Deletion","Duplication","Inversion", "Translotation")) +
  theme_pubr() + 
  theme(axis.text.y  = element_text(size = 10), axis.text.x = element_text(angle = 60, hjust = 1, size = 10), axis.title.x = element_blank(), legend.title = element_blank())

p_SV_chemo


PD_order_normal <- c("PD40521 (29M)","PD41048 (48M)","PD49237 (60F)","PD49236 (63M)","PD45534 (77F)")
p_SV_normal <- SV_sum[SV_sum$Status == "Normal",] %>% 
  pivot_longer(cols = matches("^F"), names_to = "group", values_to = "proportion") %>%
  mutate(PDid = factor(PDid, level= PD_order_normal)) %>%
  ggplot(aes(x= PDid, y = proportion, fill = group)) + ylim (0, 0.08) + 
  geom_bar(position= "stack",stat = "identity", width = 0.9, color = "black") +
  #geom_text(aes(label = SV_number, hjust = ifelse(SV_number < 0, 1.5, -1),vjust = 0.5), size = 3) + 
  ylab("Indipendently acquired SV per cell") + 
  scale_fill_manual(values = c("#7a4988","#c3b1e1","#bf40bf","#66023c"), labels = c("Deletion","Duplication","Inversion", "Translotation")) +
  theme_pubr() + 
  theme(axis.text.y  = element_text(size = 10), axis.text.x = element_text(angle = 60, hjust = 1, size = 10), axis.title.x = element_blank(), legend.title = element_blank())

p_SV_normal

pSV_all <- ggarrange(p_SV_chemo, p_SV_normal)


### Copy Number Plots

PD_order_chemo <- c("PD50308 (27F)", "PD50307 (40F)", "PD37580 (43F)", "PD47703 (48F)","PD47537 (61F)", "PD44579 (63F)")
p_CN_chemo <- CN_sum[CN_sum$Status == "Chemotherapy",] %>% 
  mutate(PDid = factor(PDid, level= PD_order_chemo)) %>%
  ggplot(aes(x= PDid, y = Copy_neutral_LOH_fraction)) +
  geom_bar(position= "stack",stat = "identity", width = 0.9, color = "black", fill = "#165caa") +
  #geom_text(aes(label = SV_number, hjust = ifelse(SV_number < 0, 1.5, -1),vjust = 0.5), size = 3) + 
  ylab("Indipendently acquired autosomal CNA per cell") + ylim(0, 0.025) + 
  theme_pubr() + 
  theme(axis.text.y  = element_text(size = 10), axis.text.x = element_text(angle = 60, hjust = 1, size = 10), axis.title.x = element_blank(), legend.title = element_blank())

p_CN_chemo


PD_order_normal <- c("PD40521 (29M)","PD41048 (48M)","PD49237 (60F)","PD49236 (63M)","PD45534 (77F)")
p_CN_normal <- CN_sum[CN_sum$Status == "Normal",] %>% 
  mutate(PDid = factor(PDid, level= PD_order_normal)) %>%
  ggplot(aes(x= PDid, y = Copy_neutral_LOH_fraction)) +
  geom_bar(position= "stack",stat = "identity", width = 0.9, color = "black", fill = "#165caa") +
  #geom_text(aes(label = SV_number, hjust = ifelse(SV_number < 0, 1.5, -1),vjust = 0.5), size = 3) + 
  ylab("Indipendently acquired autosomal CNA per cell") + ylim(0, 0.025) + 
  theme_pubr() + 
  theme(axis.text.y  = element_text(size = 10), axis.text.x = element_text(angle = 60, hjust = 1, size = 10), axis.title.x = element_blank(), legend.title = element_blank())

p_CN_normal

pCN_all <- ggarrange(p_CN_chemo, p_CN_normal)


