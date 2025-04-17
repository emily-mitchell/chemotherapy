### Signature summary - all patients

### Chemo patients

sig_hdp_all <- read.table("burden_all.txt", header=T,check.names =F, sep="\t",quote = "", row.names=1)


pdOrder = c("PD47538","PD47537","PD47702","PD47536","PD47698","PD47701","PD47697","PD47700","PD50307","PD50306","PD47540","PD44579","PD47539","PD47695","PD47696","PD60009","PD47541","PD60011","PD60010","PD47699","PD50308","PD47703-1", "PD47703-2","PD37580")

sigOrder <- c("SBSUnassigned","SBSA","SBSB","SBSC","SBSD","SBSE","SBSF","SBSG","SBSH","SBS7a","SBS9","SBSBlood","SBS1+SBS5","SBSNA")


# Chemo samples
pSigs_chemo_raw <- sig_hdp_all[sig_hdp_all$Exposure =="Chemo" & sig_hdp_all$PDid !="PD47703-2",] %>% 
  pivot_longer(cols = matches("^SBS"), names_to = "sig", values_to = "mutCount") %>% 
  mutate(PDid = factor(PDid, level= pdOrder)) %>% 
  mutate(sig = factor(sig, level = sigOrder)) %>% 
  mutate(cell_type = factor(cell_type, level= c("HSPCs","Monocytes","B cells","Naive T cells","Memory T cells"))) %>% 
  ggplot(aes(x= PDid, y= mutCount, fill = sig)) +
  geom_bar(position= "stack",stat = "identity", width = 1) +
  ylab("Number of mutations")  +
  theme_pubr() +
  facet_grid(cell_type ~ PDid, switch = "y",
             scales = "free") + 
  scale_y_continuous(limits = c(0, 15000), breaks = seq(0, 15000, by = 2000)) +
  scale_fill_manual(values = c("ivory4","#ffac3c","darkorange4","springgreen3","#0b6623","darkolivegreen2","mediumorchid3","pink","red","darkslateblue","maroon3","powderblue","wheat2","#ebebeb")) + 
  
  theme(
    legend.title = element_blank(),
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 20), 
    strip.background.x = element_blank(),
    strip.placement  = "outside",
    strip.text.x.top = element_text(angle = 90, size = 18),
    strip.text.y.right = element_text(angle = 90, size = 22),
    axis.text.y.left = element_text(size = 20),
    axis.text.x = element_blank(), 
    axis.title.y = element_blank(), axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(0,"cm"))


pSigs_chemo_raw


pSigs_chemo_proportion <- sig_hdp_all[sig_hdp_all$Exposure =="Chemo" & sig_hdp_all$PDid !="PD47703-2",] %>% 
  pivot_longer(cols = matches("^SBS"), names_to = "sig", values_to = "mutCount") %>% 
  mutate(PDid = factor(PDid, level= pdOrder)) %>% 
  mutate(sig = factor(sig, level = sigOrder)) %>% 
  mutate(cell_type = factor(cell_type, level= c("HSPCs","Monocytes","B cells","Naive T cells","Memory T cells"))) %>% 
  ggplot(aes(x= PDid, y = mutCount, fill = sig)) +
  geom_bar(position = "fill", stat = "identity", width = 1) +
  ylab("Proportion of mutations")  +
  theme_pubr() +
  facet_grid(cell_type ~ PDid, scales ="free") + 
  scale_fill_manual(values = c("ivory4","#ffac3c","darkorange4","springgreen3","#0b6623","darkolivegreen2","mediumorchid3","pink","red","darkslateblue","maroon3","powderblue","wheat2","#ebebeb")) + 
  theme(
    legend.title = element_blank(),
    legend.position = "none",
    legend.text = element_text(size = 15),
    strip.text = element_text(face = "bold", size = 20), 
    strip.background.x = element_blank(),
    strip.text.x.top = element_text(angle = 90, size = 18),
    strip.text.y.right = element_blank(),
    axis.text.x = element_blank(), 
    axis.text.y.right = element_text(size = 20),
    axis.title.y = element_blank(), axis.title.x = element_blank(),
    axis.ticks.x.bottom  = element_blank(),
    axis.line.x = element_blank(), 
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(0,"cm")) + 
  scale_y_continuous(position = "right") 


pSigs_chemo_proportion

pSigs <- ggarrange(pSigs_chemo_raw, pSigs_chemo_proportion)


### Normal patients

pdOrder_norm = c("PD40315","PD40521","PD41048","PD49237","PD43976","PD49236","PD47738","PD48402","PD45534")


sig_hdp_all$cell_type <- gsub("Naive CD4 cells", "Naive T cells", sig_hdp_all$cell_type)
sig_hdp_all$cell_type <- gsub("Memory CD4 cells", "Memory T cells", sig_hdp_all$cell_type)


pSigs_norm_raw <- sig_hdp_all[sig_hdp_all$Exposure =="Normal" & sig_hdp_all$PDid !="PD43974" & sig_hdp_all$cell_type %in% c("HSPCs","Monocytes","B cells","Naive T cells","Memory T cells"),] %>%
  pivot_longer(cols = matches("^SBS"), names_to = "sig", values_to = "mutCount") %>% 
  mutate(PDid = factor(PDid, level= pdOrder_norm)) %>% 
  mutate(sig = factor(sig, level = sigOrder)) %>% 
  mutate(cell_type = factor(cell_type, level= c("HSPCs","Monocytes","B cells","Naive T cells","Memory T cells"))) %>% 
  ggplot(aes(x= PDid, y= mutCount, fill = sig)) +
  geom_bar(position= "stack",stat = "identity", width = 1) +
  ylab("Number of mutations")  +
  theme_pubr() +
  facet_grid(cell_type ~ PDid, switch = "y",
             scales = "free") +   scale_y_continuous(limits = c(0, 15000), breaks = seq(0, 15000, by = 2000)) +
  scale_fill_manual(values = c("ivory4","#ffac3c","darkorange4","springgreen3","#0b6623","darkolivegreen2","mediumorchid3","pink","red","darkslateblue","maroon3","powderblue","wheat1","#ebebeb")) + 
  
  theme(
    legend.title = element_blank(),
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 15), 
    strip.background.x = element_blank(),
    strip.placement  = "outside",
    strip.text.x.top = element_text(angle = 90, size = 18),
    strip.text.y.right = element_text(angle = 90, size = 22),
    axis.text.y.left = element_text(size = 20),
    axis.text.x = element_blank(), 
    axis.title.y = element_blank(), axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(0,"cm"))



pSigs_norm_raw

pSigs_norm_proportion <- sig_hdp_all[sig_hdp_all$Exposure =="Normal" & sig_hdp_all$PDid !="PD43974" & sig_hdp_all$cell_type %in% c("HSPCs","Monocytes","B cells","Naive T cells","Memory T cells"),] %>% 
  pivot_longer(cols = matches("^SBS"), names_to = "sig", values_to = "mutCount") %>% 
  mutate(PDid = factor(PDid, level= pdOrder_norm)) %>%
  mutate(sig = factor(sig, level = sigOrder)) %>% 
  mutate(cell_type = factor(cell_type, level= c("HSPCs","Monocytes","B cells","Naive T cells","Memory T cells"))) %>%
  ggplot(aes(x= PDid, y = mutCount, fill = sig)) +
  geom_bar(position = "fill", stat = "identity", width = 1) +
  ylab("Proportion of mutations")  +
  theme_pubr() +
  facet_grid(cell_type ~ PDid, scales ="free") + 
  scale_fill_manual(values = c("ivory4","#ffac3c","darkorange4","springgreen3","#0b6623","darkolivegreen2","mediumorchid3","pink","red","darkslateblue","maroon3","powderblue","wheat2","#ebebeb")) + 
  theme(
    legend.title = element_blank(),
    legend.position = "none",
    legend.text = element_text(size = 15),
    strip.text = element_text(face = "bold", size = 15), 
    strip.background.x = element_blank(),
    strip.text.x.top = element_text(angle = 90, size = 18),
    strip.text.y.right = element_blank(),
    axis.text.x = element_blank(), 
    axis.text.y.right = element_text(size = 20),
    axis.title.y = element_blank(), axis.title.x = element_blank(),
    axis.ticks.x.bottom  = element_blank(),
    axis.line.x = element_blank(), 
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(0,"cm")) + 
  scale_y_continuous(position = "right") 

pSigs_norm_proportion

pSigs_norm <- ggarrange(pSigs_norm_raw, pSigs_norm_proportion)

