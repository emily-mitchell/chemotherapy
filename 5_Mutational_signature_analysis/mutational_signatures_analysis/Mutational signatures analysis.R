
library(ape)
library(ggtree)
library(RColorBrewer)
library(hdp)
library(tidyverse)
options(stringsAsFactors = F)

mutations=read.table("trinuc_mut_mat.txt")
key_table=read.table("key_table.txt")

colnames(key_table) <- key_table[1,]
key_table <- key_table[-1,]

colnames(mutations) <- mutations[1,]
mutations <- mutations[-1,]
rownames(mutations) <- mutations[,1]
mutations <- mutations[,-1]

mutations <- mutations %>% 
  mutate_all(as.numeric)


exposures <- read.table("exposures.txt", header=T,check.names =F, sep="\t",quote = "")

rownames(exposures)=rownames(mutations)

#### create a new columns for unassigned signatures (sum of all signatures with contribution < 5%)

exposures$SBSunassigned <- NA

#### oop through each row (sample) in the dataframe
for (i in 1:nrow(exposures)) {
  # Select values less than 0.05 from columns 1 to 12 for the current sample
  values_less_than_0.05 <- exposures[i, c("SBSA", "SBSB", "SBSC", "SBSD", "SBSE", "SBSF", "SBSG", "SBSH", "SBSBlood","SBS1+SBS5", "SBS7a", "SBS9")] < 0.05
  
  # Sum the selected values and store the result in the last column for the current sample
  exposures[i, "SBSunassigned"] <- sum(exposures[i, c("SBSA", "SBSB", "SBSC", "SBSD", "SBSE", "SBSF", "SBSG","SBSH","SBSBlood","SBS1+SBS5", "SBS7a", "SBS9")][values_less_than_0.05])
}

#remove signatures with <5% distribution
exposures[,c(1:12)] <- as.data.frame(exposures[,c(1:12)]) %>%
  mutate(across(where(is.numeric), ~ ifelse(. < 0.05, 0, .)))

#select colonies samples only
exposures_nanoseq = exposures[2662:2740,]

exposures_colonies = exposures[1:2661,]

exposures_colonies <- t(exposures_colonies)

patients=read.table("patients.txt")[,1]

sigs=rownames(exposures_colonies)


#### adding mutational signatures onto the phylogenetic trees 
br_lengths <- NULL
for (patient in patients){
  if(file.exists(paste0("trees/",patient,".tree"))){
    tree = read.tree(paste0("trees/",patient,".tree"))
    tree_df=fortify(tree)
    
    cols= c("wheat2","powderblue","darkslateblue","maroon3","#ffac3c","darkorange4","springgreen3","#0b6623","darkolivegreen2","mediumorchid3","pink","red2", "ivory4")
    
    names(cols)=sigs
    samples=colnames(exposures_colonies)[grepl(patient,colnames(exposures_colonies))]
    branches=sub('.*_', '', samples)
    
    pdf(paste0("trees/",patient,"_tree_signatures_small.pdf"), height=4)
    plot(tree,show.tip.label=F, cex= 0.4)
    for (k in 1:length(samples)){
      n=as.numeric(branches[k])
      x_end=tree_df$x[n]
      x_start=tree_df$x[tree_df$parent[n]]
      x_intv=x_end-x_start
      y=node.height(tree)[n]
      l=tree_df$branch.length[n]
      t=samples[k]
      br_lengths <- rbind(br_lengths,cbind(t, l))
      tipnum=sum(tree_df$isTip)
      title(main=patient)
      for (s in sigs){
        x_end=x_start+exposures_colonies[s,samples[k]]*x_intv
        rect(ybottom=y-min(0.02*tipnum,0.4),ytop=y+min(0.02*tipnum,0.4),xleft=x_start,xright=x_end,col=cols[s], border = cols[s] , lwd = 0.5)
        x_start=x_end
        
      }
    }
    axisPhylo(side = 1,backward=F)
    dev.off()
    
  }
}

#### calculate the mutation burden of colonies samples

x <- as.data.frame(t(exposures_colonies)) %>% mutate_all(as.numeric)
brl <- as.data.frame(br_lengths)
brl$l <- as.numeric(brl$l)

x <- x[order(rownames(x)),]
brl <- brl[order(brl$t),]

x1 <- sweep(x, MARGIN = 1, brl[,2], "*")

key <- as.data.frame(read_xlsx("PDID_colonies.xlsx", col_names = T))

key <- key[-5,]

out <- NULL
for (patient in patients) {
  add <- colSums(x1[grep(patient, rownames(x1)),])
  add <- add/key[key$PDID %in% patient,2]
  out <- rbind(out, add)
}

rownames(out) <- patients

out <- as.data.frame(out) %>%
  mutate(across(where(is.numeric), ~ ifelse(. < 1, 0, .)))

### mutational signatures attribution

### Chemo patients


sig_hdp_all <- read.table("/Users/mp29/volumes/mp29_lustre/Emily_blood/signatures/HDP_w_contaminated/hdp_PD_method_woPD50306/run/burden_all_summary.txt", header=T,check.names =F, sep="\t",quote = "", row.names=1)


pdOrder = c("PD47538","PD47537","PD47702","PD47536","PD47698","PD47701","PD47697","PD47700","PD50307","PD50306","PD47540","PD44579","PD47539","PD47695","PD47696","PD60009","PD47541","PD60011","PD60010","PD47699","PD50308","PD47703-1", "PD47703-2","PD37580")

sigOrder <- c("SBSunassigned","SBSA","SBSB","SBSC","SBSD","SBSE","SBSF","SBSG","SBSH","SBS7a","SBS9","SBSBlood","SBS1+SBS5","SBSNA")


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
  scale_fill_manual(values = c("ivory4","#ffac3c","darkorange4","springgreen3","#0b6623","darkolivegreen2","mediumorchid3","pink","red2","darkslateblue","maroon3","powderblue","wheat2","#ebebeb")) + 
  
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
  scale_fill_manual(values = c("ivory4","#ffac3c","darkorange4","springgreen3","#0b6623","darkolivegreen2","mediumorchid3","pink","red2","darkslateblue","maroon3","powderblue","wheat2","#ebebeb")) + 
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
ggsave(pSigs,filename = "/Users/mp29/Documents/Main PhD/Blood/Plots/HDPSigs_chemo_clean.pdf",dpi=300,height=20, width=22)




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
  scale_fill_manual(values = c("ivory4","#ffac3c","darkorange4","springgreen3","#0b6623","darkolivegreen2","mediumorchid3","pink","red2","darkslateblue","maroon3","powderblue","wheat2","#ebebeb")) + 
  
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
  scale_fill_manual(values = c("ivory4","#ffac3c","darkorange4","springgreen3","#0b6623","darkolivegreen2","mediumorchid3","pink","red2","darkslateblue","maroon3","powderblue","wheat2","#ebebeb")) + 
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



