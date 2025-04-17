
indel_hdp <- read.table("inferred_exposures.txt", header=T,check.names =F, sep="\t",quote = "", row.names=1)

pdOrder = c("PD47538","PD47537","PD47702","PD47536","PD47698","PD47701","PD47697","PD47700","PD50307","PD50306","PD47540","PD44579","PD47539","PD47695","PD47696","PD47541","PD60010","PD60009","PD60011","PD47699","PD50308","PD47703","PD37580","PD40315","PD40521","PD41048","PD49237","PD43976","PD49236","PD47738","PD48402","PD45534")

idOrder = c('ID1/2','ID3/5/9',"IDA","ID_unassigned")


pIndel_sigs_raw <- indel_hdp_burden %>% 
  pivot_longer(cols = matches("^ID"), names_to = "sig", values_to = "count") %>% 
  mutate(PDid = factor(PDid, level= pdOrder)) %>% 
  mutate(sig = factor(sig, level= idOrder)) %>% 
  ggplot(aes(x= PDid, y= count, fill = sig)) +
  geom_bar(position= "stack",stat = "identity") +
  ylab("Number of indels")  +
  theme_pubr() +
  scale_fill_manual(values = c('wheat2','cadetblue3','orange','ivory4')) +
  theme(axis.text.x = element_text(angle = 90), legend.position = "right")

pIndel_sigs_pro <- indel_hdp_burden %>% 
  pivot_longer(cols = matches("^ID"), names_to = "sig", values_to = "percent") %>% 
  mutate(PDid = factor(PDid, level= pdOrder)) %>% 
  mutate(sig = factor(sig, level= idOrder)) %>% 
  ggplot(aes(x= PDid, y= percent, fill = sig)) +
  geom_bar(position= "fill",stat = "identity") +
  ylab("Proportion of indel signatures")  +
  theme_pubr() +
  scale_fill_manual(values =  c('wheat2','cadetblue3','orange','ivory4')) +
  theme(axis.text.x = element_text(angle = 90), legend.position = "right")


pIndel_sigs <- ggarrange(pIndel_sigs_raw, pIndel_sigs_pro, nrow = 2)


