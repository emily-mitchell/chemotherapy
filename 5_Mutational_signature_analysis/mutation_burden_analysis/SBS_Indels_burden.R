
Summary_All <- read.delim("SBS_indel_burden_summary.txt")


### SBS analysis - high burden cases only

pburden_high <- ggplot(Summary_All)+
  theme_bw()+
  labs(x="Age", y="HSPC single base substitutions")+
  scale_x_continuous(limits = c(-5,90), breaks = c(0,20,40,60,80))+
  scale_y_continuous(limits = c(0,10000), breaks = c(1000,2000,3000,4000,5000,6000,7000,8000,9000,10000))+
  stat_smooth(data = (Summary_All[Summary_All$Exposure == "Normal",]), method = "lm",aes(x = Age, y =Number_mutations, colour = PDCODE), size = 0.5, fullrange= TRUE)+
  scale_colour_manual(values = c("black","firebrick","darkorchid4","gold2","deeppink"))+
  theme(text=element_text(size=18))+
  geom_boxplot(data = (Summary_All[Summary_All$Cancer_diagnosis %in% c("AML","Hodgkin lymphoma", "Neuroblastoma", "Normal"),]),aes(group = PDID, x = Age, y =Number_mutations), outlier.alpha = 0)+
  geom_point(data = (Summary_All[Summary_All$Cancer_diagnosis %in% c("AML","Hodgkin lymphoma", "Neuroblastoma", "Normal"),]), position = position_dodge(width = 0.75),  aes(group = PDID, x = Age, y =Number_mutations, colour = PDCODE))

pburden_high

### SBS analysis - Low burden cases only



my_colour = c("black","maroon","cadetblue2","orange","darkslateblue", "pink","azure4","skyblue3","tomato","hotpink","yellowgreen","thistle3","darkblue","wheat","darkgoldenrod","yellow2","chocolate4","darkgreen", "purple", "green")

pburden_low <- ggplot(Summary_All)+
  theme_bw()+
  labs(x="Age (years)", y="HSPC single base substitutions")+
  scale_x_continuous(limits = c(35,85), breaks = c(40,60,80))+
  scale_y_continuous(limits = c(0,2000), breaks = c(200,400,600,800,1000,1200,1400,1600,1800,2000))+
  stat_smooth(data = (Summary_All[Summary_All$Exposure == "Normal",]), method = "lm",aes(x = Age, y =Number_mutations, colour = PDCODE), size = 0.5, fullrange= TRUE)+
  scale_colour_manual(values = my_colour)+
  theme(text=element_text(size=18))+
  geom_boxplot(data = (Summary_All[Summary_All$Cancer_diagnosis %in% c("Follicular lymphoma","Marginal zone lymphoma", "Lymphoplasmacytic lymphoma", "DLBCL", "Myeloma and Colon cancer", "Colon cancer", "Lung cancer", "Normal"),]), aes(group = PDID, x = Age, y =Number_mutations), outlier.alpha = 0)+
  geom_point(data = (Summary_All[Summary_All$Cancer_diagnosis %in% c("Follicular lymphoma","Marginal zone lymphoma", "Lymphoplasmacytic lymphoma", "DLBCL", "Myeloma and Colon cancer", "Colon cancer", "Lung cancer", "Normal"),]), position = position_dodge(width = 0.75), aes(group = PDID, x = Age, y =Number_mutations, colour = PDCODE))

pburden_low


### Indel analysis - high burden cases only

ggplot(Summary_All)+
  theme_bw()+
  labs(x="Age", y="Indel burden")+
  scale_x_continuous(limits = c(-5,90), breaks = c(0,20,40,60,80))+
  scale_y_continuous(limits = c(0,200), breaks = c(0,20,40,60,80,100,120,140,160,180,200))+
  theme(text=element_text(size=18))+
  stat_smooth(data = (Summary_All[Summary_All$Exposure == "Normal",]), method = "lm",aes(x = Age, y =Number_indels), size = 0.5, fullrange= TRUE)+
  scale_colour_manual(values = c(brewer.pal(12,"Paired")[c(2,6)]))+
  geom_boxplot(data = (Summary_All[Summary_All$Cancer_diagnosis %in% c("AML","Hodgkin lymphoma", "Neuroblastoma", "Normal"),]),aes(group = PDID, x = Age, y =Number_indels), outlier.alpha = 0)+
  geom_point(data = (Summary_All[Summary_All$Cancer_diagnosis %in% c("AML","Hodgkin lymphoma", "Neuroblastoma", "Normal"),]), position = position_dodge(width = 0.75),  aes(group = PDID, x = Age, y =Number_indels, colour = Exposure))



### Indel analysis - low burden cases only


ggplot(Summary_All)+
  theme_bw()+
  labs(x="Age", y="Indel burden")+
  scale_x_continuous(limits = c(35,85), breaks = c(40,60,80))+
  scale_y_continuous(limits = c(0,120), breaks = c(0,20,40,60,80,100,120))+
  theme(text=element_text(size=18))+
  stat_smooth(data = (Summary_All[Summary_All$Exposure == "Normal",]), method = "lm",aes(x = Age, y =Number_indels), size = 0.5, fullrange= TRUE)+
  scale_colour_manual(values = c(brewer.pal(12,"Paired")[c(2,6)]))+
  geom_boxplot(data = (Summary_All[Summary_All$Cancer_diagnosis %in% c("Follicular lymphoma","Marginal zone lymphoma", "Lymphoplasmacytic lymphoma", "DLBCL", "Myeloma and Colon cancer", "Colon cancer", "Lung cancer", "Normal"),]), aes(group = PDID, x = Age, y =Number_indels), outlier.alpha = 0)+
  geom_point(data = (Summary_All[Summary_All$Cancer_diagnosis %in% c("Follicular lymphoma","Marginal zone lymphoma", "Lymphoplasmacytic lymphoma", "DLBCL", "Myeloma and Colon cancer", "Colon cancer", "Lung cancer", "Normal"),]), position = position_dodge(width = 0.75), aes(group = PDID, x = Age, y =Number_indels, colour = Exposure))

