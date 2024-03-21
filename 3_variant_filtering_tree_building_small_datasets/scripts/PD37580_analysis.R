#This script is designed to be interactive & quick to allow rapid exploring of filtering parameters
#run interactively on farm
#bsub -Is -G team78-grp -q yesterday -R 'select[mem>=10000] rusage[mem=10000]' -M10000 R

#module load vagrent

#Amend scratch folder and team folder as appropriate
#scp -r /Users/em16/Documents/PhD/Sequencing_results/Mutographs/Alkylating_agents_2178/PD47539/subs/info/PD37580_samples_retained.txt em16@farm5-login:/lustre/scratch126/casm/team154pc/em16/Mutographs/ALK/PD47539/analysis/


#Specify options to parse from the command line
#args = commandArgs(TRUE)

##SET RUN_ID AND FILEPATHS
#Run_ID = as.character(args[1])
#output_directory = toString(args[2])

##Define objects
PDID = "PD37580" #enter PDID
number_samples = "10" #enter number of samples for this individual
gender = "female"  #gender of individual 
project = "ALK"
samples_to_keep = "yes"
mpboot_tree_file=paste0("/lustre/scratch126/casm/team154pc/em16/Mutographs/",project,"/",PDID,"/analysis/tree_files/Mutations_MPBoot_",PDID,".fa.treefile")
dna_string_file = paste0("/lustre/scratch126/casm/team154pc/em16/Mutographs/",project,"/",PDID,"/analysis/tree_files/Mutations_MPBoot_",PDID,".fa")
xlim_tree = c(0,12000)

system("mkdir mutation_burden")
system("mkdir dnds")
system("mkdir driver_info")
system("mkdir hdp")
system("mkdir mutation_spectrums")
system("mkdir tree_files")
system("mkdir tree_pdfs")
system("mkdir data")
system("mkdir vcf")

##Set working directory
my_working_directory = (paste0("/lustre/scratch126/casm/team154pc/em16/Mutographs/",project,"/",PDID,"/analysis/"))
setwd(my_working_directory)

##Open libraries
library(GenomicRanges)
library(dplyr)
library(tidyr)
library(stringr)
library(ape)
library(seqinr)
library(ggtree)
library(data.table)
library(ggplot2)
library(BiocGenerics)
library(Rsamtools)
library(MASS)
source("/lustre/scratch126/casm/team154pc/em16/R_scripts/functions/filters_parallel_functions.R")
source("/lustre/scratch126/casm/team154pc/em16/R_scripts/functions/tree_functions.R")

samples <- read.table(paste0(PDID,"_samples_retained.txt"), header = F, stringsAsFactors = F)
samp <- samples$V1

##Read in snv file and remove columns for samples not used in this analysis
mut <- read.csv(paste0("/lustre/scratch126/casm/team154pc/em16/Mutographs/",project,"/",PDID,"/filtered_output_3/output/output/PDv37is/snp/", PDID, "_filtered_3_merged_cut.tsv"), sep="\t", header= T, stringsAsFactors=F)
mut_cut <- mut[,colnames(mut)%in%paste0(samp, "_MTR")
                 |colnames(mut)%in%paste0(samp, "_DEP")]
muts <- as.data.frame(cbind(mut[,c("Chrom", "Pos", "Ref", "Alt")], mut_cut))
muts$mut_type <- "snv"

##Read in indel file and remove columns for samples not used in this analysis
indel <- read.csv(paste0("/lustre/scratch126/casm/team154pc/em16/Mutographs/",project,"/",PDID,"/indels/output/output/PDv37is/indel/", PDID, "_indel_filtered_merged_cut.tsv"), sep="\t", header= T, stringsAsFactors=F)
indel_cut <- indel[,colnames(indel)%in%paste0(samp, "_MTR")
                 |colnames(indel)%in%paste0(samp, "_DEP")]
indels <- as.data.frame(cbind(indel[,c("Chrom", "Pos", "Ref", "Alt")], indel_cut))
indels$mut_type <- "indel"

all <- rbind(muts,indels)

## Read in useful files
chip_drivers = read.csv("/lustre/scratch126/casm/team154pc/em16/R_scripts/chip_drivers.csv", stringsAsFactors = FALSE, header = FALSE)
chip_drivers = chip_drivers$V1
cancer_drivers = read.csv("/lustre/scratch126/casm/team154pc/em16/R_scripts/cancer_drivers.csv", stringsAsFactors = FALSE, header = FALSE)
cancer_drivers = cancer_drivers$V1
biobank_drivers = read.csv("/lustre/scratch126/casm/team154pc/em16/R_scripts/biobank_drivers.csv", stringsAsFactors = FALSE, header = FALSE)
biobank_drivers = biobank_drivers$V1


###############################################################################################################################################################
#Filtering to remove germline and  in vitro mutations
###############################################################################################################################################################

##Generate list of column names containing "MTR" and list of sample IDs
NV <- as.matrix(all[,grep("MTR", colnames(all), value=T)])
NR <- as.matrix(all[,grep("DEP", colnames(all), value=T)])
samps <- gsub("_MTR", "", colnames(NV))
NV_cols <- grep("MTR", colnames(all), value=T) #gives list of column names containing "MTR"
NR_cols <- grep("DEP", colnames(all), value=T) #gives list of column names containing "MTR"

##Remove rows that now have no positive samples
if(samples_to_keep == "yes") {
  min_variant_reads_auto = 3
  min_variant_reads_xy = 2
  null_remove = rowSums(NV >= min_variant_reads_auto |(NV >= min_variant_reads_xy & all$Chrom %in% c("X","Y") & gender == "male")) == 0
  cat(sum(null_remove),"mutations removed as no positives in any remaining samples.\n")
  all = all[!null_remove,]
  NV = NV[!null_remove,]
  NR = NR[!null_remove,]
}

##Count number of samples with mutant reads at each site
all$count <- apply(all[,NV_cols], 1, function(row) length(which(row>0))) #apply returns list of values obtained by applying function to 1(rows) - length of vector counting NV > 0

##Identify likely germline variants
all$germline <- "som"  #create column called germline - filled with 'som' as default
all$germline[all$count>(length(NV_cols)/2)] <- "germ"  # removes subs present in more than 1/2 samples from an individual
table(all$germline == "som") 
all_germline_removed<- all[all$germline=="som",1:(ncol(all)-1)]

##Call mutations on at least 2 reads and a vaf >=0.2
NV <- all_germline_removed[,grep("MTR", colnames(all_germline_removed), value=T)]
NR <- all_germline_removed[,grep("DEP", colnames(all_germline_removed), value=T)]
NR_nonzero=NR
NR_nonzero[NR_nonzero==0]=1 
VAF <- NV/NR_nonzero
colnames(VAF) <- gsub("_MTR", "", colnames(VAF))
filtered_all <- NV
header <- paste(all_germline_removed$Chrom, all_germline_removed$Pos, all_germline_removed$Ref, all_germline_removed$Alt,sep="_")
rownames(NV)=rownames(NR)=rownames(VAF)=rownames(filtered_all)=header

# mtr cutoff
filtered_all[NV>=2] <- 1
filtered_all[NV<2] <- 0

XY_chromosomal=grepl("X|Y",rownames(NR))
autosomal=!XY_chromosomal

# vaf cutoff
if(gender=="female"){
  filtered_all[VAF<=0.2] <- 0
  colnames(filtered_all) <- gsub("_MTR", "", colnames(filtered_all))
}#use if female

if(gender=="male"){
  filtered_all[VAF<=0.2&autosomal] <- 0
  filtered_all[VAF<=0.4&XY_chromosomal] <- 0
  colnames(filtered_all) <- gsub("_MTR", "", colnames(filtered_all))
}#use if male

#remove high and low depth sites
if(gender=="female"){
  filtered_all[rowMeans(NR)<8&rowMeans(NR)>50] <- 0 
}#use if female

if(gender=="male"){
  filtered_all[rowMeans(NR)<8&rowMeans(NR)>50&autosomal|rowMeans(NR)<4&rowMeans(NR)>25&XY_chromosomal] <- 0 
}#use if male

# bind back on the chrom pos ref alt and mut_type info
filtered_all <- as.data.frame(cbind(all_germline_removed[,c("Chrom", "Pos", "Ref", "Alt", "mut_type")], filtered_all))
filtered_all_clean <- filtered_all[(rowSums(filtered_all[6:ncol(filtered_all)])>0),]
NR <- as.data.frame(cbind(all_germline_removed[,c("Chrom", "Pos", "Ref", "Alt", "mut_type")], NR))
NR <- NR[(rowSums(filtered_all[6:ncol(filtered_all)])>0),]
NV <- as.data.frame(cbind(all_germline_removed[,c("Chrom", "Pos", "Ref", "Alt", "mut_type")], NV))
NV <- NV[(rowSums(filtered_all[6:ncol(filtered_all)])>0),]
VAF <- VAF[(rowSums(filtered_all[6:ncol(filtered_all)])>0),]


write.table(filtered_all_clean, paste0("/lustre/scratch126/casm/team154pc/em16/Mutographs/",project,"/",PDID,"/analysis/data/",PDID,"_mutations.txt"), sep="\t", col.names = T, row.names = F, quote = F)

##Sum columns for each sample and calculate mean depth for each sample
mutsums <- as.data.frame(colSums(filtered_all_clean[filtered_all_clean$mut_type == "snv",grep("PD", colnames(filtered_all_clean))]))
indelsums <- as.data.frame(colSums(filtered_all_clean[filtered_all_clean$mut_type == "indel",grep("PD", colnames(filtered_all_clean))]))
meandep <- as.data.frame((colSums(NR[grep("PD", colnames(NR))]))/(nrow(NR)))

#create column with sample names  
mutsums$Sample <- rownames(mutsums)
colnames(mutsums)[1] <- "Number_mutations"
colnames(indelsums)[1] <- "Number_indels"
colnames(meandep)[1] <- "Mean_depth"
summary <- as.data.frame(cbind(indelsums[,c("Number_indels")], mutsums))
summary <- as.data.frame(cbind(meandep[,c("Mean_depth")], summary))
sum <- as.data.frame(cbind(summary[,c(4,3,2,1)]))
colnames(sum)[3] <- "Number_indels"
colnames(sum)[4] <- "Mean_depth"

write.table(sum, paste0("/lustre/scratch126/casm/team154pc/em16/Mutographs/",project,"/",PDID,"/analysis/mutation_burden/",PDID,"_mut_dep.txt"), sep="\t", col.names = T, row.names = F, quote=F)

##Graph of snv mutation number vs depth

pdf(paste0("/lustre/scratch126/casm/team154pc/em16/Mutographs/",project,"/",PDID,"/analysis/mutation_burden/",PDID, "_mutation_depth.pdf"))
ggplot(sum, aes(x = Mean_depth, y = Number_mutations))+
  scale_x_continuous(limits = c(0,40))+
  scale_y_continuous(limits = c(0,12000), breaks = c(0,1000, 2000,3000,4000,5000,6000,7000,8000,9000,10000,11000,12000))+
  theme_bw() +
  labs(x="Mean depth", y="SNV number", title=paste0(PDID, " number of snvs per sample"))+
  geom_point()
dev.off()

##Graph of indel mutation number vs depth

pdf(paste0("/lustre/scratch126/casm/team154pc/em16/Mutographs/",project,"/",PDID,"/analysis/mutation_burden/",PDID, "_indel_depth.pdf"))
ggplot(sum, aes(x = Mean_depth, y = Number_indels))+
  scale_x_continuous(limits = c(0,40))+
  scale_y_continuous(limits = c(0,200), breaks = c(0,20,40,60,80,100,120,140,160,180,200))+
  theme_bw() +
  labs(x="Mean depth", y="Indel number", title=paste0(PDID, " number of indels per sample"))+
  geom_point()
dev.off()


###############################################################################################################################################################
#SNV signature analysis
###############################################################################################################################################################

##Create data.frame for signature analysis - samples pooled
genomeFile = "/nfs/cancer_ref01/Homo_sapiens/37/genome.fa"

filtered_muts_clean = filtered_all_clean[filtered_all_clean$mut_type == "snv",]
filtered_indels_clean = filtered_all_clean[filtered_all_clean$mut_type == "indel",]

samp = colnames(filtered_muts_clean[6:ncol(filtered_muts_clean)])

subs_only = filtered_muts_clean[1:4]
nrow(subs_only)
colnames(subs_only) = c("chr", "pos", "ref", "mut")
subs_only = subs_only[(subs_only$ref %in% c("A","C","G","T")) & (subs_only$mut %in% c("A","C","G","T")) & subs_only$chr %in% c(1:22,"X","Y"),]
subs_only$trinuc_ref = as.vector(scanFa(genomeFile, GRanges(subs_only$chr, IRanges(subs_only$pos-1, subs_only$pos+1))))

# 2. Annotating the mutation from the pyrimidine base
ntcomp = c(T="A",G="C",C="G",A="T")
subs_only$sub = paste(subs_only$ref,subs_only$mut,sep=">")
subs_only$trinuc_ref_py = subs_only$trinuc_ref
for (j in 1:nrow(subs_only)) {
  if (subs_only$ref[j] %in% c("A","G")) { # Purine base
    subs_only$sub[j] = paste(ntcomp[subs_only$ref[j]],ntcomp[subs_only$mut[j]],sep=">")
    subs_only$trinuc_ref_py[j] = paste(ntcomp[rev(strsplit(subs_only$trinuc_ref[j],split="")[[1]])],collapse="")
  }
}

# 3. Counting subs
freqs = table(paste(subs_only$sub,paste(substr(subs_only$trinuc_ref_py,1,1),substr(subs_only$trinuc_ref_py,3,3),sep="-"),sep=","))
sub_vec = c("C>A","C>G","C>T","T>A","T>C","T>G")
ctx_vec = paste(rep(c("A","C","G","T"),each=4),rep(c("A","C","G","T"),times=4),sep="-")
full_vec = paste(rep(sub_vec,each=16),rep(ctx_vec,times=6),sep=",")
freqs_full = freqs[full_vec]
freqs_full[is.na(freqs_full)] = 0
names(freqs_full) = full_vec

xstr = paste(substr(full_vec,5,5), substr(full_vec,1,1), substr(full_vec,7,7), sep="")

pdf(paste0("/lustre/scratch126/casm/team154pc/em16/Mutographs/",project,"/",PDID,"/analysis/mutation_spectrums/",PDID,"_mutation_spectrum.pdf"),width=10,height=4, useDingbats = F)
colvec = rep(c("dodgerblue","black","red","grey70","olivedrab3","plum2"),each=16)
y = freqs_full
maxy = max(y)

h=barplot(y, las=2, col=colvec, border=NA, ylim=c(0,maxy*1.5), space=1, cex.names=0.6, names.arg=xstr, ylab="# SNVs")
mtext(side=4, text= PDID)
for (j in 1:length(sub_vec)) {
  xpos = h[c((j-1)*16+1,j*16)]
  rect(xpos[1]-0.5, maxy*1.2, xpos[2]+0.5, maxy*1.3, border=NA, col=colvec[j*16])
  text(x=mean(xpos), y=maxy*1.3, pos=3, label=sub_vec[j])
} 
dev.off()

##Run mutation spectrum analysis - individual samples
samp = colnames(filtered_muts_clean[6:ncol(filtered_muts_clean)])

for (ts in samp) {
  
  subs_only = filtered_muts_clean[(filtered_muts_clean[,ts]==1),c(1:4)]
  nrow(subs_only)
  colnames(subs_only) = c("chr", "pos", "ref", "mut")
  subs_only = subs_only[(subs_only$ref %in% c("A","C","G","T")) & (subs_only$mut %in% c("A","C","G","T")) & subs_only$chr %in% c(1:22,"X","Y"),]
  subs_only$trinuc_ref = as.vector(scanFa(genomeFile, GRanges(subs_only$chr, IRanges(subs_only$pos-1, subs_only$pos+1))))
  
  # 2. Annotating the mutation from the pyrimidine base
  ntcomp = c(T="A",G="C",C="G",A="T")
  subs_only$sub = paste(subs_only$ref,subs_only$mut,sep=">")
  subs_only$trinuc_ref_py = subs_only$trinuc_ref
  for (j in 1:nrow(subs_only)) {
    if (subs_only$ref[j] %in% c("A","G")) { # Purine base
      subs_only$sub[j] = paste(ntcomp[subs_only$ref[j]],ntcomp[subs_only$mut[j]],sep=">")
      subs_only$trinuc_ref_py[j] = paste(ntcomp[rev(strsplit(subs_only$trinuc_ref[j],split="")[[1]])],collapse="")
    }
  }
  
  # 3. Counting subs
  freqs = table(paste(subs_only$sub,paste(substr(subs_only$trinuc_ref_py,1,1),substr(subs_only$trinuc_ref_py,3,3),sep="-"),sep=","))
  sub_vec = c("C>A","C>G","C>T","T>A","T>C","T>G")
  ctx_vec = paste(rep(c("A","C","G","T"),each=4),rep(c("A","C","G","T"),times=4),sep="-")
  full_vec = paste(rep(sub_vec,each=16),rep(ctx_vec,times=6),sep=",")
  freqs_full = freqs[full_vec]
  freqs_full[is.na(freqs_full)] = 0
  names(freqs_full) = full_vec
  
  xstr = paste(substr(full_vec,5,5), substr(full_vec,1,1), substr(full_vec,7,7), sep="")
  
  pdf(paste0("/lustre/scratch126/casm/team154pc/em16/Mutographs/",project,"/",PDID,"/analysis/mutation_spectrums/",ts,"_mutation_spectrum.pdf"),width=10,height=4, useDingbats = F)
  colvec = rep(c("dodgerblue","black","red","grey70","olivedrab3","plum2"),each=16)
  y = freqs_full
  maxy = max(y)
  
  h=barplot(y, las=2, col=colvec, border=NA, ylim=c(0,maxy*1.5), space=1, cex.names=0.6, names.arg=xstr, ylab="# SNVs")
  mtext(side=4, text= ts)
  for (j in 1:length(sub_vec)) {
    xpos = h[c((j-1)*16+1,j*16)]
    rect(xpos[1]-0.5, maxy*1.2, xpos[2]+0.5, maxy*1.3, border=NA, col=colvec[j*16])
    text(x=mean(xpos), y=maxy*1.3, pos=3, label=sub_vec[j])
  } 
  dev.off()
}


###############################################################################################################################################################
#Creation SNV vcf files
###############################################################################################################################################################

##Creation vcf files for pooled data
#Write vcf files for pooled data
vcf_file = filtered_muts_clean[,1:4]

names(vcf_file) = c("#CHROM", "POS", "REF", "ALT")

#Add the required null columns
vcf_file$ID = vcf_file$QUAL = vcf_file$FILTER = vcf_file$INFO = "."

#Rearrange into correct order
vcf_file = vcf_file[,c(1,2,5,3,4,6,7,8)]

#Write files
write.table(vcf_file, sep = "\t", quote = FALSE, file = paste0("/lustre/scratch126/casm/team154pc/em16/Mutographs/",project,"/",PDID,"/analysis/vcf/",PDID,"_snvs.vcf"), row.names = FALSE)

##Creation vcf files per sample
for (ts in samp) {
  #Write vcf files for VariantCaller analysis 
  vcf_file = filtered_muts_clean[(filtered_muts_clean[,ts] >0),c(1:4)]
  names(vcf_file) = c("#CHROM", "POS", "REF", "ALT")
  
  #Add the required null columns
  vcf_file$ID = vcf_file$QUAL = vcf_file$FILTER = vcf_file$INFO = "."
  
  #Rearrange into correct order
  vcf_file = vcf_file[,c(1,2,5,3,4,6,7,8)]
  
  #Write files
  write.table(vcf_file, sep = "\t", quote = FALSE, file = paste0("/lustre/scratch126/casm/team154pc/em16/Mutographs/",project,"/",PDID,"/analysis/vcf/",ts,"_snvs.vcf"), row.names = FALSE)
}


###############################################################################################################################################################
#Creation INDEL vcf files
###############################################################################################################################################################

##Creation vcf files for pooled data
#Write vcf files for pooled data
vcf_file = filtered_indels_clean[,1:4]

names(vcf_file) = c("#CHROM", "POS", "REF", "ALT")

#Add the required null columns
vcf_file$ID = vcf_file$QUAL = vcf_file$FILTER = vcf_file$INFO = "."

#Rearrange into correct order
vcf_file = vcf_file[,c(1,2,5,3,4,6,7,8)]

#Write files
write.table(vcf_file, sep = "\t", quote = FALSE, file = paste0("/lustre/scratch126/casm/team154pc/em16/Mutographs/",project,"/",PDID,"/analysis/vcf/",PDID,"_indels.vcf"), row.names = FALSE)

##Creation vcf files per sample
for (ts in samp) {
  #Write vcf files for VariantCaller analysis 
  vcf_file = filtered_indels_clean[(filtered_indels_clean[,ts] >0),c(1:4)]
  names(vcf_file) = c("#CHROM", "POS", "REF", "ALT")
  
  #Add the required null columns
  vcf_file$ID = vcf_file$QUAL = vcf_file$FILTER = vcf_file$INFO = "."
  
  #Rearrange into correct order
  vcf_file = vcf_file[,c(1,2,5,3,4,6,7,8)]
  
  #Write files
  write.table(vcf_file, sep = "\t", quote = FALSE, file = paste0("/lustre/scratch126/casm/team154pc/em16/Mutographs/",project,"/",PDID,"/analysis/vcf/",ts,"_indels.vcf"), row.names = FALSE)
}


###############################################################################################################################################################
#Creation ALL vcf files
###############################################################################################################################################################
##Creation vcf files for pooled data
#Write vcf files for pooled data
vcf_file = filtered_all_clean[,1:4]

names(vcf_file) = c("#CHROM", "POS", "REF", "ALT")

#Add the required null columns
vcf_file$ID = vcf_file$QUAL = vcf_file$FILTER = vcf_file$INFO = "."

#Rearrange into correct order
vcf_file = vcf_file[,c(1,2,5,3,4,6,7,8)]

#Write files
write.table(vcf_file, sep = "\t", quote = FALSE, file = paste0("/lustre/scratch126/casm/team154pc/em16/Mutographs/",project,"/",PDID,"/analysis/vcf/",PDID,"_mutations_all.vcf"), row.names = FALSE)

##Creation vcf files per sample
for (ts in samp) {
  #Write vcf files for VariantCaller analysis 
  vcf_file = filtered_muts_clean[(filtered_muts_clean[,ts] >0),c(1:4)]
  names(vcf_file) = c("#CHROM", "POS", "REF", "ALT")
  
  #Add the required null columns
  vcf_file$ID = vcf_file$QUAL = vcf_file$FILTER = vcf_file$INFO = "."
  
  #Rearrange into correct order
  vcf_file = vcf_file[,c(1,2,5,3,4,6,7,8)]
  
  #Write files
  write.table(vcf_file, sep = "\t", quote = FALSE, file = paste0("/lustre/scratch126/casm/team154pc/em16/Mutographs/",project,"/",PDID,"/analysis/vcf/",ts,"_mutations_all.vcf"), row.names = FALSE)
}


###############################################################################################################################################################
#Tree building
###############################################################################################################################################################

##Create genotype_bin for tree building (subs only)
genotype_bin=as.matrix(VAF[filtered_all_clean$mut_type == "snv",])

XY_chromosomal=grepl("X|Y",rownames(filtered_muts_clean))
autosomal=!XY_chromosomal

if(gender=="female"){
  genotype_bin[genotype_bin<0.15]=0
  genotype_bin[genotype_bin>=0.15&genotype_bin<0.3]=0.5
  genotype_bin[genotype_bin>=0.3]=1
}#use if female

if(gender=="male"){
  genotype_bin[genotype_bin<0.15&autosomal]=0
  genotype_bin[genotype_bin>=0.15&genotype_bin<0.3&autosomal]=0.5
  genotype_bin[genotype_bin>=0.3&autosomal]=1
  genotype_bin[genotype_bin<0.3&XY_chromosomal]=0
  genotype_bin[genotype_bin>=0.3&genotype_bin<0.6&XY_chromosomal]=0.5
  genotype_bin[genotype_bin>=0.6&XY_chromosomal]=1
}#use if male

write.table(genotype_bin, paste0("/lustre/scratch126/casm/team154pc/em16/Mutographs/",project,"/",PDID,"/analysis/data/",PDID, "_Genotype.tsv"), sep= "\t")

##Create fasta file for tree building input
dna_strings = list()

muts <- read.table(paste0("/lustre/scratch126/casm/team154pc/em16/Mutographs/",project,"/",PDID,"/analysis/data/",PDID,"_Genotype.tsv"), sep= "\t")

muts$counts <- rowSums(muts[,])
shared_muts <- muts[muts$counts>1,]

shared_muts <- shared_muts[,1:(ncol(shared_muts)-1)]

Muts=rownames(shared_muts)
Ref = substring(Muts, nchar(Muts)-2,nchar(Muts)-2)
Alt = substring(Muts, nchar(Muts))

dna_strings[1]=paste(Ref,sep="",collapse="")

for (k in 1:ncol(shared_muts)){
  Mutations = Ref
  Mutations[shared_muts[,k]==1] = Alt[shared_muts[,k]==1]
  dna_string = paste(Mutations,sep="",collapse="")
  dna_strings[k+1]=dna_string
}
names(dna_strings)=c("Ancestral",colnames(shared_muts))
require(seqinr)
write.fasta(dna_strings, names=names(dna_strings),paste0("/lustre/scratch126/casm/team154pc/em16/Mutographs/",project,"/",PDID,"/analysis/tree_files/Mutations_MPBoot_",PDID,".fa"))

system(paste0("/lustre/scratch117/casm/team154/tc16/Programs/mpboot-sse-1.1.0-Linux/bin/mpboot -s ", dna_string_file," -bb 1000"))

#Import the tree into R using ape
tree <- read.tree(mpboot_tree_file)
tree <- drop.tip(tree,"Ancestral")
tree <- multi2di(tree)
tree$edge.length = rep(1, nrow(tree$edge)) #Initially need to assign edge lengths of 1 for the tree_muts package to work


###############################################################################################################################################################
#Assign SNV and INDEL mutations back to the tree
###############################################################################################################################################################

#Assign mutations back to the tree
setwd("/lustre/scratch119/casm/team154pc/ms56/fetal_HSC/treemut"); source("treemut.R"); setwd(paste0(my_working_directory)) #R scripts have to be sourced within containing directory or outputs an error.
df = reconstruct_genotype_summary(tree) #Define df (data frame) for treeshape

#Find columns with MTR only
mtr <- as.matrix(NV[,grep("MTR", colnames(NV), value=T)])

#Find columns with DEP only
depth <- as.matrix(NR[,grep("DEP", colnames(NR), value=T)])

Muts <- paste(NV$Chrom, NV$Pos, NV$Ref, NV$Alt,sep="_")
rownames(mtr)=rownames(depth)=Muts
colnames(mtr) <- gsub("_MTR","",colnames(mtr))
colnames(depth) <- gsub("_DEP","",colnames(depth))

#Get matrices in order, and run the main assignment functions
p.error = c(rep(0.01, ncol(mtr)))
res = assign_to_tree(mtr[,df$samples], depth[,df$samples], df, error_rate = p.error) #Get res (results!) object
treefit_pval_cutoff = 1e-3
poor_fit = res$summary$pval < treefit_pval_cutoff  #See how many mutations don't have read counts that fit the tree very well
sum(poor_fit)

#Assign edge lengths from the res object
tree$edge.length <- res$df$df$edge_length
#Add node information to the muts object
mat <- NV[,1:5]
mat$node <- tree$edge[res$summary$edge_ml,2]

#Save the res file and the tree file

save(res, file = (paste0("/lustre/scratch126/casm/team154pc/em16/Mutographs/",project,"/",PDID,"/analysis/data/",PDID,"_res_file")))
#load(res_file)

write.tree(tree, file = (paste0("/lustre/scratch126/casm/team154pc/em16/Mutographs/",project,"/",PDID,"/analysis/tree_files/",PDID,".tree")))


###############################################################################################################################################################
#Generate input for hdp (SNVs only)
###############################################################################################################################################################

all_muts <- mat[mat$mut_type == "snv",]
all_muts$SampleID <- paste0(PDID, "_", all_muts$node)
all_muts <- all_muts[,c(1:4,7)]
colnames(all_muts) <- c("Chr","Pos","Ref","Alt","SampleID")
write.table(all_muts, paste0("/lustre/scratch126/casm/team154pc/em16/Mutographs/",project,"/",PDID,"/analysis/hdp/",PDID, "_hdp_snv_input.txt"), sep="\t", col.names = T, row.names = F, quote=F)


###############################################################################################################################################################
#Generate input for hdp (indels only)
###############################################################################################################################################################

all_indels <- mat[mat$mut_type == "indel",]
all_indels$SampleID <- paste0(PDID, "_", all_indels$node)
all_indels <- all_indels[,c(1:4,7)]
colnames(all_indels) <- c("Chr","Pos","Ref","Alt","SampleID")
write.table(all_indels, paste0("/lustre/scratch126/casm/team154pc/em16/Mutographs/",project,"/",PDID,"/analysis/hdp/",PDID, "_hdp_indel_input.txt"), sep="\t", col.names = T, row.names = F, quote=F)


###############################################################################################################################################################
#Generate input for dnds (all mutations)
###############################################################################################################################################################

dnds_input <- mat[,1:4]
dnds_input$SampleID <- PDID
colnames(dnds_input) <- c("Chr","Pos","Ref","Alt","SampleID")
dnds_input <- dnds_input[,c("SampleID","Chr","Pos","Ref","Alt")]
write.table(dnds_input, paste0("/lustre/scratch126/casm/team154pc/em16/Mutographs/",project,"/",PDID,"/analysis/dnds/",PDID, "_dnds_all.txt"), sep="\t", col.names = T, row.names = F, quote=F)


###############################################################################################################################################################
#Run VAGRENT on full mutation set to annotate the filtered mutations
###############################################################################################################################################################

##create mutations all file ###

vcf_header_path = "/lustre/scratch126/casm/team154pc/em16/R_scripts/VCF_header_for_VaGrent.txt"
vcf_path = paste0("/lustre/scratch126/casm/team154pc/em16/Mutographs/",project,"/",PDID,"/analysis/vcf/",PDID,"_mutations_all.vcf")
vagrent_input_path = paste0("/lustre/scratch126/casm/team154pc/em16/Mutographs/",project,"/",PDID,"/analysis/vcf/",PDID,"_mutations_all_header.vcf")
vagrent_output_path = paste0(vagrent_input_path,".annot")

#1. paste vcf file to a dummy header file
system(paste0("cat ",vcf_header_path," ",vcf_path," > ", vagrent_input_path))
#2. commands to run vagrent
system(paste0("AnnotateVcf.pl -i ",vagrent_input_path," -o ",vagrent_output_path," -sp Human -as NCBI37 -c /lustre/scratch126/casm/team154pc/em16/R_scripts/vagrent.cache.gz"))
#3. import vagrent output
vagrent_output = fread(vagrent_output_path,skip = "#CHROM")
annot_info = as.data.frame(str_split(vagrent_output$INFO, pattern = ";",simplify = TRUE), stringsAsFactors = FALSE)
colnames(annot_info) <- c("VT","VD","VC","VW")

annot_info$VC <- gsub(x=annot_info$VC, pattern = "VC=", replacement = "")
annot_info$VT <- gsub(x=annot_info$VT, pattern = "VT=", replacement = "")
annot_info$VW <- gsub(x=annot_info$VW, pattern = "VW=", replacement = "")
annot_info$VD <- gsub(x=annot_info$VD, pattern = "VD=", replacement = "")


#Attempt to functionalize neatly
split_vagrent_output = function(df,split_col,col_IDs = c("Gene","Transcript","RNA","CDS","Protein","Type","SO_codes")) {
  col = df[[split_col]]
  output = matrix(nrow = nrow(df), ncol = length(col_IDs))
  for(i in 1:length(col_IDs)) {
    output[,i] = str_split(col, pattern = "\\|", simplify = TRUE)[,i]
  }
  colnames(output) = col_IDs
  return(as.data.frame(output))
}

mat <- cbind(mat,split_vagrent_output(df = annot_info,split_col = "VD"))

###############################################################################################################################################################
#Annotate mutation file with additional categories
###############################################################################################################################################################

##Re-label the key objects for the plot_tree functions
mat$Type <- as.character(mat$Type)
mat$Gene <- as.character(mat$Gene)
mat$variant_ID <- paste(mat$Gene, mat$Protein, sep = " ")
mat$Type[mat$Type == ""] <- "no_annotation"

##Define different classes of mutations
mat$protein_coding_mutation <- ifelse(mat$Type %in% c("protein_coding:exon:CDS:substitution:codon_variant:non_synonymous_codon",
                                                      "protein_coding:exon:CDS:insertion:frameshift_variant",
                                                      "protein_coding:exon:CDS:substitution:codon_variant:stop_gained",
                                                      "protein_coding:exon:CDS:substitution:codon_variant:initiator_codon_change",
                                                      "protein_coding:exon:CDS:deletion:frameshift_variant",
                                                      "protein_coding:exon:CDS:deletion:inframe_variant:inframe_codon_loss"),
                                          "Protein coding mutation", "no")
mat$splice_variant <- ifelse(grepl(pattern = "splice", x = mat$Type),
                                 "Splice variant", "no")
mat$exon_UTR <- ifelse(grepl(pattern = "UTR", x = mat$Type) & grepl(pattern = "exon", x = mat$Type),
                           "UTR exon", "no")
mat$intron_UTR <- ifelse(grepl(pattern = "UTR", x = mat$Type) & grepl(pattern = "intron", x = mat$Type),
                             "UTR intron", "no")

mat$protein_coding_chip_variant = ifelse(mat$Gene %in% chip_drivers & mat$protein_coding_mutation == "Protein coding mutation",
                                             "Protein coding variant in driver gene", "no")


mat$protein_coding_cancer_variant = ifelse(mat$Gene %in% cancer_drivers & mat$protein_coding_mutation == "Protein coding mutation",
                                         "Protein coding variant in cancer gene", "no")

mat$non_coding = ifelse(mat$Gene %in% c("ZFP36L2","HNRNPA3","MET","NUS1","KRT23"),
                            "Non coding driver", "no")

mat$biobank_driver = ifelse(mat$Gene %in% biobank_drivers & mat$protein_coding_mutation == "Protein coding mutation",
                        "Biobank driver mutation", "no")

#Save the annotated filtered_muts files (post tree filtering)
write.csv(mat, file = paste0("/lustre/scratch126/casm/team154pc/em16/Mutographs/",project,"/",PDID,"/analysis/data/",PDID,"_annotated_muts.csv"), row.names = F) 


###############################################################################################################################################################
#Tree pdfs with various annotations
###############################################################################################################################################################

##Plot tree with labels
pdf(paste0("/lustre/scratch126/casm/team154pc/em16/Mutographs/",project,"/",PDID,"/analysis/tree_pdfs/",PDID,"_thin.pdf"),width=8,height=6)
tree_em =plot_tree_em(tree, cex.label = 0)
dev.off()

##Plot tree no labels
pdf(paste0("/lustre/scratch126/casm/team154pc/em16/Mutographs/",project,"/",PDID,"/analysis/tree_pdfs/",PDID,"_labels.pdf"),width=8,height=6)
tree_em =plot_tree_em(tree, cex.label = 0.5)
dev.off()

pdf(paste0("/lustre/scratch126/casm/team154pc/em16/Mutographs/",project,"/",PDID,"/analysis/tree_pdfs/",PDID,"_drivers.pdf"),w=8,h=6)
tree_em =plot_tree_em(tree, cex.label = 0)
plot_tree_labels(tree_em,
                 details = mat,
                 type = "label",
                 query.field = "protein_coding_chip_variant",
                 data.frame(value="Protein coding variant in driver gene",col="red",pch = 17,stringsAsFactors = FALSE),
                 label.field = "variant_ID",
                 cex.label = 1.0)
dev.off()

pdf(paste0("/lustre/scratch126/casm/team154pc/em16/Mutographs/",project,"/",PDID,"/analysis/tree_pdfs/",PDID,"_drivers_labels.pdf"),w=8,h=6)
tree_em=plot_tree_em(tree, cex.label = 0.5)
plot_tree_labels(tree_em,
                 details = mat,
                 type = "label",
                 query.field = "protein_coding_chip_variant",
                 data.frame(value="Protein coding variant in driver gene",col="red",pch = 17,stringsAsFactors = FALSE),
                 label.field = "variant_ID",
                 cex.label = 1.0
)
dev.off()

pdf(paste0("/lustre/scratch126/casm/team154pc/em16/Mutographs/",project,"/",PDID,"/analysis/tree_pdfs/",PDID,"_cancer.pdf"),w=8,h=6)
tree_em =plot_tree_em(tree, cex.label = 0)
plot_tree_labels(tree_em,
                 details = mat,
                 type = "label",
                 query.field = "protein_coding_cancer_variant",
                 data.frame(value="Protein coding variant in cancer gene",col="red",pch = 17,stringsAsFactors = FALSE),
                 label.field = "variant_ID",
                 cex.label = 1.0)
dev.off()

pdf(paste0("/lustre/scratch126/casm/team154pc/em16/Mutographs/",project,"/",PDID,"/analysis/tree_pdfs/",PDID,"_cancer_labels.pdf"),w=8,h=6)
tree_em=plot_tree_em(tree, cex.label = 0.5)
plot_tree_labels(tree_em,
                 details = mat,
                 type = "label",
                 query.field = "protein_coding_cancer_variant",
                 data.frame(value="Protein coding variant in cancer gene",col="red",pch = 17,stringsAsFactors = FALSE),
                 label.field = "variant_ID",
                 cex.label = 1.0
)
dev.off()

pdf(paste0("/lustre/scratch126/casm/team154pc/em16/Mutographs/",project,"/",PDID,"/analysis/tree_pdfs/",PDID,"_non_coding.pdf"),w=8,h=6)
tree_em =plot_tree_em(tree, cex.label = 0)
plot_tree_labels(tree_em,
                 details = mat,
                 type = "label",
                 query.field = "non_coding",
                 data.frame(value="Non coding driver",col="red",pch = 17,stringsAsFactors = FALSE),
                 label.field = "variant_ID",
                 cex.label = 1.0)
dev.off()

pdf(paste0("/lustre/scratch126/casm/team154pc/em16/Mutographs/",project,"/",PDID,"/analysis/tree_pdfs/",PDID,"_non_coding_labels.pdf"),w=8,h=6)
tree_em=plot_tree_em(tree, cex.label = 0.5)
plot_tree_labels(tree_em,
                 details = mat,
                 type = "label",
                 query.field = "non_coding",
                 data.frame(value="Non coding driver",col="red",pch = 17,stringsAsFactors = FALSE),
                 label.field = "variant_ID",
                 cex.label = 1.0
)
dev.off()

pdf(paste0("/lustre/scratch126/casm/team154pc/em16/Mutographs/",project,"/",PDID,"/analysis/tree_pdfs/",PDID,"_coding_labels.pdf"),w=8,h=6)
tree_em=plot_tree_em(tree, cex.label = 0.5)
plot_tree_labels(tree_em,
                 details = mat,
                 type = "label",
                 query.field = "protein_coding_mutation",
                 data.frame(value="Protein coding mutation",col="red",pch = 17,stringsAsFactors = FALSE),
                 label.field = "variant_ID",
                 cex.label = 1.0
)
dev.off()

pdf(paste0("/lustre/scratch126/casm/team154pc/em16/Mutographs/",project,"/",PDID,"/analysis/tree_pdfs/",PDID,"_biobank.pdf"),w=8,h=6)
tree_em=plot_tree_em(tree, cex.label = 0)
plot_tree_labels(tree_em,
                 details = mat,
                 type = "label",
                 query.field = "non_coding",
                 data.frame(value="Non coding driver",col="red",pch = 17,stringsAsFactors = FALSE),
                 label.field = "variant_ID",
                 cex.label = 1.0
)
dev.off()

pdf(paste0("/lustre/scratch126/casm/team154pc/em16/Mutographs/",project,"/",PDID,"/analysis/tree_pdfs/",PDID,"_biobank_labels.pdf"),w=8,h=6)
tree_em=plot_tree_em(tree, cex.label = 0.5)
plot_tree_labels(tree_em,
                 details = mat,
                 type = "label",
                 query.field = "biobank_driver",
                 data.frame(value="Biobank driver mutation",col="red",pch = 17,stringsAsFactors = FALSE),
                 label.field = "variant_ID",
                 cex.label = 1.0
)
dev.off()


###############################################################################################################################################################
#Save table of mutations of interest
###############################################################################################################################################################

drivers <- mat[mat$biobank_driver == "Biobank driver mutation" | mat$non_coding == "Non coding driver" | mat$protein_coding_cancer_variant == "Protein coding variant in cancer gene" | mat$protein_coding_chip_variant == "Protein coding variant in driver gene",]
drivers$Donor <- PDID
write.csv(drivers, paste0("/lustre/scratch126/casm/team154pc/em16/Mutographs/",project,"/",PDID,"/analysis/driver_info/",PDID,"_drivers.csv"), row.names = FALSE)

#scp -r em16@farm5-login:/lustre/scratch126/casm/team154pc/em16/Mutographs/ALK/PD37580/analysis/* /Users/em16/Documents/PhD/Sequencing_results/Mutographs/Alkylating_agents_2178/PD37580/analysis/
  
#