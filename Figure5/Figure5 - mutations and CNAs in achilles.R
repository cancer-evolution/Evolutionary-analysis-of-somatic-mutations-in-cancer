library(ggplot2)
library(reshape2)
library(readr)
library(gridExtra)
library(ggrepel)
library(ggsci)

load("regulator_classification.Rdata")

tumours <- c("LUAD", "LUSC", "BRCA", "PRAD", "LIHC", "COAD", "STAD")

mixed_regulators <- regulator_classification[regulator_classification$Regulator_class == "Mixed_downstream",]
UC_regulators <- regulator_classification[regulator_classification$Regulator_class == "UC_downstream",]
EM_regulators <- regulator_classification[regulator_classification$Regulator_class == "EM_downstream",]


genes_phy <- read.csv("gene_phylostrata.txt")  ###File with gene name, entrez and phylostratum as column
genes_phy$Age <- ifelse(genes_phy$Phylostrata %in% 1:3, "UC",
                        ifelse(genes_phy$Phylostrata %in% 4:9, "EM",
                               ifelse(genes_phy$Phylostrata %in% 10:16, "MM", NA)))
EM_genes <- as.character(genes_phy[genes_phy$Age == "EM", "GeneID"])

dependency <- read.csv("gene_dependency avana_public_18Q1, v3.csv", sep=",")  #File from Achilles

dependency_t <- t(dependency)
colnames(dependency_t) <- dependency_t[1,]
dependency_t <- dependency_t[-1,]
rownames(dependency_t) <- gsub("\\...*","",rownames(dependency_t))

dependency_melt <- melt(dependency_t)
colnames(dependency_melt) <- c("Gene", "Cell_line", "Probability")
dependency_melt$Probability <- as.numeric(as.character(dependency_melt$Probability))

blood_cell_lines <- dependency_melt[grep("HAEMATOPOIETIC_AND_LYMPHOID_TISSUE", dependency_melt$Cell_line),"Cell_line"]
blood_cell_lines <- unique(blood_cell_lines)

dependency_melt <- dependency_melt[!(dependency_melt$Cell_line %in% blood_cell_lines),]

#all_cell_lines <- unique(dependency_melt$Cell_line)
#cell_lines <- c("LUNG", "BREAST", "LIVER", "LARGE_INTESTINE", "STOMACH")
#cell_lines_interest <- grep(paste(cell_lines, collapse="|"), all_cell_lines, value=TRUE)
#dependency_melt2 <- dependency_melt[dependency_melt$Cell_line %in% cell_lines_interest,]


mutated_cell_lines_o <- read.table("CCLE_DepMap_18Q1_maf_20180207.txt",  #File from CCLE database
                                   header=TRUE, sep="\t")

#mutated_cell_lines <- mutated_cell_lines_o[mutated_cell_lines_o$isDeleterious == TRUE,]
mutated_cell_lines <- mutated_cell_lines_o[mutated_cell_lines_o$Variant_Classification != "Silent",]

mutated_cell_lines <- mutated_cell_lines[,c("Hugo_Symbol", "Tumor_Sample_Barcode", "Variant_Classification")]

colnames(mutated_cell_lines) <- c("Regulator", "Cell_line", "Variant")

mutated_cell_lines$Regulator_class <- ifelse(mutated_cell_lines$Regulator %in% UC_regulators$Regulator,
                                             "UC_downstream",
                                             ifelse(mutated_cell_lines$Regulator %in% mixed_regulators$Regulator,
                                                    "Mixed_downstream",
                                                    ifelse(mutated_cell_lines$Regulator %in% EM_regulators$Regulator,
                                                           "EM_downstream", "NOT")))

mutated_cell_lines <- mutated_cell_lines[mutated_cell_lines$Regulator_class != "NOT",]

mutated_cell_lines$Cell_line_original <- mutated_cell_lines$Cell_line

for(remove in c("_LIVER", "_LARGE_INTESTINE", "_LUNG", "_STOMACH", "_BREAST",
                "_SKIN", "_OESOPHAGUS", "_SOFT_TISSUE", "_UPPER_AERODIGESTIVE_TRACT",
                "_ENDOMETRIUM", "_KIDNEY", "_URINARY_TRACT", "_CENTRAL_NERVOUS_SYSTEM",
                "_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE", "_CERVIX", "_THYROID", "_TESTIS",
                "_BONE", "_OVARY", "_PLEURA", "_AUTONOMIC_GANGLIA", "_SMALL_INTESTINE",
                "_MATCHED_NORMAL_TISSUE", "_PANCREAS", "_PLACENTA", "_PROSTATE",
                "_SALIVARY_GLAND", "_BILIARY_TRACT", "FIBROBLAST", "_ADRENAL_CORTEX")){
  mutated_cell_lines$Cell_line <- gsub(remove, "", mutated_cell_lines$Cell_line)
}

blood_cell_lines <- mutated_cell_lines[grep("HAEMATOPOIETIC_AND_LYMPHOID_TISSUE", mutated_cell_lines$Cell_line_original),"Cell_line"]

mutated_cell_lines <- mutated_cell_lines[!(mutated_cell_lines$Cell_line %in% blood_cell_lines),]


##Load CNV data
CNV_cell_lines_o <- read.table("CCLE_copynumber_byGene_2013-12-03.txt",   #File from CCLE
                                   header=TRUE, sep="\t")
CNV_cell_lines <- CNV_cell_lines_o
CNV_cell_lines$EGID <- NULL
CNV_cell_lines$CHR <- NULL
CNV_cell_lines$CHRLOC <- NULL
CNV_cell_lines$CHRLOCEND <- NULL


CNV_cell_lines <- CNV_cell_lines[CNV_cell_lines$SYMBOL %in% c(UC_regulators$Regulator, mixed_regulators$Regulator, EM_regulators$Regulator),]

CNV_cell_lines <- melt(CNV_cell_lines)

colnames(CNV_cell_lines) <- c("Regulator", "Cell_line", "CNV")


CNV_cell_lines$Regulator_class <- ifelse(CNV_cell_lines$Regulator %in% UC_regulators$Regulator,
                                             "UC_downstream",
                                             ifelse(CNV_cell_lines$Regulator %in% mixed_regulators$Regulator,
                                                    "Mixed_downstream",
                                                    ifelse(CNV_cell_lines$Regulator %in% EM_regulators$Regulator,
                                                           "EM_downstream", "NOT")))



CNV_cell_lines$Cell_line_original <- CNV_cell_lines$Cell_line

for(remove in c("_LIVER", "_LARGE_INTESTINE", "_LUNG", "_STOMACH", "_BREAST",
                "_SKIN", "_OESOPHAGUS", "_SOFT_TISSUE", "_UPPER_AERODIGESTIVE_TRACT",
                "_ENDOMETRIUM", "_KIDNEY", "_URINARY_TRACT", "_CENTRAL_NERVOUS_SYSTEM",
                "_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE", "_CERVIX", "_THYROID", "_TESTIS",
                "_BONE", "_OVARY", "_PLEURA", "_AUTONOMIC_GANGLIA", "_SMALL_INTESTINE",
                "_MATCHED_NORMAL_TISSUE", "_PANCREAS", "_PLACENTA", "_PROSTATE",
                "_SALIVARY_GLAND", "_BILIARY_TRACT", "FIBROBLAST", "_ADRENAL_CORTEX")){
  CNV_cell_lines$Cell_line <- gsub(remove, "", CNV_cell_lines$Cell_line)
}

blood_cell_lines <- CNV_cell_lines[grep("HAEMATOPOIETIC_AND_LYMPHOID_TISSUE", CNV_cell_lines$Cell_line_original),"Cell_line"]

CNV_cell_lines <- CNV_cell_lines[!(CNV_cell_lines$Cell_line %in% blood_cell_lines),]


del_cut <- quantile(CNV_cell_lines$CNV[CNV_cell_lines$CNV < 0], 0.10)
amp_cut <- quantile(CNV_cell_lines$CNV[CNV_cell_lines$CNV > 0], 0.90)

to_keep <- ifelse(CNV_cell_lines$CNV >= amp_cut | CNV_cell_lines$CNV <= del_cut, "Keep", "Dont")

CNV_cell_lines$Keep <- to_keep
CNV_cell_lines_top <- CNV_cell_lines[CNV_cell_lines$Keep == "Keep",]


CNV_cell_lines_top$CNV_type <- ifelse(CNV_cell_lines_top$CNV > 0, "Amp", "Del")

Amp_cell_lines <- CNV_cell_lines_top[CNV_cell_lines_top$CNV_type == "Amp",c("Regulator", "Cell_line", "Cell_line_original", "Regulator_class")]

table(CNV_cell_lines_top$Regulator_class)
table(CNV_cell_lines$Regulator_class)

lines_with_mut_data <- unique(mutated_cell_lines$Cell_line_original)
lines_with_CNA_data <- unique(CNV_cell_lines_top$Cell_line_original)
#lines_with_dependency_data <- unique(dependency_melt2$Cell_line)
lines_with_dependency_data <- unique(dependency_melt$Cell_line)

lines_common <- lines_with_mut_data[lines_with_mut_data %in% lines_with_CNA_data]
lines_common <- lines_common[lines_common %in% lines_with_dependency_data]

mutated_cell_lines <- mutated_cell_lines[mutated_cell_lines$Cell_line_original %in% lines_common,]
CNV_cell_lines_top <- CNV_cell_lines_top[CNV_cell_lines_top$Cell_line_original %in% lines_common,]
dependency_melt <- dependency_melt[dependency_melt$Cell_line %in% lines_common,]
#dependency_melt2 <- dependency_melt2[dependency_melt2$Cell_line %in% lines_common,]


##Cell lines dependent because gene mutated or CNV?

dependency_p_mut <- vector()
dependency_p_amp <- vector()

all_reg <- sort(unique(c(as.character(mutated_cell_lines$Regulator), as.character(CNV_cell_lines_top$Regulator))))

for(reg in all_reg){
  local_depend <- dependency_melt[dependency_melt$Gene == reg,]
  #local_depend <- dependency_melt2[dependency_melt2$Gene == reg,]
  
  local_mut_lines <- mutated_cell_lines[mutated_cell_lines$Regulator == reg,]
  
  local_CNV_lines <- CNV_cell_lines_top[CNV_cell_lines_top$Regulator == reg,]
  
  local_Amp_lines <- local_CNV_lines[local_CNV_lines$CNV_type == "Amp",]
  #local_Del_lines <- local_CNV_lines[local_CNV_lines$CNV_type == "Del",]
  
  local_depend$Mutated <- ifelse(local_depend$Cell_line %in% local_mut_lines$Cell_line_original, "Mutated", "WT")
  local_depend$Amp <- ifelse(local_depend$Cell_line %in% local_Amp_lines$Cell_line_original, "Amp", "WT")
  #local_depend$Del <- ifelse(local_depend$Cell_line %in% local_Del_lines$Cell_line_original, "Del", "WT")
  
  if(sum(local_depend$Mutated == "Mutated") >= 3){
    p <- wilcox.test(local_depend$Probability[local_depend$Mutated == "Mutated"],
                     local_depend$Probability[local_depend$Mutated != "Mutated"])$p.val
    dependency_p_mut <- rbind(dependency_p_mut, c(reg, p, local_mut_lines$Regulator_class[1]))
  
  }
  if(sum(local_depend$Amp == "Amp") >= 3){
    p <- wilcox.test(local_depend$Probability[local_depend$Amp == "Amp"],
                     local_depend$Probability[local_depend$Amp != "Amp"])$p.val
    dependency_p_amp <- rbind(dependency_p_amp, c(reg, p, local_Amp_lines$Regulator_class[1]))
    
  }
  
  print(reg)
}

#save(dependency_p_mut, file="dependency_p_mut.Rdata")
#save(dependency_p_amp, file="dependency_p_amp.Rdata")

load("dependency_p_mut.Rdata")
load("dependency_p_amp.Rdata")


dependency_p_mut <- as.data.frame(dependency_p_mut)
colnames(dependency_p_mut) <- c("Regulator", "p", "Regulator_class")
dependency_p_mut$p <- as.numeric(as.character(dependency_p_mut$p))
dependency_p_mut$adj.p <- p.adjust(dependency_p_mut$p, method="BH")

dependency_p_amp <- as.data.frame(dependency_p_amp)
colnames(dependency_p_amp) <- c("Regulator", "p", "Regulator_class")
dependency_p_amp$p <- as.numeric(as.character(dependency_p_amp$p))
dependency_p_amp$adj.p <- p.adjust(dependency_p_amp$p, method="BH")


dependency_p_sig_mut <- dependency_p_mut[dependency_p_mut$p < 0.05,]
sum(dependency_p_sig_mut$Regulator_class == "Mixed_downstream")
sum(dependency_p_mut$Regulator_class == "Mixed_downstream")


dependency_p_sig_amp <- dependency_p_amp[dependency_p_amp$p < 0.05,]
sum(dependency_p_sig_amp$Regulator_class == "Mixed_downstream")
sum(dependency_p_amp$Regulator_class == "Mixed_downstream")



n_mixed_mut <- aggregate(Regulator ~ Regulator_class, dependency_p_sig_mut, length)

n_mixed_mut$Regulator_class <- factor(n_mixed_mut$Regulator_class,
                                  levels=c("UC_downstream", "Mixed_downstream", "EM_downstream"))

n_mixed_amp <- aggregate(Regulator ~ Regulator_class, dependency_p_sig_amp, length)

n_mixed_amp$Regulator_class <- factor(n_mixed_amp$Regulator_class,
                                      levels=c("UC_downstream", "Mixed_downstream", "EM_downstream"))

n_mixed_amp$Regulator/sum(n_mixed_amp$Regulator)

##Difference in median dependency
difference_in_median_mut <- vector()
for(reg in sort(unique(mutated_cell_lines$Regulator))){
  local_mut_lines <- mutated_cell_lines[mutated_cell_lines$Regulator == reg,]
  
  local_depend <- dependency_melt[dependency_melt$Gene == reg,]
  #local_depend <- dependency_melt2[dependency_melt2$Gene == reg,]
  
  local_depend$Mutated <- ifelse(local_depend$Cell_line %in% local_mut_lines$Cell_line_original, "Mutated", "WT")
  
  if(sum(local_depend$Mutated == "Mutated") >= 3){
    m1 <- median(local_depend$Probability[local_depend$Mutated == "Mutated"])
    m2 <- median(local_depend$Probability[local_depend$Mutated != "Mutated"])
    
    difference_in_median_mut <- rbind(difference_in_median_mut,
                                        c(reg, m1, m2, m1-m2))
  }
  print(reg)
}

difference_in_median_mut <- as.data.frame(difference_in_median_mut)

colnames(difference_in_median_mut) <- c("Regulator", "Median_when_mutated", "Median_when_WT", "Difference_in_median")
difference_in_median_mut$Median_when_mutated <- as.numeric(as.character(difference_in_median_mut$Median_when_mutated))
difference_in_median_mut$Median_when_WT <- as.numeric(as.character(difference_in_median_mut$Median_when_WT))
difference_in_median_mut$Difference_in_median <- as.numeric(as.character(difference_in_median_mut$Difference_in_median))

#save(difference_in_median_mut, file="difference_in_median_mut.Rdata")

load("difference_in_median_mut.Rdata")


dependency_p_sig_mut$Difference_in_median_all <- difference_in_median_mut[match(dependency_p_sig_mut$Regulator, difference_in_median_mut$Regulator), "Difference_in_median"]

dependency_p_sig_mut$Regulator_class <- factor(dependency_p_sig_mut$Regulator_class,
                                               levels=c("UC_downstream", "Mixed_downstream", "EM_downstream"))

dependency_p_sig_mut$Regulator <- factor(dependency_p_sig_mut$Regulator,
                                         levels=dependency_p_sig_mut$Regulator[order(dependency_p_sig_mut$Difference_in_median_all)])

##Difference in median dependency - amp
difference_in_median_amp <- vector()
for(reg in sort(unique(Amp_cell_lines$Regulator))){
  local_mut_lines <- Amp_cell_lines[Amp_cell_lines$Regulator == reg,]
  
  local_depend <- dependency_melt[dependency_melt$Gene == reg,]
  #local_depend <- dependency_melt2[dependency_melt2$Gene == reg,]
  
  local_depend$Amplified <- ifelse(local_depend$Cell_line %in% local_mut_lines$Cell_line_original, "Amplified", "WT")
  
  if(sum(local_depend$Amplified == "Amplified") >= 3){
    m1 <- median(local_depend$Probability[local_depend$Amplified == "Amplified"])
    m2 <- median(local_depend$Probability[local_depend$Amplified != "Amplified"])
    
    difference_in_median_amp <- rbind(difference_in_median_amp,
                                      c(reg, m1, m2, m1-m2))
  }
  print(reg)
}

difference_in_median_amp <- as.data.frame(difference_in_median_amp)

colnames(difference_in_median_amp) <- c("Regulator", "Median_when_Amplified", "Median_when_WT", "Difference_in_median")
difference_in_median_amp$Median_when_Amplified <- as.numeric(as.character(difference_in_median_amp$Median_when_Amplified))
difference_in_median_amp$Median_when_WT <- as.numeric(as.character(difference_in_median_amp$Median_when_WT))
difference_in_median_amp$Difference_in_median <- as.numeric(as.character(difference_in_median_amp$Difference_in_median))


#save(difference_in_median_amp, file="difference_in_median_amp.Rdata")

load("difference_in_median_amp.Rdata")

dependency_p_sig_amp$Difference_in_median_all <- difference_in_median_amp[match(dependency_p_sig_amp$Regulator, difference_in_median_amp$Regulator), "Difference_in_median"]

dependency_p_sig_amp$Regulator_class <- factor(dependency_p_sig_amp$Regulator_class,
                                           levels=c("UC_downstream", "Mixed_downstream", "EM_downstream"))

dependency_p_sig_amp$Regulator <- factor(dependency_p_sig_amp$Regulator,
                                     levels=dependency_p_sig_amp$Regulator[order(dependency_p_sig_amp$Difference_in_median_all)])

dependency_p_sig_mut$Sig_by <- "Mutation"
dependency_p_sig_amp$Sig_by <- "Amplification"

dependency_p_both <- rbind(dependency_p_sig_mut, dependency_p_sig_amp)


dependency_p_both$Regulator <- as.character(dependency_p_both$Regulator)

dependency_p_both$Regulator_star <- ifelse(dependency_p_both$Sig_by == "Mutation",
                                           paste("*", dependency_p_both$Regulator, sep=""), as.character(dependency_p_both$Regulator))

dependency_p_both$Regulator_star <- factor(dependency_p_both$Regulator_star,
                                           levels=dependency_p_both$Regulator_star[order(dependency_p_both$Difference_in_median_all)])

dependency_p_both_print <- dependency_p_both[,c("Regulator","Regulator_class", "Difference_in_median_all","p","Sig_by")]
colnames(dependency_p_both_print) <- c("Regulator", "Regulator_class", "Difference_in_median_dependency", "p_value", "Significant_by")

write.table(dependency_p_both_print, file="Supp_table_Difference_in_dependency.txt",
            sep="\t", quote=FALSE, row.names=FALSE)

dependency_p_both_top <- dependency_p_both[abs(dependency_p_both$Difference_in_median_all) >= 0.20,]
top_mut <- dependency_p_both_top[dependency_p_both_top$Sig_by == "Mutation",]
top_amp <- dependency_p_both_top[dependency_p_both_top$Sig_by == "Amplification",]


sum(top_mut$Difference_in_median_all > 0)/nrow(top_mut)
sum(top_amp$Difference_in_median_all > 0)/nrow(top_amp)


pdf("Figure5_Mutations_and_dependency.pdf",
    width=6, height=6.5)
g <- ggplot(dependency_p_both_top, aes(x=Regulator_star, y=Difference_in_median_all))+
  geom_bar(stat='identity', aes(fill=Regulator_class, color=Sig_by), size=0.75)+
  geom_hline(yintercept=0)+
  scale_fill_jco()+
  scale_colour_manual(values=c(Mutation="black", Amplification=NA))+
  coord_flip()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
print(g)
dev.off()


temp <- aggregate(Regulator ~ Regulator_class, dependency_p_both_top, length)
total <- sum(temp$Regulator)

temp$Percentage <- NA
temp$Percentage[temp$Regulator_class == "UC_downstream"] <- (temp$Regulator[temp$Regulator_class == "UC_downstream"]/total)*100
temp$Percentage[temp$Regulator_class == "Mixed_downstream"] <- (temp$Regulator[temp$Regulator_class == "Mixed_downstream"]/total)*100
temp$Percentage[temp$Regulator_class == "EM_downstream"] <- (temp$Regulator[temp$Regulator_class == "EM_downstream"]/total)*100


