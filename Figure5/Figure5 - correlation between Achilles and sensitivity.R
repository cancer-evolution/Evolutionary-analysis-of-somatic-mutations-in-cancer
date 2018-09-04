library(ggplot2)
library(reshape2)
library(readr)
library(gridExtra)
library(ggrepel)
library(ggsci)

load("regulator_classification.Rdata")

tumours <- c("LUAD", "LUSC", "BRCA", "PRAD", "LIHC", "COAD", "STAD")

genes_phy <- read.csv("gene_phylostrata.txt")  ###File with gene name, entrez and phylostratum as column
genes_phy$Age <- ifelse(genes_phy$Phylostrata %in% 1:3, "UC",
                        ifelse(genes_phy$Phylostrata %in% 4:9, "EM",
                               ifelse(genes_phy$Phylostrata %in% 10:16, "MM", NA)))

regulator_classification$Age <- genes_phy[match(regulator_classification$Regulator, genes_phy$GeneID), "Age"]


dependency <- read.csv("gene_dependency avana_public_18Q1, v3.csv", sep=",")  #Gene dependency data from Achilles

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

dependency_melt$Regulator_class <- regulator_classification[match(dependency_melt$Gene, regulator_classification$Regulator), "Regulator_class"]
dependency_melt$Regulator <- ifelse(is.na(dependency_melt$Regulator_class), "Non_regulator", "Regulator")
dependency_melt$Age <- genes_phy[match(dependency_melt$Gene, genes_phy$GeneID), "Age"]

dependency_melt$Age <- factor(dependency_melt$Age,
                              levels=c("UC", "EM", "MM"))

dependency_melt_regulator <- dependency_melt[dependency_melt$Regulator == "Regulator",]

dependency_melt_regulator <- dependency_melt_regulator[dependency_melt_regulator$Regulator_class != "None",]
dependency_melt_regulator <- dependency_melt_regulator[dependency_melt_regulator$Regulator_class != "MM_downstream",]


dependency_melt_regulator$Regulator_class <- factor(dependency_melt_regulator$Regulator_class,
                                                    levels=c("UC_downstream",
                                                             "Mixed_downstream",
                                                             "EM_downstream"))


dependency_melt_regulator$Cell_line_short <- dependency_melt_regulator$Cell_line

for(remove in c("_LIVER", "_LARGE_INTESTINE", "_LUNG", "_STOMACH", "_BREAST",
                "_SKIN", "_OESOPHAGUS", "_SOFT_TISSUE", "_UPPER_AERODIGESTIVE_TRACT",
                "_ENDOMETRIUM", "_KIDNEY", "_URINARY_TRACT", "_CENTRAL_NERVOUS_SYSTEM",
                "_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE", "_CERVIX", "_THYROID", "_TESTIS",
                "_BONE", "_OVARY", "_PLEURA", "_AUTONOMIC_GANGLIA", "_SMALL_INTESTINE",
                "_MATCHED_NORMAL_TISSUE", "_PANCREAS", "_PLACENTA", "_PROSTATE",
                "_SALIVARY_GLAND", "_BILIARY_TRACT", "FIBROBLAST", "_ADRENAL_CORTEX")){
  dependency_melt_regulator$Cell_line_short <- gsub(remove, "", dependency_melt_regulator$Cell_line_short)
}


###Load files form the Genomics of Drug Sensitivity in Cancer site
gds_cell_lines <- read.table("Cell_Lines_Details.txt", sep="\t", header=TRUE)

gds_compounds <- read.table("Screened_compounds.txt", sep="\t", header=TRUE)

gds_response <- read.table("v17_fitted_dose_response.txt", sep="\t", header=TRUE)

gds_response$Drug <- gds_compounds[match(gds_response$DRUG_ID, gds_compounds$Drug_ID), "Drug_Name"]

gds_response$Cell_line <- gds_cell_lines[match(gds_response$COSMIC_ID, gds_cell_lines$COSMIC), "Sample"]
gds_response$Cell_line_tissue <- gds_cell_lines[match(gds_response$COSMIC_ID, gds_cell_lines$COSMIC), "Tissue_1"]
gds_response$Cell_line_TCGA <- gds_cell_lines[match(gds_response$COSMIC_ID, gds_cell_lines$COSMIC), "TCGA_label"]

gds_response$Cell_line <- gsub("-", "", gds_response$Cell_line)

overlapping_cell_lines <- unique(dependency_melt_regulator$Cell_line_short)[unique(dependency_melt_regulator$Cell_line_short) %in% unique(gds_response$Cell_line)]
##214 overlapping cell lines


gds_response_overlap <- gds_response[gds_response$Cell_line %in% overlapping_cell_lines,]
dependency_overlap <- dependency_melt_regulator[dependency_melt_regulator$Cell_line_short %in% overlapping_cell_lines,]


for(drug in unique(gds_response_overlap$Drug)){
  corr_drug_depend <- vector()
  local_response <- gds_response_overlap[gds_response_overlap$Drug == drug,]

  for(reg in unique(dependency_overlap$Gene)){
    local_dependency <- dependency_overlap[dependency_overlap$Gene == reg,]
    reg_class <- as.character(local_dependency[1,4])
    local_dependency$IC50 <- local_response[match(local_dependency$Cell_line_short, local_response$Cell_line), "LN_IC50"]
    c1 <- cor(local_dependency$IC50, local_dependency$Probability, use="pairwise.complete.obs", method="sp")
    p1 <- cor.test(local_dependency$IC50, local_dependency$Probability, use="pairwise.complete.obs", method="sp")$p.value

    corr_drug_depend <- rbind(corr_drug_depend, c(drug, reg, reg_class, c1, p1))
  }
  if(drug == "VNLG/124"){
    drug <- "VNLG-124"
  }
  save(corr_drug_depend, file=paste("corr_drug_depend_sp_", drug, ".Rdata", sep=""))
  print(drug)
}



corr_drug_depend_all <- vector()
for(drug in unique(gds_response_overlap$Drug)){
  if(drug == "VNLG/124"){
    drug <- "VNLG-124"
  }
  load(paste("corr_drug_depend_sp_", drug, ".Rdata", sep=""))
  corr_drug_depend_all <- rbind(corr_drug_depend_all, corr_drug_depend)
}

corr_drug_depend <- corr_drug_depend_all

corr_drug_depend <- as.data.frame(corr_drug_depend)
colnames(corr_drug_depend) <- c("Drug", "Regulator", "Regulator_class", "Spear_corr", "P")
corr_drug_depend$Spear_corr <- as.numeric(as.character(corr_drug_depend$Spear_corr))
corr_drug_depend$P <- as.numeric(as.character(corr_drug_depend$P))
corr_drug_depend$adj.p <- p.adjust(corr_drug_depend$P, method="BH")


corr_drug_depend$Regulator_class <- factor(corr_drug_depend$Regulator_class,
                                           levels=c("UC_downstream", "Mixed_downstream", "EM_downstream"))

pdf("Supp5_Correlation_IC50_dependency.pdf",
    width=8, height=3)
g <- ggplot(corr_drug_depend, aes(x=Spear_corr))+
  geom_density(aes(fill=Regulator_class))+
  scale_fill_jco()+
  facet_grid(.~Regulator_class)+
  xlab("Correlation between dependency and IC50")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
print(g)
dev.off()



corr_drug_depend_sig <- corr_drug_depend[corr_drug_depend$adj.p < 0.05,]
corr_drug_depend_sig <- corr_drug_depend_sig[corr_drug_depend_sig$Spear_corr < -0.25,]
table(corr_drug_depend_sig$Regulator_class)  #38 UC, 11 Mixed, 23, EM

corr_drug_depend_sig$Regulator_age <- genes_phy[match(corr_drug_depend_sig$Regulator, genes_phy$GeneID), "Age"]


##Enrichment, percentage, number of correlated drugs

table(corr_drug_depend$Regulator_class)

table(corr_drug_depend_sig$Regulator_class)

View(corr_drug_depend_sig[corr_drug_depend_sig$Regulator_class == "Mixed_downstream",])


regulator_drugs <- corr_drug_depend_sig[corr_drug_depend_sig$Regulator_class == "Mixed_downstream",c("Drug", "Regulator")]

relevant_dependency_mixed <- vector()
for(i in 1:nrow(regulator_drugs)){
  drug <- as.character(regulator_drugs[i,"Drug"])
  reg <- as.character(regulator_drugs[i, "Regulator"])
  
  local_response <- gds_response_overlap[gds_response_overlap$Drug == drug,]
  local_dependency <- dependency_overlap[dependency_overlap$Gene == reg,]
  
  local_dependency$IC50 <- local_response[match(local_dependency$Cell_line_short, local_response$Cell_line), "LN_IC50"]
  
  local_dependency <- local_dependency[!is.na(local_dependency$IC50),]
  
  local_dependency$Name <- paste(drug, reg, sep=" - ")
  local_dependency$Drug <- drug
  relevant_dependency_mixed <- rbind(relevant_dependency_mixed, local_dependency)
}

pdf("Figure5_Correlation_IC50.pdf",
    width=7, height=8)
g <- ggplot(relevant_dependency_mixed, aes(x=Probability, y =IC50))+
  geom_point()+
  geom_smooth(method = "lm", se=FALSE)+
  facet_wrap(~Name, scale='free_y', nrow=5)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank())
print(g)
dev.off()


##Distribution of correlations per drug
drug_relevant <- unique(relevant_dependency_mixed$Drug)

pdf("Supp5_Distribution_correlations_per_drug.pdf",
    height=2.5, width=3)
for(drug in drug_relevant){
  corr_drug_depend_local <- corr_drug_depend[corr_drug_depend$Drug == drug,]
  highlight <- relevant_dependency_mixed[relevant_dependency_mixed$Drug == drug,]
  if(length(unique(highlight$Gene)) == 1){
    gene <- unique(highlight$Gene)
    highlight_cor <- cor(highlight$IC50, highlight$Probability, use="pairwise.complete.obs", method="sp")
  }
  
  g <- ggplot(corr_drug_depend_local, aes(x = Spear_corr))+
    geom_density(fill='grey')+
    geom_vline(xintercept=highlight_cor, color="red")+
    ggtitle(paste(drug, gene, sep=" - "))+
    xlab("Correlation")+ylab("Density")+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          strip.background = element_blank())
  print(g)
}
dev.off()



pdf("Figure5_Correlation_IC50_individual.pdf",
    width=3.5, height=3)
for(local_name in unique(relevant_dependency_mixed$Name)){
  g <- ggplot(subset(relevant_dependency_mixed, Name == local_name), aes(x=Probability, y =IC50))+
    geom_point()+
    geom_smooth(method = "lm", se=FALSE, color="firebrick2")+
    ggtitle(local_name)+
    coord_cartesian(xlim = c(0, 1))+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          strip.background = element_blank())
  print(g)
}
dev.off()



##Calculating correlations based on AUC instead of IC50

corr_drug_depend_AUC <- vector()
for(drug in unique(relevant_dependency_mixed$Drug)){
  local_response <- gds_response_overlap[gds_response_overlap$Drug == drug,]
  
  for(reg in unique(relevant_dependency_mixed$Gene)){
    local_dependency <- dependency_overlap[dependency_overlap$Gene == reg,]
    reg_class <- as.character(local_dependency[1,4])
    local_dependency$AUC <- local_response[match(local_dependency$Cell_line_short, local_response$Cell_line), "AUC"]
    c1 <- cor(local_dependency$AUC, local_dependency$Probability, use="pairwise.complete.obs", method="sp")
    p1 <- cor.test(local_dependency$AUC, local_dependency$Probability, use="pairwise.complete.obs", method="sp")$p.value
    
    corr_drug_depend_AUC <- rbind(corr_drug_depend_AUC, c(drug, reg, reg_class, c1, p1))
  }
  print(drug)
}

corr_drug_depend_AUC <- as.data.frame(corr_drug_depend_AUC)
colnames(corr_drug_depend_AUC) <- c("Drug", "Regulator", "Regulator_class", "Spear_corr", "P")
corr_drug_depend_AUC$Spear_corr <- as.numeric(as.character(corr_drug_depend_AUC$Spear_corr))
corr_drug_depend_AUC$P <- as.numeric(as.character(corr_drug_depend_AUC$P))
corr_drug_depend_AUC$adj.p <- p.adjust(corr_drug_depend_AUC$P, method="BH")




##Cell lines with high dependency and low IC50
#drug <- "Temsirolimus"

cor_tissues <- vector()
reg <- "PPRC1"
for(drug in c("Temsirolimus", "Dactolisib", "YK-4-279", "Docetaxel")){
  local_response <- gds_response_overlap[gds_response_overlap$Drug == drug,]
  
  local_dependency <- dependency_overlap[dependency_overlap$Gene == reg,]
  local_dependency$IC50 <- local_response[match(local_dependency$Cell_line_short, local_response$Cell_line), "LN_IC50"]
  
  l1 <- strsplit(as.character(local_dependency$Cell_line), "_")
  l2 <- lapply(l1, function(x) x[-1])
  l3 <- sapply(l2, function(x){ paste(x, collapse="_")} )
  
  local_dependency$Cell_line_tissue <- l3
  
  for(tissue in unique(local_dependency$Cell_line_tissue)){
    temp <- local_dependency[local_dependency$Cell_line_tissue == tissue,]
    local_cor <- cor(temp$Probability, temp$IC50, use="pairwise.complete.obs", method="sp")
    number <- nrow(temp)
    cor_tissues <- rbind(cor_tissues, c(reg, drug, tissue, local_cor, number))
  }
}

cor_tissues <- as.data.frame(cor_tissues)
colnames(cor_tissues) <- c("Regulator", "Drug", "Cell_line_tissue", "Correlation", "Number_of_cell_lines")
cor_tissues$Correlation <- as.numeric(as.character(cor_tissues$Correlation))

cor_tissues <- cor_tissues[cor_tissues$Cell_line_tissue != "PLEURA",]
cor_tissues <- cor_tissues[cor_tissues$Cell_line_tissue != "THYROID",]

pdf("Supp5_Distribution_correlations_with PPRC1.pdf",
    height=6, width=7)
g <- ggplot(cor_tissues, aes(x=Cell_line_tissue, y=Correlation))+
  geom_bar(stat='identity')+
  geom_hline(yintercept=-0.25, color="red")+
  #scale_y_reverse()+
  ylim(0, -1)+
  ylab("Correlation between drug sensitivity\nand PPRC1 dependency")+
  xlab("Cell line tissue type")+
  facet_wrap(~Drug, nrow=2)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
print(g)
dev.off()
