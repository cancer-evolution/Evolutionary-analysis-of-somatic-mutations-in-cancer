##Expression of targets of regulators with point mutations

library(igraph)
library(readr)
library(ggplot2)
library(ggsci)

source("helper_functions.R")

genes_phy <- read.csv("gene_phylostrata.txt")  ###File with gene name, entrez and phylostratum as column
genes_phy$Age <- ifelse(genes_phy$Phylostrata %in% 1:3, "UC",
                        ifelse(genes_phy$Phylostrata %in% 4:9, "EM",
                               ifelse(genes_phy$Phylostrata %in% 10:16, "MM", NA)))
UC_genes <- as.character(genes_phy[genes_phy$Age == "UC", "GeneID"])
EM_genes <- as.character(genes_phy[genes_phy$Age == "EM", "GeneID"])
MM_genes <- as.character(genes_phy[genes_phy$Age == "MM", "GeneID"])

##Determine the downstream targets of each regulator
GRN_network <- load_PathwayCommons_GRN()

regulators <- unique(GRN_network$Gene1)

mut_to_exclude <- read_csv("Genes to exclude.txt", 
                           col_names = FALSE)
mut_to_exclude <- unname(unlist(mut_to_exclude[,1]))

regulators <- regulators[!(regulators %in% mut_to_exclude)]

regulator_targets <- lapply(regulators, function(reg){
  local_edges <- GRN_network[GRN_network$Gene1 == reg,]
  return(unique(local_edges$Gene2))
})

names(regulator_targets) <- regulators

#Find in which samples the regulators are altered

load("variants_LoF.Rdata")
load("variants_missense.Rdata")

tumours <- c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "ESCA", "GBM", 
             "HNSC", "KICH", "KIRC", "KIRP", "LGG",  "LIHC", "LUAD", "LUSC", 
             "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", 
             "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM")

load("patients_with_mut_info.Rdata")

##Mutations in patients
patients_with_regulator_mutations <- vector()
for(tumour in tumours){
  local_variants_LoF <- variants_LoF[[tumour]]
  local_variants_LoF <- local_variants_LoF[,c("Hugo_Symbol", "Tumor_Sample_Barcode")]
  colnames(local_variants_LoF) <- c("Gene", "Patient")
  local_variants_LoF$Patient <- substr(local_variants_LoF$Patient, 9,12)
  local_variants_LoF$Mut <- "LoF"
  local_variants_LoF$Tumour <- tumour
  regulator_variants_LoF <- local_variants_LoF[local_variants_LoF$Gene %in% regulators,]
  
  local_variants_missense <- variants_missense[[tumour]]
  local_variants_missense <- local_variants_missense[,c("Hugo_Symbol", "Tumor_Sample_Barcode")]
  colnames(local_variants_missense) <- c("Gene", "Patient")
  local_variants_missense$Patient <- substr(local_variants_missense$Patient, 9,12)
  local_variants_missense$Mut <- "Miss"
  local_variants_missense$Tumour <- tumour
  regulator_variants_missense <- local_variants_missense[local_variants_missense$Gene %in% regulators,]
  
  patients_with_regulator_mutations <- rbind(patients_with_regulator_mutations,
                                             regulator_variants_LoF, regulator_variants_missense)
}

##Load expression data, keeping only samples with mutation data

expression <- list()
for(tumour in tumours){
  path <- paste(tumour, "/gdac.broadinstitute.org_", tumour, 
                ".Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2016012800.0.0/", 
                tumour, 
                ".rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt",
                sep="")
  
  exp <- read.delim(path) ##Load normalized gene expression data from TCGA
  exp <- exp[-1,]
  
  exp_normal <- exp[,substr(colnames(exp), 14, 16) %in% c("11A", "11B", "11C")]
  exp_tumour <- exp[,substr(colnames(exp), 14, 16) %in% c("01A", "01B", "01C", "01D")]
  
  if(class(exp_normal) != "factor"){  ##If there is only 1 sample, ignore
    colnames(exp_normal) <- substr(colnames(exp_normal),9,12)
  }
  colnames(exp_tumour) <- substr(colnames(exp_tumour),9,12)
  
  all_genes <- as.character(exp[,1])
    
  all_genes <- unlist(strsplit(all_genes, "\\|"))[seq(1, length(all_genes)*2, 2)]
    
  genes_remove <- c("?", "SLC35E2")
  genes_remove_index <- which(all_genes %in% genes_remove)
  all_genes <- all_genes[-genes_remove_index]
    
  if(class(exp_normal) != "factor"){ 
    exp_normal <- exp_normal[-genes_remove_index,]
  }
  
  exp_tumour <- exp_tumour[-genes_remove_index,]
  
  if(!is.null(nrow(exp_normal))){
    rownames(exp_normal) <- all_genes
  }
  if(!is.null(nrow(exp_tumour))){
    rownames(exp_tumour) <- all_genes
  }
  
  expression[[tumour]][["normal"]] <- exp_normal
  expression[[tumour]][["tumour"]] <- exp_tumour
  print(tumour)
  
}


expression_pat <- list()
number_pat_to_keep <- vector()
for(tumour in tumours){
  patients_to_keep <- patients_with_mut_info[[tumour]]
  patients_normal <- colnames(expression[[tumour]][["normal"]])
  patients_tumour <- colnames(expression[[tumour]][["tumour"]])
  
  patients_to_keep <- patients_to_keep[patients_to_keep %in% patients_normal]
  patients_to_keep <- patients_to_keep[patients_to_keep %in% patients_tumour]
  
  number_pat_to_keep <- c(number_pat_to_keep, length(patients_to_keep))
  
  if(length(patients_to_keep) != 0){
    expression_pat[[tumour]][["tumour"]] <- expression[[tumour]][["tumour"]][,colnames(expression[[tumour]][["tumour"]]) %in% patients_to_keep]
    expression_pat[[tumour]][["normal"]] <- expression[[tumour]][["normal"]][,colnames(expression[[tumour]][["normal"]]) %in% patients_to_keep]
    print(tumour)
  }
}

names(number_pat_to_keep) <- tumours



#Are the patient names unique across tumour types? YES
all_patients <- vector()
for(tumour in tumours){
  patients <- colnames(expression_pat[[tumour]][["tumour"]])
  all_patients <- c(all_patients, patients)
}


##For each regulator, separate samples where the regulator is altered vs. those where it is not

patients_with_regulator_mutations <- patients_with_regulator_mutations[patients_with_regulator_mutations$Patient %in% all_patients,]

#Will do miss and LoF together
differential_expression_targets <- vector()
for(regulator in regulators){
  
  patients_with_mut <- patients_with_regulator_mutations[patients_with_regulator_mutations$Gene == regulator,]
  
  if(nrow(patients_with_mut) >= 3){
    tumours_where_mut_found <- unique(patients_with_mut$Tumour)
    expression_of_relevant_tumours <- rep(NA, 20500) ##Number of genes
    for(tumour in tumours_where_mut_found){
      expression_of_relevant_tumours <- cbind(expression_of_relevant_tumours, expression_pat[[tumour]][["tumour"]])
    }
    expression_of_relevant_tumours <- expression_of_relevant_tumours[,-1]
    
    local_targets <- regulator_targets[[regulator]]
    
    
    expression_of_relevant_tumours <- expression_of_relevant_tumours[rownames(expression_of_relevant_tumours) %in% local_targets,]
    
    expression_targets_with_mut <- expression_of_relevant_tumours[,colnames(expression_of_relevant_tumours) %in% patients_with_mut$Patient]
    
    expression_targets_without_mut <- expression_of_relevant_tumours[,!(colnames(expression_of_relevant_tumours) %in% patients_with_mut$Patient)]
    
    local_targets <- local_targets[local_targets %in% rownames(expression_targets_with_mut)]
    for(target in local_targets){
      expression_with_reg_mut <- expression_targets_with_mut[rownames(expression_targets_with_mut) == target,]
      expression_without_reg_mut <- expression_targets_without_mut[rownames(expression_targets_without_mut) == target,]
      
      p <- wilcox.test(as.numeric(as.character(unlist(expression_with_reg_mut))), 
                       as.numeric(as.character(unlist(expression_without_reg_mut))))$p.val
      
      differential_expression_targets <- rbind(differential_expression_targets, c(regulator, target, p))
    }
  }
  print(regulator)
}

#save(differential_expression_targets, file="differential_expression_targets.Rdata")

load("differential_expression_targets.Rdata")
differential_expression_targets <- as.data.frame(differential_expression_targets)
colnames(differential_expression_targets) <- c("Regulator", "Target", "Diff_exp")
differential_expression_targets$Diff_exp <- as.numeric(as.character(differential_expression_targets$Diff_exp))

##Calculate percentage of downstream targets that are DE for each regulator

per_targets_diff_exp <- vector()
for(regulator in unique(differential_expression_targets$Regulator)){
  local_diff_exp <- differential_expression_targets[differential_expression_targets$Regulator == regulator,]  
  local_diff_exp <- local_diff_exp[!is.na(local_diff_exp$Diff_exp),]
  
  per_targets_diff_exp <- rbind(per_targets_diff_exp,
                               c(regulator, sum(local_diff_exp$Diff_exp < 0.05)/nrow(local_diff_exp)*100))
}

per_targets_diff_exp <- as.data.frame(per_targets_diff_exp)
colnames(per_targets_diff_exp) <- c("Regulator", "Percentage_of_targets_DE")

per_targets_diff_exp$Percentage_of_targets_DE <- as.numeric(as.character(per_targets_diff_exp$Percentage_of_targets_DE))

genes_phy <- read.csv("gene_phylostrata.txt")  ###File with gene name, entrez and phylostratum as column
genes_phy$Age <- ifelse(genes_phy$Phylostrata %in% 1:3, "UC",
                        ifelse(genes_phy$Phylostrata %in% 4:9, "EM",
                               ifelse(genes_phy$Phylostrata %in% 10:16, "MM", NA)))

per_targets_diff_exp$Age <- genes_phy[match(per_targets_diff_exp$Regulator, genes_phy$GeneID), "Age"]
per_targets_diff_exp <- per_targets_diff_exp[!is.na(per_targets_diff_exp$Age),]

per_targets_diff_exp_fraction <- per_targets_diff_exp

per_targets_diff_exp_fraction$Fraction <- ifelse(per_targets_diff_exp_fraction$Percentage_of_targets_DE < 5, "<5%", ">5%")

counts <- aggregate(Regulator ~ Age+Fraction, per_targets_diff_exp_fraction, length)

counts$Fraction <- factor(counts$Fraction, levels=c("<5%", ">5%"))

counts$Age <- factor(counts$Age, levels=c("UC", "EM", "MM"))


ratio_of_fractions <- vector()
for(fraction in unique(counts$Fraction)){
  local_fraction <- counts[counts$Fraction == fraction,]
  ratio <- local_fraction[local_fraction$Age == "EM", "Regulator"]/local_fraction[local_fraction$Age == "UC", "Regulator"]
  ratio_of_fractions <- rbind(ratio_of_fractions, ratio)
}
rownames(ratio_of_fractions) <- NULL
ratio_of_fractions <- as.data.frame(ratio_of_fractions)
ratio_of_fractions$Fractions <- unique(counts$Fraction)
colnames(ratio_of_fractions)[1] <- "Ratio_of_fractions"

pdf("Supp3_Prevelance_EM.pdf",
    height=3.5, width=4)
g <- ggplot(ratio_of_fractions, aes(x=Fractions, y =log10(Ratio_of_fractions)))+
  geom_bar(stat='identity', fill="blueviolet")+
  geom_hline(yintercept = 0)+
  xlab("Percentage of differentially expressed\ndownstream targets")+
  ylab("Prevelance of EM regulators")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
print(g)
dev.off()

##Divide regulators by whether they are disrupting UC-EM, UC or EM

load("regulator_classification.Rdata")

per_targets_diff_exp_fraction$Regulator_class <- regulator_classification[match(per_targets_diff_exp_fraction$Regulator, regulator_classification$Regulator), "Regulator_class"]

per_targets_diff_exp_fraction <- per_targets_diff_exp_fraction[per_targets_diff_exp_fraction$Regulator_class != "None",]

counts_class <- aggregate(Regulator ~ Regulator_class+Fraction, per_targets_diff_exp_fraction, length)

counts_class$Fraction <- factor(counts_class$Fraction, levels=c("<5%", ">5%"))

counts_class$Regulator_class <- factor(counts_class$Regulator_class, levels=c("UC_downstream", "Mixed_downstream",
                                                                  "EM_downstream", "MM_downstream", "None"))

##Divide EM genes in those that regulate EM and mixed neighbourhoods

relevant_genes <- per_targets_diff_exp_fraction[per_targets_diff_exp_fraction$Regulator_class == "Mixed_downstream",]
relevant_genes <- relevant_genes[relevant_genes$Age == "EM",]

per_targets_diff_exp_fraction$Age2 <- paste(per_targets_diff_exp_fraction$Age, per_targets_diff_exp_fraction$Regulator_class, sep="_")

counts_class2 <- aggregate(Regulator ~ Age2+Fraction+Regulator_class+Age, per_targets_diff_exp_fraction, length)

counts_class2$Fraction <- factor(counts_class2$Fraction, levels=c("<5%", ">5%"))

counts_class2$Age2 <- factor(counts_class2$Age2, levels=c("UC_UC_downstream",
                                                          "UC_Mixed_downstream",
                                                          "UC_EM_downstream", 
                                                          "EM_UC_downstream",
                                                          "EM_Mixed_downstream",
                                                          "EM_EM_downstream",
                                                          "MM_EM_downstream",
                                                          "UC_MM_downstream",
                                                          "UC_None"))

counts_class2 <- counts_class2[counts_class2$Regulator_class != "None",]
counts_class2 <- counts_class2[counts_class2$Regulator_class != "MM_downstream",]


counts_class2$Regulator_class <- factor(counts_class2$Regulator_class,
                                        levels=c("UC_downstream", "Mixed_downstream", "EM_downstream"))
counts_class2$Age <- factor(counts_class2$Age, levels=c("MM", "EM", "UC"))


pdf("Figure3_Pie_chart_more_5.pdf",
    height=4, width=4)
g <- ggplot(subset(counts_class2, Fraction == ">5%"), aes(x=1, y=Regulator))+
  geom_bar(aes(fill=Regulator_class), stat='identity', position='stack')+
  scale_fill_manual(values=c(UC_downstream="#0073C2FF",
                             Mixed_downstream="#EFC000FF",
                             EM_downstream="#868686FF"))+
  ggtitle(">5%")+
  coord_polar("y")+
  theme_void()
print(g)
dev.off()


pdf("Figure3_Pie_chart_less_5.pdf",
    height=4, width=4)
g <- ggplot(subset(counts_class2, Fraction == "<5%"), aes(x=1, y=Regulator))+
  geom_bar(aes(fill=Regulator_class), stat='identity', position='stack')+
  scale_fill_manual(values=c(UC_downstream="#0073C2FF",
                             Mixed_downstream="#EFC000FF",
                             EM_downstream="#868686FF"))+
  ggtitle("<5%")+
  coord_polar("y")+
  theme_void()
print(g)
dev.off()

sum(counts_class2$Regulator[counts_class2$Fraction == "<5%"]) ##51
sum(counts_class2$Regulator[counts_class2$Fraction == ">5%"]) ##16
##Mixed < 5: 22
##Mixed > 5: 12
fisher.test(cbind(c(12, 16), c(22,51)), alternative="greater")

#Of all the mixed ones: 34
###12 are high impact, 22 low impact.

fisher.test(cbind(c(12, 34), c(16,67)), alternative="greater")

prop.test(x=c(12, 22), n=c(16, 51), alternative="greater")



pdf("Figure3_Pie_chart_Mixed_more_5.pdf",
    height=4, width=4)
g <- ggplot(subset(counts_class2, Regulator_class=="Mixed_downstream" & Fraction == ">5%"), aes(x=1, y=Regulator))+
  geom_bar(aes(fill=Age), stat='identity', position='stack')+
  scale_fill_manual(values=c(UC="#F8766DFF",
                             EM="#00BA38FF"))+
  coord_polar("y")+
  theme_void()
print(g)
dev.off()


pdf("Figure3_Pie_chart_Mixed_less_5.pdf",
    height=4, width=4)
g <- ggplot(subset(counts_class2, Regulator_class=="Mixed_downstream" & Fraction == "<5%"), aes(x=1, y=Regulator))+
  geom_bar(aes(fill=Age), stat='identity', position='stack')+
  scale_fill_manual(values=c(UC="#F8766DFF",
                             EM="#00BA38FF"))+
  coord_polar("y")+
  theme_void()
print(g)
dev.off()


subset(counts_class2, Regulator_class=="Mixed_downstream" & Fraction %in% c(">5%"))
subset(counts_class2, Regulator_class=="Mixed_downstream" & Fraction %in% c("<5%"))
fisher.test(cbind(c(6,3), c(3,7)), alternative="greater")
