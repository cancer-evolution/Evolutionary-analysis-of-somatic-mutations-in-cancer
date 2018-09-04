library(ggplot2)
library(reshape2)
library(gridExtra)
library(Biobase)
library(clinfun)
library(igraph)
library(ggsci)

source("helper_functions.R")

load("CNVs_curated2.Rdata")
load("CNV_classification.Rdata")
load("patients_with_CNV_info.Rdata")

genes_phy <- read.csv("gene_phylostrata.txt")  ###File with gene name, entrez and phylostratum as column
genes_phy_categorical <- genes_phy
genes_phy_categorical$Phylostrata <- ifelse(genes_phy_categorical$Phylostrata %in% 1:3, "UC",
                                            ifelse(genes_phy_categorical$Phylostrata %in% 4:9, "EM",
                                                   ifelse(genes_phy_categorical$Phylostrata %in% 10:16, "MM", NA)))

UC_genes <- as.character(genes_phy_categorical[genes_phy_categorical$Phylostrata == "UC", "GeneID"])
EM_genes <- as.character(genes_phy_categorical[genes_phy_categorical$Phylostrata == "EM", "GeneID"])
MM_genes <- as.character(genes_phy_categorical[genes_phy_categorical$Phylostrata == "MM", "GeneID"])

tumours <- c("LUAD", "LUSC", "BRCA", "PRAD", "LIHC", "COAD", "STAD")

##Load expression data
expression <- list()
for(tumour in tumours){
  for(sample_type in c("normal", "tumour")){
    exp <- read.delim() ##Load normalized gene expression data from TCGA
    if (sample_type == "normal"){
      colnames(exp)[2:ncol(exp)] <- substr(colnames(exp)[2:ncol(exp)],9,12)
    } else{
      colnames(exp)[2:ncol(exp)] <- substr(colnames(exp)[2:ncol(exp)],9,12)
    }
    samples_low_correlation <- readLines() #Load list of samples with low correlation with other samples
    if (substr(samples_low_correlation[1],9,12) != "")  #if there are samples to be removed
    {
      exp <- exp[,!(substr(colnames(exp),1,4) %in% substr(samples_low_correlation,9,12))]
    }
    
    all_genes <- as.character(exp[,1])
    
    all_genes <- unlist(strsplit(all_genes, "\\|"))[seq(1, length(all_genes)*2, 2)]
    
    genes_remove <- c("?", "SLC35E2")
    genes_remove_index <- which(all_genes %in% genes_remove)
    all_genes <- all_genes[-genes_remove_index]
    
    
    temp_exp <- exp[,2:ncol(exp)]
    temp_exp <- temp_exp[-genes_remove_index, ]
    rownames(temp_exp) <- all_genes
 
    expression[[tumour]][[sample_type]] <- temp_exp
  }
}

expression_pat <- list()
for(tumour in tumours){
  patients_to_keep <- patients_with_CNV_info[[tumour]]
  patients_normal <- colnames(expression[[tumour]][["normal"]])
  patients_tumour <- colnames(expression[[tumour]][["tumour"]])
  
  patients_to_keep <- patients_to_keep[patients_to_keep %in% patients_normal]
  patients_to_keep <- patients_to_keep[patients_to_keep %in% patients_tumour]
  
  expression_pat[[tumour]][["tumour"]] <- expression[[tumour]][["tumour"]][,colnames(expression[[tumour]][["tumour"]]) %in% patients_to_keep]
  expression_pat[[tumour]][["normal"]] <- expression[[tumour]][["normal"]][,colnames(expression[[tumour]][["normal"]]) %in% patients_to_keep]
}


patient_ratios <- list()
for(tumour in tumours){
  local_normal_exp <- expression_pat[[tumour]][["normal"]]
  local_tumour_exp <- expression_pat[[tumour]][["tumour"]]
  
  local_normal_exp <- local_normal_exp[rownames(local_normal_exp) %in% c(CNVs_curated[[tumour]]$amplifications$Gene,
                                                                         CNVs_curated[[tumour]]$deletions$Gene), ]

  local_tumour_exp <- local_tumour_exp[rownames(local_tumour_exp) %in% c(CNVs_curated[[tumour]]$amplifications$Gene,
                                                                         CNVs_curated[[tumour]]$deletions$Gene), ]
  
  genes_in_common <- rownames(local_normal_exp)[rownames(local_normal_exp) %in% rownames(local_tumour_exp)]
  
  local_normal_exp <- local_normal_exp[match(genes_in_common, rownames(local_normal_exp)),]
  local_tumour_exp <- local_tumour_exp[match(genes_in_common, rownames(local_tumour_exp)),]
  
  samples <- colnames(local_normal_exp)
  local_tumour_exp <- local_tumour_exp[,match(samples,colnames(local_tumour_exp))]
  
  patient_ratios[[tumour]] <- vector()
  for(i in 1:length(samples)){
    patient_ratio <- local_tumour_exp[,i]/local_normal_exp[,i]
    patient_ratios[[tumour]] <- cbind(patient_ratios[[tumour]], patient_ratio)
  }
 colnames(patient_ratios[[tumour]]) <- samples
 patient_ratios[[tumour]] <- as.data.frame(patient_ratios[[tumour]])
 rownames(patient_ratios[[tumour]]) <- genes_in_common
}
#save(patient_ratios, file="patient_ratios.Rdata")

p_values_df <- vector()
for(tumour in tumours){
  local_amp <- CNVs_curated[[tumour]]$amplifications
  local_del <- CNVs_curated[[tumour]]$deletions
  
  local_amp$Sample <- substr(local_amp$Sample, 9,12)
  local_del$Sample <- substr(local_del$Sample, 9,12)
  
  amp_genes <- unique(local_amp$Gene)
  del_genes <- unique(local_del$Gene)
  
  genes_in_samples <- rownames(patient_ratios[[tumour]])
  
  amp_genes <- amp_genes[amp_genes %in% genes_in_samples]
  del_genes <- del_genes[del_genes %in% genes_in_samples]
  
  amp_pvalues <- calculate_diff_exp_CNV2(amp_genes, local_amp, "Amp")
  del_pvalues <- calculate_diff_exp_CNV2(del_genes, local_del, "Del")
  
  amp_pvalues <- amp_pvalues[!is.na(amp_pvalues$Gene_age),]
  del_pvalues <- del_pvalues[!is.na(del_pvalues$Gene_age),]
  
  amp_pvalues$adj_p <- p.adjust(amp_pvalues$p_val, method="BH")
  del_pvalues$adj_p <- p.adjust(del_pvalues$p_val, method="BH")
  
  p_values_df <- rbind(p_values_df, amp_pvalues, del_pvalues)
  print(tumour)
}

#save(p_values_df, file="CNV_exp_p_values_df.Rdata")

load("CNV_exp_p_values_df.Rdata")
load("patient_ratios.Rdata")

##Use only targets

network <- load_PathwayCommons_GRN()
regulators <- unique(network$Gene1)
targets <- unique(network$Gene2)

p_values_df <- p_values_df[p_values_df$Gene %in% targets,]

for(tumour in tumours){
  patient_ratios[[tumour]] <- patient_ratios[[tumour]][rownames(patient_ratios[[tumour]]) %in% targets,]
  
}


p_values_df_sig <- p_values_df[p_values_df$adj_p < 0.05,]

observed <- aggregate(Gene_age ~ Tumour+CNV, p_values_df_sig, table)

observed$UC <- observed$Gene_age[,1]
observed$EM <- observed$Gene_age[,2]
observed$MM <- observed$Gene_age[,3]

observed$Gene_age <- NULL

observed_melt <- melt(observed)
colnames(observed_melt)[3:4] <- c("Age", "Number_sig")

#Calculate expected
expected <- vector()
for(tumour in tumours){
  local_amp <- CNVs_curated[[tumour]]$amplifications
  local_del <- CNVs_curated[[tumour]]$deletions
  
  local_amp$Sample <- substr(local_amp$Sample, 9,12)
  local_del$Sample <- substr(local_del$Sample, 9,12)
  
  amp_genes <- unique(local_amp$Gene)
  del_genes <- unique(local_del$Gene)
  
  genes_in_samples <- rownames(patient_ratios[[tumour]])
  
  amp_genes <- amp_genes[amp_genes %in% genes_in_samples]
  del_genes <- del_genes[del_genes %in% genes_in_samples]
  expected <- rbind(expected, c(tumour, sum(amp_genes %in% UC_genes), sum(amp_genes %in% EM_genes), sum(amp_genes %in% MM_genes),
                                sum(del_genes %in% UC_genes), sum(del_genes %in% EM_genes), sum(del_genes %in% MM_genes)))
}
expected <- as.data.frame(expected)
colnames(expected) <- c("Tumour", "N_Amp_UC", "N_Amp_EM", "N_Amp_MM", "N_Del_UC", "N_Del_EM", "N_Del_MM")

expected[,2:7] <- apply(expected[,2:7], 2, function(x){as.numeric(as.character(x))})

expected_melt <- melt(expected)
expected_melt$Age <- substr(expected_melt$variable, 7,8)
expected_melt$Type <- substr(expected_melt$variable, 3,5)

observed_melt$Name <- paste(observed_melt$Tumour, observed_melt$CNV, observed_melt$Age)
expected_melt$Name <- paste(expected_melt$Tumour, expected_melt$Type, expected_melt$Age)

counts_obs_exp <- data.frame(observed_melt,
                             Expected = expected_melt[match(observed_melt$Name, expected_melt$Name), "value"])

counts_obs_exp$Ratio <- counts_obs_exp$Number_sig/counts_obs_exp$Expected

counts_obs_exp$Percentage <- (counts_obs_exp$Ratio)*100

pdf("Figure4_Amp_expression_by_age.pdf",
    height=3.2, width=4)
g <- ggplot(subset(counts_obs_exp, CNV=="Amp"), aes(x=Age, y =Percentage))+
  geom_point(aes(colour=Tumour))+
  geom_boxplot()+
  geom_line(aes(group=Tumour), size=0.25, colour="grey")+
  geom_point(aes(colour=Tumour), size=2.5)+
  ylab("Percentage of differentially expressed genes")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())
print(g)
dev.off()


pdf("Supp4_CNA_expression_by_age.pdf",
    height=3.5)
g <- ggplot(counts_obs_exp, aes(x=Age, y =Percentage))+
  geom_point(aes(colour=Tumour))+
  geom_boxplot()+
  geom_line(aes(group=Tumour), size=0.25, colour="grey")+
  geom_point(aes(colour=Tumour), size=2.5)+
  ylab("Percentage of differentially expressed genes")+
  facet_grid(.~CNV)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())
print(g)
dev.off()

for(CNV in c("Amp", "Del")){
  temp <- counts_obs_exp[counts_obs_exp$CNV == CNV, ]
  groups <- temp$Age
  groups <- factor(groups, levels=c("UC", "EM", "MM"), ordered=TRUE)
  
  print(jonckheere.test(g=groups, x=temp$Percentage, alternative="decreasing"))$p.value
}

aggregate(Percentage ~ CNV+Age, counts_obs_exp, median)



#Are targets more likely to be DE than regulators? By age, regulator type

##Determine the downstream targets of each regulator
regulators <- unique(network$Gene1)
targets <- unique(network$Gene2)

regulator_targets <- lapply(regulators, function(reg){
  local_edges <- network[network$Gene1 == reg,]
  return(unique(local_edges$Gene2))
})

names(regulator_targets) <- regulators

load("regulator_classification.Rdata")

DE_reg_tar <- data.frame(p_values_df,
                         Class = ifelse(p_values_df$Gene %in% regulators, "Regulator", "Target"))

expected_reg_tar <- aggregate(Gene ~ Tumour+CNV+Class, DE_reg_tar, length)
expected_reg_tar$Label <- paste(expected_reg_tar$Tumour, expected_reg_tar$CNV, expected_reg_tar$Class, sep="_")

DE_reg_tar_sig <- DE_reg_tar[DE_reg_tar$adj_p < 0.05,]

significant_reg_tar <- aggregate(Gene ~ Tumour+CNV+Class, DE_reg_tar_sig, length)
significant_reg_tar$Label <- paste(significant_reg_tar$Tumour, significant_reg_tar$CNV, significant_reg_tar$Class, sep="_")

all <- expected_reg_tar$Label

per_sig_reg_tar <- data.frame(Label=all, 
                              expected_reg_tar[match(all, expected_reg_tar$Label),],
                              Significant=significant_reg_tar[match(all, significant_reg_tar$Label), "Gene"])


per_sig_reg_tar$Per_sig <- (per_sig_reg_tar$Significant/per_sig_reg_tar$Gene)*100

per_sig_reg_tar$Per_sig[is.na(per_sig_reg_tar$Per_sig)] <- 0

per_sig_reg_tar_cast <- vector()
tar_greater_all <- vector()
for(tumour in tumours){
  for(CNV in c("Amp", "Del")){
    temp <- per_sig_reg_tar[per_sig_reg_tar$Tumour == tumour,]
    temp <- temp[temp$CNV == CNV,]
    
    tar_greater <- temp$Per_sig[temp$Class == "Target"] > temp$Per_sig[temp$Class == "Regulator"]
    tar_greater_all <- rbind(tar_greater_all, c(tumour, CNV, tar_greater))
    per_sig_reg_tar_cast <- rbind(per_sig_reg_tar_cast,
                                  c(tumour, CNV, temp$Per_sig[temp$Class == "Regulator"], temp$Per_sig[temp$Class == "Target"]))
  }
}


colnames(tar_greater_all) <- c("Tumour", "CNV", "Target_greater")
tar_greater_all <- as.data.frame(tar_greater_all)


colnames(per_sig_reg_tar_cast) <- c("Tumour", "CNV", "Per_reg", "Per_tar")
per_sig_reg_tar_cast <- as.data.frame(per_sig_reg_tar_cast)
per_sig_reg_tar_cast$Per_reg <- as.numeric(as.character(per_sig_reg_tar_cast$Per_reg))
per_sig_reg_tar_cast$Per_tar <- as.numeric(as.character(per_sig_reg_tar_cast$Per_tar))

per_sig_reg_tar_cast$Diff <- per_sig_reg_tar_cast$Per_tar-per_sig_reg_tar_cast$Per_reg

per_sig_reg_tar_cast$Tumour <- factor(per_sig_reg_tar_cast$Tumour, levels=tumours)

pdf("Supp4_Diff_expressed_targets_regulators.pdf",
    height=4, width=8)
g <- ggplot(per_sig_reg_tar_cast, aes(x=Tumour, y=Diff))+
  geom_bar(aes(fill=Tumour), stat='identity')+
  geom_hline(yintercept=0)+
  ylab("Difference in the percentage of differentially expressed\ntargets and regulators")+
  facet_grid(.~CNV)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
print(g)
dev.off()


##Percentage of targets DE based on their regulator class
per_tar_de <- vector()
for(reg in regulators){
  reg_class <- regulator_classification[regulator_classification$Regulator == reg, "Regulator_class"]
  
  targets <- regulator_targets[[reg]]
  
  tested <- p_values_df[p_values_df$Gene %in% targets,]
  
  for(tumour in unique(tested$Tumour)){
    for(CNV in unique(tested$CNV)){
      local_tested <- tested[tested$Tumour == tumour,]
      local_tested <- local_tested[local_tested$CNV == CNV,]
      if(nrow(local_tested) > 0){
        #per_sig <- (sum(local_tested$adj_p < 0.05)/length(targets))*100
        per_sig <- (sum(local_tested$adj_p < 0.05)/nrow(local_tested))*100
        
        per_tar_de <- rbind(per_tar_de, 
                            c(reg, tumour, CNV, per_sig))
      }
    }
  }  
  print(reg)
}

per_tar_de <- as.data.frame(per_tar_de)
colnames(per_tar_de) <- c("Regulator", "Tumour", "CNV", "Per_de")
per_tar_de$Per_de <- as.numeric(as.character(per_tar_de$Per_de))

per_tar_de$Regulator_class <- regulator_classification[match(per_tar_de$Regulator, regulator_classification$Regulator), "Regulator_class"]

per_tar_de_relevant <- per_tar_de[per_tar_de$Regulator_class %in% c("UC_downstream", "Mixed_downstream", "EM_downstream"),]

per_tar_de_relevant <- per_tar_de_relevant[per_tar_de_relevant$Per_de != 0,]

per_tar_de_relevant <- per_tar_de_relevant[!is.na(per_tar_de_relevant$Per_de),]

##Average across tumours

per_tar_de_relevant_median <- aggregate(Per_de ~ Regulator+CNV+Regulator_class, per_tar_de_relevant, median)

per_tar_de_relevant_median$Regulator_class <- factor(per_tar_de_relevant_median$Regulator_class,
                                                     levels=c("UC_downstream", "Mixed_downstream", "EM_downstream"))

#Average across regulators

per_tar_de_relevant_median_reg <- aggregate(Per_de ~ Tumour+CNV+Regulator_class, per_tar_de_relevant, median)

per_tar_de_relevant_median_reg$Regulator_class <- factor(per_tar_de_relevant_median_reg$Regulator_class,
                                                     levels=c("UC_downstream", "Mixed_downstream", "EM_downstream"))

per_tar_de_relevant_median_reg$Regulator_class <- factor(per_tar_de_relevant_median_reg$Regulator_class,
                                                         levels=c("EM_downstream", "Mixed_downstream", "UC_downstream"))

per_tar_de_relevant_median_reg$Regulator_class <- factor(per_tar_de_relevant_median_reg$Regulator_class,
                                                         levels=c("UC_downstream", "Mixed_downstream", "EM_downstream"))

pdf("Figure4_DE_targets_by_regulator.pdf",
    height=3.5, width=5.5)
g <- ggplot(per_tar_de_relevant_median_reg, aes(x=CNV, y =Per_de))+
  geom_boxplot(aes(fill=Regulator_class, alpha=CNV))+
  scale_alpha_manual(values=c(1,0.25))+
  scale_fill_jco()+
  facet_grid(.~Regulator_class)+
  ylab("Percentage of differentially expressed targets\nwith a CNV")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), strip.background = element_blank())
print(g)
dev.off()


per_tar_de_relevant_median_UC <- per_tar_de_relevant_median[per_tar_de_relevant_median$Regulator_class == "UC_downstream",]
wilcox.test(per_tar_de_relevant_median_UC[per_tar_de_relevant_median_UC$CNV == "Amp","Per_de"],
            per_tar_de_relevant_median_UC[per_tar_de_relevant_median_UC$CNV == "Del","Per_de"], alternative="greater")

median(per_tar_de_relevant_median_UC[per_tar_de_relevant_median_UC$CNV == "Amp","Per_de"])
median(per_tar_de_relevant_median_UC[per_tar_de_relevant_median_UC$CNV == "Del","Per_de"])

per_tar_de_relevant_median_EM <- per_tar_de_relevant_median[per_tar_de_relevant_median$Regulator_class == "EM_downstream",]
wilcox.test(per_tar_de_relevant_median_EM[per_tar_de_relevant_median_EM$CNV == "Del","Per_de"],
            per_tar_de_relevant_median_EM[per_tar_de_relevant_median_EM$CNV == "Amp","Per_de"], alternative="greater")

per_tar_de_relevant_median_mixed <- per_tar_de_relevant_median[per_tar_de_relevant_median$Regulator_class == "Mixed_downstream",]
wilcox.test(per_tar_de_relevant_median_mixed[per_tar_de_relevant_median_mixed$CNV == "Del","Per_de"],
            per_tar_de_relevant_median_mixed[per_tar_de_relevant_median_mixed$CNV == "Amp","Per_de"])
