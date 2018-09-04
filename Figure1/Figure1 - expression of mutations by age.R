###Gene expression of genes with mutations
#Does the expression of genes with miss/LoF depend on the ages of the genes?

##Will do with normal/tumour pairs first, but likely to be too few patients

library(ggplot2)
library(reshape2)
library(readr)
library(limma)
library(edgeR)

source('herlper_functions.R')
load("gene_mut_properties.Rdata")
load("patients_with_mut_info.Rdata")
load("variants_LoF.Rdata")
load("variants_missense.Rdata")

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
expression_tumour <- list()
for(tumour in tumours){
  for(sample_type in c("tumour")){
    exp <- read.delim()  ##Normalized tumour expression
    colnames(exp)[2:ncol(exp)] <- substr(colnames(exp)[2:ncol(exp)],9,12)
    
    samples_low_correlation <- ##select samples with low correlation from TCGA
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
    
    patients_to_keep <- patients_with_mut_info[[tumour]]
    patients_tumour <- colnames(temp_exp)
    
    patients_to_keep <- patients_to_keep[patients_to_keep %in% patients_tumour]
    
    expression_tumour[[tumour]] <- temp_exp[,colnames(temp_exp) %in% patients_to_keep]
    print(dim(expression_tumour[[tumour]]))
    
  }
}

mutations_df <- load_mutations()

mutations_df <- mutations_df[mutations_df$Number_mutations >= 3,]
mutations_df <- mutations_df[mutations_df$Syn_ratio > 1,]


mut_curated <- list()
for(tumour in tumours){
  local_mutations_df <- mutations_df[mutations_df$Tumour == tumour,]
  mut_curated[[tumour]] <- local_mutations_df
}


p_lof_tumour <- list()
p_miss_tumour <- list()

for(tumour in tumours){
  
  expression_tumour_melt <- melt(as.matrix(expression_tumour[[tumour]]))
  colnames(expression_tumour_melt) <- c("Gene", "Patient", "Expression_tumour")
  
  local_LoF_variants <- variants_LoF[[tumour]][,c("Hugo_Symbol", "Entrez_Gene_Id", "Tumor_Sample_Barcode", "Consequence")]
  local_LoF_variants$Tumor_Sample_Barcode <- substr(local_LoF_variants$Tumor_Sample_Barcode, 9, 12)
  
  local_LoF_variants$Expression_tumour <- expression_tumour_melt[match(paste(as.character(local_LoF_variants$Hugo_Symbol), local_LoF_variants$Tumor_Sample_Barcode),
                                                                   paste(as.character(expression_tumour_melt$Gene), as.character(expression_tumour_melt$Patient))), "Expression_tumour"]
  local_LoF_variants <- local_LoF_variants[!is.na(local_LoF_variants$Expression_tumour),]
  
  local_miss_variants <- variants_missense[[tumour]][,c("Hugo_Symbol", "Entrez_Gene_Id", "Tumor_Sample_Barcode", "Consequence")]
  local_miss_variants$Tumor_Sample_Barcode <- substr(local_miss_variants$Tumor_Sample_Barcode, 9, 12)
  local_miss_variants$Expression_tumour <- expression_tumour_melt[match(paste(as.character(local_miss_variants$Hugo_Symbol), local_miss_variants$Tumor_Sample_Barcode),
                                                                    paste(as.character(expression_tumour_melt$Gene), as.character(expression_tumour_melt$Patient))), "Expression_tumour"]
  local_miss_variants <- local_miss_variants[!is.na(local_miss_variants$Expression_tumour),]
  
  
  patients_genes_mut <- c(paste(as.character(local_LoF_variants$Hugo_Symbol), local_LoF_variants$Tumor_Sample_Barcode),
                          paste(as.character(local_miss_variants$Hugo_Symbol), local_miss_variants$Tumor_Sample_Barcode))
  
  patients_no_mut <- expression_tumour_melt[!(paste(as.character(expression_tumour_melt$Gene), 
                                                 as.character(expression_tumour_melt$Patient)) %in% patients_genes_mut),]
  
  p_lof_tumour[[tumour]] <- vector()
  for(gene in unique(as.character(local_LoF_variants$Hugo_Symbol))){
    local_mut <- local_LoF_variants[local_LoF_variants$Hugo_Symbol == gene,]
    if(nrow(local_mut) >= 2){
      #print(gene)
      exp_mut <- local_mut$Expression_tumour
      local_not <- patients_no_mut[patients_no_mut$Gene == gene,]
      exp_not <- local_not$Expression_tumour
      p1 <- wilcox.test(exp_mut, exp_not)$p.val
      p2 <- wilcox.test(exp_mut, exp_not, alternative="greater")$p.val
      p3 <- wilcox.test(exp_mut, exp_not, alternative="less")$p.val
      
      p_lof_tumour[[tumour]] <- rbind(p_lof_tumour[[tumour]], c(gene, p1, p2, p3))
      
    }
  }
  p_lof_tumour[[tumour]] <- as.data.frame(p_lof_tumour[[tumour]])
  if(nrow(p_lof_tumour[[tumour]]) > 1){
    colnames(p_lof_tumour[[tumour]]) <- c("Gene", "p_two_sided", "p_greater", "p_less")
  }
  
  
  p_miss_tumour[[tumour]] <- vector()
  for(gene in unique(as.character(local_miss_variants$Hugo_Symbol))){
    local_mut <- local_miss_variants[local_miss_variants$Hugo_Symbol == gene,]
    if(nrow(local_mut) >= 3){
      exp_mut <- local_mut$Expression_tumour
      local_not <- patients_no_mut[patients_no_mut$Gene == gene,]
      exp_not <- local_not$Expression_tumour
      p1 <- wilcox.test(exp_mut, exp_not)$p.val
      p2 <- wilcox.test(exp_mut, exp_not, alternative="greater")$p.val
      p3 <- wilcox.test(exp_mut, exp_not, alternative="less")$p.val
      
      p_miss_tumour[[tumour]] <- rbind(p_miss_tumour[[tumour]], c(gene, p1, p2, p3))
      
    }
  }
  
  p_miss_tumour[[tumour]] <- as.data.frame(p_miss_tumour[[tumour]])
  colnames(p_miss_tumour[[tumour]]) <- c("Gene", "p_two_sided", "p_greater", "p_less")
  
  print(tumour)
}


DE_genes <- vector()
for(tumour in tumours){
    p_lof_tumour[[tumour]]$p_two_sided <- as.numeric(as.character(p_lof_tumour[[tumour]]$p_two_sided))
    p_lof_tumour[[tumour]]$p_greater <- as.numeric(as.character(p_lof_tumour[[tumour]]$p_greater))
    p_lof_tumour[[tumour]]$p_less <- as.numeric(as.character(p_lof_tumour[[tumour]]$p_less))
    
    p_lof_tumour[[tumour]]$adj_two_sided <- p.adjust(p_lof_tumour[[tumour]]$p_two_sided, method="BH")
    p_lof_tumour[[tumour]]$adj_greater <- p.adjust(p_lof_tumour[[tumour]]$p_greater, method="BH")
    p_lof_tumour[[tumour]]$adj_less <- p.adjust(p_lof_tumour[[tumour]]$p_less, method="BH")
  
  p_miss_tumour[[tumour]]$p_two_sided <- as.numeric(as.character(p_miss_tumour[[tumour]]$p_two_sided))
  p_miss_tumour[[tumour]]$p_greater <- as.numeric(as.character(p_miss_tumour[[tumour]]$p_greater))
  p_miss_tumour[[tumour]]$p_less <- as.numeric(as.character(p_miss_tumour[[tumour]]$p_less))
  
  p_miss_tumour[[tumour]]$adj_two_sided <- p.adjust(p_miss_tumour[[tumour]]$p_two_sided, method="BH")
  p_miss_tumour[[tumour]]$adj_greater <- p.adjust(p_miss_tumour[[tumour]]$p_greater, method="BH")
  p_miss_tumour[[tumour]]$adj_less <- p.adjust(p_miss_tumour[[tumour]]$p_less, method="BH")
  
  DE_genes <- rbind(DE_genes,
                      cbind(Tumour = tumour, Mut="LoF", p_lof_tumour[[tumour]]),
                      cbind(Tumour = tumour, Mut="Miss", p_miss_tumour[[tumour]]))
}

DE_genes$Age <- genes_phy_categorical[match(DE_genes$Gene, genes_phy_categorical$GeneID), "Phylostrata"]
DE_genes <- DE_genes[!is.na(DE_genes$Age),]
DE_genes$Age <- factor(DE_genes$Age, levels=c("UC", "EM", "MM"))

##Using nominal
DE_genes_greater <- DE_genes[DE_genes$p_greater < 0.05,]
DE_genes_less <- DE_genes[DE_genes$p_less < 0.05,]

observed_greater <- aggregate(p_greater ~ Tumour+Mut+Age, DE_genes_greater, length)
observed_less <- aggregate(p_less ~ Tumour+Mut+Age, DE_genes_less, length)

observed_greater$Name <- paste(observed_greater$Tumour, observed_greater$Mut, observed_greater$Age)
observed_less$Name <- paste(observed_less$Tumour, observed_less$Mut, observed_less$Age)

observed <- data.frame(Tumour = rep(tumours, each=3), Mut = "LoF", Age = rep(c("UC", "EM", "MM"), 7))
observed <- data.frame(rbind(observed,
                             data.frame(Tumour = rep(tumours, each=3), Mut = "Miss", Age = rep(c("UC", "EM", "MM"), 7))))
observed$Name <- paste(observed$Tumour, observed$Mut, observed$Age)
observed$N_greater <- observed_greater[match(observed$Name, observed_greater$Name), "p_greater"]
observed$N_less <- observed_less[match(observed$Name, observed_less$Name), "p_less"]
observed[is.na(observed)] <- 0

observed$Age <- factor(observed$Age, levels=c("UC", "EM", "MM"))
observed$Tumour <- factor(observed$Tumour, levels=tumours)

observed$N_total <- observed$N_greater+observed$N_less

observed$Ratio_greater_less <- observed$N_greater/observed$N_less

observed_no_inf <- observed
observed_no_inf$Ratio_greater_less[is.infinite(observed_no_inf$Ratio_greater_less)] <- NA

observed_mean <- aggregate(Ratio_greater_less ~ Mut+Age, observed_no_inf, mean, na.omit=TRUE)

observed_non_nan <- observed[!is.nan(observed$Ratio_greater_less),]


observed_non_nan$Group <- paste(observed_non_nan$Age, observed_non_nan$Mut, sep="_")

observed_non_nan_non_inf <- observed_non_nan[is.finite(observed_non_nan$Ratio_greater_less),]

observed_non_nan_non_inf$Mut <- factor(observed_non_nan_non_inf$Mut, levels=c("Miss", "LoF"))

observed_non_nan_non_inf$Group <- factor(observed_non_nan_non_inf$Group,
                                         levels=c("UC_Miss", "UC_LoF","EM_Miss","EM_LoF","MM_Miss","MM_LoF"))

observed_mean$Mut <- factor(observed_mean$Mut, levels=c("Miss", "LoF"))

pdf("Figure1_DE_of_mutated_genes.pdf",
    width=5, height=3)
g <- ggplot(observed_non_nan_non_inf, aes(x=Mut, y =Ratio_greater_less, group=Mut))+
  geom_point(aes(colour=Tumour, shape=Mut), 
             position = position_dodge(width = 0.75), size=4)+
  geom_hline(yintercept=1, colour="grey")+
  geom_line(aes(group=Group), 
            position = position_dodge(width = 0.75))+
  geom_point(aes(colour=Tumour, shape=Mut), 
             position = position_dodge(width = 0.75), size=4)+
  #geom_point(data=observed_mean, pch=4, size=3.5, position= position_dodge(width = 0.75))+
  ylab("Ratio of up and downregulated genes")+
  facet_grid(.~Age, scale='free_x')+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank())
print(g)
dev.off()

