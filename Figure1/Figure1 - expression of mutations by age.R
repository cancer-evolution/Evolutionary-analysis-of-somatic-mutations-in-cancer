###Gene expression of genes with mutations
#Does the expression of genes with miss/LoF depend on the ages of the genes?

##Will do with normal/tumour pairs first, but likely to be too few patients

library(ggplot2)
library(reshape2)
library(readr)
library(limma)
library(edgeR)

source("helper_functions.R")
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

tumours <- c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "ESCA", "GBM", 
             "HNSC", "KICH", "KIRC", "KIRP", "LGG",  "LIHC", "LUAD", "LUSC", 
             "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", 
             "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM")  

##Load expression data
expression_tumour <- list()
for(tumour in tumours){
  for(sample_type in c("tumour")){
    path <- paste(tumour, "/gdac.broadinstitute.org_", tumour, 
                  ".Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2016012800.0.0/", 
                  tumour, 
                  ".rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt",
                  sep="")
    
    exp <- read.delim(path) ##Load normalized gene expression data from TCGA
    exp <- exp[-1,]
    
    exp_tumour <- exp[,substr(colnames(exp), 14, 16) %in% c("01A", "01B", "01C", "01D")]
    
    colnames(exp_tumour) <- substr(colnames(exp_tumour),9,12)
    
    all_genes <- as.character(exp[,1])
    
    all_genes <- unlist(strsplit(all_genes, "\\|"))[seq(1, length(all_genes)*2, 2)]
    
    genes_remove <- c("?", "SLC35E2")
    genes_remove_index <- which(all_genes %in% genes_remove)
    all_genes <- all_genes[-genes_remove_index]
    
    exp_tumour <- exp_tumour[-genes_remove_index,]
    
    if(!is.null(nrow(exp_tumour))){
      rownames(exp_tumour) <- all_genes
    }
    
    patients_to_keep <- patients_with_mut_info[[tumour]]
    patients_tumour <- colnames(exp_tumour)
    
    patients_to_keep <- patients_to_keep[patients_to_keep %in% patients_tumour]
    
    expression_tumour[[tumour]] <- exp_tumour[,colnames(exp_tumour) %in% patients_to_keep]
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

for(tumour in tumours[6:length(tumours)]){
  
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
    if(nrow(local_mut) >= 3){###changed to 3 from 2 in the paper
      #print(gene)
      exp_mut <- local_mut$Expression_tumour
      local_not <- patients_no_mut[patients_no_mut$Gene == gene,]
      exp_not <- local_not$Expression_tumour
      p1 <- wilcox.test(as.numeric(as.character(exp_mut)), as.numeric(as.character(exp_not)))$p.val
      p2 <- wilcox.test(as.numeric(as.character(exp_mut)), as.numeric(as.character(exp_not)), alternative="greater")$p.val
      p3 <- wilcox.test(as.numeric(as.character(exp_mut)), as.numeric(as.character(exp_not)), alternative="less")$p.val
      
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
      p1 <- wilcox.test(as.numeric(as.character(exp_mut)), as.numeric(as.character(exp_not)))$p.val
      p2 <- wilcox.test(as.numeric(as.character(exp_mut)), as.numeric(as.character(exp_not)), alternative="greater")$p.val
      p3 <- wilcox.test(as.numeric(as.character(exp_mut)), as.numeric(as.character(exp_not)), alternative="less")$p.val
      
      p_miss_tumour[[tumour]] <- rbind(p_miss_tumour[[tumour]], c(gene, p1, p2, p3))
      
    }
  }
  
  p_miss_tumour[[tumour]] <- as.data.frame(p_miss_tumour[[tumour]])
  if(nrow(p_miss_tumour[[tumour]])>0){
    colnames(p_miss_tumour[[tumour]]) <- c("Gene", "p_two_sided", "p_greater", "p_less")
  }
  
  print(tumour)
  save(p_miss_tumour, file="p_miss_tumour.Rdata")
  save(p_lof_tumour, file="p_lof_tumour.Rdata")
}


load("p_miss_tumour.Rdata")
load("p_lof_tumour.Rdata")

sum(sapply(p_miss_tumour, nrow))
sum(sapply(p_lof_tumour, nrow))

DE_genes <- vector()
for(tumour in tumours){
  if(nrow(p_lof_tumour[[tumour]]) != 0){
    colnames(p_lof_tumour[[tumour]]) <- c("Gene", "p_two_sided", "p_greater", "p_less")
    p_lof_tumour[[tumour]]$p_two_sided <- as.numeric(as.character(p_lof_tumour[[tumour]]$p_two_sided))
    p_lof_tumour[[tumour]]$p_greater <- as.numeric(as.character(p_lof_tumour[[tumour]]$p_greater))
    p_lof_tumour[[tumour]]$p_less <- as.numeric(as.character(p_lof_tumour[[tumour]]$p_less))
    
    p_lof_tumour[[tumour]]$adj_two_sided <- p.adjust(p_lof_tumour[[tumour]]$p_two_sided, method="BH")
    p_lof_tumour[[tumour]]$adj_greater <- p.adjust(p_lof_tumour[[tumour]]$p_greater, method="BH")
    p_lof_tumour[[tumour]]$adj_less <- p.adjust(p_lof_tumour[[tumour]]$p_less, method="BH")
    
    DE_genes <- rbind(DE_genes,
                      cbind(Tumour = tumour, Mut="LoF", p_lof_tumour[[tumour]]))
  }
  if(nrow(p_miss_tumour[[tumour]]) != 0){
    colnames(p_miss_tumour[[tumour]]) <- c("Gene", "p_two_sided", "p_greater", "p_less")
    p_miss_tumour[[tumour]]$p_two_sided <- as.numeric(as.character(p_miss_tumour[[tumour]]$p_two_sided))
    p_miss_tumour[[tumour]]$p_greater <- as.numeric(as.character(p_miss_tumour[[tumour]]$p_greater))
    p_miss_tumour[[tumour]]$p_less <- as.numeric(as.character(p_miss_tumour[[tumour]]$p_less))
    
    p_miss_tumour[[tumour]]$adj_two_sided <- p.adjust(p_miss_tumour[[tumour]]$p_two_sided, method="BH")
    p_miss_tumour[[tumour]]$adj_greater <- p.adjust(p_miss_tumour[[tumour]]$p_greater, method="BH")
    p_miss_tumour[[tumour]]$adj_less <- p.adjust(p_miss_tumour[[tumour]]$p_less, method="BH")
    
    DE_genes <- rbind(DE_genes,
                      cbind(Tumour = tumour, Mut="Miss", p_miss_tumour[[tumour]]))
    
  }
  print(tumour)
}

DE_genes$Age <- genes_phy_categorical[match(DE_genes$Gene, genes_phy_categorical$GeneID), "Phylostrata"]
DE_genes <- DE_genes[!is.na(DE_genes$Age),]
DE_genes$Age <- factor(DE_genes$Age, levels=c("UC", "EM", "MM"))

##Using nominal
DE_genes_greater <- DE_genes[DE_genes$p_greater < 0.05,]
DE_genes_less <- DE_genes[DE_genes$p_less < 0.05,]


head(sort(table(DE_genes_greater[DE_genes_greater$Age == "EM","Gene"]), decreasing=TRUE), 10)
head(sort(table(DE_genes_less[DE_genes_less$Age == "EM","Gene"]), decreasing=TRUE), 10)

observed_greater <- aggregate(p_greater ~ Tumour+Mut+Age, DE_genes_greater, length)
observed_less <- aggregate(p_less ~ Tumour+Mut+Age, DE_genes_less, length)

observed_greater$Name <- paste(observed_greater$Tumour, observed_greater$Mut, observed_greater$Age)
observed_less$Name <- paste(observed_less$Tumour, observed_less$Mut, observed_less$Age)

observed <- data.frame(Tumour = rep(tumours, each=3), Mut = "LoF", Age = rep(c("UC", "EM", "MM"), 30))
observed <- data.frame(rbind(observed,
                             data.frame(Tumour = rep(tumours, each=3), Mut = "Miss", Age = rep(c("UC", "EM", "MM"), 30))))
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


ggplot(observed_non_nan_non_inf, aes(x=Mut, y =Ratio_greater_less))+
  geom_boxplot()+
  facet_grid(.~Age)
  

library(dplyr)
counts <- vector()
for(mut in c("Miss", "LoF")){
  for(age in c("UC", "EM", "MM")){
    observed_non_nan_non_inf %>%
        filter(Age == age & Mut == mut) -> temp
        counts <- rbind(counts, c(mut, age, sum(temp$Ratio_greater_less >= 1), sum(temp$Ratio_greater_less < 1)))
  }
}

counts <- as.data.frame(counts)
colnames(counts) <- c("Mut", "Age", "Number of tumours where ratio more than 1", "Number of tumours where ratio less than 1")

counts$Proportion <- as.numeric(as.character(counts$`Number of tumours where ratio more than 1`))/as.numeric(as.character(counts$`Number of tumours where ratio less than 1`))


##Just calculate more likely to be DE
proportion_DE <- vector()
for(tumour in tumours){
  for(alt in c("Miss", "LoF")){
    temp <- DE_genes[DE_genes$Tumour == tumour,]
    temp <- temp[temp$Mut == alt,]  
    
    proportion_DE <- rbind(proportion_DE, 
                           c(tumour, alt, as.data.frame(table(temp$Age[temp$p_two_sided < 0.05])/table(temp$Age))$Freq))
  }
  print(tumour)
}
proportion_DE <- as.data.frame(proportion_DE)
colnames(proportion_DE) <- c("Tumour", "Alt", "Prop_UC", "Prop_EM", "Prop_MM")
proportion_DE$Prop_UC <- as.numeric(as.character(proportion_DE$Prop_UC))
proportion_DE$Prop_EM <- as.numeric(as.character(proportion_DE$Prop_EM))
proportion_DE$Prop_MM <- as.numeric(as.character(proportion_DE$Prop_MM))

proportion_DE <- melt(proportion_DE)

pdf("Figure1_Proportion_DE_point_mutations.pdf",
    width=4, height=2)
g <- ggplot(proportion_DE, aes(x=variable, y=value))+
  geom_boxplot(outlier.shape=NA)+
  #geom_jitter(data=subset(proportion_DE, Alt=="Miss"), 
  #            aes(shape=Alt, colour=Tumour), width = 0.25, size=2)+
  geom_jitter(data=subset(proportion_DE, Alt=="Miss"), 
              aes(x=as.numeric(variable)-.15, shape=Alt, colour=Tumour), width = 0.075, size=1)+
  geom_jitter(data=subset(proportion_DE, Alt=="LoF"), 
             aes(x=as.numeric(variable)+.15, shape=Alt, colour=Tumour), width = 0.075, size=1)+
  theme_bw()
print(g)
dev.off()

pdf("Figure1_Proportion_DE_point_mutations_legend.pdf",
    width=4, height=5)
g <- ggplot(proportion_DE, aes(x=variable, y=value))+
  geom_boxplot(outlier.shape=NA)+
  #geom_jitter(data=subset(proportion_DE, Alt=="Miss"), 
  #            aes(shape=Alt, colour=Tumour), width = 0.25, size=2)+
  geom_jitter(data=subset(proportion_DE, Alt=="Miss"), 
              aes(x=as.numeric(variable)-.15, shape=Alt, colour=Tumour), width = 0.075, size=1)+
  geom_jitter(data=subset(proportion_DE, Alt=="LoF"), 
              aes(x=as.numeric(variable)+.15, shape=Alt, colour=Tumour), width = 0.075, size=1)+
  theme_bw()
print(g)
dev.off()


wilcox.test(proportion_DE$value[proportion_DE$variable == "Prop_EM"],
            proportion_DE$value[proportion_DE$variable == "Prop_UC"], alternative="greater")

wilcox.test(proportion_DE$value[proportion_DE$variable == "Prop_EM"],
            proportion_DE$value[proportion_DE$variable == "Prop_MM"], alternative="greater")


proportion_DE_miss <- proportion_DE[proportion_DE$Alt == "Miss",]
proportion_DE_lof <- proportion_DE[proportion_DE$Alt == "LoF",]


wilcox.test(proportion_DE_miss$value[proportion_DE_miss$variable == "Prop_EM"],
            proportion_DE_miss$value[proportion_DE_miss$variable == "Prop_UC"], alternative="greater")

wilcox.test(proportion_DE_miss$value[proportion_DE_miss$variable == "Prop_EM"],
            proportion_DE_miss$value[proportion_DE_miss$variable == "Prop_MM"], alternative="greater")

wilcox.test(proportion_DE_lof$value[proportion_DE_lof$variable == "Prop_EM"],
            proportion_DE_lof$value[proportion_DE_lof$variable == "Prop_UC"], alternative="greater")

wilcox.test(proportion_DE_lof$value[proportion_DE_lof$variable == "Prop_EM"],
            proportion_DE_lof$value[proportion_DE_lof$variable == "Prop_MM"], alternative="greater")







ggplot(proportion_DE, aes(x=variable, y=value))+
  geom_boxplot(aes(fill=Alt), position='dodge')+
  theme_bw()



ggplot(proportion_DE, aes(x=variable, y=value))+
  geom_boxplot()+
  geom_jitter(aes(shape=Alt, colour=Tumour), width = 0.25)+
  theme_bw()+
  facet_grid(.~Alt)


ggplot(proportion_DE, aes(x=Tumour, y=value))+
  geom_bar(stat='identity', aes(fill=Alt), position='dodge')+
  facet_grid(.~variable)


##Waterfall plot

observed_plot <- observed
observed_plot$N_less <- -observed_plot$N_less
observed_plot_melt <- melt(observed_plot[,c("Tumour", "Mut", "Age", "N_greater", "N_less")])


ggplot(observed_plot_melt, aes(x=Tumour, y=value))+
  geom_bar(stat='identity')+
  facet_grid(Mut~Age)+
  coord_flip()

for(age in c("UC", "EM", "MM")){
  for(alt in c("Miss", "LoF")){
    local_temp <- observed_plot[intersect(which(observed_plot$Age == age),
                                               which(observed_plot$Mut == alt)),]
    local_temp <- local_temp[local_temp$N_total != 0,]
    local_order <- local_temp$Tumour[order(local_temp$N_total)]
    local_temp_melt <- melt(local_temp[,c("Tumour", "Mut", "Age", "N_greater", "N_less")])
    local_temp_melt$Tumour <- factor(local_temp_melt$Tumour, levels=rev(as.character(local_order)))
    
    g <- ggplot(local_temp_melt, aes(x=Tumour, y=value))+
      geom_bar(stat='identity', aes(fill=variable))+
      facet_grid(Mut~Age)+
      coord_flip()
    print(g)
  }
}

##Add not only up and down, but what is DE and what is not
table_counts <- vector()
for(tumour in tumours){
  for(mut in c("Miss", "LoF")){
    local_temp <- DE_genes[intersect(which(DE_genes$Tumour == tumour),
                                     which(DE_genes$Mut == mut)),]
    t1 <- table(local_temp$Age[which(local_temp$p_greater < 0.05)])
    t2 <- table(local_temp$Age[which(local_temp$p_less < 0.05)])
    t3 <- table(local_temp$Age[which(local_temp$p_less > 0.05 & local_temp$p_greater > 0.05)])
    
    table_counts <- rbind(table_counts, 
    rbind(cbind(tumour, mut, "Greater", t1, c("UC", "EM", "MM")),
          cbind(tumour, mut, "Less", t2, c("UC", "EM", "MM")),
          cbind(tumour, mut, "Not_DE", t3, c("UC", "EM", "MM"))))
  
  }
}
table_counts <- as.data.frame(table_counts)
colnames(table_counts) <- c("Tumour", "Mut", "Direction", "Count", "Age")

table_counts$Count <- as.numeric(as.character(table_counts$Count))

ggplot(table_counts, aes(x=Tumour, y=Count))+
  geom_bar(stat='identity', aes(fill=Direction))+
  facet_grid(Mut ~ Age)


##Just calculate more likely to be DE (but breakdown into up and down)
proportion_DE <- vector()
for(tumour in tumours){
  for(alt in c("Miss", "LoF")){
    temp <- DE_genes[DE_genes$Tumour == tumour,]
    temp <- temp[temp$Mut == alt,]  
    
    proportion_DE <- rbind(proportion_DE, 
                           c(tumour, alt, as.data.frame(table(temp$Age[temp$p_two_sided < 0.05])/table(temp$Age))$Freq,
                             as.data.frame(table(temp$Age[temp$p_greater < 0.05])/table(temp$Age))$Freq,
                             as.data.frame(table(temp$Age[temp$p_less < 0.05])/table(temp$Age))$Freq,
                             as.data.frame(table(temp$Age[temp$p_greater < 0.05 | temp$p_less < 0.05])/table(temp$Age))$Freq))
  }
  print(tumour)
}
proportion_DE <- as.data.frame(proportion_DE)
colnames(proportion_DE) <- c("Tumour", "Alt", "Prop_UC_DE", "Prop_EM_DE", "Prop_MM_DE",
                             "Prop_UC_up", "Prop_EM_up", "Prop_MM_up",
                             "Prop_UC_down", "Prop_EM_down", "Prop_MM_down",
                             "Prop_UC_both", "Prop_EM_both", "Prop_MM_both")
proportion_DE[,3:ncol(proportion_DE)] <- apply(proportion_DE[,3:ncol(proportion_DE)], 2, function(x){as.numeric(as.character(x))})

proportion_DE <- melt(proportion_DE)

ggplot(proportion_DE, aes(x=variable, y=value))+
  geom_boxplot()+
  geom_jitter(aes(shape=Alt, colour=Tumour), width = 0.25)+
  theme_bw()

proportion_DE_up_down <- proportion_DE[grep(paste(c("up", "down"), collapse="|"), 
                                            proportion_DE$variable),]

proportion_DE_up_down$Age <- substr(proportion_DE_up_down$variable, 6,7)
proportion_DE_up_down$Direction <- substr(proportion_DE_up_down$variable, 9, 10)


ggplot(proportion_DE_up_down, aes(x=Tumour, y=value))+
  geom_bar(aes(fill=Direction), stat='identity')+
  facet_grid(Alt~Age)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
