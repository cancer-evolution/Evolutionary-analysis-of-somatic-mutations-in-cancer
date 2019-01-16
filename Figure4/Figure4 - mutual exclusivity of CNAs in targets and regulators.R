library(igraph)
library(readr)
library(ggplot2)
library(reshape2)
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
RN_network <- load_PathwayCommons_GRN()

regulators <- unique(RN_network$Gene1)

regulator_targets <- lapply(regulators, function(reg){
  local_edges <- RN_network[RN_network$Gene1 == reg,]
  return(unique(local_edges$Gene2))
})

names(regulator_targets) <- regulators


##Load CNVs per patient
load("CNVs_curated2.Rdata")

tumours <- c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "ESCA", "GBM", 
             "HNSC", "KICH", "KIRC", "KIRP", "LGG",  "LIHC", "LUAD", "LUSC", 
             "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", 
             "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM")

patient_CNVs <- vector()
for(tumour in tumours[1:3]){
  tumour_amp <- CNVs_curated[[tumour]]$amplifications
  tumour_del <- CNVs_curated[[tumour]]$deletions
  
  tumour_amp$Sample <- substr(tumour_amp$Sample, 9, 12)
  tumour_del$Sample <- substr(tumour_del$Sample, 9, 12)
  
  tumour_amp$CNV <- "Amp"
  tumour_del$CNV <- "Del"
  tumour_amp$Tumour <- tumour
  tumour_del$Tumour <- tumour
  tumour_amp$Segment_Mean <- NULL
  tumour_del$Segment_Mean <- NULL
  
  patient_CNVs <- rbind(patient_CNVs, tumour_amp, tumour_del)
  print(tumour)
}


me_p <- vector()
for(regulator in regulators){
  regulator_targets_CNVs <- vector()
  
  local_targets <- regulator_targets[[regulator]]
  
  targets_with_CNVs <- patient_CNVs[patient_CNVs$Gene %in% local_targets,]
  targets_with_CNVs$Target_age <- genes_phy[match(targets_with_CNVs$Gene, genes_phy$GeneID), "Age"]
  
  CNV_in_regulator <- patient_CNVs[patient_CNVs$Gene %in% regulator,]
  
  if(nrow(CNV_in_regulator) >= 3){
    all_patients_involved <- unique(c(targets_with_CNVs$Sample, CNV_in_regulator$Sample))
    
    for(patient in all_patients_involved){
      
      patient_targets_with_CNVs <- targets_with_CNVs[targets_with_CNVs$Sample == patient,]
      
      patient_targets_with_CNVs_UC <- patient_targets_with_CNVs[patient_targets_with_CNVs$Target_age == "UC",]
      patient_targets_with_CNVs_EM <- patient_targets_with_CNVs[patient_targets_with_CNVs$Target_age == "EM",]
      patient_targets_with_CNVs_MM <- patient_targets_with_CNVs[patient_targets_with_CNVs$Target_age == "MM",]
      
      patient_targets_with_CNVs_UC <- patient_targets_with_CNVs_UC[!is.na(patient_targets_with_CNVs_UC$Gene),]
      patient_targets_with_CNVs_EM <- patient_targets_with_CNVs_EM[!is.na(patient_targets_with_CNVs_EM$Gene),]
      patient_targets_with_CNVs_MM <- patient_targets_with_CNVs_MM[!is.na(patient_targets_with_CNVs_MM$Gene),]
      
      patient_regulator_with_CNV <- patient %in% CNV_in_regulator$Sample
      
      regulator_targets_CNVs <- rbind(regulator_targets_CNVs,
                                      c(patient, patient_regulator_with_CNV, nrow(patient_targets_with_CNVs)/length(local_targets),
                                        nrow(patient_targets_with_CNVs_UC), nrow(patient_targets_with_CNVs_EM), nrow(patient_targets_with_CNVs_MM)))
    }
    
    fraction_targets_when_regulator_CNV <- as.numeric(regulator_targets_CNVs[regulator_targets_CNVs[,2] == TRUE, 3])
    fraction_targets_when_regulator_not_CNV <- as.numeric(regulator_targets_CNVs[regulator_targets_CNVs[,2] == FALSE, 3])
    n_UC_targets_affected <- as.numeric(regulator_targets_CNVs[, 4])
    n_EM_targets_affected <- as.numeric(regulator_targets_CNVs[, 5])
    n_MM_targets_affected <- as.numeric(regulator_targets_CNVs[, 6])
    
    if((length(fraction_targets_when_regulator_not_CNV) >= 3) && (length(fraction_targets_when_regulator_CNV) >= 3)){
      p1 <- wilcox.test(fraction_targets_when_regulator_not_CNV, fraction_targets_when_regulator_CNV, alternative="greater")$p.val
      p2 <- wilcox.test(fraction_targets_when_regulator_CNV, fraction_targets_when_regulator_not_CNV, alternative="greater")$p.val
      me_p <- rbind(me_p, c(regulator, median(fraction_targets_when_regulator_CNV), median(fraction_targets_when_regulator_not_CNV), p1, p2,
                            median(n_UC_targets_affected), median(n_EM_targets_affected), median(n_MM_targets_affected)))
    }  
  }
  print(regulator)
}

#save(me_p, file="me_p.Rdata")
load("me_p.Rdata")

me_p <- as.data.frame(me_p)
colnames(me_p) <- c("Regulator", "Fraction_targets_regulator_CNV",
                    "Fraction_targets_regulator_not_CNV", "Preference_CNV_targets", "Preference_CNV_regulator", 
                    "Median_target_UC_genes_CNVs", "Median_target_EM_genes_CNVs", "Median_target_MM_genes_CNVs")
me_p$Fraction_targets_regulator_CNV <- as.numeric(as.character(me_p$Fraction_targets_regulator_CNV))
me_p$Fraction_targets_regulator_not_CNV <- as.numeric(as.character(me_p$Fraction_targets_regulator_not_CNV))
me_p$Preference_CNV_targets <- as.numeric(as.character(me_p$Preference_CNV_targets))
me_p$Preference_CNV_regulator <- as.numeric(as.character(me_p$Preference_CNV_regulator))

me_p$Preference_CNV_targets_adj <- p.adjust(me_p$Preference_CNV_targets, method="BH")
me_p$Preference_CNV_regulator_adj <- p.adjust(me_p$Preference_CNV_regulator, method="BH")

n_targets_per_regulator <- sapply(regulator_targets, length)
n_targets_per_regulator <- data.frame(Regulator=names(n_targets_per_regulator),
                                      N_targets = n_targets_per_regulator)

me_p$N_targets <- n_targets_per_regulator[match(me_p$Regulator, n_targets_per_regulator$Regulator), "N_targets"]

me_p <- me_p[me_p$N_targets >= 2,]
me_p <- me_p[!is.na(me_p$Regulator),]

me_p$Age <- genes_phy[match(me_p$Regulator, genes_phy$GeneID), "Age"]

me_p$Age <- factor(me_p$Age, levels=c("UC", "EM", "MM"))

me_p_fractions <- me_p[,c("Regulator", "Fraction_targets_regulator_CNV", "Fraction_targets_regulator_not_CNV")]
me_p_fractions <- melt(me_p_fractions, id="Regulator")

colnames(me_p_fractions) <- c("Regulator", "Regulator_status", "Fraction_targets_with_CNAs")

me_p_fractions$Regulator_status <- ifelse(me_p_fractions$Regulator_status == "Fraction_targets_regulator_CNV", "CNA", "CNN")

###Supplementary
pdf("Supp4_Fraction_targets_CNA_by_regulator_status.pdf",
    width=4, height=3)
g <- ggplot(me_p_fractions, aes(x=Regulator_status, y=Fraction_targets_with_CNAs))+
  geom_boxplot(aes(fill=Regulator_status))+
  ylab("Median fraction of targets with CNAs")+
  xlab("Regulator status")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())
print(g)
dev.off()

wilcox.test(me_p_fractions$Fraction_targets_with_CNAs[me_p_fractions$Regulator_status == "CNN"],
            me_p_fractions$Fraction_targets_with_CNAs[me_p_fractions$Regulator_status == "CNA"], alternative="greater")$p.value



load("regulator_classification.Rdata")

me_p_fractions$Quantiles1 <- regulator_classification[match(me_p_fractions$Regulator, regulator_classification$Regulator), "Regulator_class"]

me_p_fractions$Quantiles1 <- factor(me_p_fractions$Quantiles1, levels=c("UC_downstream", "Mixed_downstream", "EM_downstream"))


me_p_fractions_no_CNVs_in_regulators <- me_p_fractions[!is.na(me_p_fractions$Quantiles1),]
me_p_fractions_no_CNVs_in_regulators <- me_p_fractions_no_CNVs_in_regulators[me_p_fractions_no_CNVs_in_regulators$Regulator_status == "CNN",]


pdf("Figure4_Targets_CNA_by_class.pdf",
    width=5.25, height=4)
g <- ggplot(me_p_fractions_no_CNVs_in_regulators, aes(x=Quantiles1, y=Fraction_targets_with_CNAs))+
  geom_boxplot(aes(fill=Quantiles1))+
  scale_fill_jco()+
  xlab("Regulator class")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())
print(g)
dev.off()

aggregate(Fraction_targets_with_CNAs ~ Quantiles1, me_p_fractions_no_CNVs_in_regulators, median)

wilcox.test(me_p_fractions_no_CNVs_in_regulators[me_p_fractions_no_CNVs_in_regulators$Quantiles1 == "Mixed_downstream", "Fraction_targets_with_CNAs"],
            me_p_fractions_no_CNVs_in_regulators[me_p_fractions_no_CNVs_in_regulators$Quantiles1 == "UC_downstream", "Fraction_targets_with_CNAs"])

wilcox.test(me_p_fractions_no_CNVs_in_regulators[me_p_fractions_no_CNVs_in_regulators$Quantiles1 == "Mixed_downstream", "Fraction_targets_with_CNAs"],
            me_p_fractions_no_CNVs_in_regulators[me_p_fractions_no_CNVs_in_regulators$Quantiles1 == "EM_downstream", "Fraction_targets_with_CNAs"])$p.value


me_p_fractions$Regulator_status <- factor(me_p_fractions$Regulator_status, levels=c("CNN", "CNA"))

pdf("Figure4_Targets_CNA_by_class_facet.pdf",
    width=8, height=4)
g <- ggplot(me_p_fractions[!is.na(me_p_fractions$Quantiles1),], aes(x=Regulator_status, y=Fraction_targets_with_CNAs))+
  geom_boxplot(aes(fill=Quantiles1, alpha=Regulator_status))+
  scale_alpha_manual(values=c(CNN=1, CNA=0.3))+
  scale_fill_jco()+
  theme_bw()+
  xlab("")+
  facet_grid(.~Quantiles1)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank())
print(g)
dev.off()


me_p_fractions_UC_t <- me_p_fractions[me_p_fractions$Quantiles1 == "UC_downstream",]
wilcox.test(me_p_fractions_UC_t$Fraction_targets_with_CNAs[me_p_fractions_UC_t$Regulator_status == "CNN"],
            me_p_fractions_UC_t$Fraction_targets_with_CNAs[me_p_fractions_UC_t$Regulator_status == "CNA"], alternative="greater")

me_p_fractions_EM_t <- me_p_fractions[me_p_fractions$Quantiles1 == "EM_downstream",]
wilcox.test(me_p_fractions_EM_t$Fraction_targets_with_CNAs[me_p_fractions_EM_t$Regulator_status == "CNN"],
            me_p_fractions_EM_t$Fraction_targets_with_CNAs[me_p_fractions_EM_t$Regulator_status == "CNA"], alternative="greater")$p.value

me_p_fractions_UC_EM_i <- me_p_fractions[me_p_fractions$Quantiles1 == "Mixed_downstream",]
wilcox.test(me_p_fractions_UC_EM_i$Fraction_targets_with_CNAs[me_p_fractions_UC_EM_i$Regulator_status == "CNN"],
            me_p_fractions_UC_EM_i$Fraction_targets_with_CNAs[me_p_fractions_UC_EM_i$Regulator_status == "CNA"], alternative="greater")

wilcox.test(me_p_fractions_UC_EM_i$Fraction_targets_with_CNAs[me_p_fractions_UC_EM_i$Regulator_status == "CNN"],
            me_p_fractions_UC_EM_i$Fraction_targets_with_CNAs[me_p_fractions_UC_EM_i$Regulator_status == "CNA"], alternative="less")


###Bar chart (waterfall plot)

me_p2 <- me_p

me_p2$Quantiles1 <- regulator_classification[match(me_p2$Regulator, regulator_classification$Regulator), "Regulator_class"]

me_p2$Diff_fraction <- (me_p2$Fraction_targets_regulator_CNV-me_p2$Fraction_targets_regulator_not_CNV)*100

me_p2 <- me_p2[order(me_p2$Diff_fraction, decreasing=FALSE),]
me_p2$Regulator <- factor(me_p2$Regulator, levels=me_p2$Regulator)

me_p2 <- me_p2[me_p2$Diff_fraction != 0,]

me_p2$Quantiles1 <- factor(me_p2$Quantiles1, levels=c("UC_downstream", "Mixed_downstream", "EM_downstream"))



me_p2_UC_EM <- me_p2[!is.na(me_p2$Quantiles1),]
me_p2_UC_EM <- me_p2_UC_EM[!is.na(me_p2_UC_EM$Age),]
me_p2_UC_EM <- me_p2_UC_EM[me_p2_UC_EM$Age != "MM",]

pdf("Supp4_Waterfall_by_age.pdf",
    height=5, width=8)
g <- ggplot(me_p2_UC_EM, aes(x=Regulator, y=Diff_fraction))+
  geom_bar(stat='identity', aes(fill=Quantiles1))+
  #scale_fill_manual(values=c(UC="#F8766D",
  #                           EM="#00BA38"))+
  scale_fill_jco()+
  ylim(-100,100)+
  coord_flip()+
  geom_hline(yintercept=0)+
  facet_wrap(Age~Quantiles1, scales='free_y')+
  theme_bw()+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        strip.text = element_blank())
print(g)
dev.off()


me_p2_UC_t <- me_p2_UC_EM[me_p2_UC_EM$Quantiles1 == "UC_downstream",]
sum(me_p2_UC_t$Diff_fraction[me_p2_UC_t$Age == "UC"] < 0)/nrow(me_p2_UC_t[me_p2_UC_t$Age == "UC",])
sum(me_p2_UC_t$Diff_fraction[me_p2_UC_t$Age == "EM"] < 0)/nrow(me_p2_UC_t[me_p2_UC_t$Age == "EM",])

me_p2_EM_t <- me_p2_UC_EM[me_p2_UC_EM$Quantiles1 == "EM_downstream",]
sum(me_p2_EM_t$Diff_fraction[me_p2_EM_t$Age == "UC"] < 0)/nrow(me_p2_EM_t[me_p2_EM_t$Age == "UC",])
sum(me_p2_EM_t$Diff_fraction[me_p2_EM_t$Age == "EM"] < 0)/nrow(me_p2_EM_t[me_p2_EM_t$Age == "EM",])

me_p2_UC_EM_i <- me_p2_UC_EM[me_p2_UC_EM$Quantiles1 == "Mixed_downstream",]
sum(me_p2_UC_EM_i$Diff_fraction[me_p2_UC_EM_i$Age == "UC"] < 0)/nrow(me_p2_UC_EM_i[me_p2_UC_EM_i$Age == "UC",])
sum(me_p2_UC_EM_i$Diff_fraction[me_p2_UC_EM_i$Age == "EM"] < 0)/nrow(me_p2_UC_EM_i[me_p2_UC_EM_i$Age == "EM",])


pdf("Figure4_Waterfall.pdf",
    height=3, width=8)
g <- ggplot(me_p2_UC_EM, aes(x=Regulator, y=Diff_fraction))+
  geom_bar(stat='identity', aes(fill=Quantiles1))+
  scale_fill_jco()+
  ylim(-100,100)+
  coord_flip()+
  geom_hline(yintercept=0)+
  facet_wrap(~Quantiles1, scales='free_y')+
  theme_bw()+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        strip.text = element_blank())
print(g)
dev.off()

me_p2_UC_t <- me_p2_UC_EM[me_p2_UC_EM$Quantiles1 == "UC_downstream",]
sum(me_p2_UC_t$Diff_fraction < 0)/nrow(me_p2_UC_t)
sum(me_p2_UC_t$Diff_fraction > 0)/nrow(me_p2_UC_t)

me_p2_EM_t <- me_p2_UC_EM[me_p2_UC_EM$Quantiles1 == "EM_downstream",]
sum(me_p2_EM_t$Diff_fraction < 0)/nrow(me_p2_EM_t)
sum(me_p2_EM_t$Diff_fraction > 0)/nrow(me_p2_EM_t)

me_p2_UC_EM_i <- me_p2_UC_EM[me_p2_UC_EM$Quantiles1 == "Mixed_downstream",]
sum(me_p2_UC_EM_i$Diff_fraction < 0)/nrow(me_p2_UC_EM_i)
sum(me_p2_UC_EM_i$Diff_fraction > 0)/nrow(me_p2_UC_EM_i)

View(me_p2_UC_EM[me_p2_UC_EM$Regulator == "MDM2",])
