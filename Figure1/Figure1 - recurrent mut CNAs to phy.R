###Recurrent Mutations and CNAs to phylostrata

library(ggplot2)
library(reshape2)
library(gridExtra)
library(Biobase)
library(clinfun)
library(readr)

source("CNV_functions.R")
source('helper_functions.R')

tumours <- c("LUAD", "LUSC", "BRCA", "PRAD", "LIHC", "COAD", "STAD")


genes_phy <- read.csv("gene_phylostrata.txt")  ###File with gene name, entrez and phylostratum as column
genes_phy_categorical <- genes_phy
genes_phy_categorical$Phylostrata <- ifelse(genes_phy_categorical$Phylostrata %in% 1:3, "UC",
                                            ifelse(genes_phy_categorical$Phylostrata %in% 4:9, "EM",
                                                   ifelse(genes_phy_categorical$Phylostrata %in% 10:16, "MM", NA)))
n_genes_phy <- table(genes_phy$Phylostrata)
n_genes_phy <- data.frame(Phy = names(n_genes_phy), Number = as.vector(n_genes_phy))

n_genes_phy$Phy <- factor(n_genes_phy$Phy, levels=1:16)

n_genes_phy$Age <- ifelse(n_genes_phy$Phy %in% 1:3, "UC",
                          ifelse(n_genes_phy$Phy %in% 4:9, "EM",
                                 ifelse(n_genes_phy$Phy %in% 10:16, "MM", NA)))

n_genes_phy$Age <- factor(n_genes_phy$Age, levels=c("UC", "EM", "MM"))

pdf("Supp1_Gene_number_in_phy.pdf",
    height=3, width=6)
g <- ggplot(n_genes_phy, aes(x=Phy, y=Number))+
  geom_bar(stat='identity', aes(fill=Age))+
  ylab("Number of genes")+
  xlab("Phylostratum")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())
print(g)
dev.off()

UC_genes <- as.character(genes_phy_categorical[genes_phy_categorical$Phylostrata == "UC", "GeneID"])
EM_genes <- as.character(genes_phy_categorical[genes_phy_categorical$Phylostrata == "EM", "GeneID"])
MM_genes <- as.character(genes_phy_categorical[genes_phy_categorical$Phylostrata == "MM", "GeneID"])


number_each_phy <- aggregate(GeneID ~ Phylostrata, genes_phy, length)



##Mutations
load("patients_with_Mut_info.Rdata")
sum(sapply(patients_with_mut_info, length))

mutations_df <- load_mutations()

n_mut_df <- vector()
for(tumour in tumours){
  local_mut_top <- mutations_df[mutations_df$Tumour == tumour,]
  
  local_mut_top <- local_mut_top[!is.na(local_mut_top$Gene_age),]
  
  local_mut_top$Frequency <- NA
  
  local_mut_top$Frequency[intersect(which(local_mut_top$Number_mutations >=3),
                                    which(local_mut_top$Syn_ratio >1))] <- "Recurrent"  
  
  non_recurrent <- local_mut_top[is.na(local_mut_top$Frequency), ]
  #non_recurrent <- non_recurrent[non_recurrent$Number_synonymous >=3,]
  
  local_mut_top$Frequency[which(local_mut_top$Syn_ratio <=1)] <- "Purified"
  
  local_mut_top_miss <- local_mut_top[local_mut_top$Variant_type == "Missense",]
  local_mut_top_lof <- local_mut_top[local_mut_top$Variant_type == "LoF",]
  
  n_enriched_miss <- table(local_mut_top_miss[,c("Gene_age", "Frequency")])
  n_enriched_lof <- table(local_mut_top_lof[,c("Gene_age", "Frequency")])
  n_syn <- table(non_recurrent[,c("Gene_age")])
  
  local_results <- data.frame(Tumour = tumour, Phy = 1:16, 
                              Frequency=c(rep("Recurrent", 16),
                                          rep("Purified", 16)))
  
  local_results$N_miss <- c(n_enriched_miss[,"Recurrent"][match(1:16, rownames(n_enriched_miss))],
                            n_enriched_miss[,"Purified"][match(1:16, rownames(n_enriched_miss))])
  local_results$N_lof <- c(n_enriched_lof[,"Recurrent"][match(1:16, rownames(n_enriched_lof))],
                           n_enriched_lof[,"Purified"][match(1:16, rownames(n_enriched_lof))])
  local_results$N_syn <- n_syn[match(1:16,names(n_syn))]
  
  n_mut_df <- rbind(n_mut_df, local_results)
}

n_mut_df$Phy <- factor(n_mut_df$Phy, levels=1:16)

n_mut_df$Total <- number_each_phy[match(n_mut_df$Phy, number_each_phy$Phylostrata),"GeneID"]

colnames(n_mut_df) <- c("Tumour", "Phylostrata", "Frequency", "N_genes_miss", "N_genes_lof", "N_genes_syn", "Total_genes")

n_mut_df$Fraction_genes_miss <- n_mut_df$N_genes_miss/n_mut_df$Total_genes
n_mut_df$Fraction_genes_lof <- n_mut_df$N_genes_lof/n_mut_df$Total_genes
n_mut_df$Fraction_genes_syn <- n_mut_df$N_genes_syn/n_mut_df$Total_genes



temp_n_mut_df <- subset(n_mut_df, Frequency=="Recurrent")

#Calculate ranks
n_mut_df2 <- vector()
for(tumour in tumours){
  local_n_mut_df <- temp_n_mut_df[temp_n_mut_df$Tumour == tumour,]
  
  local_n_mut_df$Rank_miss1 <- rank(-local_n_mut_df$Fraction_genes_miss)
  local_n_mut_df$Rank_lof1 <- rank(-local_n_mut_df$Fraction_genes_lof)
  
  n_mut_df2 <- rbind(n_mut_df2, local_n_mut_df)
}

enrichment_of_recurrent_mutations <- n_mut_df2 
#save(enrichment_of_recurrent_mutations, file="enrichment_of_recurrent_mutations.Rdata")

pdf("Figure1_freq_miss.pdf",
    height=2.5, width=6)
g <- ggplot(n_mut_df2, aes(x=Phylostrata, y = Rank_miss1))+
  geom_point(aes(colour=Tumour))+
  geom_path(aes(group=Phylostrata), colour="grey")+
  geom_hline(yintercept=8, colour="darkgrey")+
  geom_point(aes(colour=Tumour))+
  scale_y_reverse(breaks = 1:16)+
  ylab("Ranked fraction of\ngenes with missense mutations")+
  xlab("")+
  geom_smooth(aes(x=as.numeric(as.character(Phylostrata)), y =Rank_miss1), method="loess", se=FALSE, colour="black", size=0.5)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())
print(g)
dev.off()

pdf("Figure1_freq_lof.pdf",
    height=2.5, width=6)
g <- ggplot(n_mut_df2, aes(x=Phylostrata, y = Rank_lof1))+
  geom_point(aes(colour=Tumour))+
  geom_path(aes(group=Phylostrata), colour="grey")+
  geom_hline(yintercept=8, colour="darkgrey")+
  geom_point(aes(colour=Tumour))+
  scale_y_reverse(breaks = 1:16)+
  ylab("Ranked fraction of\ngenes with lof mutations")+
  xlab("")+
  geom_smooth(aes(x=as.numeric(as.character(Phylostrata)), y =Rank_lof1), method="loess", se=FALSE, colour="black", size=0.5)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())
print(g)
dev.off()

aggregate(Rank_miss1 ~ Phylostrata, n_mut_df2, median)
aggregate(Rank_lof1 ~ Phylostrata, n_mut_df2, median)

View(n_mut_df2[n_mut_df2$Phylostrata %in% 4:9, c("Tumour", "Phylostrata", "Rank_miss1")])
View(n_mut_df2[n_mut_df2$Phylostrata %in% 4:9, c("Tumour", "Phylostrata", "Rank_lof1")])

p_values_fraction <- vector()
for(tumour in tumours){
  local_n_mut_df2 <- n_mut_df2[n_mut_df2$Tumour == tumour,]
  
  groups <- local_n_mut_df2$Phylostrata
  groups <- factor(groups,levels=1:16, ordered=TRUE) 
  p_miss <- jonckheere.test(g=groups, x=local_n_mut_df2$Fraction_genes_miss, alternative="decreasing")$p.value
  p_lof <- jonckheere.test(g=groups, x=local_n_mut_df2$Fraction_genes_lof, alternative="decreasing")$p.value
  p_values_fraction <- rbind(p_values_fraction, c(tumour, p_miss, p_lof))
}

p_values_fraction <- data.frame(p_values_fraction)
colnames(p_values_fraction) <- c("Tumour", "p_miss", "p_lof")

p_values_fraction[,2:3] <- apply(p_values_fraction[,2:3], 2, p.adjust, method="BH")




###CNAs

##Load in CNVs
#Options: any, focal, broad, >10, <10

###Fraction of patients with CNV
load("CNVs_curated2.Rdata")
load("patients_with_CNV_info.Rdata")
sum(sapply(patients_with_CNV_info, length))

CNVs_df <- load_CNVs(only_focal="ANY")
CNVs_df[CNVs_df$Genes == "TNFRSF17",]  ##and for all other genes to get their amp/del frequency
CNVs_df[CNVs_df$Genes == "EGFR",]
mean(CNVs_df[CNVs_df$Genes == "TP53","Patients_del"])

patient_frequency_CNV_df <- calculate_fraction_of_patients_with_CNV(CNVs_df)

n_CNV_df <- subset(patient_frequency_CNV_df, Frequency==">0.10")

n_CNV_df2 <- vector()
for(tumour in tumours){
  local_n_CNV_df <- n_CNV_df[n_CNV_df$Tumour == tumour,]
  
  local_n_CNV_df$Rank_amp1 <- rank(-local_n_CNV_df$Fraction_genes_amp)
  local_n_CNV_df$Rank_del1 <- rank(-local_n_CNV_df$Fraction_genes_del)
  
  n_CNV_df2 <- rbind(n_CNV_df2, local_n_CNV_df)
}

View(n_CNV_df2[n_CNV_df2$Phylostrata %in% 4:9, c("Tumour", "Phylostrata", "Rank_amp1")])
View(n_CNV_df2[n_CNV_df2$Phylostrata %in% 4:9, c("Tumour", "Phylostrata", "Rank_del1")])

enrichment_of_recurrent_CNVs <- n_CNV_df2 
#save(enrichment_of_recurrent_CNVs, file="enrichment_of_recurrent_CNVs.Rdata")


p_values_fraction <- vector()
for(tumour in tumours){
  local_n_CNV_df2 <- n_CNV_df2[n_CNV_df2$Tumour == tumour,]
  
  groups <- local_n_CNV_df2$Phylostrata
  groups <- factor(groups,levels=1:16, ordered=TRUE) 
  p_amp <- jonckheere.test(g=groups, x=local_n_CNV_df2$Fraction_genes_amp, alternative="decreasing")$p.value
  p_del <- jonckheere.test(g=groups, x=local_n_CNV_df2$Fraction_genes_del, alternative="decreasing")$p.value
  p_values_fraction <- rbind(p_values_fraction, c(tumour, p_amp, p_del))
}

p_values_fraction <- data.frame(p_values_fraction)
colnames(p_values_fraction) <- c("Tumour", "p_amp", "p_del")

p_values_fraction[,2:3] <- apply(p_values_fraction[,2:3], 2, p.adjust, method="BH")

pdf("Figure1_freq_amp.pdf",
    height=2.5, width=6)
g <- ggplot(n_CNV_df2, aes(x=Phylostrata, y = Rank_amp1))+
  geom_point(aes(colour=Tumour))+
  geom_path(aes(group=Phylostrata), colour="grey")+
  geom_hline(yintercept=8, colour="darkgrey")+
  geom_point(aes(colour=Tumour))+
  scale_y_reverse(breaks = 1:16)+
  ylab("Ranked fraction of\namplified genes")+
  xlab("")+
  geom_smooth(aes(x=as.numeric(as.character(Phylostrata)), y =Rank_amp1), method="loess", se=FALSE, colour="black", size=0.5)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())
print(g)
dev.off()

pdf("Figure1_freq_del.pdf",
    height=2.5, width=6)
g <- ggplot(n_CNV_df2, aes(x=Phylostrata, y = Rank_del1))+
  geom_point(aes(colour=Tumour))+
  geom_path(aes(group=Phylostrata), colour="grey")+
  geom_hline(yintercept=8, colour="darkgrey")+
  geom_point(aes(colour=Tumour))+
  scale_y_reverse(breaks = 1:16)+
  ylab("Ranked fraction of\ndeleted genes")+
  xlab("")+
  geom_smooth(aes(x=as.numeric(as.character(Phylostrata)), y =Rank_del1), method="loess", se=FALSE, colour="black", size=0.5)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())
print(g)
dev.off()



##Non-recurrent alterations

##Of all altered genes, calculate the fraction that is recurrent, and the 
#fraction that is non-recurrent

total_amp_genes <- aggregate(N_genes_amp ~ Tumour+Phylostrata, patient_frequency_CNV_df, sum)
colnames(total_amp_genes)[3] <- "Total_amp_genes"

total_del_genes <- aggregate(N_genes_del ~ Tumour+Phylostrata, patient_frequency_CNV_df, sum)
colnames(total_del_genes)[3] <- "Total_del_genes"

fraction_rec_non_CNV <- patient_frequency_CNV_df[,c("Tumour", "Frequency", "Phylostrata",
                                                    "N_genes_amp", "N_genes_del")]

fraction_rec_non_CNV$Total_amp_genes <- total_amp_genes[match(paste(fraction_rec_non_CNV$Tumour, fraction_rec_non_CNV$Phylostrata),
                                                              paste(total_amp_genes$Tumour, total_amp_genes$Phylostrata)),"Total_amp_genes"]

fraction_rec_non_CNV$Total_del_genes <- total_del_genes[match(paste(fraction_rec_non_CNV$Tumour, fraction_rec_non_CNV$Phylostrata),
                                                              paste(total_del_genes$Tumour, total_del_genes$Phylostrata)),"Total_del_genes"]
fraction_rec_non_CNV$Frac_amp <- fraction_rec_non_CNV$N_genes_amp/fraction_rec_non_CNV$Total_amp_genes

fraction_rec_non_CNV$Frac_del <- fraction_rec_non_CNV$N_genes_del/fraction_rec_non_CNV$Total_del_genes

pdf("Supp1_Non_recurrent_amp.pdf",
    height=2.5, width=16)
g <- ggplot(subset(fraction_rec_non_CNV, Frequency=="<0.10"), aes(x=Phylostrata, y=Frac_amp))+
  #geom_bar(aes(fill=Frequency), position='dodge', stat='identity')+
  geom_point()+
  geom_line(aes(group=Tumour))+
  facet_grid(.~Tumour)+
  geom_rect(aes(xmin=0.5, xmax=3.5, ymin=0.45, ymax=Inf), fill="#F8766D", alpha=0.02)+
  geom_rect(aes(xmin=3.5, xmax=9.5, ymin=0.45, ymax=Inf), fill="#00BA38", alpha=0.02)+
  geom_rect(aes(xmin=9.5, xmax=16.5, ymin=0.45, ymax=Inf), fill="#619CFF", alpha=0.02)+
  geom_point(size=2)+
  geom_line(aes(group=Tumour))+
  geom_smooth(aes(x=as.numeric(as.character(Phylostrata)), y =Frac_amp), 
              method="loess", se=FALSE, colour="red", size=0.5)+
  theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
print(g)
dev.off()



pdf("Supp1_Non_recurrent_del.pdf",
    height=2.5, width=16)
g <- ggplot(subset(fraction_rec_non_CNV, Frequency=="<0.10"), aes(x=Phylostrata, y=Frac_del))+
  #geom_bar(aes(fill=Frequency), position='dodge', stat='identity')+
  geom_point()+
  geom_line(aes(group=Tumour))+
  facet_grid(.~Tumour)+
  geom_rect(aes(xmin=0.5, xmax=3.5, ymin=0.7, ymax=Inf), fill="#F8766D", alpha=0.02)+
  geom_rect(aes(xmin=3.5, xmax=9.5, ymin=0.7, ymax=Inf), fill="#00BA38", alpha=0.02)+
  geom_rect(aes(xmin=9.5, xmax=16.5, ymin=0.7, ymax=Inf), fill="#619CFF", alpha=0.02)+
  geom_line(aes(group=Tumour))+
  geom_point(size=2)+
  geom_smooth(aes(x=as.numeric(as.character(Phylostrata)), y =Frac_del), 
              method="loess", se=FALSE, colour="red", size=0.5)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
print(g)
dev.off()



p_values_fraction <- vector()
for(tumour in tumours){
  local_n_CNV_df2 <- subset(fraction_rec_non_CNV, Frequency=="<0.10")[subset(fraction_rec_non_CNV, Frequency=="<0.10")$Tumour == tumour,]
  
  groups <- local_n_CNV_df2$Phylostrata
  groups <- factor(groups,levels=1:16, ordered=TRUE) 
  p_amp <- jonckheere.test(g=groups, x=local_n_CNV_df2$Frac_amp, alternative="increasing")$p.value
  p_del <- jonckheere.test(g=groups, x=local_n_CNV_df2$Frac_del, alternative="increasing")$p.value
  p_values_fraction <- rbind(p_values_fraction, c(tumour, p_amp, p_del))
}

p_values_fraction <- data.frame(p_values_fraction)
colnames(p_values_fraction) <- c("Tumour", "p_amp", "p_del")

p_values_fraction[,2:3] <- apply(p_values_fraction[,2:3], 2, p.adjust, method="BH")


total_miss_genes <- aggregate(N_genes_miss ~ Tumour+Phylostrata, n_mut_df, sum)
colnames(total_miss_genes)[3] <- "Total_miss_genes"

total_lof_genes <- aggregate(N_genes_lof ~ Tumour+Phylostrata, n_mut_df, sum)
colnames(total_lof_genes)[3] <- "Total_lof_genes"

fraction_rec_non_mut <- n_mut_df[,c("Tumour", "Frequency", "Phylostrata",
                                                    "N_genes_miss", "N_genes_lof")]

fraction_rec_non_mut$Total_miss_genes <- total_miss_genes[match(paste(fraction_rec_non_mut$Tumour, fraction_rec_non_mut$Phylostrata),
                                                                paste(total_miss_genes$Tumour, total_miss_genes$Phylostrata)),"Total_miss_genes"]

fraction_rec_non_mut$Total_lof_genes <- total_lof_genes[match(paste(fraction_rec_non_mut$Tumour, fraction_rec_non_mut$Phylostrata),
                                                              paste(total_lof_genes$Tumour, total_lof_genes$Phylostrata)),"Total_lof_genes"]
fraction_rec_non_mut$Frac_miss <- fraction_rec_non_mut$N_genes_miss/fraction_rec_non_mut$Total_miss_genes

fraction_rec_non_mut$Frac_lof <- fraction_rec_non_mut$N_genes_lof/fraction_rec_non_mut$Total_lof_genes


pdf("Supp1_Non_recurrent_miss.pdf",
    height=2.5, width=16)
g <- ggplot(subset(fraction_rec_non_mut, Frequency=="Purified"), aes(x=Phylostrata, y=Frac_miss))+
  #geom_bar(aes(fill=Frequency), position='dodge', stat='identity')+
  geom_point()+
  geom_line(aes(group=Tumour))+
  facet_grid(.~Tumour)+
  geom_rect(aes(xmin=0.5, xmax=3.5, ymin=0.45, ymax=Inf), fill="#F8766D", alpha=0.02)+
  geom_rect(aes(xmin=3.5, xmax=9.5, ymin=0.45, ymax=Inf), fill="#00BA38", alpha=0.02)+
  geom_rect(aes(xmin=9.5, xmax=16.5, ymin=0.45, ymax=Inf), fill="#619CFF", alpha=0.02)+
  geom_point(size=2)+
  geom_line(aes(group=Tumour))+
  geom_smooth(aes(x=as.numeric(as.character(Phylostrata)), y =Frac_miss), 
              method="loess", se=FALSE, colour="red", size=0.5)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
print(g)
dev.off()

pdf("Supp1_Non_recurrent_lof.pdf",
    height=2.5, width=16)
g <- ggplot(subset(fraction_rec_non_mut, Frequency=="Purified"), aes(x=Phylostrata, y=Frac_lof))+
  #geom_bar(aes(fill=Frequency), position='dodge', stat='identity')+
  geom_point()+
  geom_line(aes(group=Tumour))+
  facet_grid(.~Tumour)+
  geom_rect(aes(xmin=0.5, xmax=3.5, ymin=0.7, ymax=Inf), fill="#F8766D", alpha=0.02)+
  geom_rect(aes(xmin=3.5, xmax=9.5, ymin=0.7, ymax=Inf), fill="#00BA38", alpha=0.02)+
  geom_rect(aes(xmin=9.5, xmax=16.5, ymin=0.7, ymax=Inf), fill="#619CFF", alpha=0.02)+
  geom_line(aes(group=Tumour))+
  geom_point(size=2)+
  geom_smooth(aes(x=as.numeric(as.character(Phylostrata)), y =Frac_lof), 
              method="loess", se=FALSE, colour="red", size=0.5)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
print(g)
dev.off()
##

p_values_fraction <- vector()
for(tumour in tumours){
  local_n_mut_df2 <- subset(fraction_rec_non_mut, Frequency=="Purified")[subset(fraction_rec_non_mut, Frequency=="Purified")$Tumour == tumour,]
  
  groups <- local_n_mut_df2$Phylostrata
  groups <- factor(groups,levels=1:16, ordered=TRUE) 
  p_miss <- jonckheere.test(g=groups, x=local_n_mut_df2$Frac_miss, alternative="increasing")$p.value
  p_lof <- jonckheere.test(g=groups, x=local_n_mut_df2$Frac_lof, alternative="increasing")$p.value
  p_values_fraction <- rbind(p_values_fraction, c(tumour, p_miss, p_lof))
}

p_values_fraction <- data.frame(p_values_fraction)
colnames(p_values_fraction) <- c("Tumour", "p_miss", "p_lof")

p_values_fraction[,2:3] <- apply(p_values_fraction[,2:3], 2, p.adjust, method="BH")

