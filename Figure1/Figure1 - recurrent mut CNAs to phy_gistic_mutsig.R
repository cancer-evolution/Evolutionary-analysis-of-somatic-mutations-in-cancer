###Recurrent Mutations and CNAs to phylostrata

library(ggplot2)
library(reshape2)
library(gridExtra)
library(Biobase)
library(clinfun)
library(readr)

source('helper_functions.R')

tumours <- c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "ESCA", "GBM", 
             "HNSC", "KICH", "KIRC", "KIRP", "LGG",  "LIHC", "LUAD", "LUSC", 
             "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", 
             "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM")

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
load("patients_with_mut_info.Rdata")
sum(sapply(patients_with_mut_info, length))

mutations_df <- load_mutations_mutsig()

n_mut_df <- vector()
for(tumour in tumours){
  local_mut_top <- mutations_df[mutations_df$Tumour == tumour,]
  
  local_mut_top <- local_mut_top[!is.na(local_mut_top$Gene_age),]
  
  local_mut_top$Frequency <- "Recurrent"  
  
  n_enriched <- table(local_mut_top[,c("Gene_age", "Frequency")])
  
  local_results <- data.frame(Tumour = tumour, Phy = 1:16, 
                              Frequency=c(rep("Recurrent", 16)))
  
  if("Recurrent" %in% local_mut_top$Frequency){
    local_results$N_mut <- c(n_enriched[,"Recurrent"][match(1:16, rownames(n_enriched))])
  }
  n_mut_df <- rbind(n_mut_df, local_results)
}

n_mut_df$Phy <- factor(n_mut_df$Phy, levels=1:16)

n_mut_df$Total <- number_each_phy[match(n_mut_df$Phy, number_each_phy$Phylostrata),"GeneID"]

colnames(n_mut_df) <- c("Tumour", "Phylostrata", "Frequency", "N_genes_mut", "Total_genes")

n_mut_df$Fraction_genes_mut <- n_mut_df$N_genes_mut/n_mut_df$Total_genes

temp_n_mut_df <- subset(n_mut_df, Frequency=="Recurrent")

#Calculate ranks
n_mut_df2 <- vector()
for(tumour in tumours){
  local_n_mut_df <- temp_n_mut_df[temp_n_mut_df$Tumour == tumour,]
  
  local_n_mut_df$Rank_mut1 <- rank(-local_n_mut_df$Fraction_genes_mut)
  local_n_mut_df[is.na(local_n_mut_df$Fraction_genes_mut), "Rank_mut1"] <- NA
  
  n_mut_df2 <- rbind(n_mut_df2, local_n_mut_df)
}

enrichment_of_recurrent_mutations <- n_mut_df2 
save(enrichment_of_recurrent_mutations, file="enrichment_of_recurrent_mutations.Rdata")

pdf("Figure1_freq_mut.pdf",
    height=2.5, width=6)
g <- ggplot(n_mut_df2, aes(x=Phylostrata, y = Rank_mut1))+
  geom_point(aes(colour=Tumour))+
  geom_path(aes(group=Phylostrata), colour="grey")+
  geom_hline(yintercept=8, colour="darkgrey")+
  geom_point(aes(colour=Tumour))+
  scale_y_reverse(breaks = 1:16)+
  ylab("Ranked fraction of\ngenes with point mutations")+
  xlab("")+
  geom_smooth(aes(x=as.numeric(as.character(Phylostrata)), y =Rank_mut1), method="loess", se=FALSE, colour="black", size=0.5)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())
print(g)
dev.off()

n_mut_df2$Age <- ifelse(n_mut_df2$Phylostrata %in% 1:3, "UC",
                        ifelse(n_mut_df2$Phylostrata %in% 4:9, "EM", "MM"))
n_mut_df2$Age <- factor(n_mut_df2$Age, levels=c("UC", "EM", "MM"))

pdf("Figure1_freq_mut_v2.pdf",
    height=2.5, width=6)
g <- ggplot(n_mut_df2, aes(x=Phylostrata, y = Rank_mut1))+
  geom_hline(yintercept=8, colour="darkgrey")+
  geom_boxplot(aes(colour=Age))+
  scale_y_reverse(breaks = 1:16)+
  ylab("Ranked fraction of\ngenes with point mutations")+
  xlab("")+
  geom_smooth(aes(x=as.numeric(as.character(Phylostrata)), y =Rank_mut1), method="loess", se=FALSE, colour="black", size=0.5)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())
print(g)
dev.off()




aggregate(Rank_mut1 ~ Phylostrata, n_mut_df2, median)

View(n_mut_df2[n_mut_df2$Phylostrata %in% 4:9, c("Tumour", "Phylostrata", "Rank_mut1")])

p_values_fraction <- vector()
for(tumour in tumours){
  local_n_mut_df2 <- n_mut_df2[n_mut_df2$Tumour == tumour,]
  
  groups <- local_n_mut_df2$Phylostrata
  groups <- factor(groups,levels=1:16, ordered=TRUE) 
  temp <- local_n_mut_df2$Fraction_genes_mut
  temp[is.na(temp)] <- 0
  p_mut <- jonckheere.test(g=groups, x=temp, alternative="decreasing")$p.value
  p_values_fraction <- rbind(p_values_fraction, c(tumour, p_mut))
}

p_values_fraction <- data.frame(p_values_fraction)
colnames(p_values_fraction) <- c("Tumour", "p_mut")

p_values_fraction[,2] <- p.adjust(p_values_fraction[,2], method="BH")


##How many phylostrata with at least 1 miss or LoF mutation?
n_phy_with_mut <- vector()
for(tumour in tumours){
  local_temp <- n_mut_df2[n_mut_df2$Tumour == tumour,]
  
  n_phy_with_mut <- rbind(n_phy_with_mut,
                          c(tumour, sum(local_temp$Fraction_genes_mut != 0, na.rm=TRUE)))
  
}

n_phy_with_mut <- as.data.frame(n_phy_with_mut)
colnames(n_phy_with_mut) <- c("Tumour", "N_phy_with_mut")


n_phy_with_mut$N_phy_with_mut <- as.numeric(as.character(n_phy_with_mut$N_phy_with_mut))
length(n_phy_with_mut[n_phy_with_mut$N_phy_with_mut >= 5,"Tumour"])

tumours_enough_mut <- as.character(n_phy_with_mut[n_phy_with_mut$N_phy_with_mut >= 5,"Tumour"])
View(n_mut_df2[intersect(which(n_mut_df2$Phylostrata %in% 4:9),
                         which(n_mut_df2$Tumour %in% tumours_enough_mut)), c("Tumour", "Phylostrata", "Rank_mut1")])

top_5_mut <- vector()
for(tumour in tumours_enough_mut){
  temp <- n_mut_df2[n_mut_df2$Tumour == tumour,]
  temp <- temp[temp$Age == "EM",]
  top_5_mut <- c(top_5_mut, sum(temp$Rank_mut1 <= 5, na.rm=TRUE))
}

sum(top_5_mut >= 3)



###CNAs

##Load in CNVs
#Options: any, focal, broad, >10, <10

###Fraction of patients with CNV
load("CNVs_curated_gistic.Rdata")
load("patients_with_CNV_info.Rdata")
sum(sapply(patients_with_CNV_info, length))

CNVs_df <- load_CNVs_gistic(only_focal="ANY")
CNVs_df[CNVs_df$Genes == "TNFRSF17",]  ##and for all other genes to get their amp/del frequency
CNVs_df[CNVs_df$Genes == "EGFR",]

n_CNV_df <- vector()
for(tumour in tumours){
  temp <- CNVs_df[CNVs_df$Tumour == tumour,]
  temp <- as.data.frame(table(temp[,c("CNA", "Gene_age")]))
  temp$Tumour <- tumour
  temp$Total <- number_each_phy[match(temp$Gene_age, number_each_phy$Phylostrata), "GeneID"]
  
  temp_amp <- temp[temp$CNA == "Amplification",]
  temp_del <- temp[temp$CNA == "Deletion",]
  temp_amp$Fraction_amp <- temp_amp$Freq/temp_amp$Total
  temp_del$Fraction_del <- temp_del$Freq/temp_del$Total
  
  temp_df <- data.frame(Phylostrata = 1:16,
                        Tumour = tumour,
                        N_amp = temp_amp[match(1:16, temp_amp$Gene_age), "Freq"],
                        N_del = temp_del[match(1:16, temp_del$Gene_age), "Freq"],
                        Fraction_amp = temp_amp[match(1:16, temp_amp$Gene_age), "Fraction_amp"],
                        Fraction_del = temp_del[match(1:16, temp_del$Gene_age), "Fraction_del"])
  
  n_CNV_df <- rbind(n_CNV_df, temp_df)
}

n_CNV_df2 <- vector()
for(tumour in tumours){
  local_n_CNV_df <- n_CNV_df[n_CNV_df$Tumour == tumour,]
  
  local_n_CNV_df$Rank_amp1 <- rank(-local_n_CNV_df$N_amp)
  local_n_CNV_df$Rank_del1 <- rank(-local_n_CNV_df$N_del)
  
  local_n_CNV_df[is.na(local_n_CNV_df$Fraction_genes_amp), "Rank_amp1"] <- NA
  local_n_CNV_df[is.na(local_n_CNV_df$Fraction_genes_del), "Rank_del1"] <- NA
  
  n_CNV_df2 <- rbind(n_CNV_df2, local_n_CNV_df)
}

enrichment_of_recurrent_CNVs <- n_CNV_df2 
#save(enrichment_of_recurrent_CNVs, file="enrichment_of_recurrent_CNVs.Rdata")


p_values_fraction <- vector()
for(tumour in tumours){
  local_n_CNV_df2 <- n_CNV_df2[n_CNV_df2$Tumour == tumour,]
  
  groups <- local_n_CNV_df2$Phylostrata
  groups <- factor(groups,levels=1:16, ordered=TRUE) 
  p_amp <- jonckheere.test(g=groups, x=local_n_CNV_df2$Fraction_amp, alternative="decreasing")$p.value
  p_del <- jonckheere.test(g=groups, x=local_n_CNV_df2$Fraction_del, alternative="decreasing")$p.value
  p_values_fraction <- rbind(p_values_fraction, c(tumour, p_amp, p_del))
}

p_values_fraction <- data.frame(p_values_fraction)
colnames(p_values_fraction) <- c("Tumour", "p_amp", "p_del")

p_values_fraction[,2:3] <- apply(p_values_fraction[,2:3], 2, p.adjust, method="BH")

pdf("Figure1_freq_amp.pdf",
    height=2.5, width=6.5)
g <- ggplot(n_CNV_df2, aes(x=Phylostrata, y = Rank_amp1))+
  geom_point(aes(colour=Tumour))+
  geom_path(aes(group=Phylostrata), colour="grey")+
  geom_hline(yintercept=8, colour="darkgrey")+
  geom_point(aes(colour=Tumour))+
  scale_y_reverse(breaks = 1:16)+
  scale_x_continuous(breaks = 1:16)+
  ylab("Ranked fraction of\namplified genes")+
  xlab("")+
  #geom_smooth(aes(x=as.numeric(as.character(Phylostrata)), y =Rank_amp1), method="loess", se=FALSE, colour="black", size=0.5)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())
print(g)
dev.off()

pdf("Figure1_freq_del.pdf",
    height=2.5, width=6.5)
g <- ggplot(n_CNV_df2, aes(x=Phylostrata, y = Rank_del1))+
  geom_point(aes(colour=Tumour))+
  geom_path(aes(group=Phylostrata), colour="grey")+
  geom_hline(yintercept=8, colour="darkgrey")+
  geom_point(aes(colour=Tumour))+
  scale_y_reverse(breaks = 1:16)+
  scale_x_continuous(breaks = 1:16)+
  ylab("Ranked fraction of\ndeleted genes")+
  xlab("")+
  #geom_smooth(aes(x=as.numeric(as.character(Phylostrata)), y =Rank_del1), method="loess", se=FALSE, colour="black", size=0.5)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())
print(g)
dev.off()

View(n_CNV_df2[n_CNV_df2$Phylostrata %in% 4:9,])


n_CNV_df2$Age <- ifelse(n_CNV_df2$Phylostrata %in% 1:3, "UC",
                        ifelse(n_CNV_df2$Phylostrata %in% 4:9, "EM", "MM"))
n_CNV_df2$Age <- factor(n_CNV_df2$Age, levels=c("UC", "EM", "MM"))
n_CNV_df2$Phylostrata <- factor(n_CNV_df2$Phylostrata, levels=1:16)

pdf("Figure1_freq_amp_v2.pdf",
    height=2.5, width=6)
g <- ggplot(n_CNV_df2, aes(x=Phylostrata, y = Rank_amp1))+
  geom_boxplot(aes(colour=Age))+
  scale_y_reverse(breaks = 1:16)+
  ylab("Ranked fraction of\namplified genes")+
  xlab("")+
  geom_smooth(aes(x=as.numeric(as.character(Phylostrata)), y =Rank_amp1), method="loess", se=FALSE, colour="black", size=0.5)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())
print(g)
dev.off()

pdf("Figure1_freq_del_v2.pdf",
    height=2.5, width=6)
g <- ggplot(n_CNV_df2, aes(x=Phylostrata, y = Rank_del1))+
  geom_boxplot(aes(colour=Age))+
  scale_y_reverse(breaks = 1:16)+
  ylab("Ranked fraction of\ndeleted genes")+
  xlab("")+
  geom_smooth(aes(x=as.numeric(as.character(Phylostrata)), y =Rank_del1), method="loess", se=FALSE, colour="black", size=0.5)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())
print(g)
dev.off()




n_phy_with_amp <- vector()
n_phy_with_del <- vector()
for(tumour in tumours){
  local_temp <- n_CNV_df2[n_CNV_df2$Tumour == tumour,]
  
  n_phy_with_amp <- rbind(n_phy_with_amp,
                          c(tumour, sum(local_temp$Fraction_amp != 0, na.rm=TRUE)))
  n_phy_with_del <- rbind(n_phy_with_del,
                          c(tumour, sum(local_temp$Fraction_del != 0, na.rm=TRUE)))
  
}

n_phy_with_amp <- as.data.frame(n_phy_with_amp)
colnames(n_phy_with_amp) <- c("Tumour", "N_phy_with_amp")
n_phy_with_del <- as.data.frame(n_phy_with_del)
colnames(n_phy_with_del) <- c("Tumour", "N_phy_with_del")

n_phy_with_amp$N_phy_with_amp <- as.numeric(as.character(n_phy_with_amp$N_phy_with_amp))
length(n_phy_with_amp[n_phy_with_amp$N_phy_with_amp >= 5,"Tumour"])

tumours_enough_amp <- as.character(n_phy_with_amp[n_phy_with_amp$N_phy_with_amp >= 5,"Tumour"])

n_phy_with_del$N_phy_with_del <- as.numeric(as.character(n_phy_with_del$N_phy_with_del))
length(n_phy_with_del[n_phy_with_del$N_phy_with_del >= 5,"Tumour"])

tumours_enough_del <- as.character(n_phy_with_del[n_phy_with_del$N_phy_with_del >= 5,"Tumour"])


top_5_amp <- vector()
for(tumour in tumours_enough_amp){
  temp <- n_CNV_df2[n_CNV_df2$Tumour == tumour,]
  temp <- temp[temp$Age == "EM",]
  top_5_amp <- c(top_5_amp, sum(temp$Rank_amp1 <= 5, na.rm=TRUE))
}

sum(top_5_amp >= 3)


top_5_del <- vector()
for(tumour in tumours_enough_del){
  temp <- n_CNV_df2[n_CNV_df2$Tumour == tumour,]
  temp <- temp[temp$Age == "EM",]
  top_5_del <- c(top_5_del, sum(temp$Rank_del1 <= 5, na.rm=TRUE))
}

sum(top_5_del >= 3)
