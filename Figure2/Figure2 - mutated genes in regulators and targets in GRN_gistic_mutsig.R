##Figure 2. Network analyses

library(ggplot2)
library(igraph)
library(reshape2)
library(gridExtra)
library(readr)
library(clinfun)
library(ggrepel)

source("helper_functions.R")

genes_phy <- read.csv("gene_phylostrata.txt")  ###File with gene name, entrez and phylostratum as column
genes_phy$Age <- ifelse(genes_phy$Phylostrata %in% 1:3, "UC",
                        ifelse(genes_phy$Phylostrata %in% 4:9, "EM",
                               ifelse(genes_phy$Phylostrata %in% 10:16, "MM", NA)))

UC_genes <- as.character(genes_phy[which(genes_phy$Phylostrata %in% 1:3),1])
EM_genes <- as.character(genes_phy[which(genes_phy$Phylostrata %in% 4:9),1])
MM_genes <- as.character(genes_phy[which(genes_phy$Phylostrata %in% 10:16),1])

tumours <- c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "ESCA", "GBM", 
             "HNSC", "KICH", "KIRC", "KIRP", "LGG",  "LIHC", "LUAD", "LUSC", 
             "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", 
             "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM")  

network_list <- list()
databases <- c("PathwayCommons", "Biogrid","WebIM", "GRN")

for(database in databases){
  if(database == "Biogrid"){
    network <- read.delim("BIOGRID-ORGANISM-Homo_sapiens-3.4.152.tab2.txt")
    network <- network[,c("Official.Symbol.Interactor.A", "Official.Symbol.Interactor.B")]
    colnames(network) <- c("Gene1", "Gene2")  
    network_graph <- graph.data.frame(network, directed=FALSE)
    network_graph <- simplify(network_graph, remove.loops=TRUE, remove.multiple=TRUE)
    
    network <- get.data.frame(network_graph)
    colnames(network) <- c("Gene1", "Gene2")
    
  }else if(database == "PathwayCommons"){
    network <- read.delim("PathwayCommons9.All.hgnc.sif", header=FALSE, dec=",")
    network <- network[!(network[,2] %in% c("controls-production-of", "controls-transport-of-chemical",
                                            "chemical-affects", "consumption-controlled-by", "reacts-with",
                                            "used-to-produce")),]  ##To remove all reactions involving chemicals
    network <- network[!(network[,2] %in% c("controls-expression-of")),]
    network <- network[,c(1,3)]
    colnames(network) <- c("Gene1", "Gene2")
    network_graph <- graph.data.frame(network, directed=FALSE)
    network_graph <- simplify(network_graph, remove.loops=TRUE, remove.multiple=TRUE)
    
    network <- get.data.frame(network_graph)
    colnames(network) <- c("Gene1", "Gene2")
    
  }else if(database == "WebIM"){
    network <- read.delim("InBio_Map_core_2016_09_12/core.psimitab",
                          header=FALSE)
    
    network <- network[, c("V5", "V6")]
    colnames(network) <- c("Gene1", "Gene2")
    network <- as.data.frame(network)
    network$Gene1 <- substr(as.character(network$Gene1), 11, nchar(as.character(network$Gene1)))
    network$Gene2 <- substr(as.character(network$Gene2), 11, nchar(as.character(network$Gene2)))
    
    network$Gene1 <- gsub("\\(.*","",network$Gene1)
    network$Gene2 <- gsub("\\(.*","",network$Gene2)
    
    network_graph <- graph.data.frame(network, directed=FALSE)
    network_graph <- simplify(network_graph, remove.loops=TRUE, remove.multiple=TRUE)
    
    network <- get.data.frame(network_graph)
    colnames(network) <- c("Gene1", "Gene2")
    
  }else if(database == "GRN"){
    network <- read.delim("PathwayCommons9.All.hgnc.sif", header=FALSE, dec=",")
    network <- network[!(network[,2] %in% c("controls-production-of", "controls-transport-of-chemical",
                                            "chemical-affects", "consumption-controlled-by", "reacts-with",
                                            "used-to-produce")),]  ##To remove all reactions involving chemicals
    network <- network[,c(1,3,2)]
    colnames(network) <- c("Gene1", "Gene2", "Edge_type")
    
    network <- network[network$Edge_type == "controls-expression-of",]
    
    network_graph <- graph.data.frame(network, directed=TRUE)
    network_graph <- simplify(network_graph, remove.loops=TRUE, remove.multiple=TRUE)
    
    network <- get.data.frame(network_graph)
    colnames(network) <- c("Gene1", "Gene2")
    
  }
  
  network$Age1 <- ifelse(network$Gene1 %in% UC_genes, "UC", 
                         ifelse(network$Gene1 %in% EM_genes, "EM", 
                                ifelse(network$Gene1 %in% MM_genes, "MM", NA)))
  network$Age2 <- ifelse(network$Gene2 %in% UC_genes, "UC", 
                         ifelse(network$Gene2 %in% EM_genes, "EM", 
                                ifelse(network$Gene2 %in% MM_genes, "MM", NA)))
  
  network <- network[!is.na(network$Age1),]
  network <- network[!is.na(network$Age2),]
  
  network$Edge_age <- paste(network$Age1, network$Age2, sep="_")
  
  network_list[[database]] <- network
}


degree_df <- vector()
for(database in databases){
  local_network <- network_list[[database]]
  network_graph <- graph.data.frame(local_network, directed=FALSE)
  degree_network <- degree(network_graph)
  local_degree_df <- data.frame(Genes=names(degree_network), Degree=degree_network, Database = database,
                                Age = genes_phy[match(names(degree_network), genes_phy$GeneID), "Age"])
  degree_df <- rbind(degree_df, local_degree_df)
}

degree_df$Age <- factor(degree_df$Age, levels=c("UC", "EM", "MM"))


##Divide by the median of each database

median_per_db <- aggregate(Degree ~ Database, degree_df, median)

degree_df$Median <- median_per_db[match(degree_df$Database, median_per_db$Database), "Degree"]
degree_df$Norm_degree <- degree_df$Degree/degree_df$Median

pdf("Supp2_Degree_of_PPI_databases.pdf",
    height=3.5)
g <- ggplot(degree_df, aes(x=Age, y =log10(Norm_degree)))+
  geom_boxplot(aes(fill=Age))+
  geom_hline(yintercept=0, color="grey")+
  geom_boxplot(aes(fill=Age))+
  xlab("Gene age")+
  ylab("Normalized Degree (log10)")+
  facet_grid(.~Database)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
print(g)
dev.off()


decreasing_degree_p <- vector()
for(database in databases){
  local_degree_df <- degree_df[degree_df$Database == database,]
  groups <- factor(local_degree_df$Age, levels=c("UC", "EM", "MM"), ordered=TRUE)
  decreasing_degree_p <- c(decreasing_degree_p, 
                           jonckheere.test(g=groups, x=local_degree_df$Degree, alternative="decreasing")$p.value)
}
names(decreasing_degree_p) <- databases

decreasing_degree_p <- p.adjust(decreasing_degree_p, method="BH")

GRN_degree_df <- degree_df[degree_df$Database == "GRN",]
degree_UC_genes <- GRN_degree_df[GRN_degree_df$Age == "UC", "Degree"]
degree_EM_genes <- GRN_degree_df[GRN_degree_df$Age == "EM", "Degree"]
degree_MM_genes <- GRN_degree_df[GRN_degree_df$Age == "MM", "Degree"]
wilcox.test(degree_EM_genes, degree_UC_genes, alternative="greater")
wilcox.test(degree_EM_genes, degree_MM_genes, alternative="greater")


##In-degree in RN
GRN_network <- network_list[["GRN"]]
in_degree <- table(GRN_network$Gene2)
in_degree_df <- data.frame(Genes=names(in_degree),
                           In_degree = as.vector(in_degree),
                           Age = genes_phy[match(names(in_degree), genes_phy$GeneID), "Age"])

in_degree_df$Age <- factor(in_degree_df$Age, levels=c("UC", "EM", "MM"))

pdf("Supp2_Indegree.pdf",
    height=3.5, width=4)
g <- ggplot(in_degree_df, aes(x=Age, y =log10(In_degree)))+
  geom_boxplot(aes(fill=Age))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(g)
dev.off()

in_degree_UC <- in_degree_df[in_degree_df$Age == "UC", "In_degree"]
in_degree_EM <- in_degree_df[in_degree_df$Age == "EM", "In_degree"]
in_degree_MM <- in_degree_df[in_degree_df$Age == "MM", "In_degree"]

wilcox.test(in_degree_EM, in_degree_UC, alternative="greater")
wilcox.test(in_degree_EM, in_degree_MM, alternative="greater")

aggregate(In_degree ~ Age, in_degree_df, mean)


summary_in_degree <- aggregate(In_degree ~ Age, in_degree_df, mean)
colnames(summary_in_degree) <- c("Age", "Mean_indegree")

pdf("Figure2_Number_incoming_edges.pdf",
    height=3, width=3)
g <- ggplot(summary_in_degree, aes(x=Age, y=Mean_indegree))+
  geom_bar(aes(fill=Age), stat='identity')+
  ylab("Average number of\nincoming edges")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(g)
dev.off()


out_degree <- table(GRN_network$Gene1)
out_degree_df <- data.frame(Genes=names(out_degree),
                            out_degree = as.vector(out_degree),
                            Age = genes_phy[match(names(out_degree), genes_phy$GeneID), "Age"])

out_degree_df$Age <- factor(out_degree_df$Age, levels=c("UC", "EM", "MM"))

quantile(out_degree_df$out_degree, 0.75)
  
out_degree_df$Group <- "GRN"

pdf("Supp2_Outdegree_distribution.pdf",
    height=2, width=4)
g <- ggplot(out_degree_df, aes(y=log10(out_degree), x=Group))+
  geom_boxplot()+
  #geom_hline(yintercept = 1, color="red")+
  geom_boxplot()+
  ylab("Out-degree (log10)")+
  xlab("")+
  coord_flip()+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
print(g)
dev.off()
  

n_out_10_UC <- sum(out_degree_df$out_degree[out_degree_df$Age == "UC"] >= 10)
n_out_10_EM <- sum(out_degree_df$out_degree[out_degree_df$Age == "EM"] >= 10)
n_out_10_MM <- sum(out_degree_df$out_degree[out_degree_df$Age == "MM"] >= 10)
total_out_10 <- sum(out_degree_df$out_degree >= 10)

n_out_10_UC/total_out_10
n_out_10_EM/total_out_10
n_out_10_MM/total_out_10

##Enrichment of EM genes are regulators
n_EM_genes <- length(EM_genes)
total_genes <- nrow(genes_phy)

total_regulators <- length(unique(GRN_network$Gene1))

n_EM_regulators <- sum(unique(GRN_network$Gene1) %in% EM_genes)
n_UC_regulators <- sum(unique(GRN_network$Gene1) %in% UC_genes)
n_MM_regulators <- sum(unique(GRN_network$Gene1) %in% MM_genes)

(n_UC_regulators/total_regulators)*100
(n_EM_regulators/total_regulators)*100
(n_MM_regulators/total_regulators)*100

fisher.test(cbind(c(n_EM_regulators, total_regulators),
                  c(n_EM_genes, total_genes)), alternative="greater")

percentage_regulators <- data.frame(Age=c("UC", "EM", "MM"),
                                    Percentage=c((n_UC_regulators/total_regulators)*100,
                                                 (n_EM_regulators/total_regulators)*100,
                                                 (n_MM_regulators/total_regulators)*100))

percentage_regulators$Age <- factor(percentage_regulators$Age, levels=c("UC", "EM", "MM"))


pdf("Figure2_Percentage_regulators_by_age.pdf",
    height=3, width=3)
g <- ggplot(percentage_regulators, aes(x=Age, y=Percentage))+
  geom_bar(aes(fill=Age), stat='identity')+
  ylab("Percentage regulators")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(g)
dev.off()

##Read in somatic mutations

mutations_df <- load_mutations_mutsig()

all_mut <- vector()
mut_df <- vector()

for(tumour in tumours){
  print(tumour)
  local_mut <- mutations_df[mutations_df$Tumour == tumour,]
  
  mut <- as.character(local_mut[, "Hugo_Symbol"])
  
  if(length(mut) != 0){
    mut_df <- rbind(mut_df, cbind(Tumour=tumour, Alt="Mut", Genes=mut))
    all_mut <- c(all_mut, mut)
  }
}

##Read in CNVs
##Using only focal ones
load("CNVs_curated_gistic.Rdata")
load("patients_with_CNV_info.Rdata")

CNV_df <- load_CNVs_gistic(only_focal="ANY")

all_amp <- vector()
all_del <- vector()
amp_df <- vector()
del_df <- vector()
for(tumour in tumours){
  n_pat_df <- CNV_df[CNV_df$Tumour == tumour,]
  amp_genes <- as.character(n_pat_df[n_pat_df$CNA == "Amplification", "Genes"])
  del_genes <- as.character(n_pat_df[n_pat_df$CNA == "Deletion", "Genes"])
  
  print(tumour)
  print(length(amp_genes))
  print(length(del_genes))
  
  all_amp <- c(all_amp, amp_genes)
  all_del <- c(all_del, del_genes)
  
  if(length(amp_genes) > 0){
    amp_df <- rbind(amp_df, cbind(tumour, "Amp", amp_genes))
    
  }
  if(length(del_genes) > 0){
    del_df <- rbind(del_df, cbind(tumour, "Del", del_genes))
  }
}

only_mut <- unique(all_mut)
all_mut_common <- only_mut

only_amp <- unique(all_amp)
only_del <- unique(all_del)
both_CNVs <- only_amp[only_amp %in% only_del]
only_amp <- only_amp[!(only_amp %in% both_CNVs)]
only_del <- only_del[!(only_del %in% both_CNVs)]


#Common alterations
all_amp_common <- names(table(all_amp))[table(all_amp)>=7]
all_del_common <- names(table(all_del))[table(all_del)>=7]
all_CNV_common_both <- all_amp_common[all_amp_common %in% all_del_common]
all_amp_common <- all_amp_common[!(all_amp_common %in% all_CNV_common_both)]
all_del_common <- all_del_common[!(all_del_common %in% all_CNV_common_both)]


all_mut_common <- names(table(all_mut))[table(all_mut)>=7]

genes_in_GRN <- unique(c(GRN_network$Gene1, GRN_network$Gene2))
all_mut_common <- all_mut_common[all_mut_common %in% genes_in_GRN]
all_amp_common <- all_amp_common[all_amp_common %in% genes_in_GRN]
all_del_common <- all_del_common[all_del_common %in% genes_in_GRN]


###Enrichment of mutations and CNVs as targets, regulators
regulators <- unique(GRN_network$Gene1)
targets <- unique(GRN_network$Gene2)
non_target_regulators <- regulators[!(regulators %in% targets)]
non_regulator_targets <- targets[!(targets %in% regulators)]


total_regulators <- length(regulators)
total_non_target_regulators <- length(non_target_regulators)

alt_genes <- rbind(mut_df, amp_df, del_df)
alt_genes <- data.frame(alt_genes)
alt_genes <- alt_genes[alt_genes$Genes %in% unique(c(GRN_network$Gene1, GRN_network$Gene2)),]


##Percentage of genes with alterations that are regulators, targets, etc
alt_genes$Alt_simplified <- ifelse(alt_genes$Alt %in% c("Mut"), "Point",
                                   ifelse(alt_genes$Alt %in% c("Amp", "Del"), "CNV", NA))

fraction_gene_class_binary <- vector()
number_genes <- vector()
for(tumour in tumours){
  local_alt <- alt_genes[alt_genes$Tumour == tumour,]
  for(alt in unique(alt_genes$Alt_simplified)){
    local_alt2 <- local_alt[local_alt$Alt_simplified == alt,]
    if(nrow(local_alt2) >= 3){
      number_genes <- c(number_genes, nrow(local_alt2))
      per_regulators <- sum(unique(local_alt2$Genes) %in% regulators)/length(unique(local_alt2$Genes))
      per_non_r_targets <- sum(unique(local_alt2$Genes) %in% non_regulator_targets)/length(unique(local_alt2$Genes))
      
      fraction_gene_class_binary <- rbind(fraction_gene_class_binary,
                                          c(Tumour=tumour, Alteration=alt, 
                                            Per_regulators=per_regulators, 
                                            Per_non_r_targets=per_non_r_targets))
    }
  }
}

fraction_gene_class_binary <- as.data.frame(fraction_gene_class_binary)
fraction_gene_class_binary$Per_regulators <- as.numeric(as.character(fraction_gene_class_binary$Per_regulators))
fraction_gene_class_binary$Per_non_r_targets <- as.numeric(as.character(fraction_gene_class_binary$Per_non_r_targets))


fraction_gene_class_binary_melt <- melt(fraction_gene_class_binary[,c("Tumour", "Alteration", "Per_regulators", "Per_non_r_targets")])

colnames(fraction_gene_class_binary_melt) <- c("Tumour", "Alteration", "Gene_role", "Fraction")

fraction_gene_class_binary_melt$Gene_role <- as.character(fraction_gene_class_binary_melt$Gene_role)
fraction_gene_class_binary_melt$Gene_role[fraction_gene_class_binary_melt$Gene_role == "Per_non_r_targets"] <- "Target"
fraction_gene_class_binary_melt$Gene_role[fraction_gene_class_binary_melt$Gene_role == "Per_regulators"] <- "Regulator"

fraction_gene_class_binary_melt$Alteration <- factor(fraction_gene_class_binary_melt$Alteration, levels=c("Point", "CNV"))

#Looking at the genes that are recurrent across tumours
sum(c(all_mut_common) %in% regulators)/length(c(all_mut_common))*100
sum(c(all_amp_common, all_del_common) %in% regulators)/length(c(all_amp_common, all_del_common))*100

sum(c(all_mut_common) %in% non_regulator_targets)/length(c(all_mut_common))*100
sum(c(all_amp_common, all_del_common) %in% non_regulator_targets)/length(c(all_amp_common, all_del_common))*100


fraction_gene_class_binary_melt$Tumour <- factor(fraction_gene_class_binary_melt$Tumour,
                                                 levels=tumours)
fraction_gene_class_binary_melt$Alteration <- factor(fraction_gene_class_binary_melt$Alteration,
                                                     levels=c("Point", "CNV"))

fraction_gene_class_binary_mean <- aggregate(Fraction ~ Alteration+Gene_role, fraction_gene_class_binary_melt, mean)

fraction_gene_class_binary_mean$Alteration <- factor(fraction_gene_class_binary_mean$Alteration,
                                                     levels=c("Point", "CNV"))

pdf("Figure2_Fraction_mutated_regulators.pdf",
    width=3.5, height=3)
g <- ggplot(subset(fraction_gene_class_binary_melt, Gene_role=="Regulator"), aes(x=Alteration, y =Fraction))+
  geom_point(aes(colour=Tumour))+
  #coord_cartesian(ylim=c(0, 0.37))+
  geom_bar(data=subset(fraction_gene_class_binary_mean, Gene_role=="Regulator"), aes(x=Alteration, y=Fraction), stat='identity', fill="grey92", color="black")+
  geom_path(aes(group=Alteration))+
  geom_point(aes(colour=Tumour), size=2)+
  ggtitle("Regulator")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())
print(g)
dev.off()

##In how many tumour types is the fraction of point mutation greater than the fraction with CNV
greater_reg_point_fraction <- vector()
for(tumour in tumours){
  local_temp <- subset(fraction_gene_class_binary_melt, Gene_role=="Regulator" & Tumour == tumour)
  print(tumour)
  print(nrow(local_temp))
  greater_reg_point_fraction <- c(greater_reg_point_fraction,
                                  local_temp[local_temp$Alteration == "Point","Fraction"] > local_temp[local_temp$Alteration == "CNV","Fraction"])
  
}
sum(greater_reg_point_fraction)/length(greater_reg_point_fraction)###80.77% of tumours


pdf("Figure2_Fraction_mutated_targets.pdf",
    width=3.5, height=3)
g <- ggplot(subset(fraction_gene_class_binary_melt, Gene_role=="Target"), aes(x=Alteration, y =Fraction))+
  geom_point(aes(colour=Tumour))+
  #coord_cartesian(ylim=c(0.6,1))+
  geom_bar(data=subset(fraction_gene_class_binary_mean, Gene_role=="Target"), 
           aes(x=Alteration, y=Fraction), stat='identity', fill="grey92", color="black")+
  geom_path(aes(group=Alteration))+
  #geom_path(aes(group=Tumour))+
  geom_point(aes(colour=Tumour), size=2)+
  #ylim(0.5, 1)+
  scale_y_continuous(position = "right")+
  ggtitle("Target")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())
print(g)
dev.off()

greater_tar_point_fraction <- vector()
for(tumour in tumours){
  local_temp <- subset(fraction_gene_class_binary_melt, Gene_role=="Target" & Tumour == tumour)
  greater_tar_point_fraction <- c(greater_tar_point_fraction,
                                  local_temp[local_temp$Alteration == "Point","Fraction"] < local_temp[local_temp$Alteration == "CNV","Fraction"])
  
}
sum(greater_tar_point_fraction)/length(greater_tar_point_fraction)###80.77% of tumours



aggregate(Fraction ~ Alteration+Gene_role, fraction_gene_class_binary_melt, mean)

fraction_gene_class_binary_melt_reg <- fraction_gene_class_binary_melt[fraction_gene_class_binary_melt$Gene_role == "Regulator",]
fraction_gene_class_binary_melt_tar <- fraction_gene_class_binary_melt[fraction_gene_class_binary_melt$Gene_role == "Target",]


wilcox.test(fraction_gene_class_binary_melt_reg$Fraction[fraction_gene_class_binary_melt_reg$Alteration == "Point"],
            fraction_gene_class_binary_melt_reg$Fraction[fraction_gene_class_binary_melt_reg$Alteration == "CNV"], alternative="greater")

wilcox.test(fraction_gene_class_binary_melt_tar$Fraction[fraction_gene_class_binary_melt_tar$Alteration == "CNV"],
            fraction_gene_class_binary_melt_tar$Fraction[fraction_gene_class_binary_melt_tar$Alteration == "Point"], alternative="greater")





##Similar to the above but keeping the classification of each alteration
fraction_gene_class_all <- vector()
for(tumour in tumours){
  local_alt <- alt_genes[alt_genes$Tumour == tumour,]
  for(alt in unique(alt_genes$Alt)){
    local_alt2 <- local_alt[local_alt$Alt == alt,]
    if(nrow(local_alt2) >=3){
      per_regulators <- sum(unique(local_alt2$Genes) %in% regulators)/length(unique(local_alt2$Genes))
      per_non_r_targets <- sum(unique(local_alt2$Genes) %in% non_regulator_targets)/length(unique(local_alt2$Genes))
      
      fraction_gene_class_all <- rbind(fraction_gene_class_all,
                                       c(Tumour=tumour, Alteration=alt, 
                                         Per_regulators=per_regulators, 
                                         Per_non_r_targets=per_non_r_targets))
      
    }
  }
}

fraction_gene_class_all <- as.data.frame(fraction_gene_class_all)
fraction_gene_class_all$Per_regulators <- as.numeric(as.character(fraction_gene_class_all$Per_regulators))
fraction_gene_class_all$Per_non_r_targets <- as.numeric(as.character(fraction_gene_class_all$Per_non_r_targets))

fraction_gene_class_all_melt <- melt(fraction_gene_class_all[,c("Tumour", "Alteration", "Per_regulators", "Per_non_r_targets")])

colnames(fraction_gene_class_all_melt) <- c("Tumour", "Alteration", "Gene_role", "Fraction")

fraction_gene_class_all_melt$Gene_role <- as.character(fraction_gene_class_all_melt$Gene_role)
fraction_gene_class_all_melt$Gene_role[fraction_gene_class_all_melt$Gene_role == "Per_non_r_targets"] <- "Target"
fraction_gene_class_all_melt$Gene_role[fraction_gene_class_all_melt$Gene_role == "Per_regulators"] <- "Regulator"

fraction_gene_class_all_melt$Alteration <- factor(fraction_gene_class_all_melt$Alteration, levels=c("Mut", "Amp", "Del"))


pdf("Supp2_Fraction_mutated_regulators.pdf",
    width=4, height=3)
g <- ggplot(fraction_gene_class_all_melt[fraction_gene_class_all_melt$Gene_role == "Regulator",], aes(x=Alteration, y =Fraction))+
  geom_boxplot(aes(fill=Alteration))+
  ggtitle("Regulators")+
  xlab("")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())
print(g)
dev.off()


pdf("Supp2_Fraction_mutated_targets.pdf",
    width=4, height=3)
g <- ggplot(fraction_gene_class_all_melt[fraction_gene_class_all_melt$Gene_role == "Target",], aes(x=Alteration, y =Fraction))+
  geom_boxplot(aes(fill=Alteration))+
  scale_y_continuous(position = "right")+
  ggtitle("Targets")+
  xlab("")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())
print(g)
dev.off()

aggregate(Fraction ~ Alteration+Gene_role, fraction_gene_class_all_melt, mean)

fraction_gene_class_all_melt_reg <- fraction_gene_class_all_melt[fraction_gene_class_all_melt$Gene_role == "Regulator",]

wilcox.test(fraction_gene_class_all_melt_reg$Fraction[fraction_gene_class_all_melt_reg$Alteration == "Mut"],
            fraction_gene_class_all_melt_reg$Fraction[fraction_gene_class_all_melt_reg$Alteration %in% c("Amp", "Del")], alternative = "greater")



fraction_gene_class_all_melt_tar <- fraction_gene_class_all_melt[fraction_gene_class_all_melt$Gene_role == "Target",]

wilcox.test(fraction_gene_class_all_melt_tar$Fraction[fraction_gene_class_all_melt_tar$Alteration == "Amp"],
            fraction_gene_class_all_melt_tar$Fraction[fraction_gene_class_all_melt_tar$Alteration %in% c("Mut")],  alternative = "greater")

wilcox.test(fraction_gene_class_all_melt_tar$Fraction[fraction_gene_class_all_melt_tar$Alteration == "Del"],
            fraction_gene_class_all_melt_tar$Fraction[fraction_gene_class_all_melt_tar$Alteration %in% c("Mut")],  alternative = "greater")





###Calculation of in and out degree ratio per tumour type
out_degree_df$Database <- "GRN"
colnames(out_degree_df)[2] <- "Degree"
in_degree_df$Database <- "GRN"
colnames(in_degree_df)[2] <- "Degree"

#First calculating the ratio for each gene, then averaging
median_ratio <- vector()
p_EM <- vector()
median_overall <- vector()
#options(warn=2)
for(tumour in tumours){
  local_mut <- mut_df[mut_df[,1] == tumour,]
  local_amp <- amp_df[amp_df[,1] == tumour,]
  local_del <- del_df[del_df[,1] == tumour,]
  
  local_out_degree_df <- add_mut_CNV_to_degree_df(out_degree_df, local_mut, local_mut, local_amp, local_del)
  local_in_degree_df <- add_mut_CNV_to_degree_df(in_degree_df, local_mut, local_mut, local_amp, local_del)
  
  local_out_degree_df$Alt_simplified <- ifelse(local_out_degree_df$Alt %in% c("Miss", "LoF"), "Point",
                                               ifelse(local_out_degree_df$Alt %in% c("Amp", "Del"), "CNV", "WT"))
  local_in_degree_df$Alt_simplified <- ifelse(local_in_degree_df$Alt %in% c("Miss", "LoF"), "Point",
                                              ifelse(local_in_degree_df$Alt %in% c("Amp", "Del"), "CNV", "WT"))
  colnames(local_out_degree_df)[2] <- "Out_degree"
  local_degree_df <- local_out_degree_df
  
  local_degree_df$Alt[local_degree_df$Alt == "Miss"] <- "Mut"
  
  #Only considering genes with dual roles
  local_degree_df <- data.frame(local_degree_df, In_degree=local_in_degree_df[match(local_degree_df$Genes, local_in_degree_df$Genes), "Degree"])
  
  local_degree_df$Ratio_out_in <- local_degree_df$Out_degree/local_degree_df$In_degree
  
  local_median_ratio <- aggregate(Ratio_out_in ~ Age+Alt_simplified, local_degree_df, median)  ##It was median
  #local_median_ratio <- local_degree_df

  local_median_ratio$Tumour <- tumour
  
  median_ratio <- rbind(median_ratio, local_median_ratio)
}


upper_quantile <- aggregate(Ratio_out_in ~ Age+Alt_simplified, median_ratio, quantile, 0.75)
lower_quantile <- aggregate(Ratio_out_in ~ Age+Alt_simplified, median_ratio, quantile, 0.25)
median_values <- aggregate(Ratio_out_in ~ Age+Alt_simplified, median_ratio, median)

colnames(median_values)[3] <- "Median_values"
median_values$Upper <- upper_quantile$Ratio_out_in
median_values$Lower <- lower_quantile$Ratio_out_in
median_values$Full_name <- paste(median_values$Age, median_values$Alt_simplified, sep="_")

median_values$Age <- factor(median_values$Age, levels=c("MM", "EM", "UC"))


pdf("Figure2_Ratio_out_in_degree.pdf",
    width=7, height=3)
g <- ggplot(median_values, aes(x=Alt_simplified, y = log2(Median_values)))+
  geom_point(aes(colour=Age, shape=Alt_simplified), size=4,
             position = position_dodge(width = 0.75))+
  geom_linerange(aes(ymin=log2(Lower), ymax=log2(Upper), colour=Age),
                 position = position_dodge(width = 0.75), size=0.75)+
  coord_flip()+
  #geom_hline(yintercept = 1, size=0.5, colour="grey")+
  geom_point(aes(colour=Age, shape=Alt_simplified), size=4,
             position = position_dodge(width = 0.75))+
  geom_linerange(aes(ymin=log2(Lower), ymax=log2(Upper), colour=Age),
                 position = position_dodge(width = 0.75), size=0.75)+
  scale_colour_manual(values=c(UC=gg_color_hue(3)[1], EM=gg_color_hue(3)[2], MM=gg_color_hue(3)[3]))+
  ylab("Ratio out-degree/in-degree")+
  xlab("")+
  facet_wrap("Alt_simplified", scales='free_y', strip.position="left", nrow=2)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(colour=NA, fill=NA),
        panel.background = element_rect(fill = NA))
print(g)
dev.off()



ggplot(median_values, aes(x=Alt_simplified, y = Median_values))+
  geom_point(aes(colour=Age, shape=Alt_simplified), size=4,
             position = position_dodge(width = 0.75))+
  geom_linerange(aes(ymin=Lower, ymax=Upper, colour=Age),
                 position = position_dodge(width = 0.75), size=0.75)+
  #ylim(0.25,1.75)+
  coord_flip(ylim=c(0.3,8.5))+
  #geom_hline(yintercept = 1, size=0.5, colour="grey")+
  geom_point(aes(colour=Age, shape=Alt_simplified), size=4,
             position = position_dodge(width = 0.75))+
  geom_linerange(aes(ymin=Lower, ymax=Upper, colour=Age),
                 position = position_dodge(width = 0.75), size=0.75)+
  scale_colour_manual(values=c(UC=gg_color_hue(3)[1], EM=gg_color_hue(3)[2], MM=gg_color_hue(3)[3]))+
  ylab("Ratio out-degree/in-degree")+
  xlab("")+
  facet_wrap("Alt_simplified", scales='free_y', strip.position="left", nrow=2)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(colour=NA, fill=NA),
        panel.background = element_rect(fill = NA))





median_ratio[intersect(which(median_ratio$Age == "EM"),
                       which(median_ratio$Alt_simplified == "Point")),]

median_ratio[intersect(which(median_ratio$Age == "EM"),
                       which(median_ratio$Alt_simplified == "CNV")),]


median_values_print <- median_values[,c("Age","Alt_simplified","Median_values","Upper","Lower")]
colnames(median_values_print) <- c("Age", "Mutation", "Median", "Upperquantile", "Lowerquantile")
write.table(median_values_print, file="Supp_table_Ratio_out_in.txt",
            sep="\t", quote=FALSE, row.names=FALSE)
