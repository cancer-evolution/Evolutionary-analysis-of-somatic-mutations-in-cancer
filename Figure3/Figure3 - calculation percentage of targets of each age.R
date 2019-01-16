library(igraph)
library(ggplot2)
library(reshape2)
library(gridExtra)

#Percentage of target of each age

source("helper_functions.R")

genes_phy <- read.csv("gene_phylostrata.txt")  ###File with gene name, entrez and phylostratum as column
genes_phy$Age <- ifelse(genes_phy$Phylostrata %in% 1:3, "UC",
                        ifelse(genes_phy$Phylostrata %in% 4:9, "EM",
                               ifelse(genes_phy$Phylostrata %in% 10:16, "MM", NA)))
UC_genes <- as.character(genes_phy[genes_phy$Age == "UC", "GeneID"])
EM_genes <- as.character(genes_phy[genes_phy$Age == "EM", "GeneID"])
MM_genes <- as.character(genes_phy[genes_phy$Age == "MM", "GeneID"])



##Determine the downstream targets of each regulator
network <- load_PathwayCommons_GRN()

network$Edge_age <- paste(network$Age1, network$Age2, sep="_")

source_nodes <- unique(as.character(network[,1]))

degree_sources <- table(as.character(network[,1]))

percentage_UC_per_gene_directed <- lapply(1:length(source_nodes), function(i){
  local_gene <- source_nodes[i]
  local_network <- network[network[,1] == local_gene,]
  
  other_genes <- as.character(unique(local_network[,2]))
  other_genes <- other_genes[other_genes != local_gene]
  
  per_UC <- (sum(other_genes %in% UC_genes)/length(other_genes))*100
  per_EM <- (sum(other_genes %in% EM_genes)/length(other_genes))*100
  per_MM <- (sum(other_genes %in% MM_genes)/length(other_genes))*100
  
  return(c(per_UC, per_EM, per_MM))
})

percentage_UC_per_gene_directed <- as.data.frame(do.call('rbind', percentage_UC_per_gene_directed))
rownames(percentage_UC_per_gene_directed) <- source_nodes
colnames(percentage_UC_per_gene_directed) <- c("Per_UC", "Per_EM", "Per_MM")

percentage_UC_per_gene_directed$Gene_age <- ifelse(rownames(percentage_UC_per_gene_directed) %in% UC_genes, "UC",
                                                   ifelse(rownames(percentage_UC_per_gene_directed) %in% EM_genes, "EM", "MM"))

percentage_UC_per_gene_directed$Gene_age <- factor(percentage_UC_per_gene_directed$Gene_age,
                                                   levels=c("UC", "EM", "MM"))

percentage_UC_per_gene_directed$Degree <- degree_sources[match(rownames(percentage_UC_per_gene_directed), 
                                                               names(degree_sources))]
percentage_UC_per_gene_directed$Degree <- as.vector(percentage_UC_per_gene_directed$Degree)


save(percentage_UC_per_gene_directed,
     file="percentage_UC_per_gene_directed_PC2.Rdata")
