##Classify regulators using data from the various databases

##Degree and neighbourhood of directed networks

library(igraph)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(pheatmap)

source("helper_functions.R")

genes_phy <- read.csv("gene_phylostrata.txt")  ###File with gene name, entrez and phylostratum as column
UC_genes <- as.character(genes_phy[which(genes_phy$Phylostrata %in% 1:3),1])
EM_genes <- as.character(genes_phy[which(genes_phy$Phylostrata %in% 4:9),1])
MM_genes <- as.character(genes_phy[which(genes_phy$Phylostrata %in% 10:16),1])


network <- load_TTRUST(); database <- "TTRUST"

network <- load_RegNetwork(); database <- "RegNetwork"

network <- load_PathwayCommons_GRN(); database <- "PathwayCommons"


source_nodes <- unique(as.character(network[,1]))

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

percentage_UC_per_gene_directed$Age <- ifelse(rownames(percentage_UC_per_gene_directed) %in% UC_genes, "UC",
                                                   ifelse(rownames(percentage_UC_per_gene_directed) %in% EM_genes, "EM", "MM"))

percentage_UC_per_gene_directed$Age <- factor(percentage_UC_per_gene_directed$Age,
                                                   levels=c("UC", "EM", "MM"))


percentage_UC_per_gene_directed$Database <- "RN"
percentage_UC_per_gene_directed$Genes <- rownames(percentage_UC_per_gene_directed)


save(percentage_UC_per_gene_directed,
     file=paste("percentage_UC_per_gene_directed_", database, ".Rdata", sep=""))


##Classify targets into quantiles based on age

##Separate by %UC and %EM

#At least >2/3 UC to be UC, at least >2/3 EM to be EM, and to be mixed at least 10% of each

per_UC_EM_downstream <- percentage_UC_per_gene_directed[,c("Genes", "Per_UC", "Per_EM", "Per_MM", "Age")]

per_UC_EM_downstream$Regulator_class <- ifelse(per_UC_EM_downstream$Per_UC >= 200/3, "UC_downstream",
                                               ifelse(per_UC_EM_downstream$Per_EM >= 200/3, "EM_downstream", 
                                                      ifelse(per_UC_EM_downstream$Per_MM >= 200/3, "MM_downstream",
                                                             ifelse(per_UC_EM_downstream$Per_UC >= 10 & per_UC_EM_downstream$Per_EM >= 10 & per_UC_EM_downstream$Per_MM < 10, "Mixed_downstream", "None"))))

regulator_classification <- per_UC_EM_downstream[,c("Per_UC", "Per_EM", "Per_MM", "Regulator_class")]
regulator_classification$Regulator <- rownames(regulator_classification)

save(regulator_classification, 
     file=paste("regulator_classification_", database, ".Rdata", sep=""))



#####All regulators
all_regulators <- vector()
for(database in c("TTRUST", "RegNetwork", "PathwayCommons")){
  load(paste("regulator_classification_", database, ".Rdata", sep=""))
  
  all_regulators <- unique(c(all_regulators, regulator_classification$Regulator))  
}

##Regulators in PathwayCommons
load("regulator_classification_PathwayCommons.Rdata")
reg_class_PC <- regulator_classification[,c("Regulator", "Regulator_class")]
  
reg_class_PC <- reg_class_PC[reg_class_PC$Regulator_class %in% c("UC_downstream", "Mixed_downstream", "EM_downstream"),]
reg_class_PC <- reg_class_PC[order(reg_class_PC$Regulator_class),]


regulators_df <- matrix(0, nrow=length(reg_class_PC$Regulator), ncol=5)
colnames(regulators_df) <- c("UC_downstream", "Mixed_downstream", "EM_downstream", "MM_downstream", "None")
rownames(regulators_df) <- reg_class_PC$Regulator

##Number databases concordant

for(reg in all_regulators){
  for(database in c("TTRUST", "RegNetwork", "PathwayCommons")){
    load(paste("regulator_classification_", database, ".Rdata", sep=""))
    temp <- regulator_classification[regulator_classification$Regulator == reg,]
    if(nrow(temp) > 0){
      regulators_df[rownames(regulators_df) == reg,
                    colnames(regulators_df) == temp$Regulator_class] <- regulators_df[rownames(regulators_df) == reg,
                                                                                      colnames(regulators_df) == temp$Regulator_class]+1
    }
    
  }
}


annotation <- reg_class_PC[match(rownames(regulators_df), reg_class_PC$Regulator), ]
annotation$Regulator <- NULL
annotation <- as.data.frame(apply(annotation, 2, rev))


ann_colors <- c("#0073C2FF", "#EFC000FF", "#868686FF")
names(ann_colors) <- c("UC_downstream", "Mixed_downstream", "EM_downstream")
ann_colors <- list(Regulator_class=ann_colors)

pdf("Comparison_regulator_class.pdf",
    width=4, height=9)
pheatmap(regulators_df[,c("UC_downstream", "Mixed_downstream", "EM_downstream")], 
         cluster_rows=FALSE, cluster_cols=FALSE,
         color=c("white", "brown", "red", "green"),
         annotation_colors=ann_colors,
         annotation_row=annotation, show_rownames=FALSE)
dev.off()


regulators_df_more_1 <- regulators_df[, c("UC_downstream", "Mixed_downstream", "EM_downstream")]
temp <- apply(regulators_df_more_1, 1, sum)
regulators_df_more_1 <- regulators_df_more_1[temp > 2,]
regulators_df_more_1 <- as.data.frame(regulators_df_more_1)
regulators_df_more_1$PC <- reg_class_PC[match(rownames(regulators_df_more_1), reg_class_PC$Regulator), "Regulator_class"]

reg_more_1_mixed <- regulators_df_more_1[regulators_df_more_1$PC == "Mixed_downstream",]
sum(reg_more_1_mixed$Mixed_downstream >= 2)/nrow(reg_more_1_mixed)

reg_more_1_UC <- regulators_df_more_1[regulators_df_more_1$PC == "UC_downstream",]
sum(reg_more_1_UC$UC_downstream >= 2)/nrow(reg_more_1_UC)

reg_more_1_EM <- regulators_df_more_1[regulators_df_more_1$PC == "EM_downstream",]
sum(reg_more_1_EM$EM_downstream >= 2)/nrow(reg_more_1_EM)

