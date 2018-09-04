library(ggplot2)
library(reshape2)
library(readr)
library(gridExtra)
library(ggrepel)
library(ggsci)

load("regulator_classification.Rdata")

tumours <- c("LUAD", "LUSC", "BRCA", "PRAD", "LIHC", "COAD", "STAD")

mixed_regulators <- regulator_classification[regulator_classification$Regulator_class == "Mixed_downstream",]
UC_regulators <- regulator_classification[regulator_classification$Regulator_class == "UC_downstream",]
EM_regulators <- regulator_classification[regulator_classification$Regulator_class == "EM_downstream",]


genes_phy <- read.csv("gene_phylostrata.txt")  ###File with gene name, entrez and phylostratum as column
genes_phy$Age <- ifelse(genes_phy$Phylostrata %in% 1:3, "UC",
                        ifelse(genes_phy$Phylostrata %in% 4:9, "EM",
                               ifelse(genes_phy$Phylostrata %in% 10:16, "MM", NA)))
EM_genes <- as.character(genes_phy[genes_phy$Age == "EM", "GeneID"])

mixed_regulators <- mixed_regulators[mixed_regulators$Regulator %in% EM_genes,]

cancer_census <- read.delim("Census_allTue Dec 12 03_33_11 2017.csv", sep=",") #File from the Cancer Census database

mixed_regulators <- mixed_regulators[!(mixed_regulators$Regulator %in% cancer_census$Gene.Symbol),]

dependency <- read.csv("gene_dependency avana_public_18Q1, v3.csv", sep=",")  #File from Achilles

dependency_t <- t(dependency)
colnames(dependency_t) <- dependency_t[1,]
dependency_t <- dependency_t[-1,]
rownames(dependency_t) <- gsub("\\...*","",rownames(dependency_t))

dependency_melt <- melt(dependency_t)
colnames(dependency_melt) <- c("Gene", "Cell_line", "Probability")
dependency_melt$Probability <- as.numeric(as.character(dependency_melt$Probability))

regulator_classification$Age <- genes_phy[match(regulator_classification$Regulator, genes_phy$GeneID), "Age"]

blood_cell_lines <- dependency_melt[grep("HAEMATOPOIETIC_AND_LYMPHOID_TISSUE", dependency_melt$Cell_line),"Cell_line"]
blood_cell_lines <- unique(blood_cell_lines)

dependency_melt <- dependency_melt[!(dependency_melt$Cell_line %in% blood_cell_lines),]

dependency_melt$Regulator_class <- regulator_classification[match(dependency_melt$Gene, regulator_classification$Regulator), "Regulator_class"]
dependency_melt$Regulator <- ifelse(is.na(dependency_melt$Regulator_class), "Non_regulator", "Regulator")
dependency_melt$Age <- genes_phy[match(dependency_melt$Gene, genes_phy$GeneID), "Age"]

dependency_melt$Age <- factor(dependency_melt$Age,
                              levels=c("UC", "EM", "MM"))

dependency_melt_regulator <- dependency_melt[dependency_melt$Regulator == "Regulator",]

dependency_melt_regulator <- dependency_melt_regulator[dependency_melt_regulator$Regulator_class != "None",]
dependency_melt_regulator <- dependency_melt_regulator[dependency_melt_regulator$Regulator_class != "MM_downstream",]


dependency_melt_regulator$Regulator_class <- factor(dependency_melt_regulator$Regulator_class,
                                                    levels=c("UC_downstream",
                                                             "Mixed_downstream",
                                                             "EM_downstream"))
pdf("Supp5_Distribution_dependency.pdf",
    height=3.2, width=9)
g <-ggplot(dependency_melt_regulator, aes(x=Probability))+
  geom_density(aes(fill=Regulator_class))+
  scale_fill_jco()+
  facet_grid(.~Regulator_class)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank())
print(g)
dev.off()

pdf("Supp5_Distribution_dependency_inset.pdf",
    height=3.2, width=9)
g <-ggplot(dependency_melt_regulator, aes(x=Probability))+
  geom_density(aes(fill=Regulator_class))+
  xlim(0.75,1)+
  scale_fill_jco()+
  facet_grid(.~Regulator_class)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank())
print(g)
dev.off()



fraction_cell_lines <- vector()
for(cutoff in c(0.75, 0.90, 0.95, 0.99)){
  or_dependency_lines <- vector()
  for(line in unique(dependency_melt_regulator$Cell_line)){
    local_dependency <- dependency_melt_regulator[dependency_melt_regulator$Cell_line == line,]
    high <- local_dependency[local_dependency$Probability >= cutoff,]
    
    high_mixed <- high[high$Regulator_class == "Mixed_downstream",]
    high_mixed_EM <- high_mixed[high_mixed$Age == "EM",]
    high_mixed_UC <- high_mixed[high_mixed$Age == "UC",]
    
    local_dependency_mixed <- local_dependency[local_dependency$Regulator_class == "Mixed_downstream",]
    local_dependency_mixed_EM <- local_dependency_mixed[local_dependency_mixed$Age == "EM",]
    local_dependency_mixed_UC <- local_dependency_mixed[local_dependency_mixed$Age == "UC",]
    
    or1 <- fisher.test(cbind(c(sum(high$Regulator_class == "Mixed_downstream"),
                               length(high$Regulator_class)),
                             c(sum(local_dependency$Regulator_class == "Mixed_downstream"),
                               length(local_dependency$Regulator_class))), alternative="greater")$estimate
    
    or2 <- fisher.test(cbind(c(sum(high$Regulator_class == "UC_downstream"),
                               length(high$Regulator_class)),
                             c(sum(local_dependency$Regulator_class == "UC_downstream"),
                               length(local_dependency$Regulator_class))), alternative="greater")$estimate
    
    or3 <- fisher.test(cbind(c(sum(high$Regulator_class == "EM_downstream"),
                               length(high$Regulator_class)),
                             c(sum(local_dependency$Regulator_class == "EM_downstream"),
                               length(local_dependency$Regulator_class))), alternative="greater")$estimate
    
    or_dependency_lines <- rbind(or_dependency_lines, c(line, or1, or2, or3))
  }
  or_dependency_lines <- as.data.frame(or_dependency_lines)
  colnames(or_dependency_lines) <- c("Cell_line", "OR_Mixed_downstream", "OR_UC_downstream", "OR_EM_downstream")
  or_dependency_lines$OR_Mixed_downstream <- as.numeric(as.character(or_dependency_lines$OR_Mixed_downstream))
  or_dependency_lines$OR_UC_downstream <- as.numeric(as.character(or_dependency_lines$OR_UC_downstream))
  or_dependency_lines$OR_EM_downstream <- as.numeric(as.character(or_dependency_lines$OR_EM_downstream))
  
  or_dependency_lines_melt <- melt(or_dependency_lines)
  
  or_dependency_lines_melt$variable <- gsub("OR_", "", or_dependency_lines_melt$variable)
  colnames(or_dependency_lines_melt)[2:3] <- c("Regulator_class", "Odds_ratio")
  or_dependency_lines_melt$Regulator_class <- factor(or_dependency_lines_melt$Regulator_class,
                                                     levels=c("UC_downstream", "Mixed_downstream", "EM_downstream"))
  
  g <- ggplot(or_dependency_lines_melt, aes(x=Cell_line, y =Odds_ratio))+
    geom_point(aes(colour=Regulator_class), size=2.5)+
    scale_colour_jco()+
    geom_hline(yintercept=1)+
    ylab("Odds ratio")+
    xlab("Cell lines")+
    theme_bw()+
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  print(g)
 
  g <- ggplot(or_dependency_lines_melt, aes(x=Cell_line, y =Odds_ratio))+
    geom_point(aes(colour=Regulator_class), size=1.5)+
    scale_colour_jco()+
    geom_hline(yintercept=1)+
    ylab("Odds ratio")+
    xlab("Cell lines")+
    theme_bw()+
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  print(g)
  
  c1 <- sum(or_dependency_lines$OR_Mixed_downstream > 1)/length(or_dependency_lines$OR_Mixed_downstream)
  c2 <- sum(or_dependency_lines$OR_UC_downstream > 1)/length(or_dependency_lines$OR_UC_downstream)
  c3 <- sum(or_dependency_lines$OR_EM_downstream > 1)/length(or_dependency_lines$OR_EM_downstream)
  
  fraction_cell_lines <- rbind(fraction_cell_lines,
                               c(cutoff, c1, c2, c3))
}

colnames(fraction_cell_lines) <- c("Probability", "OR_Mixed_downstream",
                                   "OR_UC_downstream", "OR_EM_downstream")
