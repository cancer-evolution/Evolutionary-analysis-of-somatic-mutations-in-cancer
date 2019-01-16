##Replicationg in other GRN

library(ggplot2)
library(igraph)
library(reshape2)
library(gridExtra)
library(readr)
library(clinfun)
library(ggrepel)
library(ggsci)

source("helper_functions.R")

tumours <- c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "ESCA", "GBM", 
             "HNSC", "KICH", "KIRC", "KIRP", "LGG",  "LIHC", "LUAD", "LUSC", 
             "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", 
             "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM")

genes_phy <- read.csv("gene_phylostrata.txt")  ###File with gene name, entrez and phylostratum as column
genes_phy$Age <- ifelse(genes_phy$Phylostrata %in% 1:3, "UC",
                        ifelse(genes_phy$Phylostrata %in% 4:9, "EM",
                               ifelse(genes_phy$Phylostrata %in% 10:16, "MM", NA)))

UC_genes <- as.character(genes_phy[which(genes_phy$Phylostrata %in% 1:3),1])
EM_genes <- as.character(genes_phy[which(genes_phy$Phylostrata %in% 4:9),1])
MM_genes <- as.character(genes_phy[which(genes_phy$Phylostrata %in% 10:16),1])

load("CNVs_curated2.Rdata")
load("patients_with_CNV_info.Rdata")

load("CNVs_above_fraction_0.25.Rdata")


#TTRUST v2
##From literature mining
#http://www.grnpedia.org/trrust/
#last updated 2017.05.25

network <- load_TRRUST()

nrow(network)
length(unique(network$Gene1))
length(unique(network$Gene2))

reproduce_results_other_databases_fig_1_2(network, "Trrust")
reproduce_results_figure_3C("Trrust")
reproduce_results_figure_4A(network, "Trrust")

##RegNetwork
#http://www.regnetworkweb.org/home.jsp

network <- load_RegNetwork()
nrow(network)
length(unique(network$Gene1))
length(unique(network$Gene2))
reproduce_results_other_databases_fig_1_2(network, "RegNetwork")
reproduce_results_figure_3C("RegNetwork")
reproduce_results_figure_4A(network, "RegNetwork")

##PathwayCommons

network <- load_PathwayCommons_GRN()
nrow(network)
length(unique(network$Gene1))
length(unique(network$Gene2))
reproduce_results_other_databases_fig_1_2(network, "PathwayCommons")
reproduce_results_figure_3C("PathwayCommons")
#reproduce_results_figure_4A(network, "PathwayCommons")

