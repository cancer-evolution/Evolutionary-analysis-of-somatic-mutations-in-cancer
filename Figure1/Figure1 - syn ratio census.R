##Syn ratio

library(ggplot2)
library(readr)

source('helper_functions.R')

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


tumours <- c("LUAD", "LUSC", "BRCA", "PRAD", "LIHC", "COAD", "STAD")


mutations_df <- load_mutations()


cancer_census <- read.delim("Census_allTue Dec 12 03_33_11 2017.csv", sep=",")  ##File downloaded from the Cancer Census database
cancer_census$Gene.Symbol <- toupper(cancer_census$Gene.Symbol)

#Using the cancer genes defined for each tumour types

cancer_census$Tumour.Types.Somatic. <- as.character(cancer_census$Tumour.Types.Somatic.)

cancer_census_per_tumour <- list()
cancer_census_per_tumour[["LUAD"]] <- cancer_census[grep("lung", cancer_census$Tumour.Types.Somatic.), "Gene.Symbol"]
cancer_census_per_tumour[["LUSC"]] <- cancer_census[grep("lung", cancer_census$Tumour.Types.Somatic.), "Gene.Symbol"]
cancer_census_per_tumour[["BRCA"]] <- cancer_census[grep("breast", cancer_census$Tumour.Types.Somatic.), "Gene.Symbol"]
cancer_census_per_tumour[["PRAD"]] <- cancer_census[grep("prostate", cancer_census$Tumour.Types.Somatic.), "Gene.Symbol"]
cancer_census_per_tumour[["LIHC"]] <- cancer_census[grep("hepatocellular", cancer_census$Tumour.Types.Somatic.), "Gene.Symbol"]
cancer_census_per_tumour[["COAD"]] <- cancer_census[grep("color", cancer_census$Tumour.Types.Somatic.), "Gene.Symbol"]
cancer_census_per_tumour[["STAD"]] <- cancer_census[c(grep("stomach", cancer_census$Tumour.Types.Somatic.),
                                                    grep("gastric", cancer_census$Tumour.Types.Somatic.)), "Gene.Symbol"]


sapply(cancer_census_per_tumour, length)

mutations_df2 <- vector()
for(tumour in tumours){
  local_cancer_census <- cancer_census_per_tumour[[tumour]]
  
  local_ratio <- mutations_df[mutations_df$Tumour == tumour,]
  
  local_ratio$Category <- "Genes_other"
  local_ratio$Category[local_ratio$Hugo_Symbol %in% local_cancer_census] <- "Known_cancer_genes"
  mutations_df2 <- rbind(mutations_df2, local_ratio)  
}



mutations_df2$Tumour <- factor(mutations_df2$Tumour, levels=tumours)

aggregate(Category ~ Tumour+Variant_type, mutations_df2, table)


mutations_df2$Variant_type <- factor(mutations_df2$Variant_type, levels=c("Missense", "LoF"))

pdf("Supp1_Syn_ratio_census.pdf",
    height=3, width=11)
g <- ggplot(mutations_df2, aes(x=Syn_ratio))+
  geom_density(aes(fill=Category), alpha=0.7)+
  facet_grid(Variant_type~Tumour)+
  ylab("Density")+
  xlab("Ratio of number of missense or LoF mutations\nand synonymous mutations")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())
print(g)

dev.off()







