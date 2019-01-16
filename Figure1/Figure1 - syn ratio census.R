##Syn ratio

library(ggplot2)
library(readr)

source("helper_functions.R")

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

tumours <- c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "ESCA", "GBM", 
             "HNSC", "KICH", "KIRC", "KIRP", "LGG",  "LIHC", "LUAD", "LUSC", 
             "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", 
             "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM")  

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

cancer_census_per_tumour[["ACC"]] <- cancer_census[grep("adrenocortical", cancer_census$Tumour.Types.Somatic.), "Gene.Symbol"]
cancer_census_per_tumour[["BLCA"]] <- cancer_census[grep("bladder", cancer_census$Tumour.Types.Somatic.), "Gene.Symbol"]
cancer_census_per_tumour[["CESC"]] <- cancer_census[grep("cervical", cancer_census$Tumour.Types.Somatic.), "Gene.Symbol"]
cancer_census_per_tumour[["CHOL"]] <- cancer_census[grep("cholangiocarcinoma", cancer_census$Tumour.Types.Somatic.), "Gene.Symbol"]
cancer_census_per_tumour[["GBM"]] <- cancer_census[grep("glioblastoma", cancer_census$Tumour.Types.Somatic.), "Gene.Symbol"]
cancer_census_per_tumour[["HNSC"]] <- cancer_census[grep("head", cancer_census$Tumour.Types.Somatic.), "Gene.Symbol"]
cancer_census_per_tumour[["KIRC"]] <- cancer_census[grep("clear cell renal", cancer_census$Tumour.Types.Somatic.), "Gene.Symbol"]
cancer_census_per_tumour[["KIRP"]] <- cancer_census[grep("papillary renal", cancer_census$Tumour.Types.Somatic.), "Gene.Symbol"]
cancer_census_per_tumour[["SARC"]] <- cancer_census[grep("sarcoma", cancer_census$Tumour.Types.Somatic.), "Gene.Symbol"]
cancer_census_per_tumour[["TGCT"]] <- cancer_census[grep("testicular", cancer_census$Tumour.Types.Somatic.), "Gene.Symbol"]
cancer_census_per_tumour[["THCA"]] <- cancer_census[grep("thyroid", cancer_census$Tumour.Types.Somatic.), "Gene.Symbol"]
cancer_census_per_tumour[["UCEC"]] <- cancer_census[grep("endometrial", cancer_census$Tumour.Types.Somatic.), "Gene.Symbol"]
cancer_census_per_tumour[["UVM"]] <- cancer_census[grep("uveal melanoma", cancer_census$Tumour.Types.Somatic.), "Gene.Symbol"]

cancer_census_per_tumour[["ESCA"]] <- cancer_census[grep(paste(c("oesophageal","oesophagus"),collapse="|"), cancer_census$Tumour.Types.Somatic.), "Gene.Symbol"]

cancer_census_per_tumour[["READ"]] <- cancer_census[grep("colorectal", cancer_census$Tumour.Types.Somatic.), "Gene.Symbol"]

cancer_census_per_tumour[["PCPG"]] <- cancer_census[grep(paste(c("pheochromocytoma", "paraganglioma"),collapse="|"), cancer_census$Tumour.Types.Somatic.), "Gene.Symbol"]

#For melanoma
exclude <- grep(paste(c("uveal melanoma","malignant melanoma of soft parts"),collapse="|"), cancer_census$Tumour.Types.Somatic.)
include <- grep("melanoma", cancer_census$Tumour.Types.Somatic.)
include <- include[!(include %in% exclude)]
cancer_census_per_tumour[["SKCM"]] <- cancer_census[include, "Gene.Symbol"]

#For ovarian
exclude <- grep("clear cell ovarian carcinoma", cancer_census$Tumour.Types.Somatic.)
include <- grep("ovarian", cancer_census$Tumour.Types.Somatic.)
include <- include[!(include %in% exclude)]
cancer_census_per_tumour[["OV"]] <- cancer_census[include, "Gene.Symbol"]

##Pancreatic
exclude <- grep(paste(c("pancreatic neuroendocrine tumours", "pancreatic acinar cell carcinoma",	
                        "pancreatic intraductal",	"pancreas acinar carcinoma"),collapse="|"), 
                cancer_census$Tumour.Types.Somatic.)
include <- grep(paste(c("pancreatic",	"pancreas"),collapse="|"), cancer_census$Tumour.Types.Somatic.)
include <- include[!(include %in% exclude)]
cancer_census_per_tumour[["PAAD"]] <- cancer_census[include, "Gene.Symbol"]

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


##ACC, CESC, CHOL - only 1 known cancer gene
##ESCA, KIRP, OV, PRAD - multiple known cancer genes, but only 1 with a Syn ratio diff from inf
##LGG, TGCT, THYM, UCS, KICH - no known cancer genes in the dataset
##PAAD, UVM, PCPG no known cancer genes with a sYn ratio diff from inf

mutations_df2 <- mutations_df2[!(mutations_df2$Tumour %in% c("ACC", "CESC", "CHOL", 
"ESCA", "KIRP", "OV", "PRAD", "LGG", "TGCT", "THYM", "UCS", "KICH",
"PAAD", "UVM", "PCPG")),]



pdf("Supp1_Syn_ratio_census.pdf",
    height=15, width=5)
g <- ggplot(mutations_df2, aes(x=Syn_ratio))+
  geom_density(aes(fill=Category), alpha=0.7)+
  facet_grid(Tumour~Variant_type, scales="free_y")+
  ylab("Density")+
  xlab("Ratio of number of missense or LoF mutations\nand synonymous mutations")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())
print(g)
dev.off()



pdf("Supp1_Syn_ratio_census_part1.pdf",
    height=9, width=5)
g <- ggplot(subset(mutations_df2, Tumour %in% c("BLCA", "BRCA",
                                                "COAD", "GBM", "HNSC", "KIRC", "LIHC", "LUAD")), 
                   aes(x=Syn_ratio))+
  geom_density(aes(fill=Category), alpha=0.7)+
  facet_grid(Tumour~Variant_type, scales="free_y")+
  ylab("Density")+
  xlab("Ratio of number of missense or LoF mutations\nand synonymous mutations")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())
print(g)
dev.off()


pdf("Supp1_Syn_ratio_census_part2.pdf",
    height=8.5, width=5)
g <- ggplot(subset(mutations_df2, Tumour %in% c("LUSC", "READ", "SARC", "SKCM", "STAD", "THCA", "UCEC")),
                   aes(x=Syn_ratio))+
  geom_density(aes(fill=Category), alpha=0.7)+
  facet_grid(Tumour~Variant_type, scales="free_y")+
  ylab("Density")+
  xlab("Ratio of number of missense or LoF mutations\nand synonymous mutations")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())
print(g)
dev.off()
