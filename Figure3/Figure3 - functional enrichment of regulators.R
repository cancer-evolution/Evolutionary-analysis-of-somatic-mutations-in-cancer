load("regulator_classification.Rdata")

temp <- regulator_classification[regulator_classification$Regulator_class == "Mixed_downstream",]

temp_point <- temp[temp$Alt_binary == "Point",]

temp_CNA <- temp[temp$Alt_binary == "CNV",]

library(gProfileR)

results <- gprofiler(temp$Regulator, organism = "hsapiens",
                     exclude_iea=T, max_set_size=1000)
results_not_go <- results[!(results$domain %in% c("BP", "MF", "CC")),]
results_not_go <- results_not_go[results_not_go$domain != "hpa",]
results_not_go <- results_not_go[results_not_go$domain != "mir",]
results_not_go <- results_not_go[results_not_go$domain != "tf",]
results_not_go <- results_not_go[results_not_go$domain != "hp",]

write.table(results_not_go, file="Supp_table_Functional_enrichment_UC-EM-i_regulators.txt",
            quote=FALSE, sep="\t")


##Calculate ssGSEA of all samples
load("patients_with_mut_info.Rdata")

tumours <- c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "ESCA", "GBM", 
             "HNSC", "KICH", "KIRC", "KIRP", "LGG",  "LIHC", "LUAD", "LUSC", 
             "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", 
             "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM")  

expression_tumour <- list()
for(tumour in tumours){
  for(sample_type in c("tumour")){
    path <- paste(tumour, "/gdac.broadinstitute.org_", tumour, 
                  ".Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2016012800.0.0/", 
                  tumour, 
                  ".rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt",
                  sep="")
    
    exp <- read.delim(path) ##Load normalized gene expression data from TCGA
    exp <- exp[-1,]
    
    exp_tumour <- exp[,substr(colnames(exp), 14, 16) %in% c("01A", "01B", "01C", "01D")]
    
    colnames(exp_tumour) <- substr(colnames(exp_tumour),9,12)
    
    all_genes <- as.character(exp[,1])
    
    all_genes <- unlist(strsplit(all_genes, "\\|"))[seq(1, length(all_genes)*2, 2)]
    
    genes_remove <- c("?", "SLC35E2")
    genes_remove_index <- which(all_genes %in% genes_remove)
    all_genes <- all_genes[-genes_remove_index]
    
    exp_tumour <- exp_tumour[-genes_remove_index,]
    
    if(!is.null(nrow(exp_tumour))){
      rownames(exp_tumour) <- all_genes
    }
    
    patients_to_keep <- patients_with_mut_info[[tumour]]
    patients_tumour <- colnames(exp_tumour)
    
    patients_to_keep <- patients_to_keep[patients_to_keep %in% patients_tumour]
    
    expression_tumour[[tumour]] <- exp_tumour[,colnames(exp_tumour) %in% patients_to_keep]
    print(dim(expression_tumour[[tumour]]))
  }
}


library(GSVA)
library(qusage)
##Get pathways from MSigDB
kegg_id <- read.gmt("c2.cp.kegg.v6.2.symbols.gmt")
kegg_id <- kegg_id[-grep("DISEASE", names(kegg_id))]
kegg_id <- kegg_id[-grep("DEPRESSION", names(kegg_id))]
kegg_id <- kegg_id[-grep("CANCER", names(kegg_id))]
kegg_id <- kegg_id[-grep("MELANOMA", names(kegg_id))]
kegg_id <- kegg_id[-grep("INFECTION", names(kegg_id))]
kegg_id <- kegg_id[-grep("CARCINOMA", names(kegg_id))]
kegg_id <- kegg_id[-grep("LEUKEMIA", names(kegg_id))]
kegg_id <- kegg_id[-grep("CARDIOMYOPATHY", names(kegg_id))]
kegg_id <- kegg_id[-grep("MYOCARDITIS", names(kegg_id))]
kegg_id <- kegg_id[-grep("ASTHMA", names(kegg_id))]
kegg_id <- kegg_id[-grep("MUSCLE", names(kegg_id))]
kegg_id <- kegg_id[-grep("DIABETES", names(kegg_id))]
kegg_id <- kegg_id[-grep("AMYOTROPHIC_LATERAL_SCLEROSIS", names(kegg_id))]
kegg_id <- kegg_id[-grep("GLIOMA", names(kegg_id))]
kegg_id <- kegg_id[-grep("SYSTEMIC_LUPUS_ERYTHEMATOSUS", names(kegg_id))]
kegg_id <- kegg_id[-grep("KEGG_ALLOGRAFT_REJECTION", names(kegg_id))]
kegg_id <- kegg_id[-grep("KEGG_PRIMARY_IMMUNODEFICIENCY", names(kegg_id))]

library(Biobase)
gene_pathways <- reverseSplit(kegg_id)

ssGSEA_local <- list()
for(tumour in tumours){
  ssGSEA_local[[tumour]] <- gsva(as.matrix(expression_tumour[[tumour]]), kegg_id, method="ssgsea", verbose=FALSE, rnaseq=FALSE, ssgsea.norm=TRUE)
}
#save(ssGSEA_local, file="~/Documents/Paper 3/Objects_second_submission/ssGSEA_local.Rdata")


##Compare the ssGSEA by mutations
library(readr)
source("helper_functions.R")
load("patients_with_mut_info.Rdata")
load("variants_LoF.Rdata")
load("variants_missense.Rdata")

genes_phy <- read.csv("gene_phylostrata.txt")  ###File with gene name, entrez and phylostratum as column

mutations_df <- load_mutations()

mutations_df <- mutations_df[mutations_df$Number_mutations >= 3,]
mutations_df <- mutations_df[mutations_df$Syn_ratio > 1,]

load("regulator_classification.Rdata")

mixed_regulators <- regulator_classification[regulator_classification$Regulator_class == "Mixed_downstream",]
mutations_df <- mutations_df[mutations_df$Hugo_Symbol %in% mixed_regulators$Regulator,]

pathways_DE <- vector()
for(tumour in tumours){
  local_recurrent_mut <- mutations_df[mutations_df$Tumour == tumour,]
  local_point_patients <- rbind(variants_LoF[[tumour]][,c("Hugo_Symbol", "Tumor_Sample_Barcode")],
                                variants_missense[[tumour]][,c("Hugo_Symbol", "Tumor_Sample_Barcode")])
  
  for(gene in unique(local_point_patients$Hugo_Symbol)){
    local_pathways <- gene_pathways[[gene]]
    local_patients <- local_point_patients[local_point_patients$Hugo_Symbol == gene,]
    local_patients <- substr(local_patients$Tumor_Sample_Barcode, 9, 12)
    
    if(!is.null(local_pathways) & sum(colnames(ssGSEA_local[[tumour]]) %in% local_patients) >= 3){
      for(path in local_pathways){
        p <- wilcox.test(ssGSEA_local[[tumour]][path,colnames(ssGSEA_local[[tumour]]) %in% local_patients],
        ssGSEA_local[[tumour]][local_pathways,!(colnames(ssGSEA_local[[tumour]]) %in% local_patients)])$p.value
        pathways_DE <- rbind(pathways_DE, c(tumour, gene, path, p))
      }
    }
    print(gene)
  }
}

pathways_DE <- as.data.frame(pathways_DE)
colnames(pathways_DE) <- c("Tumour", "Gene", "Pathway", "P_value")
pathways_DE$P_value <- as.numeric(as.character(pathways_DE$P_value))
pathways_DE$Adj_p <- p.adjust(pathways_DE$P_value, method="BH")

write.table(pathways_DE, file="Supp_ssGSEA_of_pathways_of_point_mutated_UC-EM-i_regulators.txt",
            quote=FALSE, sep="\t", row.names=FALSE)


View(pathways_DE[pathways_DE$Adj_p < 0.05,])

View(pathways_DE[pathways_DE$Adj_p < 0.05 & pathways_DE$Gene == "PIK3R1",])
length(unique(pathways_DE[pathways_DE$Adj_p < 0.05 & pathways_DE$Gene == "TP53","Tumour"]))

