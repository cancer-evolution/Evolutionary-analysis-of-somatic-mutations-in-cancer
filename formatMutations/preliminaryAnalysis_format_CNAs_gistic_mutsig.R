library(GenomicRanges)
library(Homo.sapiens)

tumours <- c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "ESCA", "GBM", 
             "HNSC", "KICH", "KIRC", "KIRP", "LGG",  "LIHC", "LUAD", "LUSC", 
             "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", 
             "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM")

#Selecting only CNVs in tumour samples, and only genes that are either deleted or amplified,
#With a segment mean in the top 10%.
CNVs_curated <- list()
CNVs_curated_focal <- list()
for(tumour in tumours){
  for(tumour in tumours){
    file <- paste("gdac.broadinstitute.org_",
                  tumour, "-TP.CopyNumber_Gistic2.Level_4.2016012800.0.0", sep="")
    if(tumour == "SKCM"){
      file <- paste("gdac.broadinstitute.org_",
                    tumour, "-TM.CopyNumber_Gistic2.Level_4.2016012800.0.0", sep="")
    }
    if(file.exists(file)){
      if(tumour != "SKCM"){
        local_gistic_amp <- read.delim(paste("gdac.broadinstitute.org_",
                                             tumour, "-TP.CopyNumber_Gistic2.Level_4.2016012800.0.0/table_amp.conf_99.txt", sep=""))
        local_gistic_del <- read.delim(paste("gdac.broadinstitute.org_",
                                             tumour, "-TP.CopyNumber_Gistic2.Level_4.2016012800.0.0/table_del.conf_99.txt", sep=""))
      }else{
        local_gistic_amp <- read.delim(paste("gdac.broadinstitute.org_",
                                             tumour, "-TM.CopyNumber_Gistic2.Level_4.2016012800.0.0/table_amp.conf_99.txt", sep=""))
        local_gistic_del <- read.delim(paste("gdac.broadinstitute.org_",
                                             tumour, "-TM.CopyNumber_Gistic2.Level_4.2016012800.0.0/table_del.conf_99.txt", sep=""))
      }
      
      if(nrow(local_gistic_amp) > 0){
        local_gistic_amp_genes_region <- unlist(strsplit(as.character(local_gistic_amp$genes_in_region), ","))
        local_gistic_amp_genes_region <- unlist(strsplit(local_gistic_amp_genes_region, "\\|"))
        local_gistic_amp_genes_peak <- unlist(strsplit(as.character(local_gistic_amp$genes_in_peak), ","))
        local_gistic_amp_genes_peak <- unlist(strsplit(local_gistic_amp_genes_peak, "\\|"))
        
        local_gistic_amp_genes_region <- local_gistic_amp_genes_region[-grep("^ENSG", local_gistic_amp_genes_region)]
        local_gistic_amp_genes_peak <- local_gistic_amp_genes_peak[-grep("^ENSG", local_gistic_amp_genes_peak)]
        
        CNVs_curated[[tumour]][["amplifications"]] <- local_gistic_amp_genes_region
        CNVs_curated_focal[[tumour]][["amplifications"]] <- local_gistic_amp_genes_peak
      }
      
      if(nrow(local_gistic_del) > 0){
        local_gistic_del_genes_region <- unlist(strsplit(as.character(local_gistic_del$genes_in_region), ","))
        local_gistic_del_genes_region <- unlist(strsplit(local_gistic_del_genes_region, "\\|"))
        
        local_gistic_del_genes_peak <- unlist(strsplit(as.character(local_gistic_del$genes_in_peak), ","))
        local_gistic_del_genes_peak <- unlist(strsplit(local_gistic_del_genes_peak, "\\|"))
        
        local_gistic_del_genes_region <- local_gistic_del_genes_region[-grep("^ENSG", local_gistic_del_genes_region)]
        local_gistic_del_genes_peak <- local_gistic_del_genes_peak[-grep("^ENSG", local_gistic_del_genes_peak)]
        
        CNVs_curated[[tumour]][["deletions"]] <- local_gistic_del_genes_region
        CNVs_curated_focal[[tumour]][["deletions"]] <- local_gistic_del_genes_peak
      }
    }else{
      print(tumour)
    }
  }
}

save(CNVs_curated, file="CNVs_curated_gistic.Rdata")
save(CNVs_curated, file="CNVs_curated_gistic_focal.Rdata")

