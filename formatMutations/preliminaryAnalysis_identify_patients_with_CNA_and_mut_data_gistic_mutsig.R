###Obtain background list of patients with data on CNVs and somatic mutations

tumours <- c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "ESCA", "GBM", 
             "HNSC", "KICH", "KIRC", "KIRP", "LGG",  "LIHC", "LUAD", "LUSC", 
             "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", 
             "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM")


patients_with_mut_info <- list()
for(tumour in tumours){
  path <- ##Path to MAF files
  filenames <- list.files(path, pattern="*.maf$", full.names=TRUE)  ##MAF files from TCGA
  maf_files <- lapply(filenames, read.delim, comment.char = "#")
  names(maf_files) <- c("muse", "mutect", "somaticsniper", "varscan")
  patients <- unique(c(as.character(maf_files$muse$Tumor_Sample_Barcode),
                       as.character(maf_files$mutec$Tumor_Sample_Barcode),
                       as.character(maf_files$somaticsniper$Tumor_Sample_Barcode),
                       as.character(maf_files$varscan$Tumor_Sample_Barcode)))
  patients <- substr(patients, 9, 12)
  patients_with_mut_info[[tumour]] <- patients
  print(tumour)
}
save(patients_with_mut_info, file="patients_with_mut_info.Rdata")


patients_with_CNV_info <- list()
for(tumour in tumours){
  file <- paste("gdac.broadinstitute.org_",
                tumour, "-TP.CopyNumber_Gistic2.Level_4.2016012800.0.0", sep="")
  if(tumour == "SKCM"){
    file <- paste("gdac.broadinstitute.org_",
                  tumour, "-TM.CopyNumber_Gistic2.Level_4.2016012800.0.0", sep="")
  }
  if(file.exists(file)){
    if(tumour != "SKCM"){
      patients <- read.delim(paste("gdac.broadinstitute.org_",
                                   tumour, "-TP.CopyNumber_Gistic2.Level_4.2016012800.0.0/all_thresholded.by_genes.txt", sep=""))
     }else{
       patients <- read.delim(paste("gdac.broadinstitute.org_",
                                    tumour, "-TM.CopyNumber_Gistic2.Level_4.2016012800.0.0/all_thresholded.by_genes.txt", sep=""))
     }
    
    patients <- colnames(patients)
    patients <- substr(patients, 9, 12)
    patients <- patients[-c(1:3)]

    patients_with_CNV_info[[tumour]] <- patients
    }else{
    print(tumour)
  }
}

save(patients_with_CNV_info, file="patients_with_CNV_info.Rdata")

