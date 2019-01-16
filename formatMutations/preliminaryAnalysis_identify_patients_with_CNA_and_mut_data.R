###Obtain background list of patients with data on CNVs and somatic mutations

tumours <- c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "ESCA", "GBM", 
             "HNSC", "KICH", "KIRC", "KIRP", "LGG",  "LIHC", "LUAD", "LUSC", 
             "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", 
             "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM")

patients_with_mut_info <- list()
for(tumour in tumours){
  path <- ###Path to maf files  
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

sum(sapply(patients_with_mut_info, length))


patients_with_CNV_info <- list()
for(tumour in tumours){
  TCGA_files <- read.table(paste(tumour, "/gdac.broadinstitute.org_", tumour, ".Merge_snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.Level_3.2016012800.0.0/", tumour, ".snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.seg.txt", sep=""),
                           header=TRUE)  #File from TCGA

  #Select tumour samples
  samples <- TCGA_files$Sample
  sample_codes <- substring(samples, 14, 16)
  
  #Sample types available: "10A" "01A" "11A" "11B" "01B" "10B" "06A"
  
  #normal_samples <- which(sample_codes %in% c("11A", "11B"))
  tumour_samples <- which(sample_codes %in% c("01A", "01B", "01C", "01D"))
  
  TCGA_files <- TCGA_files[c(tumour_samples),]
  patients <- unique(substr(TCGA_files$Sample, 9, 12))
  patients_with_CNV_info[[tumour]] <- patients
  print(tumour)
}
save(patients_with_CNV_info, file="patients_with_CNV_info.Rdata")

load("patients_with_CNV_info.Rdata")
sum(sapply(patients_with_CNV_info, length))
