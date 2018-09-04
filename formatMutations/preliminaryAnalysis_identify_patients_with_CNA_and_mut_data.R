###Obtain background list of patients with data on CNVs and somatic mutations

tumours <- c("LUAD", "LUSC", "BRCA", "PRAD", "LIHC", "COAD", "STAD")

patients_with_mut_info <- list()
for(tumour in tumours){
  filenames <- list.files(pattern="*.maf", full.names=TRUE)  ##MAF files from TCGA
  maf_files <- lapply(filenames, read.delim, comment.char = "#")
  names(maf_files) <- c("muse", "mutect", "somaticsniper", "varscan")
  patients <- unique(c(as.character(maf_files$muse$Tumor_Sample_Barcode),
                       as.character(maf_files$mutec$Tumor_Sample_Barcode),
                       as.character(maf_files$somaticsniper$Tumor_Sample_Barcode),
                       as.character(maf_files$varscan$Tumor_Sample_Barcode)))
  patients <- substr(patients, 9, 12)
  patients_with_mut_info[[tumour]] <- patients
}
save(patients_with_mut_info, file="patients_with_mut_info.Rdata")


patients_with_CNV_info <- list()
for(tumour in tumours){
  TCGA_files <- read.table(paste("gdac.broadinstitute.org_", tumour, ".Merge_snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.Level_3.2016012800.0.0/", tumour, ".snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.seg.txt", sep=""),
                           header=TRUE)
  
  #Select tumour samples
  samples <- TCGA_files$Sample
  sample_codes <- substring(samples, 14, 16)
  
  #Sample types available: "10A" "01A" "11A" "11B" "01B" "10B" "06A"
  
  normal_samples <- which(sample_codes %in% c("11A", "11B"))
  tumour_samples <- which(sample_codes %in% c("01A", "01B"))
  
  TCGA_files <- TCGA_files[c(tumour_samples, normal_samples),]
  patients <- unique(substr(TCGA_files$Sample, 9, 12))
  patients_with_CNV_info[[tumour]] <- patients
}
save(patients_with_CNV_info, file="patients_with_CNV_info.Rdata")

