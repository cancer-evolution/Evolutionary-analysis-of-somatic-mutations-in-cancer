###Filtering of variants
#Script to reduce variants to the relevant ones per patient

library(stringr)
library(readr)

#Calculating overlap between somatic mutations obtained by different methods in TCGA

tumours <- c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "ESCA", "GBM", 
             "HNSC", "KICH", "KIRC", "KIRP", "LGG",  "LIHC", "LUAD", "LUSC", 
             "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", 
             "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM")

maf_files <- list()
variants <- list()

pdf("Concordance_between_variant_callers.pdf")
for(tumour in tumours){
  path <-  ###path to maf files
  filenames <- list.files(path, pattern="*.maf$", full.names=TRUE)  ##MAF files from TCGA
  maf_files[[tumour]] <- lapply(filenames, read.delim, comment.char = "#")
  names(maf_files[[tumour]]) <- c("muse", "mutect", "somaticsniper", "varscan")
  patients <- unique(c(as.character(maf_files[[tumour]]$muse$Tumor_Sample_Barcode),
                       as.character(maf_files[[tumour]]$mutec$Tumor_Sample_Barcode),
                       as.character(maf_files[[tumour]]$somaticsniper$Tumor_Sample_Barcode),
                       as.character(maf_files[[tumour]]$varscan$Tumor_Sample_Barcode)))
  number_of_patients <- length(patients)
  maf_files[[tumour]]$muse$Variant_summary <- apply(maf_files[[tumour]]$muse[,c(1, 6, 7, 11:13, 16:17)], 1, paste, collapse="_")
  maf_files[[tumour]]$mutect$Variant_summary <- apply(maf_files[[tumour]]$mutect[,c(1, 6, 7, 11:13, 16:17)], 1, paste, collapse="_")
  maf_files[[tumour]]$somaticsniper$Variant_summary <- apply(maf_files[[tumour]]$somaticsniper[,c(1, 6, 7, 11:13, 16:17)], 1, paste, collapse="_")
  maf_files[[tumour]]$varscan$Variant_summary <- apply(maf_files[[tumour]]$varscan[,c(1, 6, 7, 11:13, 16:17)], 1, paste, collapse="_")
  all_variants <- c(maf_files[[tumour]]$muse$Variant_summary,
                    maf_files[[tumour]]$mutect$Variant_summary,
                    maf_files[[tumour]]$somaticsniper$Variant_summary,
                    maf_files[[tumour]]$varscan$Variant_summary)
  variant_table <- table(all_variants)
  hist(variant_table, main=tumour, xlab="Number of methods in which variants were identified",
       ylab="Number of variants")
  variants_to_choose <- names(variant_table)[variant_table==4]
  
  columns_to_remove <- c("t_depth", "t_ref_count", "t_alt_count", "n_depth", "src_vcf_id", "PICK", "MINIMISED", "HGNC_ID")
  variants_to_choose_muse <- maf_files[[tumour]]$muse[maf_files[[tumour]]$muse$Variant_summary %in% variants_to_choose, 
                                                      which(!(colnames(maf_files[[tumour]]$muse) %in% columns_to_remove))]
  
  variants_to_choose_mutect <- maf_files[[tumour]]$mutect[maf_files[[tumour]]$mutect$Variant_summary %in% variants_to_choose,
                                                          which(!(colnames(maf_files[[tumour]]$mutect) %in% columns_to_remove))]
  
  variants_to_choose_somaticsniper <- maf_files[[tumour]]$somaticsniper[maf_files[[tumour]]$somaticsniper$Variant_summary %in% variants_to_choose,
                                                                        which(!(colnames(maf_files[[tumour]]$somaticsniper) %in% columns_to_remove))]
  
  variants_to_choose_varscan <- maf_files[[tumour]]$varscan[maf_files[[tumour]]$varscan$Variant_summary %in% variants_to_choose,
                                                            which(!(colnames(maf_files[[tumour]]$varscan) %in% columns_to_remove))]
  variants_to_choose_muse[is.na(variants_to_choose_muse)] <- ""
  variants_to_choose_mutect[is.na(variants_to_choose_mutect)] <- ""
  variants_to_choose_somaticsniper[is.na(variants_to_choose_somaticsniper)] <- ""
  variants_to_choose_varscan[is.na(variants_to_choose_varscan)] <- ""
  
  variants_to_choose_muse[,"all_effects"] <- as.character(variants_to_choose_muse[,"all_effects"])
  variants_to_choose_mutect[,"all_effects"] <- as.character(variants_to_choose_mutect[,"all_effects"])
  variants_to_choose_somaticsniper[,"all_effects"] <- as.character(variants_to_choose_somaticsniper[,"all_effects"])
  variants_to_choose_varscan[,"all_effects"] <- as.character(variants_to_choose_varscan[,"all_effects"])
  
  variants_to_choose_muse[,"all_effects"] <- sapply(variants_to_choose_muse[,"all_effects"], function(x){ paste(sort(unlist(strsplit(x, ";"))), collapse=";") })
  variants_to_choose_mutect[,"all_effects"] <- sapply(variants_to_choose_mutect[,"all_effects"], function(x){ paste(sort(unlist(strsplit(x, ";"))), collapse=";") })
  variants_to_choose_somaticsniper[,"all_effects"] <- sapply(variants_to_choose_somaticsniper[,"all_effects"], function(x){ paste(sort(unlist(strsplit(x, ";"))), collapse=";") })
  variants_to_choose_varscan[,"all_effects"] <- sapply(variants_to_choose_varscan[,"all_effects"], function(x){ paste(sort(unlist(strsplit(x, ";"))), collapse=";") })
  
  variants[[tumour]] <- unique(rbind(variants_to_choose_muse,
                                     variants_to_choose_mutect,
                                     variants_to_choose_somaticsniper,
                                     variants_to_choose_varscan))
  if(nrow(variants[[tumour]]) != length(variants_to_choose)){
    print(tumour)
    print("CHECK!")
  }
  print(tumour)
  number_of_patients <- c(number_of_patients, length(unique(variants[[tumour]]$Tumor_Sample_Barcode)))
  print(number_of_patients)    
}
dev.off()

save(variants, file="variants_in_common.Rdata")

source("helper_functions.R")

#Filtering
variants_no_exac <- list()
variants_no_common <- list()
variants_no_common_canonical <- list()
variants_filtered <- list()
for(tumour in tumours){
  variants_no_exac[[tumour]] <- exclude_ExAC(variants[[tumour]])
  variants_no_common[[tumour]] <- exclude_1000_genomes(variants_no_exac[[tumour]])
  #Use data only from canonical transcripts
  variants_no_common_canonical[[tumour]] <- variants_no_common[[tumour]][which(variants_no_common[[tumour]]$CANONICAL == "YES"),]
  
  variants_filtered[[tumour]] <- variants_no_common_canonical[[tumour]]
}

save(variants_filtered, file="variants_filtered.Rdata")

#load("variants_filtered.Rdata")


##Identify all possible types of variants

variant_types <- vector()
for(tumour in tumours){
  variant_types <- c(variant_types, unique(as.character(variants_filtered[[tumour]]$One_Consequence)))
}
unique(variant_types)



#Do filtering based on loss-of-function
loss_variants <- c("frameshift_variant", "stop_gained", 
                   "splice_acceptor_variant", "splice_donor_variant")

missense_variants <- c("missense_variant", "inframe_deletion", "inframe_insertion")



#View(variants_no_common_canonical[[tumour]][which(variants_no_common_canonical[[tumour]]$Exon_Number == ""),])
#Use variants that are not in the last exon
#Note that those with empty exon numbers are the ones that occur in splice sites
# exons <- as.character(variants_no_common_canonical[[tumour]]$Exon_Number)
# split_exons <- as.data.frame(str_split_fixed(exons, "/", n=2))
# split_exons <- apply(split_exons, 2, as.numeric)
# to_keep <- which(apply(split_exons, 1, function(x){ x[1] != x[2] }))  #variants not in last exon
# to_keep <- sort(c(to_keep, which(exons == "")))




variants_LoF <- list()
variants_missense <- list()
for(tumour in tumours){
  #Only keep variants with loss-of-function properties
  variants_LoF[[tumour]] <- variants_filtered[[tumour]]
  matches_LoF <- unique(grep(paste(loss_variants,collapse="|"), 
                         variants_LoF[[tumour]]$One_Consequence, value=FALSE))
  
  variants_LoF[[tumour]] <- variants_LoF[[tumour]][matches_LoF,]
  ###LoF variants do not have a SIFT, Polyphen score

  #Only keep variants with missense/possible oncogene properties
  variants_missense[[tumour]] <- variants_filtered[[tumour]]
  matches_missense <- unique(grep(paste(missense_variants,collapse="|"), 
                             variants_missense[[tumour]]$Consequence, value=FALSE))
  
  variants_missense[[tumour]] <- variants_missense[[tumour]][matches_missense,]
  #Missense variants do have a SIFT, Polyphen score
  #Keep those marked as "deleterious" from SIFT AND Polyphen probably_damaging (more than possibly_damaging)
  SIFT_deleterious <- grep("deleterious\\(", variants_missense[[tumour]]$SIFT)  #Need to add '(' because if not will pick up deleterious_low_confidence 
  PolyPhen_damaging <- grep("probably_damaging", variants_missense[[tumour]]$PolyPhen)  
  
  missense_to_keep <- intersect(SIFT_deleterious, PolyPhen_damaging)
  variants_missense[[tumour]] <- variants_missense[[tumour]][missense_to_keep,]
}

##Include only data from primary tumour samples

for(tumour in tumours){
  print(tumour)
  samples_miss <- variants_missense[[tumour]]$Tumor_Sample_Barcode
  temp_miss <- substr(samples_miss, 14, 16)
  variants_missense[[tumour]] <- variants_missense[[tumour]][which(temp_miss %in% c("01A", "01B", "01C", "01D")),]
  
  samples_lof <- variants_LoF[[tumour]]$Tumor_Sample_Barcode
  temp_lof <- substr(samples_lof, 14, 16)
  variants_LoF[[tumour]] <- variants_LoF[[tumour]][which(temp_lof %in% c("01A", "01B", "01C", "01D")),]
  
}


save(variants_LoF, file="variants_LoF.Rdata")
save(variants_missense, file="variants_missense.Rdata")





###Use only mutations in genes enriched in LoF or synonymous mutations based on the LoF/Syn ratio
genes_to_exclude <- read_csv("Genes to exclude.txt", 
                             col_names = FALSE)
genes_to_exclude <- unname(unlist(genes_to_exclude[,1]))

GC_length <- read.table("GC_length_manual.txt", 
                        header=TRUE)

gene_mut_properties <- list()
for(tumour in tumours){
  
  #Analyse missense and LoF separately
  
  local_LoF_variants <- variants_LoF[[tumour]][,c("Hugo_Symbol", "Entrez_Gene_Id", "Tumor_Sample_Barcode", "One_Consequence", "Start_Position", "End_Position", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2")]
  
  local_LoF_variants_ag <- aggregate( Tumor_Sample_Barcode ~ Hugo_Symbol,local_LoF_variants, length)
  print(sum(local_LoF_variants_ag$Hugo_Symbol %in% genes_to_exclude))
  local_LoF_variants_ag <- local_LoF_variants_ag[!(local_LoF_variants_ag$Hugo_Symbol %in% genes_to_exclude),]
  colnames(local_LoF_variants_ag) <- c("Hugo_Symbol", "Number_mutations")
  
  local_miss_variants <- variants_missense[[tumour]][,c("Hugo_Symbol", "Entrez_Gene_Id", "Tumor_Sample_Barcode", "One_Consequence", "Start_Position", "End_Position", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2")]
  local_miss_variants_ag <- aggregate( Tumor_Sample_Barcode ~ Hugo_Symbol,local_miss_variants, length)
  print(sum(local_miss_variants_ag$Hugo_Symbol %in% genes_to_exclude))
  local_miss_variants_ag <- local_miss_variants_ag[!(local_miss_variants_ag$Hugo_Symbol %in% genes_to_exclude),]
  colnames(local_miss_variants_ag) <- c("Hugo_Symbol", "Number_mutations")
    
  local_variants_ag <- rbind(cbind(local_LoF_variants_ag, Variant_type = "LoF"),
                             cbind(local_miss_variants_ag, Variant_type = "Missense"))
  
  local_variants_ag$Gene_length <- GC_length[match(local_variants_ag$Hugo_Symbol, GC_length$Gene_name), "Length_coding"]
  local_variants_ag$GC_percentage <- GC_length[match(local_variants_ag$Hugo_Symbol, GC_length$Gene_name), "GC_percentage"]
  local_variants_ag$LoF_MB <- local_variants_ag$Number_mutations/(local_variants_ag$Gene_length/1000000)
  
  temp <- variants_filtered[[tumour]][,c("Hugo_Symbol", "Tumor_Sample_Barcode", "Consequence")]
  temp <- temp[temp$Consequence == "synonymous_variant",]
  temp_ag <- aggregate(Tumor_Sample_Barcode ~ Hugo_Symbol,temp, length)
  
  local_variants_ag$Number_synonymous <- temp_ag[match(local_variants_ag$Hugo_Symbol, temp_ag$Hugo_Symbol), "Tumor_Sample_Barcode"]
  local_variants_ag$Number_synonymous[is.na(local_variants_ag$Number_synonymous)] <- 0 
  
  local_variants_ag$Syn_MB <- local_variants_ag$Number_synonymous/(local_variants_ag$Gene_length/1000000)
  
  #LoF over synonymous ratio
  local_variants_ag$Syn_ratio <- local_variants_ag$Number_mutations/local_variants_ag$Number_synonymous
  local_variants_ag$Syn_ratio <- log2(local_variants_ag$Syn_ratio)
  
  gene_mut_properties[[tumour]] <- local_variants_ag
}

save(gene_mut_properties, file="gene_mut_properties.Rdata")
