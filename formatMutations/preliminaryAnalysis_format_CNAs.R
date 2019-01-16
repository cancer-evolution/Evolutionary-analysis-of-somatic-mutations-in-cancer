library(GenomicRanges)
library(Homo.sapiens)

tumours <- c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "ESCA", "GBM", 
             "HNSC", "KICH", "KIRC", "KIRP", "LGG",  "LIHC", "LUAD", "LUSC", 
             "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", 
             "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM")

#Select the maximum segment mean

for(tumour in tumours){
  if(!(tumour %in% c("ACC","BLCA","BRCA","CESC","CHOL","COAD","ESCA","GBM","HNSC","KICH",
                   "PAAR", "PCPG", "READ", "SARC", "TGCT","THYM",
                   "UCS","UVM", "MESO"))){
    CNV <- list()
    pdf(paste("Distribution of percentage of gene covered",
              tumour, sep=" - ", ".pdf"))
    
    TCGA_files <- read.table(paste(tumour, "/gdac.broadinstitute.org_", tumour, ".Merge_snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.Level_3.2016012800.0.0/", tumour, ".snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.seg.txt", sep=""),
                             header=TRUE)  #File from TCGA
    
    #Select tumour samples
    samples <- TCGA_files$Sample
    sample_codes <- substring(samples, 14, 16)
    
    #Sample types available: "10A" "01A" "11A" "11B" "01B" "10B" "06A"
    
    normal_samples <- which(sample_codes %in% c("11A", "11B"))
    tumour_samples <- which(sample_codes %in% c("01A", "01B", "01C", "01D"))
    
    TCGA_files <- TCGA_files[c(tumour_samples, normal_samples),]
    TCGA_files$Chromosome <- paste("chr", TCGA_files$Chromosome, sep="")
    
    ranges_TCGA <- makeGRangesFromDataFrame(TCGA_files, keep.extra.columns = TRUE)
    
    g <- genes(Homo.sapiens, columns="SYMBOL")
    
    olaps <- findOverlaps(g, ranges_TCGA) #query is the first one, subject is the second
    
    g_df <- as.data.frame(g)  #note that there are chromosomes here not found in TCGA, and it is giving a warning
    
    olaps_df <- as.data.frame(olaps)
    
    ranges_TCGA_df <- as.data.frame(ranges_TCGA)
    
    olaps_df <- data.frame(olaps_df,
                           g_df[olaps_df$queryHits,])
    
    colnames(olaps_df)[3:8] <- paste(colnames(olaps_df)[3:8], "query", sep=".")
    
    olaps_df <- data.frame(olaps_df,
                           ranges_TCGA_df[olaps_df$subjectHits,])
    colnames(olaps_df)[9:13] <- paste(colnames(olaps_df)[9:13], "subject", sep=".")
    
    ##The above is using "any", so the segment can include only the start, the end, or be within
    
    ##Amplification should be within probes, including the start
    #Deletions just any part
    
    olaps_df$Gene <- unname(unlist(olaps_df$SYMBOL.query))
    
    #At least 10 probes
    olaps_df <- olaps_df[olaps_df$Num_Probes > 10,]
    
    
    #will use just positive and negative here. I will exclude genes that have both deletions
    #and amplification
    olaps_df_amp <- olaps_df[olaps_df$Segment_Mean > 0,]
    olaps_df_del <- olaps_df[olaps_df$Segment_Mean < 0,]
    
    ##start.subject needs to be smaller than start.query in both cases
    olaps_df_amp <- olaps_df_amp[which(olaps_df_amp$start.subject <= olaps_df_amp$start.query),]
    olaps_df_del <- olaps_df_del[which(olaps_df_del$start.subject <= olaps_df_del$start.query),]
    
    ##Calculate percentage of gene covered
    
    olaps_df_amp$Percentage_gene <- unname(apply(olaps_df_amp, 1, function(row){
      portion_covered <- min(row$end.query, row$end.subject)-
        max(row$start.query, row$start.subject)+1
      return(portion_covered/row$width.query*100)
    }))
    hist(olaps_df_amp$Percentage_gene, main=paste("Amplifications", tumour), breaks=100)
    
    olaps_df_del$Percentage_gene <- unname(apply(olaps_df_del, 1, function(row){
      portion_covered <- min(row$end.query, row$end.subject)-
        max(row$start.query, row$start.subject)+1
      return(portion_covered/row$width.query*100)
    }))
    hist(olaps_df_del$Percentage_gene, main=paste("Deletions", tumour), breaks=100)
    
    
    print(tumour)
    print(sum(olaps_df_del$Percentage_gene != 100))  
    print(sum(olaps_df_amp$Percentage_gene != 100))  
    
    olaps_df_amp <- olaps_df_amp[olaps_df_amp$Percentage_gene == 100,]
    olaps_df_del <- olaps_df_del[olaps_df_del$Percentage_gene == 100,]
    
    
    #Use the maximum Segment_Mean for amplification, and the minimum for deletions
    
    gene_deletions_min <- aggregate(Segment_Mean ~ Gene+Sample, olaps_df_del, min)
    gene_amplifications_max <- aggregate(Segment_Mean ~ Gene+Sample, olaps_df_amp, max)
    
    CNV[[tumour]] <- list(amplifications = gene_amplifications_max,
                          deletions = gene_deletions_min)
    print(tumour)
    dev.off()
    save(CNV, file=paste("CNV2_", tumour, ".Rdata", sep=""))  #now excluding amplification that do not cover entire gene
    
  }
  
}

##Merging into a single object
for(tumour in tumours){
  CNV <- list()
  CNV[[tumour]] <- CNV_o[[tumour]]
}

save(CNV, file=paste("CNV2_", tumour, ".Rdata", sep=""))  #now excluding amplification that do not cover entire gene

###Comparing CNVs of normal and tumour samples
#And checking that genes are either only amplified or only deleted, not both

library(ggplot2)

tumours <- c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "ESCA", "GBM", 
             "HNSC", "KICH", "KIRC", "KIRP", "LGG",  "LIHC", "LUAD", "LUSC", 
             "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", 
             "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM")

#Selecting only CNVs in tumour samples, and only genes that are either deleted or amplified,
#With a segment mean in the top 10%.
CNVs_curated <- list()

for(tumour in tumours){
  load(paste(tumour, ".Rdata", sep=""))
  
  for(type in c("amplifications", "deletions")){
    local_CNV <- CNV[[tumour]][[type]]
    
    samples <- local_CNV$Sample
    sample_codes <- substring(samples, 14, 16)
    
    #Sample types available: "10A" "01A" "11A" "11B" "01B" "10B" "06A"
    
    tumour_samples <- which(sample_codes %in% c("01A", "01B", "01C", "01D"))
    
    local_CNV_tumour <- local_CNV[tumour_samples,]
    
    cutoff <- quantile(abs(local_CNV_tumour$Segment_Mean), 0.90)
    
    if(type == "amplifications"){
      local_CNV_tumour_cut <- local_CNV_tumour[local_CNV_tumour$Segment_Mean >= cutoff,]
    }else{
      local_CNV_tumour_cut <- local_CNV_tumour[local_CNV_tumour$Segment_Mean <= -cutoff,]
    }
    print(tumour)
    CNVs_curated[[tumour]][[type]] <- local_CNV_tumour_cut
  }
}


for(tumour in tumours){
  local_CNV_amp <- CNVs_curated[[tumour]][["amplifications"]]
  local_CNV_del <- CNVs_curated[[tumour]][["deletions"]]
  
  genes_samples_amp <- paste(local_CNV_amp$Gene, local_CNV_amp$Sample, sep=".")
  genes_samples_del <- paste(local_CNV_del$Gene, local_CNV_del$Sample, sep=".")
  
  amp_to_exclude <- which(genes_samples_amp %in% genes_samples_del)
  del_to_exclude <- which(genes_samples_del %in% genes_samples_amp)
  
  print(tumour)
  print(length(amp_to_exclude))
  print(length(del_to_exclude))
  
  if(length(amp_to_exclude) != 0){
    local_CNV_amp <- local_CNV_amp[-amp_to_exclude,]
  }
  if(length(del_to_exclude) != 0){
    local_CNV_del <- local_CNV_del[-del_to_exclude,]  
  }
  CNVs_curated[[tumour]][["amplifications"]] <- local_CNV_amp
  CNVs_curated[[tumour]][["deletions"]] <- local_CNV_del
}


save(CNVs_curated, file="CNVs_curated2.Rdata")

