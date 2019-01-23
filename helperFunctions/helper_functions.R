get_gene_lengths <- function(){
  #To get entrez IDs
  load("gene_lengths_human.Rdata")
  genes_with_length <- names(exonic.gene.sizes)
  
  #mart <- useMart(biomart = 'ensembl', dataset = 'hsapiens_gene_ensembl' )
  #ensg_hsap <- getBM(mart = mart, attributes=c('ensembl_gene_id', 'external_gene_name'), filter='ensembl_gene_id', values = genes_with_length)
  #save( ensg_hsap, file='ensg_hsap.Rdata')
  
  load(file='ensg_hsap.Rdata')
  
  names(exonic.gene.sizes) <- ensg_hsap[match(names(exonic.gene.sizes),ensg_hsap[,1]), 2]
  
  genes_length <- data.frame("Gene" = names(exonic.gene.sizes),
                             "Length" = unname(exonic.gene.sizes))
  return(genes_length)
}


exclude_ExAC <- function(local_vep_results){
  #Exclude variants with a frequency in ExAC different to 0
  ExAC_values <- unique(c(local_vep_results$ExAC_AF_EAS, local_vep_results$ExAC_AF_AMR,
                          local_vep_results$ExAC_AF_NFE, local_vep_results$ExAC_AF_AFR,
                          local_vep_results$ExAC_AF, local_vep_results$ExAC_AF_FIN,
                          local_vep_results$ExAC_AF_OTH, local_vep_results$ExAC_AF_SAS,
                          local_vep_results$ExAC_AF_Adj))
  
  #Overall, frequencies are small. 
  
  #Exclude variants that are not equal to 0 and not NA
  ExAC_names <- c("ExAC_AF_EAS", "ExAC_AF_AMR", "ExAC_AF_NFE", "ExAC_AF_AFR",
                  "ExAC_AF", "ExAC_AF_FIN", "ExAC_AF_OTH", "ExAC_AF_SAS", "ExAC_AF_Adj")
  
  #Keeping variants with any ExAC equal to zero or NA
  local_vep_no_exac <- local_vep_results
  for(ExAC in ExAC_names){
    print(ExAC)
    print(dim(local_vep_no_exac))
    local_vep_no_exac <- local_vep_no_exac[c(which(is.na(local_vep_no_exac[,ExAC])),
                                             which(local_vep_no_exac[,ExAC] == 0),
                                             which(local_vep_no_exac[,ExAC] == "")),]
    print(length(c(which(is.na(local_vep_no_exac[,ExAC])),
                   which(local_vep_no_exac[,ExAC] == 0),
                   which(local_vep_no_exac[,ExAC] == ""))))
    #print(dim(local_vep_no_exac))
  }
  
  ## ExAC_EAS_AF	: ExAC East Asian Allele frequency
  ## ExAC_AMR_AF	: ExAC American Allele frequency
  ## ExAC_NFE_AF	: ExAC Non-Finnish European Allele frequency
  ## ExAC_AFR_AF	: ExAC African/African American Allele frequency
  ## ExAC_AF	: ExAC Allele frequency in genotypes, for each ALT allele, in the same order as listed
  ## ExAC_FIN_AF	: ExAC Finnish Allele frequency
  ## ExAC_OTH_AF	: ExAC Other Allele frequency
  ## ExAC_SAS_AF	: ExAC South Asian Allele frequency
  return(local_vep_no_exac)
}

exclude_1000_genomes <- function(local_vep_no_exac){
  #Exclude variants with a maf in 1000 Genomes different to zero and NA
  
  ## GMAF : Minor allele and frequency of existing variant in 1000 Genomes Phase 1 combined population
  ## AFR_MAF : Frequency of existing variant in 1000 Genomes Phase 1 combined African population
  ## AMR_MAF : Frequency of existing variant in 1000 Genomes Phase 1 combined American population
  ## ASN_MAF : Frequency of existing variant in 1000 Genomes Phase 1 combined Asian population
  ## EUR_MAF : Frequency of existing variant in 1000 Genomes Phase 1 combined European population
  ## AA_MAF : Frequency of existing variant in NHLBI-ESP African American population
  ## EA_MAF : Frequency of existing variant in NHLBI-ESP European American population
  
  thousand_genomes <- c("GMAF", "EUR_MAF","AFR_MAF","ASN_MAF","AMR_MAF","AA_MAF","EA_MAF","EAS_MAF","SAS_MAF")
  
  local_vep_no_common <- local_vep_no_exac
  for(thousand in thousand_genomes){
    print(thousand)
    print(dim(local_vep_no_common))
    local_vep_no_common <- local_vep_no_common[c(which(is.na(local_vep_no_common[,thousand])),
                                                 which(local_vep_no_common[,thousand] == 0),
                                                 which(local_vep_no_common[,thousand] == "")),]
    print(length(c(which(is.na(local_vep_no_common[,thousand])),
                   which(local_vep_no_common[,thousand] == 0),
                   which(local_vep_no_common[,thousand] == ""))))
    print(dim(local_vep_no_common))
  }
  
  #33456 variants remain. 
  #ExAc filtering removed 3757 variants, 
  #thousand genomes and NHLBI-ESP removed another 102 variants
  return(local_vep_no_common)
}


load_mutations <- function(){
  load("gene_mut_properties.Rdata")
  mut_to_exclude <- read_csv("Genes to exclude.txt", 
                             col_names = FALSE)
  mut_to_exclude <- unname(unlist(mut_to_exclude[,1]))
  
  mutations_df <- vector()
  for(tumour in tumours){
    local_mut <- gene_mut_properties[[tumour]]
    local_mut <- local_mut[!(local_mut$Hugo_Symbol %in% mut_to_exclude),]
    
    ##At least 3 patients with actual mutation
    #local_mut_top <- local_mut[local_mut$Number_mutations >=3,]
    #local_mut_top <- local_mut_top[local_mut_top$Syn_ratio >= 1,]
    
    local_mut$Replication_time <- NULL
    local_mut$HiC <- NULL
    local_mut$Gene_length <- NULL
    local_mut$GC_percentage <- NULL
    local_mut$LoF_MB <- NULL
    local_mut$Syn_MB <- NULL
    
    local_mut$Gene_age <- genes_phy[match(local_mut$Hugo_Symbol, genes_phy$GeneID), "Phylostrata"]
    local_mut$Tumour <- tumour
    mutations_df <- rbind(mutations_df, local_mut)
    
  }
  return(mutations_df)
}

calculate_rank_of_fraction_phy_mut <- function(n_mut_melt){
  n_mut_melt_rank <- vector()
  for(tumour in tumours){
    local_df <- n_mut_melt[n_mut_melt$Tumour == tumour,]
    local_df$Number_ranked <- rank(-(local_df$Number))
    local_df$Phy_class <- paste(local_df$Phylostrata, local_df$Mutations, sep="-")
    n_mut_melt_rank <- rbind(n_mut_melt_rank, local_df)
  }
  return(n_mut_melt_rank)
}

median_fraction_of_phy_mut <- function(n_mut_melt){
  n_mut_melt$Phy_class <- paste(n_mut_melt$Phylostrata, n_mut_melt$Mutations, sep="_")
  n_mut_melt$Phy_class <- as.factor(n_mut_melt$Phy_class)
  
  median_norm <- aggregate(Number_ranked ~ Phylostrata+Mutations+Phy_class, n_mut_melt, median)
  median_norm$Median <- "Median"
  median_norm$Mutations <- factor(median_norm$Mutations, levels=c("N_enriched_lof", "N_enriched_miss"))
  return(median_norm)  
}

median_fraction_of_phy_CNVs <- function(local_frequency){
  median_amp <- aggregate(Frac_amp2_rank ~ Phylostrata, local_frequency, median)  
  
  median_del <- aggregate(Frac_del2_rank ~ Phylostrata, local_frequency, median)  
  
  median_df <- data.frame(Phylostrata=median_amp$Phylostrata,
                          Median_amp = median_amp$Frac_amp2_rank,
                          Median_del = median_del$Frac_del2_rank,
                          Group="ALL")
  return(median_df)
}



load_CNVs <- function(only_focal){
  CNVs_df <- vector()
  for(tumour in tumours){
    local_amp <- CNVs_curated[[tumour]]$amplifications
    local_del <- CNVs_curated[[tumour]]$deletions
    
    local_amp$Sample <- substr(local_amp$Sample, 9, 12)
    local_del$Sample <- substr(local_del$Sample, 9, 12)
    
    if(only_focal == "FOCAL"){
      focal_amp <- as.character(CNVs_above_fraction_0.25[[tumour]][["Amp"]])
      focal_del <- as.character(CNVs_above_fraction_0.25[[tumour]][["Del"]])
      
      local_amp <- local_amp[local_amp$Gene %in% focal_amp,]
      local_del <- local_del[local_del$Gene %in% focal_del,]
    }else if(only_focal == "GLOBAL"){
      
      focal_amp <- as.character(CNVs_above_fraction_0.25[[tumour]][["Amp"]])
      focal_del <- as.character(CNVs_above_fraction_0.25[[tumour]][["Del"]])
      
      local_amp <- local_amp[!(local_amp$Gene %in% focal_amp),]
      local_del <- local_del[!(local_del$Gene %in% focal_del),]
    }else if(only_focal=="ANY"){
      #No subsetting
    }
    
    n_pat_amp <- aggregate(Segment_Mean ~ Gene, local_amp, length)
    colnames(n_pat_amp) <- c("Gene", "N_patients_amp")
    
    n_pat_del <- aggregate(Segment_Mean ~ Gene, local_del, length)
    colnames(n_pat_del) <- c("Gene", "N_patients_del")
    
    all_genes <- unique(c(n_pat_amp$Gene, n_pat_del$Gene))
    
    n_pat_df <- data.frame(Genes = all_genes,
                           Patients_amp =n_pat_amp[match(all_genes, n_pat_amp$Gene), "N_patients_amp"],
                           Patients_del =n_pat_del[match(all_genes, n_pat_del$Gene), "N_patients_del"])
    
    local_patients <- length(patients_with_CNV_info[[tumour]])
    n_pat_df$Patients_amp <- n_pat_df$Patients_amp/local_patients
    n_pat_df$Patients_del <- n_pat_df$Patients_del/local_patients
    
    n_pat_df$Patients_amp[is.na(n_pat_df$Patients_amp)] <- 0
    n_pat_df$Patients_del[is.na(n_pat_df$Patients_del)] <- 0
    
    #hist(n_pat_df$Patients_amp)
    #hist(n_pat_df$Patients_del)
    
    
    n_pat_df$Gene_age <- genes_phy[match(n_pat_df$Genes, genes_phy$GeneID), "Phylostrata"]
    #n_pat_df <- n_pat_df[!is.na(n_pat_df$Gene_age),]
    n_pat_df$Gene_age <- factor(n_pat_df$Gene_age, levels=1:16)
    n_pat_df$Tumour <- tumour
    CNVs_df <- rbind(CNVs_df, n_pat_df)
  }
  return(CNVs_df)
}

calculate_fraction_of_patients_with_CNV <- function(CNVs_df){
  
  patient_frequency_CNV_df <- vector()
  for(tumour in tumours){
    n_pat_df <- CNVs_df[CNVs_df$Tumour == tumour,]
      
    #3 categories: <0.10, 0.10-0.20, >0.20
    n_pat_df$Pat_amp_cat <- ifelse(n_pat_df$Patients_amp > 0 & n_pat_df$Patients_amp <= 0.1, "<0.10",
                                   ifelse(n_pat_df$Patients_amp > 0.10 & n_pat_df$Patients_amp <= 0.20, "0.10<x<0.20",
                                          ifelse(n_pat_df$Patients_amp > 0.20, ">0.20", NA)))
    n_pat_df$Pat_del_cat <- ifelse(n_pat_df$Patients_del > 0 & n_pat_df$Patients_del <= 0.1, "<0.10",
                                   ifelse(n_pat_df$Patients_del > 0.10 & n_pat_df$Patients_del <= 0.20, "0.10<x<0.20",
                                          ifelse(n_pat_df$Patients_del > 0.20, ">0.20", NA)))
    
    ##2 categories: < 0.10, >0.10
    n_pat_df$Pat_amp_cat <- ifelse(n_pat_df$Patients_amp > 0 & n_pat_df$Patients_amp <= 0.1, "<0.10",
                                   ifelse(n_pat_df$Patients_amp > 0.10, ">0.10", NA))
    n_pat_df$Pat_del_cat <- ifelse(n_pat_df$Patients_del > 0 & n_pat_df$Patients_del <= 0.1, "<0.10",
                                   ifelse(n_pat_df$Patients_del > 0.10, ">0.10", NA))
    
    
    n_pat_cat_amp_df <- aggregate(Gene_age ~ Pat_amp_cat, n_pat_df, table)
    n_pat_cat_amp_df <- data.frame(Frequency = n_pat_cat_amp_df$Pat_amp_cat, n_pat_cat_amp_df$Gene_age)
    colnames(n_pat_cat_amp_df)[2:17] <- paste("Phy", 1:16, sep="_")
    n_pat_cat_amp_df_melt <- melt(n_pat_cat_amp_df)
    
    n_pat_cat_del_df <- aggregate(Gene_age ~ Pat_del_cat, n_pat_df, table)
    n_pat_cat_del_df <- data.frame(Frequency = n_pat_cat_del_df$Pat_del_cat, n_pat_cat_del_df$Gene_age)
    colnames(n_pat_cat_del_df)[2:17] <- paste("Phy", 1:16, sep="_")
    n_pat_cat_del_df_melt <- melt(n_pat_cat_del_df)
    
    #frequencies <- c("<0.05", "0.05<x<0.15", "0.15<x<0.25", ">0.25")
    #frequencies <- c("<0.10", "0.10<x<0.20", ">0.20")
    frequencies <- c("<0.10", ">0.10")
    
    patient_frequency_CNV <- data.frame(Tumour =tumour, Frequency = rep(frequencies, 16), Phylostrata = rep(paste("Phy", 1:16, sep="_"), each=length(frequencies)))
    
    patient_frequency_CNV$N_genes_amp <- n_pat_cat_amp_df_melt[match(paste(patient_frequency_CNV$Frequency, patient_frequency_CNV$Phylostrata),
                                                                     paste(n_pat_cat_amp_df_melt$Frequency, n_pat_cat_amp_df_melt$variable)),"value"]
    patient_frequency_CNV$N_genes_del <- n_pat_cat_del_df_melt[match(paste(patient_frequency_CNV$Frequency, patient_frequency_CNV$Phylostrata),
                                                                     paste(n_pat_cat_del_df_melt$Frequency, n_pat_cat_del_df_melt$variable)),"value"]
    patient_frequency_CNV$N_genes_amp[is.na(patient_frequency_CNV$N_genes_amp)] <- 0
    patient_frequency_CNV$N_genes_del[is.na(patient_frequency_CNV$N_genes_del)] <- 0
    
    patient_frequency_CNV$Phylostrata <- substr(patient_frequency_CNV$Phylostrata, 5,6)
    patient_frequency_CNV$Phylostrata <- factor(patient_frequency_CNV$Phylostrata, levels=1:16)
    
    patient_frequency_CNV$Frequency <- factor(patient_frequency_CNV$Frequency, levels=rev(frequencies))
    patient_frequency_CNV$Total_genes <- n_genes_phy[match(patient_frequency_CNV$Phylostrata, n_genes_phy$Phy), "Number"]
    
    patient_frequency_CNV$Fraction_genes_amp <- patient_frequency_CNV$N_genes_amp/patient_frequency_CNV$Total_genes
    patient_frequency_CNV$Fraction_genes_del <- patient_frequency_CNV$N_genes_del/patient_frequency_CNV$Total_genes
    
    patient_frequency_CNV_df <- rbind(patient_frequency_CNV_df, patient_frequency_CNV)
  } 
  return(patient_frequency_CNV_df)
}


add_mutation_information <- function(gene_mut_properties, interest){
  mutations_int <- vector()
  for(tumour in tumours){
    local_mut <- gene_mut_properties[[tumour]]
    
    local_mut_interest <- local_mut[local_mut$Hugo_Symbol %in% interest,]
    local_mut_interest <- local_mut_interest[,c("Hugo_Symbol", "Variant_type", "Number_mutations", 
                                                "Number_synonymous", "Syn_ratio")]
    local_mut_interest$Tumour <- tumour
    mutations_int <- rbind(mutations_int, local_mut_interest)
  }
  mutations_int_long <- vector()

  for(gene in unique(mutations_int$Hugo_Symbol)){
    temp <- mutations_int[mutations_int$Hugo_Symbol == gene,] 
    
    for(tumour in tumours){
      temp2 <- temp[temp$Tumour == tumour,]
      if(nrow(temp2) > 1){
        temp_df1 <- c(gene, "Missense", temp2[temp2$Variant_type=="Missense","Number_mutations"],
                      temp2[temp2$Variant_type=="Missense","Number_synonymous"],
                      temp2[temp2$Variant_type=="Missense","Syn_ratio"],
                      "LoF", temp2[temp2$Variant_type=="LoF","Number_mutations"],
                      temp2[temp2$Variant_type=="LoF","Number_synonymous"],
                      temp2[temp2$Variant_type=="LoF","Syn_ratio"],
                      temp2[temp2$Variant_type=="LoF","Tumour"])
        mutations_int_long <- rbind(mutations_int_long,
                                    temp_df1)
        
      }else if(nrow(temp2) ==1){
        temp_df2 <- c(gene, as.character(temp2$Variant_type), temp2$Number_mutations,
                      temp2$Number_synonymous, temp2$Syn_ratio, NA, NA, NA, NA, temp2$Tumour)
        mutations_int_long <- rbind(mutations_int_long,
                                    temp_df2)
      }
    }
  }
  colnames(mutations_int_long) <- c("Gene", "Variant_type_1", "Number_mutations_1",
                                    "Number_synonymous_1", "Syn_ratio_1",
                                    "Variant_type_2", "Number_mutations_2",
                                    "Number_synonymous_2", "Syn_ratio_2",
                                    "Tumour")
  
  mutations_int_long <- as.data.frame(mutations_int_long)
  mutations_int_long$Name <- paste(mutations_int_long$Gene,
                                   mutations_int_long$Tumour, sep="_")
  return(mutations_int_long)
}

add_CNV_information <- function(CNVs_curated, patients_with_CNV_info, use_selected_CNVs){
  CNVs_int <- list()
  for(tumour in tumours){
    local_amp <- CNVs_curated[[tumour]]$amplifications
    local_del <- CNVs_curated[[tumour]]$deletions
    
    if(use_selected_CNVs == TRUE){
      amp_above_context_cutoff <- as.character(CNVs_above_fraction[[tumour]][["Amp"]])
      del_above_context_cutoff <- as.character(CNVs_above_fraction[[tumour]][["Del"]])
      
      local_amp <- local_amp[local_amp$Gene %in% amp_above_context_cutoff,]
      local_del <- local_del[local_del$Gene %in% del_above_context_cutoff,]
    }

    
    local_amp$Sample <- substr(local_amp$Sample, 9, 12)
    local_del$Sample <- substr(local_del$Sample, 9, 12)
    
    n_pat_amp <- aggregate(Segment_Mean ~ Gene, local_amp, length)
    colnames(n_pat_amp) <- c("Gene", "N_patients_amp")
    
    n_pat_del <- aggregate(Segment_Mean ~ Gene, local_del, length)
    colnames(n_pat_del) <- c("Gene", "N_patients_del")
    
    all_genes <- unique(c(n_pat_amp$Gene, n_pat_del$Gene))
    
    n_pat_df <- data.frame(Genes = all_genes,
                           Patients_amp =n_pat_amp[match(all_genes, n_pat_amp$Gene), "N_patients_amp"],
                           Patients_del =n_pat_del[match(all_genes, n_pat_del$Gene), "N_patients_del"])
    
    local_patients <- length(patients_with_CNV_info[[tumour]])
    n_pat_df$Patients_amp <- n_pat_df$Patients_amp/local_patients
    n_pat_df$Patients_del <- n_pat_df$Patients_del/local_patients
    
    n_pat_df$Patients_amp[is.na(n_pat_df$Patients_amp)] <- 0
    n_pat_df$Patients_del[is.na(n_pat_df$Patients_del)] <- 0
    
    n_pat_df$Tumour <- tumour
    n_pat_df_int <- n_pat_df[n_pat_df$Genes %in% interest,]
    CNVs_int <- rbind(CNVs_int, n_pat_df_int)
  }
  CNVs_int$Name <- paste(CNVs_int$Genes, CNVs_int$Tumour, sep="_")
  return(CNVs_int)
}

add_neighbourhood <- function(interest_df){
  #In all of PathwayCommons
  load("percentage_UC_per_gene_PathwayCommons.Rdata")
  interest_df$Neigh_PC <- percentage_UC_per_gene[match(interest_df$Genes, names(percentage_UC_per_gene))]
  
  #In Biogrid
  load("percentage_UC_per_gene_BioGRID.Rdata")
  interest_df$Neigh_BG <- percentage_UC_per_gene[match(interest_df$Genes, names(percentage_UC_per_gene))]
  
  #Of all in the regulatory network
  load("percentage_UC_directed_all.Rdata")
  interest_df$Neigh_reg_PC <- percentage_UC_directed_all[match(interest_df$Genes, rownames(percentage_UC_directed_all)), "Per_UC"]
  
  #Of sources
  load("percentage_UC_per_gene_directed_PC.Rdata")
  interest_df$Neigh_sources_PC <- percentage_UC_per_gene_directed[match(interest_df$Genes, rownames(percentage_UC_per_gene_directed)), "Per_UC"]
  return(interest_df)
}

add_degree <- function(interest_df){
  #All of PathwayCommons
  network <- read.delim("PathwayCommons9.All.hgnc.sif", header=FALSE, dec=",")
  network <- network[!(network[,2] %in% c("controls-production-of", "controls-transport-of-chemical",
                                          "chemical-affects", "consumption-controlled-by", "reacts-with",
                                          "used-to-produce")),]  ##To remove all reactions involving chemicals
  network <- network[,c(1,3,2)]
  colnames(network) <- c("Gene1", "Gene2", "Edge_type")
  network_graph <- graph.data.frame(network, directed=FALSE)
  network_graph <- simplify(network_graph, remove.loops=TRUE, remove.multiple=TRUE)
  degree_network <- degree(network_graph)
  
  interest_df$Degree_PC <- degree_network[match(interest_df$Genes, names(degree_network))]
  
  ##Regulatory of PC
  network_reg <- network[network$Edge_type == "controls-expression-of",]
  network_graph_reg <- graph.data.frame(network_reg, directed=FALSE)
  network_graph_reg <- simplify(network_graph_reg, remove.loops=TRUE, remove.multiple=TRUE)
  degree_network_reg <- degree(network_graph_reg)
  
  interest_df$Degree_PC_regulatory <- degree_network_reg[match(interest_df$Genes, names(degree_network_reg))]
  
  #Sources of PC
  degree_sources <- table(network_reg$Gene1)
  interest_df$Degree_PC_sources <- degree_sources[match(interest_df$Genes, names(degree_sources))]
  
  #Biogrid
  network <- read.delim("BIOGRID-ORGANISM-Homo_sapiens-3.4.152.tab2.txt")
  network <- network[,c("Official.Symbol.Interactor.A", "Official.Symbol.Interactor.B")]
  colnames(network) <- c("Gene1", "Gene2")  
  
  network_graph <- graph.data.frame(network, directed=FALSE)
  network_graph <- simplify(network_graph, remove.loops=TRUE, remove.multiple=TRUE)
  
  degree_network <- degree(network_graph)
  
  interest_df$Degree_BG <- degree_network[match(interest_df$Genes, names(degree_network))]
  
  ##Edges ages that are broken up
  network_reg$Age1 <- ifelse(network_reg$Gene1 %in% UC_genes, "UC", 
                             ifelse(network_reg$Gene1 %in% EM_genes, "EM", 
                                    ifelse(network_reg$Gene1 %in% MM_genes, "MM", NA)))
  network_reg$Age2 <- ifelse(network_reg$Gene2 %in% UC_genes, "UC", 
                             ifelse(network_reg$Gene2 %in% EM_genes, "EM", 
                                    ifelse(network_reg$Gene2 %in% MM_genes, "MM", NA)))
  
  network_reg <- network_reg[!is.na(network_reg$Age1),]
  network_reg <- network_reg[!is.na(network_reg$Age2),]
  
  network_reg$Edge_age <- paste(network_reg$Age1, network_reg$Age2, sep="_")
  
  percentage_edges_affected <- vector()
  for(gene in as.character(unique(interest_df$Genes))){
    
    local_network_reg <- rbind(network_reg[network_reg$Gene1 == gene,],
                               network_reg[network_reg$Gene2 == gene,])
    edge_ages <- local_network_reg$Edge_age
    
    percentage_edges <- c(gene, (sum(edge_ages == "UC_UC")/length(edge_ages))*100,
                          (sum(edge_ages == "EM_EM")/length(edge_ages))*100,
                          (sum(edge_ages == "MM_MM")/length(edge_ages))*100,
                          (sum(edge_ages == "UC_EM")/length(edge_ages))*100,
                          (sum(edge_ages == "EM_UC")/length(edge_ages))*100,
                          (sum(edge_ages == "UC_MM")/length(edge_ages))*100,
                          (sum(edge_ages == "MM_UC")/length(edge_ages))*100,
                          (sum(edge_ages == "EM_MM")/length(edge_ages))*100,
                          (sum(edge_ages == "MM_EM")/length(edge_ages))*100)
    
    percentage_edges_affected <- rbind(percentage_edges_affected, 
                                       percentage_edges)
  }
  
  colnames(percentage_edges_affected) <- c("Gene", "UC_UC","EM_EM","MM_MM",
                                           "UC_EM","EM_UC","UC_MM","MM_UC",
                                           "EM_MM","MM_EM")
  rownames(percentage_edges_affected) <- NULL
  percentage_edges_affected <- as.data.frame(percentage_edges_affected)
  gene_names <- unname(as.character(percentage_edges_affected$Gene))
  percentage_edges_affected <- percentage_edges_affected[,-1]
  percentage_edges_affected <- as.data.frame(apply(percentage_edges_affected, 2, function(x){ as.numeric(as.character(x)) }))
  
  if(sum(is.na(gene_names)) > 0){
    percentage_edges_affected <- percentage_edges_affected[-which(is.na(gene_names)),]
  }
  
  gene_names <- gene_names[!is.na(gene_names)]
  rownames(percentage_edges_affected) <- gene_names
  
  percentage_edges_affected$UC_EM_EM_UC <- percentage_edges_affected$UC_EM/percentage_edges_affected$EM_UC
  percentage_edges_affected$UC_MM_MM_UC <- percentage_edges_affected$UC_MM/percentage_edges_affected$MM_UC
  percentage_edges_affected$EM_MM_MM_EM <- percentage_edges_affected$EM_MM/percentage_edges_affected$MM_EM
  
  interest_df <- data.frame(interest_df,
                            percentage_edges_affected[match(interest_df$Genes, rownames(percentage_edges_affected)),])
  
  return(interest_df)  
}

add_CNV_expression <- function(CNV_p_values_df, interest_df){
  CNV_p_values_df$Name <- paste(CNV_p_values_df$Gene, CNV_p_values_df$Tumour, sep="_")
  CNV_p_values_df_amp <- CNV_p_values_df[CNV_p_values_df$CNV == "Amp",]
  CNV_p_values_df_del <- CNV_p_values_df[CNV_p_values_df$CNV == "Del",]
  colnames(CNV_p_values_df_amp)[6] <- "Amp_p"
  colnames(CNV_p_values_df_del)[6] <- "Del_p"
  
  interest_df$Amp_p <- CNV_p_values_df_amp[match(interest_df$Name, CNV_p_values_df_amp$Name), "Amp_p"]
  interest_df$Del_p <- CNV_p_values_df_del[match(interest_df$Name, CNV_p_values_df_del$Name), "Del_p"]
  return(interest_df)  
}

add_mutation_expression <- function(p_lof_tumour, p_miss_tumour, interest_df){
  LoF_p_values_df <- vector()
  for(tumour in tumours){
    LoF_p_values_df <- rbind(LoF_p_values_df,
                             cbind(Tumour=tumour, Mut_type="LoF", 
                                   p_lof_tumour[[tumour]][,c("Gene", "p_two_sided")]))
    
  }
  colnames(LoF_p_values_df)[4] <- "LoF_p"
  
  miss_p_values_df <- vector()
  for(tumour in tumours){
    miss_p_values_df <- rbind(miss_p_values_df,
                              cbind(Tumour=tumour, Mut_type="Miss", 
                                    p_miss_tumour[[tumour]][,c("Gene", "p_two_sided")]))
  }
  colnames(miss_p_values_df)[4] <- "Miss_p"
  
  mut_p_values_df <- data.frame(miss_p_values_df, LoF_p=NA)
  mut_p_values_df <- rbind(mut_p_values_df, data.frame(LoF_p_values_df, Miss_p=NA))
  
  mut_p_values_df$Name <- paste(mut_p_values_df$Gene, mut_p_values_df$Tumour, sep="_")
  
  interest_df$Mut_Miss_exp <- mut_p_values_df[match(interest_df$Name, mut_p_values_df$Name), "Miss_p"]
  interest_df$Mut_LoF_exp <- mut_p_values_df[match(interest_df$Name, mut_p_values_df$Name), "LoF_p"]
  
  return(interest_df)  
}

prepare_trends_of_genes <- function(interest, use_selected_CNVs){
  interest_df <- data.frame(Gene=interest,
                            Age=genes_phy[match(interest, genes_phy$GeneID), "Phylostrata"])
  
  mutations_int_long <- add_mutation_information(gene_mut_properties, interest)
  print("Done adding mutations")
  
  ##Read in CNVs
  
  CNVs_int <- add_CNV_information(CNVs_curated, patients_with_CNV_info, use_selected_CNVs)
  print("Done adding CNVs")
  
  all_names <- unique(sort(c(mutations_int_long$Name, CNVs_int$Name)))
  
  alt_df <- data.frame(Name = all_names)
  alt_df <- data.frame(alt_df, CNVs_int[match(alt_df$Name, CNVs_int$Name),])
  
  alt_df <- data.frame(alt_df,
                       mutations_int_long[match(alt_df$Name, mutations_int_long$Name),])
  
  alt_df <- alt_df[,c("Genes","Patients_amp","Patients_del",
                      "Variant_type_1","Number_mutations_1",
                      "Number_synonymous_1","Syn_ratio_1",
                      "Variant_type_2","Number_mutations_2",
                      "Number_synonymous_2","Syn_ratio_2","Tumour", "Name")]
  
  alt_df$Patients_amp <- alt_df$Patients_amp*100
  alt_df$Patients_del <- alt_df$Patients_del*100
  
  interest_df <- data.frame(alt_df, Age=interest_df[match(alt_df$Genes, interest_df$Gene), "Age"])
  
  interest_df <- add_degree(interest_df)
  print("Done adding degree")
  
  interest_df <- add_neighbourhood(interest_df)
  print("Done adding neighbourhood")
  
  load("CNV_exp_p_values_df.Rdata")
  CNV_p_values_df <- p_values_df
  
  interest_df <- add_CNV_expression(CNV_p_values_df, interest_df)
  print("Done adding CNV expression data")
  
  load("p_lof_tumour.Rdata")
  load("p_miss_tumour.Rdata")
  
  interest_df <- add_mutation_expression(p_lof_tumour, p_miss_tumour, interest_df)  
  print("Done adding mutation expression data")
  
  return(interest_df)
} 

correct_ratios_with_baseline <- function(ratios_alterations, ratios_baseline){
  ratios_melt <- melt(ratios_alterations)
  colnames(ratios_melt)[2:3] <- c("Comparison", "Ratio")
  
  ratios_melt$Dominated <- ifelse(ratios_melt$Ratio > 1, "Dominated", "Not")
  
  ratios_melt$Corrected_ratios <- sapply(1:nrow(ratios_melt), function(i){
    comp <- ratios_melt$Comparison[i]
    ratio <- ratios_melt$Ratio[i]
    baseline <- ratios_baseline$Values[ratios_baseline$Ratios==comp]  
    return(ratio/baseline)
  })  
  return(ratios_melt)
}

add_mut_CNV_to_degree_df <- function(alt_degree_df, only_miss, only_lof, only_amp, only_del){
  alt_degree_df$Mut <- NA
  alt_degree_df$Mut <- ifelse(alt_degree_df$Genes %in% only_miss, "Miss",
                              ifelse(alt_degree_df$Genes %in% only_lof, "LoF", NA))
  
  alt_degree_df$CNV <- NA
  alt_degree_df$CNV <- ifelse(alt_degree_df$Genes %in% only_amp, "Amp",
                              ifelse(alt_degree_df$Genes %in% only_del, "Del", NA))
  
  alt_degree_df$WT <- NA
  alt_degree_df$WT <- apply(alt_degree_df, 1, function(x){
    if(is.na(x[5]) & is.na(x[6])){
      return("WT")
    }else{
      return(NA)
    }
  })
  
  both_mutations <- only_miss[only_miss %in% only_lof]
  both_CNVs <- only_amp[only_amp %in% only_del]
  
  alt_degree_df$WT[alt_degree_df$Genes %in% both_mutations] <- NA
  alt_degree_df$WT[alt_degree_df$Genes %in% both_CNVs] <- NA
  
  alt_degree_df1 <- alt_degree_df[,c("Genes","Degree","Database","Age","Mut")]
  alt_degree_df2 <- alt_degree_df[,c("Genes","Degree","Database","Age","CNV")]
  alt_degree_df3 <- alt_degree_df[,c("Genes","Degree","Database","Age","WT")]
  colnames(alt_degree_df1)[5] <- "Alt"
  colnames(alt_degree_df2)[5] <- "Alt"
  colnames(alt_degree_df3)[5] <- "Alt"
  
  alt_degree_df_melt <- rbind(alt_degree_df1, alt_degree_df2, alt_degree_df3)
  
  alt_degree_df_melt <- alt_degree_df_melt[!is.na(alt_degree_df_melt$Alt),]
  return(alt_degree_df_melt)
}

count_edges_affected_by_alt <- function(tumour_network,
                                        local_all_miss, local_all_lof, local_all_mut_both,
                                        local_all_amp, local_all_del, local_all_CNV_both){
  tumour_network$Mut1 <- NA
  tumour_network$Mut1[tumour_network$Gene1 %in% local_all_miss] <- "Miss"
  tumour_network$Mut1[tumour_network$Gene1 %in% local_all_lof] <- "LoF"
  tumour_network$Mut1[tumour_network$Gene1 %in% local_all_mut_both] <- "Miss&LoF"
  
  tumour_network$Mut2 <- NA
  tumour_network$Mut2[tumour_network$Gene2 %in% local_all_miss] <- "Miss"
  tumour_network$Mut2[tumour_network$Gene2 %in% local_all_lof] <- "LoF"
  tumour_network$Mut2[tumour_network$Gene2 %in% local_all_mut_both] <- "Miss&LoF"
  
  tumour_network$Mut <- paste(tumour_network$Mut1, tumour_network$Mut2, sep="_")
  
  tumour_network$CNV1 <- NA
  tumour_network$CNV1[tumour_network$Gene1 %in% local_all_amp] <- "Amp"
  tumour_network$CNV1[tumour_network$Gene1 %in% local_all_del] <- "Del"
  tumour_network$CNV1[tumour_network$Gene1 %in% local_all_CNV_both] <- "Amp&Del"
  
  tumour_network$CNV2 <- NA
  tumour_network$CNV2[tumour_network$Gene2 %in% local_all_amp] <- "Amp"
  tumour_network$CNV2[tumour_network$Gene2 %in% local_all_del] <- "Del"
  tumour_network$CNV2[tumour_network$Gene2 %in% local_all_CNV_both] <- "Amp&Del"  
  tumour_network$CNV <- paste(tumour_network$CNV1, tumour_network$CNV2, sep="_")
  
  
  tumour_network_sources <- tumour_network[, c("Edge_age", "Mut", "CNV")]
  tumour_network_sources$Edge_mut <- do.call('cbind', strsplit(tumour_network_sources$Mut, "_"))[1,]
  tumour_network_sources$Edge_CNV <- do.call('cbind', strsplit(tumour_network_sources$CNV, "_"))[1,]
  
  tumour_network_targets <- tumour_network[, c("Edge_age", "Mut", "CNV")]
  tumour_network_targets$Edge_mut <- do.call('cbind', strsplit(tumour_network_targets$Mut, "_"))[2,]
  tumour_network_targets$Edge_CNV <- do.call('cbind', strsplit(tumour_network_targets$CNV, "_"))[2,]
  
  tumour_network_either <- tumour_network[, c("Edge_age", "Mut", "CNV")]
  tumour_network_either$Edge_mut <- tumour_network_either$Mut
  tumour_network_either$Edge_mut <- gsub("_", ",", tumour_network_either$Edge_mut)
  tumour_network_either$Edge_mut <- gsub("NA,", "", tumour_network_either$Edge_mut)
  tumour_network_either$Edge_mut <- gsub(",NA", "", tumour_network_either$Edge_mut)
  
  tumour_network_either$Edge_mut[tumour_network_either$Edge_mut == "Miss,Miss"] <- "Miss"
  tumour_network_either$Edge_mut[tumour_network_either$Edge_mut == "LoF,LoF"] <- "LoF"
  
  tumour_network_either$Edge_mut[tumour_network_either$Edge_mut %in% c("Miss,Miss&LoF","Miss&LoF","Miss,LoF","Miss&LoF,Miss&LoF","Miss&LoF,Miss",
                                                                       "Miss&LoF,LoF","LoF,Miss&LoF","LoF,Miss")] <- "Miss_LoF"
  tumour_network_either$Edge_CNV <- tumour_network_either$CNV
  tumour_network_either$Edge_CNV <- gsub("_", ",", tumour_network_either$Edge_CNV)
  tumour_network_either$Edge_CNV <- gsub("NA,", "", tumour_network_either$Edge_CNV)
  tumour_network_either$Edge_CNV <- gsub(",NA", "", tumour_network_either$Edge_CNV)
  
  tumour_network_either$Edge_CNV[tumour_network_either$Edge_CNV == "Amp,Amp"] <- "Amp"
  tumour_network_either$Edge_CNV[tumour_network_either$Edge_CNV == "Del,Del"] <- "Del"
  tumour_network_either$Edge_CNV[tumour_network_either$Edge_CNV %in% c("Amp,Del","Del,Amp","Amp&Del","Amp,Amp&Del",
                                                                       "Del,Amp&Del","Amp&Del,Del","Amp&Del,Amp&Del","Amp&Del,Amp")] <- "Amp_Del"
  return(list(sources=tumour_network_sources, 
              targets=tumour_network_targets, 
              either=tumour_network_either))
}

count_edges_affected_by_alt_binary <- function(tumour_network,
                                               recurrent_point, recurrent_CNV){
  tumour_network$Mut1 <- NA
  tumour_network$Mut1[tumour_network$Gene1 %in% recurrent_point] <- "Point"
  
  tumour_network$Mut2 <- NA
  tumour_network$Mut2[tumour_network$Gene2 %in% recurrent_point] <- "Point"
  
  tumour_network$Mut <- paste(tumour_network$Mut1, tumour_network$Mut2, sep="_")
  
  tumour_network$CNV1 <- NA
  tumour_network$CNV1[tumour_network$Gene1 %in% recurrent_CNV] <- "CNV"
  
  tumour_network$CNV2 <- NA
  tumour_network$CNV2[tumour_network$Gene2 %in% recurrent_CNV] <- "CNV"
  
  tumour_network$CNV <- paste(tumour_network$CNV1, tumour_network$CNV2, sep="_")
  
  tumour_network_sources <- tumour_network[, c("Edge_age", "Mut", "CNV")]
  tumour_network_sources$Edge_mut <- do.call('cbind', strsplit(tumour_network_sources$Mut, "_"))[1,]
  tumour_network_sources$Edge_CNV <- do.call('cbind', strsplit(tumour_network_sources$CNV, "_"))[1,]
  
  tumour_network_targets <- tumour_network[, c("Edge_age", "Mut", "CNV")]
  tumour_network_targets$Edge_mut <- do.call('cbind', strsplit(tumour_network_targets$Mut, "_"))[2,]
  tumour_network_targets$Edge_CNV <- do.call('cbind', strsplit(tumour_network_targets$CNV, "_"))[2,]
  
  tumour_network_either <- tumour_network[, c("Edge_age", "Mut", "CNV")]
  tumour_network_either$Edge_mut <- NA
  tumour_network_either$Edge_mut[grep("Point", tumour_network_either$Mut)] <- "Point"
  
  tumour_network_either$Edge_CNV <- NA
  tumour_network_either$Edge_CNV[grep("CNV", tumour_network_either$CNV)] <- "CNV"
  
  return(list(sources=tumour_network_sources, 
              targets=tumour_network_targets, 
              either=tumour_network_either))
}

load_RN_network <- function(){
  network <- read.delim("PathwayCommons9.All.hgnc.sif", header=FALSE, dec=",")
  network <- network[!(network[,2] %in% c("controls-production-of", "controls-transport-of-chemical",
                                          "chemical-affects", "consumption-controlled-by", "reacts-with",
                                          "used-to-produce")),]  ##To remove all reactions involving chemicals
  network <- network[,c(1,3,2)]
  colnames(network) <- c("Gene1", "Gene2", "Edge_type")
  
  network <- network[network$Edge_type == "controls-expression-of",]
  
  network_graph <- graph.data.frame(network, directed=TRUE)
  network_graph <- simplify(network_graph, remove.loops=TRUE, remove.multiple=TRUE)
  
  network <- get.data.frame(network_graph)
  colnames(network) <- c("Gene1", "Gene2")
  network$Age1 <- ifelse(network$Gene1 %in% UC_genes, "UC", 
                         ifelse(network$Gene1 %in% EM_genes, "EM", 
                                ifelse(network$Gene1 %in% MM_genes, "MM", NA)))
  network$Age2 <- ifelse(network$Gene2 %in% UC_genes, "UC", 
                         ifelse(network$Gene2 %in% EM_genes, "EM", 
                                ifelse(network$Gene2 %in% MM_genes, "MM", NA)))
  
  network <- network[!is.na(network$Age1),]
  network <- network[!is.na(network$Age2),]
  
  network$Edge_age <- paste(network$Age1, network$Age2, sep="_")
  return(network)
}

calculate_percentage_downstream_genes_of_each_age <- function(RN_network){
  regulators <- unique(RN_network$Gene1)
  
  percentage_connections_each_age <- vector()
  for(regulator in regulators){
    local_RN_network <- RN_network[RN_network$Gene1 == regulator,]
    regulator_age <- unique(local_RN_network$Age1)
    result <- table(local_RN_network$Age2)
    result <- result[c("UC", "EM", "MM")]
    result <- c(result, Total=sum(result, na.rm=TRUE), c((result/sum(result, na.rm=TRUE)))*100)
    result <- c(Regulator=regulator, Age=regulator_age, result)
    percentage_connections_each_age <- rbind(percentage_connections_each_age, result)
    
  }
  percentage_connections_each_age <- as.data.frame(percentage_connections_each_age)
  colnames(percentage_connections_each_age) <- c("Regulator", "Age", "UC", "EM", "MM", "Total", "Per_UC", "Per_EM", "Per_MM")
  percentage_connections_each_age[,3:9] <- apply(percentage_connections_each_age[,3:9], 2,function(x){
    as.numeric(as.character(x)) } )
  
  percentage_connections_each_age <- percentage_connections_each_age[percentage_connections_each_age$Total > 1, ]
  
  percentage_connections_each_age_melt <- melt(percentage_connections_each_age[,c("Regulator", "Age", "Per_UC", "Per_EM", "Per_MM")])
  return(percentage_connections_each_age_melt)
}


divide_downstream_EM_into_quantiles <- function(per_EM_downstream_alt){
  temp <- per_EM_downstream_alt
  temp$Quantiles1 <- ifelse(per_EM_downstream_alt$Per_EM < 100/3, "<33%",
                 #ifelse(per_EM_downstream_alt$Per_EM < 50, "25-50%",
                 ifelse(per_EM_downstream_alt$Per_EM < (100/3)*2, "33-66%", ">66%"))
  
  quantiles1 <- aggregate(Genes ~ Alt_binary+Quantiles1, temp, length)
  colnames(quantiles1)[3] <- "Number_of_regulators"
  
  #Out of the total number of regulators
  quantiles1$Percentage_regulators <- (quantiles1$Number_of_regulators/sum(quantiles1$Number_of_regulators))*100
  
  #Out of the total number of regulators with point mutations
  quantiles1$Per_regulators_by_type <- NA
  total_regulators_point <- sum(quantiles1[quantiles1$Alt_binary == "Point", "Number_of_regulators"])
  quantiles1$Per_regulators_by_type[quantiles1$Alt_binary == "Point"] <- quantiles1[quantiles1$Alt_binary == "Point", "Number_of_regulators"]/total_regulators_point
  
  total_regulators_CNVs <- sum(quantiles1[quantiles1$Alt_binary == "CNV", "Number_of_regulators"])
  quantiles1$Per_regulators_by_type[quantiles1$Alt_binary == "CNV"] <- quantiles1[quantiles1$Alt_binary == "CNV", "Number_of_regulators"]/total_regulators_CNVs
  
  total_regulators_WT <- sum(quantiles1[quantiles1$Alt_binary == "WT", "Number_of_regulators"])
  quantiles1$Per_regulators_by_type[quantiles1$Alt_binary == "WT"] <- quantiles1[quantiles1$Alt_binary == "WT", "Number_of_regulators"]/total_regulators_WT
  
  quantiles1$Alt_binary <- factor(quantiles1$Alt_binary, levels=c("Point", "CNV", "WT"))
  
  quantiles1$Quantiles1 <- factor(quantiles1$Quantiles1, levels=c("<33%", "33-66%", ">66%"))
  
  return(quantiles1)  
}

divide_downstream_UC_EM_into_quantiles <- function(per_UC_EM_downstream_alt){
  temp <- per_UC_EM_downstream_alt
  temp$Quantiles1 <- ifelse(per_UC_EM_downstream_alt$Per_UC >= 200/3, "UC_downstream",
                            ifelse(per_UC_EM_downstream_alt$Per_EM >= 200/3, "EM_downstream", 
                                   ifelse(per_UC_EM_downstream_alt$Per_MM >= 200/3, "MM_downstream",
                                          ifelse(per_UC_EM_downstream_alt$Per_UC >= 10 & per_UC_EM_downstream_alt$Per_EM >= 10, "Mixed_downstream", "None"))))

  quantiles1 <- aggregate(Genes ~ Alt_binary+Quantiles1, temp, length)
  colnames(quantiles1)[3] <- "Number_of_regulators"
  
  #Out of the total number of regulators
  quantiles1$Percentage_regulators <- (quantiles1$Number_of_regulators/sum(quantiles1$Number_of_regulators))*100
  
  #Out of the total number of regulators with point mutations
  quantiles1$Per_regulators_by_type <- NA
  total_regulators_point <- sum(quantiles1[quantiles1$Alt_binary == "Point", "Number_of_regulators"])
  quantiles1$Per_regulators_by_type[quantiles1$Alt_binary == "Point"] <- quantiles1[quantiles1$Alt_binary == "Point", "Number_of_regulators"]/total_regulators_point
  
  total_regulators_CNVs <- sum(quantiles1[quantiles1$Alt_binary == "CNV", "Number_of_regulators"])
  quantiles1$Per_regulators_by_type[quantiles1$Alt_binary == "CNV"] <- quantiles1[quantiles1$Alt_binary == "CNV", "Number_of_regulators"]/total_regulators_CNVs
  
  total_regulators_WT <- sum(quantiles1[quantiles1$Alt_binary == "WT", "Number_of_regulators"])
  quantiles1$Per_regulators_by_type[quantiles1$Alt_binary == "WT"] <- quantiles1[quantiles1$Alt_binary == "WT", "Number_of_regulators"]/total_regulators_WT
  
  quantiles1$Alt_binary <- factor(quantiles1$Alt_binary, levels=c("Point", "CNV", "WT"))
  
  quantiles1$Quantiles1 <- factor(quantiles1$Quantiles1, levels=c("UC_downstream", "Mixed_downstream", "EM_downstream"))
  
  return(quantiles1)  
}


count_regulators_each_class <- function(per_UC_EM_downstream){
  quantiles1 <- aggregate(Genes ~ Alt_binary+Regulator_class, per_UC_EM_downstream, length)
  colnames(quantiles1)[3] <- "Number_of_regulators"
  
  quantiles1 <- quantiles1[quantiles1$Regulator_class %in% c("UC_downstream", "Mixed_downstream", "EM_downstream"),]
  
  #Out of the total number of regulators
  quantiles1$Percentage_regulators <- (quantiles1$Number_of_regulators/sum(quantiles1$Number_of_regulators))*100
  
  #Out of the total number of regulators with point mutations
  quantiles1$Per_regulators_by_type <- NA
  total_regulators_point <- sum(quantiles1[quantiles1$Alt_binary == "Point", "Number_of_regulators"])
  quantiles1$Per_regulators_by_type[quantiles1$Alt_binary == "Point"] <- quantiles1[quantiles1$Alt_binary == "Point", "Number_of_regulators"]/total_regulators_point
  
  total_regulators_CNVs <- sum(quantiles1[quantiles1$Alt_binary == "CNV", "Number_of_regulators"])
  quantiles1$Per_regulators_by_type[quantiles1$Alt_binary == "CNV"] <- quantiles1[quantiles1$Alt_binary == "CNV", "Number_of_regulators"]/total_regulators_CNVs
  
  total_regulators_WT <- sum(quantiles1[quantiles1$Alt_binary == "WT", "Number_of_regulators"])
  quantiles1$Per_regulators_by_type[quantiles1$Alt_binary == "WT"] <- quantiles1[quantiles1$Alt_binary == "WT", "Number_of_regulators"]/total_regulators_WT
  
  quantiles1$Alt_binary <- factor(quantiles1$Alt_binary, levels=c("Point", "CNV", "WT"))
  
  quantiles1$Regulator_class <- factor(quantiles1$Regulator_class, levels=c("UC_downstream", "Mixed_downstream", "EM_downstream", "MM_downstream", "None"))
  
  return(quantiles1)  
}

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

prepare_network <- function(network){
  network_graph <- graph.data.frame(network, directed=TRUE)
  network_graph <- simplify(network_graph, remove.loops=TRUE, remove.multiple=TRUE)
  
  network <- get.data.frame(network_graph)
  colnames(network) <- c("Gene1", "Gene2")
  
  network$Age1 <- ifelse(network$Gene1 %in% UC_genes, "UC", 
                         ifelse(network$Gene1 %in% EM_genes, "EM", 
                                ifelse(network$Gene1 %in% MM_genes, "MM", NA)))
  network$Age2 <- ifelse(network$Gene2 %in% UC_genes, "UC", 
                         ifelse(network$Gene2 %in% EM_genes, "EM", 
                                ifelse(network$Gene2 %in% MM_genes, "MM", NA)))
  
  network <- network[!is.na(network$Age1),]
  network <- network[!is.na(network$Age2),]
  return(network)
}

reproduce_results_other_databases_fig_1_2 <- function(network, database){
  
  sink(paste(database, "_Summary.txt", sep=""))
  
  network_graph <- graph.data.frame(network, directed=FALSE)
  degree_network <- degree(network_graph)
  degree_df <- data.frame(Genes=names(degree_network), Degree=degree_network,
                          Age = genes_phy[match(names(degree_network), genes_phy$GeneID), "Age"])
  
  degree_df$Age <- factor(degree_df$Age, levels=c("UC", "EM", "MM"))
  
  ##Divide by the median of each database
  median_per_db <- median(degree_df$Degree)
  
  degree_df$Median <- median_per_db
  degree_df$Norm_degree <- degree_df$Degree/degree_df$Median
  
  pdf(paste(database, "_Fig1.pdf", sep=""))
  g <- ggplot(degree_df, aes(x=Age, y =log10(Norm_degree)))+
    geom_boxplot(aes(fill=Age))+
    geom_hline(yintercept=0, color="grey")+
    geom_boxplot(aes(fill=Age))+
    xlab("Gene age")+
    ylab("Normalized Degree (log10)")+
    #facet_grid(.~Database)+
    theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  print(g)
  dev.off()
  
  
  degree_UC_genes <- degree_df[degree_df$Age == "UC", "Degree"]
  degree_EM_genes <- degree_df[degree_df$Age == "EM", "Degree"]
  degree_MM_genes <- degree_df[degree_df$Age == "MM", "Degree"]
  
  print("Difference in degree")
  print(wilcox.test(degree_EM_genes, degree_UC_genes, alternative="greater")$p.value)
  print(wilcox.test(degree_EM_genes, degree_MM_genes, alternative="greater")$p.value)
  
  
  ##In-degree in RN
  in_degree <- table(as.character(network$Gene2))
  in_degree_df <- data.frame(Genes=names(in_degree),
                             In_degree = as.vector(in_degree),
                             Age = genes_phy[match(names(in_degree), genes_phy$GeneID), "Age"])
  
  in_degree_df$Age <- factor(in_degree_df$Age, levels=c("UC", "EM", "MM"))
  
  pdf(paste(database, "_Fig2.pdf", sep=""))
  g <- ggplot(in_degree_df, aes(x=Age, y =log10(In_degree)))+
    geom_boxplot(aes(fill=Age))+
    theme_bw()
  print(g)
  dev.off()
  
  in_degree_UC <- in_degree_df[in_degree_df$Age == "UC", "In_degree"]
  in_degree_EM <- in_degree_df[in_degree_df$Age == "EM", "In_degree"]
  in_degree_MM <- in_degree_df[in_degree_df$Age == "MM", "In_degree"]
  
  print("Difference in in-degree")
  print(wilcox.test(in_degree_EM, in_degree_UC, alternative="greater")$p.value)
  print(wilcox.test(in_degree_EM, in_degree_MM, alternative="greater")$p.value)
  
  print(aggregate(In_degree ~ Age, in_degree_df, mean))
  
  
  summary_in_degree <- aggregate(In_degree ~ Age, in_degree_df, mean)
  colnames(summary_in_degree) <- c("Age", "Mean_indegree")
  sd_in_degree <- aggregate(In_degree ~ Age, in_degree_df, sd)
  summary_in_degree$SD <- sd_in_degree$In_degree
  
  pdf(paste(database, "_Fig3.pdf", sep=""),
      height=3, width=3)
  g <- ggplot(summary_in_degree, aes(x=Age, y=Mean_indegree))+
    geom_bar(aes(fill=Age), stat='identity')+
    ylab("Average number of\nincoming edges")+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(g)
  dev.off()
  
  
  out_degree <- table(as.character(network$Gene1))
  out_degree_df <- data.frame(Genes=names(out_degree),
                              out_degree = as.vector(out_degree),
                              Age = genes_phy[match(names(out_degree), genes_phy$GeneID), "Age"])
  
  out_degree_df$Age <- factor(out_degree_df$Age, levels=c("UC", "EM", "MM"))
  
  pdf(paste(database, "_Fig4.pdf", sep=""))
  g <- ggplot(out_degree_df, aes(x=Age, y =log10(out_degree)))+
    geom_boxplot(aes(fill=Age))+
    theme_bw()
  print(g)
  dev.off()
  
  out_degree_UC <- out_degree_df[out_degree_df$Age == "UC", "out_degree"]
  out_degree_EM <- out_degree_df[out_degree_df$Age == "EM", "out_degree"]
  out_degree_MM <- out_degree_df[out_degree_df$Age == "MM", "out_degree"]
  
  print("Difference in out-degree")
  
  print(wilcox.test(out_degree_EM, out_degree_UC, alternative="greater")$p.value)
  print(wilcox.test(out_degree_EM, out_degree_MM, alternative="greater")$p.value)
  
  print(aggregate(out_degree ~ Age, out_degree_df, mean))
  
  out_degree_df$Group <- "A"
  
  pdf(paste(database, "_Fig5.pdf", sep=""))
  g <- ggplot(out_degree_df, aes(y=log10(out_degree), x=Group))+
    geom_boxplot()+
    ylab("Out-degree (log10)")+
    xlab("")+
    coord_flip()+
    theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  print(g)
  dev.off()
  
  upper <- quantile(out_degree_df$out_degree, 0.75)
  
  n_out_10_UC <- sum(out_degree_df$out_degree[out_degree_df$Age == "UC"] >= upper, na.rm=TRUE)
  n_out_10_EM <- sum(out_degree_df$out_degree[out_degree_df$Age == "EM"] >= upper, na.rm=TRUE)
  n_out_10_MM <- sum(out_degree_df$out_degree[out_degree_df$Age == "MM"] >= upper, na.rm=TRUE)
  total_out_10 <- sum(out_degree_df$out_degree >= upper)
  
  print("Percentage upperquantile out degree UC over total over upperquantile")
  print(n_out_10_UC/total_out_10)
  print("Percentage upperquantile out degree EM over total over upperquantile")
  print(n_out_10_EM/total_out_10)
  print("Percentage upperquantile out degree MM over total over upperquantile")
  print(n_out_10_MM/total_out_10)
  
  ##Enrichment of EM genes are regulators
  genes_in_RN <- unique(c(as.character(network$Gene1), as.character(network$Gene2)))
  
  n_EM_genes <- length(EM_genes)
  total_genes <- nrow(genes_phy)
  
  total_regulators <- length(unique(network$Gene1))
  
  n_EM_regulators <- sum(unique(network$Gene1) %in% EM_genes)
  n_UC_regulators <- sum(unique(network$Gene1) %in% UC_genes)
  n_MM_regulators <- sum(unique(network$Gene1) %in% MM_genes)
  
  print("Percentage UC, EM and MM regulators")
  print((n_UC_regulators/total_regulators)*100)
  print((n_EM_regulators/total_regulators)*100)
  print((n_MM_regulators/total_regulators)*100)
  
  print("Enrichment EM genes as regulators")
  print(fisher.test(cbind(c(n_EM_regulators, total_regulators),
                          c(n_EM_genes, total_genes)), alternative="greater")$p.value)
  
  percentage_regulators <- data.frame(Age=c("UC", "EM", "MM"),
                                      Percentage=c((n_UC_regulators/total_regulators)*100,
                                                   (n_EM_regulators/total_regulators)*100,
                                                   (n_MM_regulators/total_regulators)*100))
  
  percentage_regulators$Age <- factor(percentage_regulators$Age, levels=c("UC", "EM", "MM"))
  
  pdf(paste(database, "_Fig6.pdf", sep=""),
      height=3, width=3)
  g <- ggplot(percentage_regulators, aes(x=Age, y=Percentage))+
    geom_bar(aes(fill=Age), stat='identity')+
    ylab("Percentage regulators")+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(g)
  dev.off()
  
  
  
  ##Read in somatic mutations
 
  mutations_df <- load_mutations()
  
  mutations_df <- mutations_df[mutations_df$Number_mutations >= 3,]
  mutations_df <- mutations_df[mutations_df$Syn_ratio > 1,]
  
  all_miss <- vector()
  all_lof <- vector()
  miss_df <- vector()
  lof_df <- vector()
  
  for(tumour in tumours){
    local_mut <- mutations_df[mutations_df$Tumour == tumour,]
    
    missense <- as.character(local_mut[local_mut$Variant_type == "Missense", "Hugo_Symbol"])
    lof <- as.character(local_mut[local_mut$Variant_type == "LoF", "Hugo_Symbol"])
    
    if(length(missense) != 0){
      miss_df <- rbind(miss_df, cbind(Tumour=tumour, Alt="Miss", Genes=missense))
      all_miss <- c(all_miss, missense)
    }
    if(length(lof) != 0){
      lof_df <- rbind(lof_df, cbind(Tumour=tumour, Alt="LoF", Genes=lof))
      all_lof <- c(all_lof, lof)
    }
  }
  
  ##Read in CNVs
  ##Using only focal ones
  
  CNV_df <- load_CNVs(only_focal="FOCAL")
  
  all_amp <- vector()
  all_del <- vector()
  amp_df <- vector()
  del_df <- vector()
  for(tumour in tumours){
    n_pat_df <- CNV_df[CNV_df$Tumour == tumour,]
    amp_genes <- as.character(n_pat_df[n_pat_df$Patients_amp > 0.10, "Genes"])
    del_genes <- as.character(n_pat_df[n_pat_df$Patients_del > 0.10, "Genes"])
    
    print(tumour)
    print(length(amp_genes))
    print(length(del_genes))
    
    all_amp <- c(all_amp, amp_genes)
    all_del <- c(all_del, del_genes)
    
    if(length(amp_genes) > 0){
      amp_df <- rbind(amp_df, cbind(tumour, "Amp", amp_genes))
      
    }
    if(length(del_genes) > 0){
      del_df <- rbind(del_df, cbind(tumour, "Del", del_genes))
    }
  }
  
  only_miss <- unique(all_miss)
  only_lof <- unique(all_lof)
  both_mutations <- only_miss[only_miss %in% only_lof]
  
  only_miss <- only_miss[!(only_miss %in% both_mutations)]
  only_lof <- only_lof[!(only_lof %in% both_mutations)]
  
  only_amp <- unique(all_amp)
  only_del <- unique(all_del)
  both_CNVs <- only_amp[only_amp %in% only_del]
  only_amp <- only_amp[!(only_amp %in% both_CNVs)]
  only_del <- only_del[!(only_del %in% both_CNVs)]
  
  
  #Common alterations
  all_amp_common <- names(table(all_amp))[table(all_amp)>=7]
  all_del_common <- names(table(all_del))[table(all_del)>=7]
  all_CNV_common_both <- all_amp_common[all_amp_common %in% all_del_common]
  all_amp_common <- all_amp_common[!(all_amp_common %in% all_CNV_common_both)]
  all_del_common <- all_del_common[!(all_del_common %in% all_CNV_common_both)]
  
  
  all_miss_common <- names(table(all_miss))[table(all_miss)>=7]
  all_lof_common <- names(table(all_lof))[table(all_lof)>=7]
  all_mut_common_both <- all_miss_common[all_miss_common %in% all_lof_common]
  all_miss_common <- all_miss_common[!(all_miss_common %in% all_mut_common_both)]
  all_lof_common <- all_lof_common[!(all_lof_common %in% all_mut_common_both)]
  
  genes_in_GRN <- unique(c(network$Gene1, network$Gene2))
  all_miss_common <- all_miss_common[all_miss_common %in% genes_in_GRN]
  all_lof_common <- all_lof_common[all_lof_common %in% genes_in_GRN]
  all_amp_common <- all_amp_common[all_amp_common %in% genes_in_GRN]
  all_del_common <- all_del_common[all_del_common %in% genes_in_GRN]
  
  
  ###Enrichment of mutations and CNVs as targets, regulators
  regulators <- unique(network$Gene1)
  targets <- unique(network$Gene2)
  non_target_regulators <- regulators[!(regulators %in% targets)]
  non_regulator_targets <- targets[!(targets %in% regulators)]
  
  
  total_regulators <- length(regulators)
  total_non_target_regulators <- length(non_target_regulators)
  
  alt_genes <- rbind(miss_df, lof_df, amp_df, del_df)
  alt_genes <- data.frame(alt_genes)
  alt_genes <- alt_genes[alt_genes$Genes %in% unique(c(network$Gene1, network$Gene2)),]
  
  
  ##Percentage of genes with alterations that are regulators, targets, etc
  alt_genes$Alt_simplified <- ifelse(alt_genes$Alt %in% c("Miss", "LoF"), "Point",
                                     ifelse(alt_genes$Alt %in% c("Amp", "Del"), "CNV", NA))
  
  fraction_gene_class_binary <- vector()
  for(tumour in tumours){
    local_alt <- alt_genes[alt_genes$Tumour == tumour,]
    for(alt in unique(alt_genes$Alt_simplified)){
      local_alt2 <- local_alt[local_alt$Alt_simplified == alt,]
      if(nrow(local_alt2) >= 3){
        per_regulators <- sum(unique(local_alt2$Genes) %in% regulators)/length(unique(local_alt2$Genes))
        per_non_r_targets <- sum(unique(local_alt2$Genes) %in% non_regulator_targets)/length(unique(local_alt2$Genes))
      
        fraction_gene_class_binary <- rbind(fraction_gene_class_binary,
                                          c(Tumour=tumour, Alteration=alt, 
                                            Per_regulators=per_regulators, 
                                            Per_non_r_targets=per_non_r_targets))
      }
    }
  }
  
  fraction_gene_class_binary <- as.data.frame(fraction_gene_class_binary)
  fraction_gene_class_binary$Per_regulators <- as.numeric(as.character(fraction_gene_class_binary$Per_regulators))
  fraction_gene_class_binary$Per_non_r_targets <- as.numeric(as.character(fraction_gene_class_binary$Per_non_r_targets))
  
  
  fraction_gene_class_binary_melt <- melt(fraction_gene_class_binary[,c("Tumour", "Alteration", "Per_regulators", "Per_non_r_targets")])
  
  colnames(fraction_gene_class_binary_melt) <- c("Tumour", "Alteration", "Gene_role", "Fraction")
  
  fraction_gene_class_binary_melt$Gene_role <- as.character(fraction_gene_class_binary_melt$Gene_role)
  fraction_gene_class_binary_melt$Gene_role[fraction_gene_class_binary_melt$Gene_role == "Per_non_r_targets"] <- "Target"
  fraction_gene_class_binary_melt$Gene_role[fraction_gene_class_binary_melt$Gene_role == "Per_regulators"] <- "Regulator"
  
  fraction_gene_class_binary_melt$Alteration <- factor(fraction_gene_class_binary_melt$Alteration, levels=c("Point", "CNV"))
  
  fraction_gene_class_binary_melt$Tumour <- factor(fraction_gene_class_binary_melt$Tumour,
                                                   levels=tumours)
  fraction_gene_class_binary_melt$Alteration <- factor(fraction_gene_class_binary_melt$Alteration,
                                                       levels=c("Point", "CNV"))
  
  fraction_gene_class_binary_mean <- aggregate(Fraction ~ Alteration+Gene_role, fraction_gene_class_binary_melt, mean)
  
  fraction_gene_class_binary_mean$Alteration <- factor(fraction_gene_class_binary_mean$Alteration,
                                                       levels=c("Point", "CNV"))
  
  pdf(paste(database, "_Fig7.pdf", sep=""),
      width=4, height=4)
  g <- ggplot(subset(fraction_gene_class_binary_melt, Gene_role=="Regulator"), aes(x=Alteration, y =Fraction))+
    geom_point(aes(colour=Tumour))+
    #coord_cartesian(ylim=c(0.05,0.4))+
    geom_bar(data=subset(fraction_gene_class_binary_mean, Gene_role=="Regulator"), aes(x=Alteration, y=Fraction), stat='identity', fill="grey92", color="black")+
    geom_path(aes(group=Alteration))+
    geom_point(aes(colour=Tumour), size=2)+
    ggtitle("Regulator")+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank())
  print(g)
  dev.off()
  
  if(database == "RegNetwork"){
    pdf(paste(database, "_Fig8.pdf", sep=""),
        width=4, height=4)
    g <- ggplot(subset(fraction_gene_class_binary_melt, Gene_role=="Target"), aes(x=Alteration, y =Fraction))+
      geom_point(aes(colour=Tumour))+
      coord_cartesian(ylim=c(0.6,1))+
      geom_bar(data=subset(fraction_gene_class_binary_mean, Gene_role=="Target"), 
               aes(x=Alteration, y=Fraction), stat='identity', fill="grey92", color="black")+
      geom_path(aes(group=Alteration))+
      geom_point(aes(colour=Tumour), size=2)+
      #ylim(0.5, 1)+
      scale_y_continuous(position = "right")+
      ggtitle("Target")+
      theme_bw()+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank())
    print(g)
    dev.off()
    
  }else if (database == "Trrust"){
    pdf(paste(database, "_Fig8.pdf", sep=""),
        width=4, height=4)
    g <- ggplot(subset(fraction_gene_class_binary_melt, Gene_role=="Target"), aes(x=Alteration, y =Fraction))+
      geom_point(aes(colour=Tumour))+
      coord_cartesian(ylim=c(0.25,0.95))+
      geom_bar(data=subset(fraction_gene_class_binary_mean, Gene_role=="Target"), 
               aes(x=Alteration, y=Fraction), stat='identity', fill="grey92", color="black")+
      geom_path(aes(group=Alteration))+
      geom_point(aes(colour=Tumour), size=2)+
      #ylim(0.5, 1)+
      scale_y_continuous(position = "right")+
      ggtitle("Target")+
      theme_bw()+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank())
    print(g)
    dev.off()
    
  }else{
    pdf(paste(database, "_Fig8.pdf", sep=""),
        width=4, height=4)
    g <- ggplot(subset(fraction_gene_class_binary_melt, Gene_role=="Target"), aes(x=Alteration, y =Fraction))+
      geom_point(aes(colour=Tumour))+
      #coord_cartesian(ylim=c(0.6,0.95))+
      geom_bar(data=subset(fraction_gene_class_binary_mean, Gene_role=="Target"), 
               aes(x=Alteration, y=Fraction), stat='identity', fill="grey92", color="black")+
      geom_path(aes(group=Alteration))+
      geom_point(aes(colour=Tumour), size=2)+
      #ylim(0.5, 1)+
      scale_y_continuous(position = "right")+
      ggtitle("Target")+
      theme_bw()+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank())
    print(g)
    dev.off()
    
  }

  fraction_gene_class_binary_melt_reg <- fraction_gene_class_binary_melt[fraction_gene_class_binary_melt$Gene_role == "Regulator",]
  fraction_gene_class_binary_melt_tar <- fraction_gene_class_binary_melt[fraction_gene_class_binary_melt$Gene_role == "Target",]
  
  print("Wilcoxon regulators are point mutated")  
  print(wilcox.test(fraction_gene_class_binary_melt_reg$Fraction[fraction_gene_class_binary_melt_reg$Alteration == "Point"],
              fraction_gene_class_binary_melt_reg$Per_regulators[fraction_gene_class_binary_melt_reg$Alteration == "CNV"], alternative="greater")$p.value)
  
  print("Wilcoxon targets are CNAed")  
  print(wilcox.test(fraction_gene_class_binary_melt_tar$Fraction[fraction_gene_class_binary_melt_tar$Alteration == "CNV"],
              fraction_gene_class_binary_melt_tar$Per_regulators[fraction_gene_class_binary_melt_tar$Alteration == "Point"], alternative="greater")$p.value)
  

  
  ###Calculation of in and out degree ratio per tumour type
  out_degree_df$Database <- "GRN"
  colnames(out_degree_df)[2] <- "Degree"
  in_degree_df$Database <- "GRN"
  colnames(in_degree_df)[2] <- "Degree"
  
  #First calculating the ratio for each gene, then averaging
  median_ratio <- vector()
  p_EM <- vector()
  median_overall <- vector()
  for(tumour in tumours){
    local_miss <- miss_df[miss_df[,1] == tumour,]
    local_lof <- lof_df[lof_df[,1] == tumour,]
    local_amp <- amp_df[amp_df[,1] == tumour,]
    local_del <- del_df[del_df[,1] == tumour,]
    
    local_out_degree_df <- add_mut_CNV_to_degree_df(out_degree_df, local_miss, local_lof, local_amp, local_del)
    local_in_degree_df <- add_mut_CNV_to_degree_df(in_degree_df, local_miss, local_lof, local_amp, local_del)
    
    local_out_degree_df$Alt_simplified <- ifelse(local_out_degree_df$Alt %in% c("Miss", "LoF"), "Point",
                                                 ifelse(local_out_degree_df$Alt %in% c("Amp", "Del"), "CNV", "WT"))
    local_in_degree_df$Alt_simplified <- ifelse(local_in_degree_df$Alt %in% c("Miss", "LoF"), "Point",
                                                ifelse(local_in_degree_df$Alt %in% c("Amp", "Del"), "CNV", "WT"))
    colnames(local_out_degree_df)[2] <- "Out_degree"
    local_degree_df <- local_out_degree_df
    
    #Only considering genes with dual roles
    local_degree_df <- data.frame(local_degree_df, In_degree=local_in_degree_df[match(local_degree_df$Genes, local_in_degree_df$Genes), "Degree"])
    
    local_degree_df$Ratio_out_in <- local_degree_df$Out_degree/local_degree_df$In_degree
    
    #local_median_overall <- median(local_degree_df$Ratio_out_in)
    #median_overall <- c(median_overall, local_median_overall)
    #local_degree_df$Ratio_out_in <- local_degree_df$Ratio_out_in/local_median_overall
    
    local_median_ratio <- aggregate(Ratio_out_in ~ Age+Alt_simplified, local_degree_df, median)
    local_median_ratio$Tumour <- tumour
    
    median_ratio <- rbind(median_ratio, local_median_ratio)
  }
  
  
  upper_quantile <- aggregate(Ratio_out_in ~ Age+Alt_simplified, median_ratio, quantile, 0.75)
  lower_quantile <- aggregate(Ratio_out_in ~ Age+Alt_simplified, median_ratio, quantile, 0.25)
  median_values <- aggregate(Ratio_out_in ~ Age+Alt_simplified, median_ratio, median)
  
  
  
  colnames(median_values)[3] <- "Median_values"
  median_values$Upper <- upper_quantile$Ratio_out_in
  median_values$Lower <- lower_quantile$Ratio_out_in
  median_values$Full_name <- paste(median_values$Age, median_values$Alt_simplified, sep="_")
  
  median_values$Age <- factor(median_values$Age, levels=c("MM", "EM", "UC"))
  
  
  pdf(paste(database, "_Fig9.pdf", sep=""),
      width=7, height=3)
  g <- ggplot(median_values, aes(x=Alt_simplified, y =log2(Median_values)))+
    geom_point(aes(colour=Age, shape=Alt_simplified), size=4,
               position = position_dodge(width = 0.75))+
    geom_linerange(aes(ymin=log2(Lower), ymax=log2(Upper), colour=Age),
                   position = position_dodge(width = 0.75), size=0.75)+
    #ylim(0.25,1.75)+
    #coord_flip(ylim=c(0.35,1.75))+
    coord_flip()+
    #geom_hline(yintercept = 1, size=0.5, colour="grey")+
    geom_point(aes(colour=Age, shape=Alt_simplified), size=4,
               position = position_dodge(width = 0.75))+
    geom_linerange(aes(ymin=log2(Lower), ymax=log2(Upper), colour=Age),
                   position = position_dodge(width = 0.75), size=0.75)+
    scale_colour_manual(values=c(UC=gg_color_hue(3)[1], EM=gg_color_hue(3)[2], MM=gg_color_hue(3)[3]))+
    ylab("Ratio out-degree/in-degree")+
    xlab("")+
    facet_wrap("Alt_simplified", scales='free_y', strip.position="left", nrow=2)+
    theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_rect(colour=NA, fill=NA),
          panel.background = element_rect(fill = NA))
  print(g)
  dev.off()
  
}

load_TRRUST <- function(){
  network <- read.delim("trrust_rawdata.human.tsv", sep="\t",
                        header=FALSE)
  
  network <- network[,1:2]
  colnames(network) <- c("Gene1", "Gene2")
  network <- unique(network)
  
  network <- prepare_network(network)
  return(network)  
}

load_RegNetwork <- function(){
  network <- read.delim("human.source", sep="\t",
                        header=FALSE)
  network <- network[,c(1,3)]
  colnames(network) <- c("Gene1", "Gene2")
  network <- unique(network)
  
  network <- prepare_network(network)
  return(network)
}

load_KEGG <- function(){
  network <- read.delim("kegg.human.reg.direction", sep="\t")
  network <- network[,c(1,3)]
  colnames(network) <- c("Gene1", "Gene2")
  network <- unique(network)
  
  network <- prepare_network(network)
  return(network)
}


load_PathwayCommons_GRN <- function(){
  network <- read.delim("PathwayCommons9.All.hgnc.sif", header=FALSE, dec=",")
  network <- network[!(network[,2] %in% c("controls-production-of", "controls-transport-of-chemical",
                                          "chemical-affects", "consumption-controlled-by", "reacts-with",
                                          "used-to-produce")),]  ##To remove all reactions involving chemicals
  network <- network[,c(1,3,2)]
  colnames(network) <- c("Gene1", "Gene2", "Edge_type")
  network <- network[network$Edge_type == "controls-expression-of",]
  network <- network[,1:2]
  
  network <- prepare_network(network)
  return(network)
}


reproduce_results_figure_3C <- function(database){
  
  mutations_df <- load_mutations()
  
  mutations_df <- mutations_df[mutations_df$Number_mutations >= 3,]
  mutations_df <- mutations_df[mutations_df$Syn_ratio > 1,]
  
  all_miss <- vector()
  all_lof <- vector()
  miss_df <- vector()
  lof_df <- vector()
  
  for(tumour in tumours){
    local_mut <- mutations_df[mutations_df$Tumour == tumour,]
    
    missense <- as.character(local_mut[local_mut$Variant_type == "Missense", "Hugo_Symbol"])
    lof <- as.character(local_mut[local_mut$Variant_type == "LoF", "Hugo_Symbol"])

    if(length(missense) != 0){
      miss_df <- rbind(miss_df, cbind(Tumour=tumour, Alt="Miss", Genes=missense))
      all_miss <- c(all_miss, missense)
    }
    if(length(lof) != 0){
      lof_df <- rbind(lof_df, cbind(Tumour=tumour, Alt="LoF", Genes=lof))
      all_lof <- c(all_lof, lof)
    }
    
  }
  
  ##Read in CNVs
  ##Using only focal ones
  CNV_df <- load_CNVs(only_focal="FOCAL")
  
  all_amp <- vector()
  all_del <- vector()
  amp_df <- vector()
  del_df <- vector()
  for(tumour in tumours){
    n_pat_df <- CNV_df[CNV_df$Tumour == tumour,]
    amp_genes <- as.character(n_pat_df[n_pat_df$Patients_amp > 0.10, "Genes"])
    del_genes <- as.character(n_pat_df[n_pat_df$Patients_del > 0.10, "Genes"])
    
    print(tumour)
    print(length(amp_genes))
    print(length(del_genes))
    
    all_amp <- c(all_amp, amp_genes)
    all_del <- c(all_del, del_genes)
    
    if(length(amp_genes) > 0){
      amp_df <- rbind(amp_df, cbind(tumour, "Amp", amp_genes))
      
    }
    if(length(del_genes) > 0){
      del_df <- rbind(del_df, cbind(tumour, "Del", del_genes))
    }
  }
  
  only_miss <- unique(all_miss)
  only_lof <- unique(all_lof)
  both_mutations <- only_miss[only_miss %in% only_lof]
  
  only_miss <- only_miss[!(only_miss %in% both_mutations)]
  only_lof <- only_lof[!(only_lof %in% both_mutations)]
  
  
  only_amp <- unique(all_amp)
  only_del <- unique(all_del)
  both_CNVs <- only_amp[only_amp %in% only_del]
  only_amp <- only_amp[!(only_amp %in% both_CNVs)]
  only_del <- only_del[!(only_del %in% both_CNVs)]
  
  #Common alterations
  all_amp_common <- names(table(all_amp))[table(all_amp)>=7]
  all_del_common <- names(table(all_del))[table(all_del)>=7]
  all_CNV_common_both <- all_amp_common[all_amp_common %in% all_del_common]
  all_amp_common <- all_amp_common[!(all_amp_common %in% all_CNV_common_both)]
  all_del_common <- all_del_common[!(all_del_common %in% all_CNV_common_both)]
  
  
  all_miss_common <- names(table(all_miss))[table(all_miss)>=7]
  all_lof_common <- names(table(all_lof))[table(all_lof)>=7]
  all_mut_common_both <- all_miss_common[all_miss_common %in% all_lof_common]
  all_miss_common <- all_miss_common[!(all_miss_common %in% all_mut_common_both)]
  all_lof_common <- all_lof_common[!(all_lof_common %in% all_mut_common_both)]
  
  load(paste("regulator_classification_", database, ".Rdata", sep=""))
  
  ##Enrichment of ages in each class
  per_each_age <- vector()
  p_enrich_regulator_ages <- vector()
  for(local_class in unique(regulator_classification$Regulator_class)){
    reg_of_class <- regulator_classification[regulator_classification$Regulator_class == local_class,"Regulator"]
    
    UC_reg_local <- sum(reg_of_class %in% UC_genes)
    EM_reg_local <- sum(reg_of_class %in% EM_genes)
    MM_reg_local <- sum(reg_of_class %in% MM_genes)
    total_reg_local <- UC_reg_local+EM_reg_local+MM_reg_local
    
    per <- c((UC_reg_local/total_reg_local)*100,
             (EM_reg_local/total_reg_local)*100,
             (MM_reg_local/total_reg_local)*100)
    
    per_each_age <- rbind(per_each_age,
                          c(local_class, per))
    
    total_UC <- sum(regulator_classification$Regulator %in% UC_genes)
    total_EM <- sum(regulator_classification$Regulator %in% EM_genes)
    total_MM <- sum(regulator_classification$Regulator %in% MM_genes)
    total_genes_ages <- length(regulator_classification$Regulator)
    
    total_UC <- length(UC_genes)
    total_EM <- length(EM_genes)
    total_MM <- length(MM_genes)
    total_genes_ages <- sum(total_UC, total_EM, total_MM)
    
    p1 <- fisher.test(cbind(c(UC_reg_local, total_reg_local),
                            c(total_UC, total_genes_ages)), alternative="greater")$p.val
    
    p2 <- fisher.test(cbind(c(EM_reg_local, total_reg_local),
                            c(total_EM, total_genes_ages)), alternative="greater")$p.val
    
    p3 <- fisher.test(cbind(c(MM_reg_local, total_reg_local),
                            c(total_MM, total_genes_ages)), alternative="greater")$p.val
    
    p_enrich_regulator_ages <- rbind(p_enrich_regulator_ages, c(local_class, p1, p2, p3))
    
  }
  
  colnames(p_enrich_regulator_ages) <- c("Regulator_class", "Enriched_UC", "Enriched_EM", "Enriched_MM")
  p_enrich_regulator_ages <- as.data.frame(p_enrich_regulator_ages)
  p_enrich_regulator_ages$Enriched_UC <- as.numeric(as.character(p_enrich_regulator_ages$Enriched_UC))
  p_enrich_regulator_ages$Enriched_EM <- as.numeric(as.character(p_enrich_regulator_ages$Enriched_EM))
  p_enrich_regulator_ages$Enriched_MM <- as.numeric(as.character(p_enrich_regulator_ages$Enriched_MM))
  
  p_enrich_regulator_ages <- melt(p_enrich_regulator_ages)
  
  colnames(per_each_age) <- c("Regulator_class", "Per_UC", "Per_EM", "Per_MM")
  
  
  sink(paste(database, "_Summary_of_Fig_3C.txt", sep=""))
  print("Enrichment of ages in regulators")
  print(p_enrich_regulator_ages)
  
  load(paste("percentage_UC_per_gene_directed_", database, ".Rdata", sep=""))

  percentage_UC_per_gene_directed$Degree <- NA
   
  percentage_UC_EM_MM_alt <- add_mut_CNV_to_degree_df(percentage_UC_per_gene_directed[,c("Genes","Degree","Age","Database")], all_miss_common, all_lof_common,
                                                      all_amp_common, all_del_common)
  
  percentage_UC_EM_MM_alt$Alt_binary <- ifelse(percentage_UC_EM_MM_alt$Alt %in% c("Miss", "LoF"), "Point",
                                               ifelse(percentage_UC_EM_MM_alt$Alt %in% c("Amp", "Del"), "CNV", "WT"))
  
  
  percentage_UC_EM_MM_alt$Regulator_class <- regulator_classification[match(percentage_UC_EM_MM_alt$Genes, regulator_classification$Regulator), "Regulator_class"]
  
  quantiles1 <- count_regulators_each_class(percentage_UC_EM_MM_alt)
  
  quantiles1$Alt_binary <- factor(quantiles1$Alt_binary, levels=c("WT", "CNV", "Point"))
  
  quantiles1$Regulator_class <- factor(quantiles1$Regulator_class,
                                       levels=c("EM_downstream", "Mixed_downstream", "UC_downstream"))
  
  pdf(paste(database, "_Fig3C.pdf", sep=""),
      width=5, height=2.4)
  g <- ggplot(quantiles1, aes(x=Alt_binary, y=Per_regulators_by_type))+
    geom_bar(aes(fill=Regulator_class), stat='identity')+
    scale_fill_manual(values=c("UC_downstream" = "#0073C2FF",
                               "Mixed_downstream" = "#EFC000FF",
                               "EM_downstream" = "#868686FF"))+
    #scale_fill_jco()+
    ylab("Fraction of regulators")+
    coord_flip()+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank())
  
  print(g)
  dev.off()
  
  
  percentage_UC_EM_MM_alt <- percentage_UC_EM_MM_alt[percentage_UC_EM_MM_alt$Regulator_class %in% c("UC_downstream", "Mixed_downstream", "EM_downstream"),]
  temp <- vector()
  for(mut in c("Point", "CNV")){
    for(class in c("UC_downstream", "Mixed_downstream", "EM_downstream")){
      mut_class <- length(intersect(which(percentage_UC_EM_MM_alt$Alt_binary == mut),
                                    which(percentage_UC_EM_MM_alt$Regulator_class == class)))
      all_mut <- sum(percentage_UC_EM_MM_alt$Alt_binary == mut)
      temp <- rbind(temp, c(mut, class, (mut_class/all_mut)*100))
    }
  }
  colnames(temp) <- c("Mutation", "Regulator_class", "Percentage")
  print(temp)
  sink()
}


reproduce_results_figure_4A <- function(RN_network, database){

  ##Determine the downstream targets of each regulator
  regulators <- unique(RN_network$Gene1)
  
  regulator_targets <- lapply(regulators, function(reg){
    local_edges <- RN_network[RN_network$Gene1 == reg,]
    return(unique(local_edges$Gene2))
  })
  
  names(regulator_targets) <- regulators
  
  patient_CNVs <- vector()
  for(tumour in tumours){
    tumour_amp <- CNVs_curated[[tumour]]$amplifications
    tumour_del <- CNVs_curated[[tumour]]$deletions
    
    tumour_amp$Sample <- substr(tumour_amp$Sample, 9, 12)
    tumour_del$Sample <- substr(tumour_del$Sample, 9, 12)
    
    tumour_amp$CNV <- "Amp"
    tumour_del$CNV <- "Del"
    tumour_amp$Tumour <- tumour
    tumour_del$Tumour <- tumour
    tumour_amp$Segment_Mean <- NULL
    tumour_del$Segment_Mean <- NULL
    
    patient_CNVs <- rbind(patient_CNVs, tumour_amp, tumour_del)
  }
  
  
  me_p <- vector()
  for(regulator in regulators){
    regulator_targets_CNVs <- vector()

    local_targets <- regulator_targets[[regulator]]

    targets_with_CNVs <- patient_CNVs[patient_CNVs$Gene %in% local_targets,]
    targets_with_CNVs$Target_age <- genes_phy[match(targets_with_CNVs$Gene, genes_phy$GeneID), "Age"]

    CNV_in_regulator <- patient_CNVs[patient_CNVs$Gene %in% regulator,]

    if(nrow(CNV_in_regulator) >= 3){
      all_patients_involved <- unique(c(targets_with_CNVs$Sample, CNV_in_regulator$Sample))

      for(patient in all_patients_involved){

        patient_targets_with_CNVs <- targets_with_CNVs[targets_with_CNVs$Sample == patient,]

        patient_targets_with_CNVs_UC <- patient_targets_with_CNVs[patient_targets_with_CNVs$Target_age == "UC",]
        patient_targets_with_CNVs_EM <- patient_targets_with_CNVs[patient_targets_with_CNVs$Target_age == "EM",]
        patient_targets_with_CNVs_MM <- patient_targets_with_CNVs[patient_targets_with_CNVs$Target_age == "MM",]

        patient_targets_with_CNVs_UC <- patient_targets_with_CNVs_UC[!is.na(patient_targets_with_CNVs_UC$Gene),]
        patient_targets_with_CNVs_EM <- patient_targets_with_CNVs_EM[!is.na(patient_targets_with_CNVs_EM$Gene),]
        patient_targets_with_CNVs_MM <- patient_targets_with_CNVs_MM[!is.na(patient_targets_with_CNVs_MM$Gene),]

        patient_regulator_with_CNV <- patient %in% CNV_in_regulator$Sample

        regulator_targets_CNVs <- rbind(regulator_targets_CNVs,
                                        c(patient, patient_regulator_with_CNV, nrow(patient_targets_with_CNVs)/length(local_targets),
                                          nrow(patient_targets_with_CNVs_UC), nrow(patient_targets_with_CNVs_EM), nrow(patient_targets_with_CNVs_MM)))
      }

      fraction_targets_when_regulator_CNV <- as.numeric(regulator_targets_CNVs[regulator_targets_CNVs[,2] == TRUE, 3])
      fraction_targets_when_regulator_not_CNV <- as.numeric(regulator_targets_CNVs[regulator_targets_CNVs[,2] == FALSE, 3])
      n_UC_targets_affected <- as.numeric(regulator_targets_CNVs[, 4])
      n_EM_targets_affected <- as.numeric(regulator_targets_CNVs[, 5])
      n_MM_targets_affected <- as.numeric(regulator_targets_CNVs[, 6])

      if((length(fraction_targets_when_regulator_not_CNV) >= 3) && (length(fraction_targets_when_regulator_CNV) >= 3)){
        p1 <- wilcox.test(fraction_targets_when_regulator_not_CNV, fraction_targets_when_regulator_CNV, alternative="greater")$p.val
        p2 <- wilcox.test(fraction_targets_when_regulator_CNV, fraction_targets_when_regulator_not_CNV, alternative="greater")$p.val
        me_p <- rbind(me_p, c(regulator, median(fraction_targets_when_regulator_CNV), median(fraction_targets_when_regulator_not_CNV), p1, p2,
                              median(n_UC_targets_affected), median(n_EM_targets_affected), median(n_MM_targets_affected)))
      }
    }
    print(regulator)
  }
   
  save(me_p, file=paste("me_p_", database, ".Rdata", sep=""))
  load(paste("me_p_", database, ".Rdata", sep=""))
  
  me_p <- as.data.frame(me_p)
  colnames(me_p) <- c("Regulator", "Fraction_targets_regulator_CNV",
                      "Fraction_targets_regulator_not_CNV", "Preference_CNV_targets", "Preference_CNV_regulator", 
                      "Median_target_UC_genes_CNVs", "Median_target_EM_genes_CNVs", "Median_target_MM_genes_CNVs")
  me_p$Fraction_targets_regulator_CNV <- as.numeric(as.character(me_p$Fraction_targets_regulator_CNV))
  me_p$Fraction_targets_regulator_not_CNV <- as.numeric(as.character(me_p$Fraction_targets_regulator_not_CNV))
  me_p$Preference_CNV_targets <- as.numeric(as.character(me_p$Preference_CNV_targets))
  me_p$Preference_CNV_regulator <- as.numeric(as.character(me_p$Preference_CNV_regulator))
  
  me_p$Preference_CNV_targets_adj <- p.adjust(me_p$Preference_CNV_targets, method="BH")
  me_p$Preference_CNV_regulator_adj <- p.adjust(me_p$Preference_CNV_regulator, method="BH")
  
  n_targets_per_regulator <- sapply(regulator_targets, length)
  n_targets_per_regulator <- data.frame(Regulator=names(n_targets_per_regulator),
                                        N_targets = n_targets_per_regulator)
  
  me_p$N_targets <- n_targets_per_regulator[match(me_p$Regulator, n_targets_per_regulator$Regulator), "N_targets"]
  
  me_p <- me_p[me_p$N_targets >= 2,]
  me_p <- me_p[!is.na(me_p$Regulator),]  
  me_p$Age <- genes_phy[match(me_p$Regulator, genes_phy$GeneID), "Age"]
  
  me_p$Age <- factor(me_p$Age, levels=c("UC", "EM", "MM"))
  
  me_p_fractions <- me_p[,c("Regulator", "Fraction_targets_regulator_CNV", "Fraction_targets_regulator_not_CNV")]
  me_p_fractions <- melt(me_p_fractions)
  
  colnames(me_p_fractions) <- c("Regulator", "Regulator_status", "Fraction_targets_with_CNAs")
  
  me_p_fractions$Regulator_status <- ifelse(me_p_fractions$Regulator_status == "Fraction_targets_regulator_CNV", "CNA", "CNN")
  
  ###Supplementary
  pdf(paste(database, "_Fig4A_supp.pdf", sep=""),
      width=4, height=3)
  g <- ggplot(me_p_fractions, aes(x=Regulator_status, y=Fraction_targets_with_CNAs))+
    geom_boxplot(aes(fill=Regulator_status))+
    ylab("Median fraction of targets with CNAs")+
    xlab("Regulator status")+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank())
  print(g)
  dev.off()
  
  sink(paste(database, "_Summary_of_Fig_4A.txt", sep=""))
  print("Fraction targets CNAs greater when reg CNN")
  print(wilcox.test(me_p_fractions$Fraction_targets_with_CNAs[me_p_fractions$Regulator_status == "CNN"],
              me_p_fractions$Fraction_targets_with_CNAs[me_p_fractions$Regulator_status == "CNA"], alternative="greater")$p.value)
  
  
  
  load(paste("regulator_classification_", database, ".Rdata", sep=""))
  
  me_p_fractions$Quantiles1 <- regulator_classification[match(me_p_fractions$Regulator, regulator_classification$Regulator), "Regulator_class"]
  
  me_p_fractions$Quantiles1 <- factor(me_p_fractions$Quantiles1, levels=c("UC_downstream", "Mixed_downstream", "EM_downstream"))
  
  
  me_p_fractions_no_CNVs_in_regulators <- me_p_fractions[!is.na(me_p_fractions$Quantiles1),]
  me_p_fractions_no_CNVs_in_regulators <- me_p_fractions_no_CNVs_in_regulators[me_p_fractions_no_CNVs_in_regulators$Regulator_status == "CNN",]
  
  pdf(paste(database, "_Fig4A.pdf", sep=""),
      width=5.25, height=4)
  g <- ggplot(me_p_fractions_no_CNVs_in_regulators, aes(x=Quantiles1, y=Fraction_targets_with_CNAs))+
    geom_boxplot(aes(fill=Quantiles1))+
    scale_fill_jco()+
    xlab("Regulator class")+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank())
  print(g)
  dev.off()
  
  aggregate(Fraction_targets_with_CNAs ~ Quantiles1, me_p_fractions_no_CNVs_in_regulators, median)

  print("Mixed greater fraction targets CNA when reg CNN than UC-t")  
  print(wilcox.test(me_p_fractions_no_CNVs_in_regulators[me_p_fractions_no_CNVs_in_regulators$Quantiles1 == "Mixed_downstream", "Fraction_targets_with_CNAs"],
              me_p_fractions_no_CNVs_in_regulators[me_p_fractions_no_CNVs_in_regulators$Quantiles1 == "UC_downstream", "Fraction_targets_with_CNAs"])$p.value)
  print("Mixed greater fraction targets CNA when reg CNN than EM-t")
  print(wilcox.test(me_p_fractions_no_CNVs_in_regulators[me_p_fractions_no_CNVs_in_regulators$Quantiles1 == "Mixed_downstream", "Fraction_targets_with_CNAs"],
              me_p_fractions_no_CNVs_in_regulators[me_p_fractions_no_CNVs_in_regulators$Quantiles1 == "EM_downstream", "Fraction_targets_with_CNAs"])$p.value)
  
}


prepare_CNV_genes <- function(CNV_genes){
  CNV_genes <- data.frame(CNV_genes[,c("Gene", "Sample")],
                          g[match(CNV_genes$Gene, g$SYMBOL), c("seqnames", "start", "end")])
  colnames(CNV_genes) <- c("Gene", "Sample", "Chromosome", "Start", "End")
  CNV_genes <- melt(CNV_genes)
  colnames(CNV_genes) <- c("Gene", "Sample", "Chromosome", "Location", "Coordinate")
  
  CNV_genes$Age <- genes_phy[match(CNV_genes$Gene, genes_phy$GeneID), "Phy_sum"]
  return(CNV_genes)
}

calculate_recurrent_genes <- function(CNV_genes, total_patients){
  recurrent_counts <- as.data.frame(table(CNV_genes$Gene))
  colnames(recurrent_counts) <- c("Gene", "Number_patients")
  recurrent_counts$Per_patients <- (recurrent_counts$Number_patients/total_patients)*100
  recurrent_counts$Age <- genes_phy[match(recurrent_counts$Gene, genes_phy$GeneID), "Phy_sum"]
  
  recurrent_threshold <- quantile(recurrent_counts$Per_patients, 0.9)  #top 10%
  recurrent_counts <- recurrent_counts[recurrent_counts$Per_patients >= recurrent_threshold,]
  recurrent_counts <- data.frame(recurrent_counts,
                                 g[match(recurrent_counts$Gene, g$SYMBOL), c("seqnames", "start", "end")])
  return(recurrent_counts)
}

prepare_CNV_df <- function(TCGA_files_CNV, CNV_genes, recurrent_counts, CNV_type, chr){
  local_TCGA_files_CNV <- TCGA_files_CNV[TCGA_files_CNV$Chromosome == chr,c("Sample", "Chromosome", "Start", "End")]
  
  local_TCGA_files_CNV <- melt(local_TCGA_files_CNV)
  colnames(local_TCGA_files_CNV) <- c("Sample", "Chromosome", "Location", "Coordinate")
  
  local_CNV_genes <- CNV_genes[CNV_genes$Chromosome == chr,c("Gene", "Chromosome", "Location", "Coordinate", "Age")]
  local_CNV_genes <- unique(local_CNV_genes)
  colnames(local_CNV_genes)[1] <- "Sample"
  if(nrow(local_CNV_genes) != 0){
    local_CNV_genes$Sample <- "Genes"
  }
  
  local_TCGA_files_CNV$Age <- "Irrelevant"
  
  local_TCGA_files_CNV <- rbind(local_TCGA_files_CNV, local_CNV_genes)
  
  local_recurrent_counts <- recurrent_counts[recurrent_counts$seqnames == chr,]
  local_recurrent_counts$Number_patients <- NULL
  local_recurrent_counts$Per_patients <- NULL
  local_recurrent_counts <- melt(local_recurrent_counts)
  colnames(local_recurrent_counts) <- c("Sample", "Age", "Chromosome", "Location", "Coordinate")
  local_recurrent_counts <- local_recurrent_counts[,c("Sample", "Chromosome", "Location", "Coordinate", "Age")]
  if(nrow(local_recurrent_counts) > 0){
    local_recurrent_counts$Sample <- paste("Recurrent", CNV_type, sep="_")
    local_TCGA_files_CNV <- rbind(local_TCGA_files_CNV, local_recurrent_counts)
  }
  return(local_TCGA_files_CNV)
}


plot_CNV <- function(TCGA_files_amp, TCGA_files_del, CNV_genes, tumour){
  CNV_genes_amp <- prepare_CNV_genes(CNV_genes$amplifications)
  CNV_genes_del <- prepare_CNV_genes(CNV_genes$deletions)
  
  total_patients <- length(unique(c(CNV_genes_amp$Sample, CNV_genes_del$Sample)))
  
  #Selecting recurrent genes
  recurrent_counts_amp <- calculate_recurrent_genes(CNV_genes_amp, total_patients)
  recurrent_counts_del <- calculate_recurrent_genes(CNV_genes_del, total_patients)
  
  
  
  chromosomes <- paste("chr", 1:23, sep="")
  
  for(chr in chromosomes){
    
    local_TCGA_files_amp <- prepare_CNV_df(TCGA_files_amp, CNV_genes_amp, recurrent_counts_amp, "Amp", chr)
    local_TCGA_files_del <- prepare_CNV_df(TCGA_files_del, CNV_genes_del, recurrent_counts_del, "Del", chr)
    
    local_TCGA_files_both <- data.frame(rbind(cbind(local_TCGA_files_amp, CNV="Amplification"),
                                              cbind(local_TCGA_files_del, CNV="Deletion")))
    
    levels <- unique(local_TCGA_files_both$Sample)
    levels <- levels[!(levels %in% c("Genes", "Recurrent_Amp", "Recurrent_Del"))]
    levels <- c("Genes", "Recurrent_Amp", "Recurrent_Del", as.character(levels))
    
    local_TCGA_files_both$Sample <- factor(local_TCGA_files_both$Sample, levels=levels )
    
    temp1 <- subset(local_TCGA_files_both, Sample %in% c("Genes", "Recurrent_Amp", "Recurrent_Del"))
    temp2 <- subset(local_TCGA_files_both, !(Sample %in% c("Genes", "Recurrent_Amp", "Recurrent_Del")))
    
    temp1$Age[temp1$Sample == "Recurrent_Amp"] <- "black"
    temp1$Age[temp1$Sample == "Recurrent_Del"] <- "black"
    
    
    g <- ggplot(local_TCGA_files_both, aes(x=Coordinate, y=Sample))+
      geom_line(data=temp1,aes(group=Sample, colour=Age), size=2)+
      geom_line(data=temp2, aes(group=Sample, colour=CNV), size=0.75)+
      scale_colour_manual(values=c(UC = "red", EM="green", MdM="yellow", MM = "blue", Irrelevant="black",
                                   Amplification="red", Deletion="blue", black="black"))+
      ggtitle(paste(tumour, chr))
    print(g)
  }
}

prepare_CNV_for_smooth <- function(local_CNV_genes){
  CNV_counts <- unlist(table(local_CNV_genes$Gene))
  CNV_counts_df <- data.frame(Genes = names(CNV_counts),
                              Counts = as.vector(CNV_counts))
  CNV_counts_df <- data.frame(CNV_counts_df,
                              g[match(CNV_counts_df$Genes, g$SYMBOL), c("seqnames", "start")])
  CNV_counts_df$Counts <- as.character(CNV_counts_df$Counts)
  CNV_counts_df_melt <- melt(CNV_counts_df)
  colnames(CNV_counts_df_melt) <- c("Genes", "N_patients", "Chromosome", "Location", "Coordinate")
  CNV_counts_df_melt$N_patients <- as.numeric(CNV_counts_df_melt$N_patients)
  return(CNV_counts_df_melt)
}


calculate_diff_exp_CNV <- function(CNV_genes, local_CNV, genes_class){
  CNV_pvalues <- vector()
  for(local_gene in CNV_genes){
    CNV_samples <- local_CNV[local_CNV$Gene == local_gene,"Sample"]  
    
    local_gene_exp <- unlist(patient_ratios[[tumour]][rownames(patient_ratios[[tumour]]) == local_gene,])
    
    local_gene_exp_df <- data.frame(Sample=names(local_gene_exp), Exp = local_gene_exp)
    
    local_gene_exp_df$Sample_type <- "NOT"
    local_gene_exp_df$Sample_type[local_gene_exp_df$Sample %in% CNV_samples] <- genes_class
 
     if(genes_class == "Amp"){
      direction <- "greater"  
    }else{
      direction <- "less"
    }
    
    if(sum(local_gene_exp_df$Sample %in% CNV_samples) > 1 & sum(is.finite(local_gene_exp_df$Exp)) != 0){
      
      if(sum(is.finite(local_gene_exp_df[local_gene_exp_df$Sample_type == genes_class, "Exp"])) != 0){
        if(sum(is.finite(local_gene_exp_df[local_gene_exp_df$Sample_type == "NOT", "Exp"])) != 0){
          p <- wilcox.test(local_gene_exp_df[local_gene_exp_df$Sample_type == genes_class, "Exp"],
                           local_gene_exp_df[local_gene_exp_df$Sample_type == "NOT", "Exp"], alternative=direction)$p.val
          CNV_pvalues <- rbind(CNV_pvalues,
                               c(tumour, genes_class, local_gene, p))
          
        }
      }
      
    }
  }
  CNV_pvalues <- as.data.frame(CNV_pvalues)
  colnames(CNV_pvalues) <- c("Tumour", "CNV", "Gene", "p_val")
  CNV_pvalues$p_val <- as.numeric(as.character(CNV_pvalues$p_val))
  
  CNV_pvalues$Gene_age <- genes_phy_categorical[match(CNV_pvalues$Gene, genes_phy_categorical$GeneID), 
                                                "Phylostrata"]
  CNV_pvalues$Gene_age <- factor(CNV_pvalues$Gene_age, levels=c("UC", "EM", "MM"))
  return(CNV_pvalues)
}

calculate_diff_exp_CNV2 <- function(CNV_genes, local_CNV, genes_class){
  CNV_pvalues <- vector()
  for(local_gene in CNV_genes){
    CNV_samples <- local_CNV[local_CNV$Gene == local_gene,"Sample"]  
    
    local_gene_exp <- unlist(patient_ratios[[tumour]][rownames(patient_ratios[[tumour]]) == local_gene,])
    
    local_gene_exp_df <- data.frame(Sample=names(local_gene_exp), Exp = local_gene_exp)
    
    local_gene_exp_df$Sample_type <- "NOT"
    local_gene_exp_df$Sample_type[local_gene_exp_df$Sample %in% CNV_samples] <- genes_class

    if(genes_class == "Amp"){
      direction <- "greater"  
    }else{
      direction <- "less"
    }
    
    if(sum(local_gene_exp_df$Sample %in% CNV_samples) >= 3 & sum(is.finite(local_gene_exp_df$Exp)) != 0){
      if(sum(is.finite(local_gene_exp_df[local_gene_exp_df$Sample_type == genes_class, "Exp"])) != 0){
        if(sum(is.finite(local_gene_exp_df[local_gene_exp_df$Sample_type == "NOT", "Exp"])) != 0){
          p <- wilcox.test(local_gene_exp_df[local_gene_exp_df$Sample_type == genes_class, "Exp"],
                           local_gene_exp_df[local_gene_exp_df$Sample_type == "NOT", "Exp"], alternative=direction)$p.val
          CNV_pvalues <- rbind(CNV_pvalues,
                               c(tumour, genes_class, local_gene, p))
          
        }
      }
      
    }
  }
  CNV_pvalues <- as.data.frame(CNV_pvalues)
  if(nrow(CNV_pvalues) != 0){
    colnames(CNV_pvalues) <- c("Tumour", "CNV", "Gene", "p_val")
    CNV_pvalues$p_val <- as.numeric(as.character(CNV_pvalues$p_val))
    
    CNV_pvalues$Gene_age <- genes_phy_categorical[match(CNV_pvalues$Gene, genes_phy_categorical$GeneID), 
                                                  "Phylostrata"]
    CNV_pvalues$Gene_age <- factor(CNV_pvalues$Gene_age, levels=c("UC", "EM", "MM"))
    
  }
  return(CNV_pvalues)
}

calculate_fraction_genes_altered_by_age <- function(genes_age){
  number_genes_altered <- vector()
  for(tumour in tumours){
    patients_with_CNVs <- patients_with_CNV_info[[tumour]]
    
    recurrent_amp <- CNV_classification$amplifications[CNV_classification$amplifications$Tumour == tumour,]
    recurrent_amp <- recurrent_amp[recurrent_amp$Genes %in% genes_age,]
    recurrent_amp_focal <- as.character(recurrent_amp[recurrent_amp$CNV_type == "FOCAL","Genes"])
    recurrent_amp_broad <- as.character(recurrent_amp[recurrent_amp$CNV_type == "BROAD","Genes"])
    
    recurrent_del <- CNV_classification$deletions[CNV_classification$deletions$Tumour == tumour,]
    recurrent_del <- recurrent_del[recurrent_del$Genes %in% genes_age,]
    recurrent_del_focal <- as.character(recurrent_del[recurrent_del$CNV_type == "FOCAL","Genes"])
    recurrent_del_broad <- as.character(recurrent_del[recurrent_del$CNV_type == "BROAD","Genes"])
    
    patients_amp <- CNVs_curated[[tumour]]$amplifications[,c("Gene", "Sample")]
    patients_amp <- patients_amp[patients_amp$Gene %in% genes_age,]
    patients_del <- CNVs_curated[[tumour]]$deletions[,c("Gene", "Sample")]
    patients_del <- patients_del[patients_del$Gene %in% genes_age,]
    
    patients_amp$Sample <- substr(patients_amp$Sample, 9, 12)
    patients_del$Sample <- substr(patients_del$Sample, 9, 12)
    
    for(patient in patients_with_CNVs){
      amp_genes <- patients_amp[patients_amp$Sample == patient,]
      del_genes <- patients_del[patients_del$Sample == patient,]
      
      if(nrow(amp_genes) != 0 && nrow(del_genes) != 0){
        amp_genes <- data.frame(amp_genes, recurrent_amp[match(amp_genes$Gene, recurrent_amp$Genes), c("Recurrent", "CNV_type"),])
        del_genes <- data.frame(del_genes, recurrent_del[match(del_genes$Gene, recurrent_del$Genes), c("Recurrent", "CNV_type"),])
        
        total_amp_genes <-  nrow(amp_genes)
        total_del_genes <- nrow(del_genes)
        
        total_recurrent_amp <- sum(amp_genes$Recurrent == "Recurrent_amp")
        total_recurrent_del <- sum(del_genes$Recurrent == "Recurrent_del")
        
        total_private_amp <- total_amp_genes - total_recurrent_amp
        total_private_del <- total_del_genes - total_recurrent_del
        
        total_recurrent <- total_recurrent_amp+total_recurrent_del
        total_private <- total_private_amp+total_private_del
        
        total_broad_amp <- sum(amp_genes$CNV_type == "BROAD")
        total_broad_del <- sum(del_genes$CNV_type == "BROAD")
        
        total_focal_amp <- sum(amp_genes$CNV_type == "FOCAL")
        total_focal_del <- sum(del_genes$CNV_type == "FOCAL")
        
        total_broad <- total_broad_amp+total_broad_del
        total_focal <- total_focal_amp+total_focal_del
        
        total_CNV <- total_amp_genes+total_del_genes
        
        temp <- c(tumour, patient, total_CNV,
                  total_amp_genes, total_del_genes,
                  total_recurrent_amp, total_recurrent_del,
                  total_private_amp, total_private_del,
                  total_recurrent, total_private,
                  total_broad_amp, total_broad_del, 
                  total_focal_amp, total_focal_del,
                  total_broad, total_focal)
        number_genes_altered <- rbind(number_genes_altered, temp)
      }else{
        temp <- c(tumour, patient, 0,
                  0, 0,
                  0, 0,
                  0, 0,
                  0, 0,
                  0, 0, 
                  0, 0,
                  0, 0)
        number_genes_altered <- rbind(number_genes_altered, temp)
      }
    }
  }
  
  number_genes_altered <- as.data.frame(number_genes_altered)
  rownames(number_genes_altered) <- NULL
  
  colnames(number_genes_altered) <- c("Tumour", "Patient", "Fraction_genes_CNV", "Fraction_amp_genes", "Fraction_del_genes",
                                      "Fraction_recurrent_amp", "Fraction_recurrent_del",
                                      "Fraction_private_amp", "Fraction_private_del",
                                      "Fraction_recurrent", "Fraction_private",
                                      "Fraction_broad_amp", "Fraction_broad_del", 
                                      "Fraction_focal_amp", "Fraction_focal_del",
                                      "Fraction_broad", "Fraction_focal")
  
  number_genes_altered[3:ncol(number_genes_altered)] <- as.data.frame(apply(number_genes_altered[3:ncol(number_genes_altered)], 2, function(x){
    return(as.numeric(as.character(x)))
  }))
  
  fraction_genes_altered <- number_genes_altered
  fraction_genes_altered[3:ncol(fraction_genes_altered)] <- fraction_genes_altered[3:ncol(fraction_genes_altered)]/length(genes_age)
  
  fraction_genes_altered_from_total <- number_genes_altered
  fraction_genes_altered_from_total[3:ncol(fraction_genes_altered_from_total)] <- fraction_genes_altered_from_total[3:ncol(fraction_genes_altered_from_total)]/length(all_genes)
  
  return(list(Per_category=fraction_genes_altered, With_respect_to_total=fraction_genes_altered_from_total))
}

load_mutations_mutsig <- function(){
  load("~/Documents/Paper 3/Objects_gistic_mutsig/mutsig2CV_sig.Rdata")

  mutations_df <- vector()
  for(tumour in tumours){
    local_mut <- mutsig2CV_sig[[tumour]]
    
    local_mut <- local_mut[,c(2, 19)]
    
    local_mut$Gene_age <- genes_phy[match(local_mut$Hugo_Symbol, genes_phy$GeneID), "Phylostrata"]
    local_mut$Tumour <- tumour
    mutations_df <- rbind(mutations_df, local_mut)
    
  }
  return(mutations_df)
}

load_CNVs_gistic <- function(only_focal){
  CNVs_df <- vector()
  for(tumour in tumours){
    local_amp <- CNVs_curated[[tumour]]$amplifications
    local_del <- CNVs_curated[[tumour]]$deletions
    
    if(only_focal == "FOCAL"){
      local_amp <- as.character(CNVs_curated_gistic_focal[[tumour]][["amplifications"]])
      local_amp <- as.character(CNVs_curated_gistic_focal[[tumour]][["deletions"]])
 
    }else if(only_focal=="ANY"){
      #No subsetting
    }
    
    if(length(local_amp) > 0 ){
      local_amp <- data.frame(Genes = local_amp,
                              CNA = "Amplification",
                              Gene_age = genes_phy[match(local_amp, genes_phy$GeneID), "Phylostrata"],
                              Tumour = tumour)
      CNVs_df <- rbind(CNVs_df, local_amp)
      
    }
    if(length(local_del) > 0 ){
      local_del <- data.frame(Genes = local_del,
                            CNA = "Deletion",
                            Gene_age = genes_phy[match(local_del, genes_phy$GeneID), "Phylostrata"],
                            Tumour = tumour)
      CNVs_df <- rbind(CNVs_df, local_del)
    }
  }
  return(CNVs_df)
}

calculate_fraction_of_patients_with_CNV_gistic <- function(CNVs_df){
  
  patient_frequency_CNV_df <- vector()
  for(tumour in tumours){
   
    n_pat_cat_amp_df <- aggregate(Gene_age ~ Pat_amp_cat, n_pat_df, table)
    n_pat_cat_amp_df <- data.frame(Frequency = n_pat_cat_amp_df$Pat_amp_cat, n_pat_cat_amp_df$Gene_age)
    colnames(n_pat_cat_amp_df)[2:17] <- paste("Phy", 1:16, sep="_")
    n_pat_cat_amp_df_melt <- melt(n_pat_cat_amp_df)
    
    n_pat_cat_del_df <- aggregate(Gene_age ~ Pat_del_cat, n_pat_df, table)
    n_pat_cat_del_df <- data.frame(Frequency = n_pat_cat_del_df$Pat_del_cat, n_pat_cat_del_df$Gene_age)
    colnames(n_pat_cat_del_df)[2:17] <- paste("Phy", 1:16, sep="_")
    n_pat_cat_del_df_melt <- melt(n_pat_cat_del_df)
    
    #frequencies <- c("<0.05", "0.05<x<0.15", "0.15<x<0.25", ">0.25")
    #frequencies <- c("<0.10", "0.10<x<0.20", ">0.20")
    frequencies <- c("<0.10", ">0.10")
    
    patient_frequency_CNV <- data.frame(Tumour =tumour, Frequency = rep(frequencies, 16), Phylostrata = rep(paste("Phy", 1:16, sep="_"), each=length(frequencies)))
    
    patient_frequency_CNV$N_genes_amp <- n_pat_cat_amp_df_melt[match(paste(patient_frequency_CNV$Frequency, patient_frequency_CNV$Phylostrata),
                                                                     paste(n_pat_cat_amp_df_melt$Frequency, n_pat_cat_amp_df_melt$variable)),"value"]
    patient_frequency_CNV$N_genes_del <- n_pat_cat_del_df_melt[match(paste(patient_frequency_CNV$Frequency, patient_frequency_CNV$Phylostrata),
                                                                     paste(n_pat_cat_del_df_melt$Frequency, n_pat_cat_del_df_melt$variable)),"value"]
    patient_frequency_CNV$N_genes_amp[is.na(patient_frequency_CNV$N_genes_amp)] <- 0
    patient_frequency_CNV$N_genes_del[is.na(patient_frequency_CNV$N_genes_del)] <- 0
    
    patient_frequency_CNV$Phylostrata <- substr(patient_frequency_CNV$Phylostrata, 5,6)
    patient_frequency_CNV$Phylostrata <- factor(patient_frequency_CNV$Phylostrata, levels=1:16)
    
    patient_frequency_CNV$Frequency <- factor(patient_frequency_CNV$Frequency, levels=rev(frequencies))
    patient_frequency_CNV$Total_genes <- n_genes_phy[match(patient_frequency_CNV$Phylostrata, n_genes_phy$Phy), "Number"]
    
    patient_frequency_CNV$Fraction_genes_amp <- patient_frequency_CNV$N_genes_amp/patient_frequency_CNV$Total_genes
    patient_frequency_CNV$Fraction_genes_del <- patient_frequency_CNV$N_genes_del/patient_frequency_CNV$Total_genes
    
    patient_frequency_CNV_df <- rbind(patient_frequency_CNV_df, patient_frequency_CNV)
  } 
  return(patient_frequency_CNV_df)
}



add_mut_CNV_to_degree_df_mutsig <- function(alt_degree_df, only_mut, only_amp, only_del){
  alt_degree_df$Mut <- NA
  alt_degree_df$Mut <- ifelse(alt_degree_df$Genes %in% only_mut, "Mut", NA)
  
  alt_degree_df$CNV <- NA
  alt_degree_df$CNV <- ifelse(alt_degree_df$Genes %in% only_amp, "Amp",
                              ifelse(alt_degree_df$Genes %in% only_del, "Del", NA))
  
  alt_degree_df$WT <- NA
  alt_degree_df$WT <- apply(alt_degree_df, 1, function(x){
    if(is.na(x[5]) & is.na(x[6])){
      return("WT")
    }else{
      return(NA)
    }
  })
  
  both_mutations <- only_mut
  both_CNVs <- only_amp[only_amp %in% only_del]
  
  alt_degree_df$WT[alt_degree_df$Genes %in% both_mutations] <- NA
  alt_degree_df$WT[alt_degree_df$Genes %in% both_CNVs] <- NA
  
  alt_degree_df1 <- alt_degree_df[,c("Genes","Degree","Database","Age","Mut")]
  alt_degree_df2 <- alt_degree_df[,c("Genes","Degree","Database","Age","CNV")]
  alt_degree_df3 <- alt_degree_df[,c("Genes","Degree","Database","Age","WT")]
  colnames(alt_degree_df1)[5] <- "Alt"
  colnames(alt_degree_df2)[5] <- "Alt"
  colnames(alt_degree_df3)[5] <- "Alt"
  
  alt_degree_df_melt <- rbind(alt_degree_df1, alt_degree_df2, alt_degree_df3)
  
  alt_degree_df_melt <- alt_degree_df_melt[!is.na(alt_degree_df_melt$Alt),]
  return(alt_degree_df_melt)
}
