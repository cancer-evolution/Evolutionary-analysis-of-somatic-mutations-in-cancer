##CNA classification by fraction of genome altered

library(Homo.sapiens)
library(reshape2)
library(ggplot2)
library(gridExtra)

load("CNVs_curated2.Rdata")
load("patients_with_CNV_info.Rdata")


genes_chromosome <- as.data.frame(genes(Homo.sapiens, columns="SYMBOL"))
genes_chromosome <- genes_chromosome[,c("seqnames", "SYMBOL")]
colnames(genes_chromosome) <- c("Chromosome", "Gene")
genes_chromosome$Gene <- as.character(genes_chromosome$Gene)

n_genes_per_chr <- aggregate(Gene ~ Chromosome, genes_chromosome, length)

genes_phy <- read.csv("gene_phylostrata.txt")  ###File with gene name, entrez and phylostratum as column
genes_phy$Age <- ifelse(genes_phy$Phylostrata %in% 1:3, "UC",
                        ifelse(genes_phy$Phylostrata %in% 4:9, "EM", "MM"))
UC_genes <- as.character(genes_phy[genes_phy$Age == "UC", "GeneID"])
EM_genes <- as.character(genes_phy[genes_phy$Age == "EM", "GeneID"])
MM_genes <- as.character(genes_phy[genes_phy$Age == "MM", "GeneID"])

tumours <- c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "ESCA", "GBM", 
             "HNSC", "KICH", "KIRC", "KIRP", "LGG",  "LIHC", "LUAD", "LUSC", 
             "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", 
             "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM")  

##Select CNVs that are in more focal regions
##Calculate percentage of chromosome affected when a gene is altered


for(tumour in tumours){
  context_of_CNVs <- vector()
  aveg_context_of_CNVs <- vector()
  
  local_amp <- CNVs_curated[[tumour]]$amplifications
  local_del <- CNVs_curated[[tumour]]$deletions
  
  local_amp$Sample <- substr(local_amp$Sample, 9, 12)
  local_del$Sample <- substr(local_del$Sample, 9, 12)
  
  local_amp$CNV <- "Amp"
  local_del$CNV <- "Del"
  
  local_CNVs <- rbind(local_amp, local_del)
  local_CNVs$Chr <- genes_chromosome[match(local_CNVs$Gene, genes_chromosome$Gene), "Chromosome"]
  local_CNVs$Gene_age <- genes_phy[match(local_CNVs$Gene, genes_phy$GeneID), "Phylostrata"]
  local_CNVs$Gene_age <- as.factor(local_CNVs$Gene_age)
  
  
  for(chr in unique(local_CNVs$Chr)){
    local_CNVs_chr <- local_CNVs[local_CNVs$Chr == chr,]
    
    for(sample in unique(local_CNVs_chr$Sample)){
      local_CNVs_sample <- local_CNVs_chr[local_CNVs_chr$Sample == sample,]
      
      for(CNV_type in unique(local_CNVs_sample$CNV)){
        
        local_CNVs_sample2 <- local_CNVs_sample[local_CNVs_sample$CNV == CNV_type,]
        
        if(sum(is.na(local_CNVs_sample2$Gene_age)) != nrow(local_CNVs_sample2)){
          genes_with_CNV <- local_CNVs_sample2$Gene
          CNV_chr_sample <- aggregate(Gene ~ Chr+CNV, local_CNVs_sample2, length)
          colnames(CNV_chr_sample)[3] <- "n_CNV"
          
          CNV_chr_sample$Total_genes <- n_genes_per_chr[match(CNV_chr_sample$Chr, n_genes_per_chr$Chromosome), "Gene"]
          
          CNV_chr_sample$Ratio <- CNV_chr_sample$n_CNV/CNV_chr_sample$Total_genes
          CNV_chr_sample$Patient <- sample
          
          context_of_CNVs <- rbind(context_of_CNVs,
                                   data.frame(tumour, CNV_chr_sample[,c("Patient", "Chr", "CNV", "Total_genes", "Ratio")],
                                              genes_with_CNV))
          
        }
      }
      
    }
    print(chr)
  }
  local_context_of_CNVs <- context_of_CNVs[context_of_CNVs$tumour == tumour,]
  local_aveg_context_of_CNVs <- aggregate(Ratio ~ genes_with_CNV+CNV+Chr+tumour, context_of_CNVs, median)
  
  aveg_context_of_CNVs <- rbind(aveg_context_of_CNVs, local_aveg_context_of_CNVs)
  
  aveg_context_of_CNVs$Gene_age <- genes_phy[match(aveg_context_of_CNVs$genes_with_CNV, genes_phy$GeneID), "Age"]
  colnames(aveg_context_of_CNVs) <- c("Genes","CNV","Chr","Tumour","Ratio","Gene_age")
  save(aveg_context_of_CNVs, file=paste("aveg_context_of_CNVs_median_",
  tumour, ".Rdata", sep=""))
  
  print(tumour)
  
}

aveg_context_of_CNVs_temp <- vector()
for(tumour in tumours){
  
  load(paste("aveg_context_of_CNVs_median_",
        tumour, ".Rdata", sep=""))
  
  aveg_context_of_CNVs_temp <- rbind(aveg_context_of_CNVs_temp,
                                     aveg_context_of_CNVs)
}

aveg_context_of_CNVs <- aveg_context_of_CNVs_temp

first_quantile_cutoff <- vector()
par(mfrow=c(2,1))
pdf("Distribution of fraction of chr altered with CNVs.pdf", height=5)
for(tumour in tumours){
  local_aveg_context_of_CNVs <- aveg_context_of_CNVs[aveg_context_of_CNVs$Tumour == tumour,]
  hist(local_aveg_context_of_CNVs$Ratio[local_aveg_context_of_CNVs$CNV == "Amp"],
       main=paste("Fraction of chromosome altered", "Amp", tumour, sep="\n"))
  abline(v=quantile(local_aveg_context_of_CNVs$Ratio[local_aveg_context_of_CNVs$CNV == "Amp"], 0.25),
         col="red")
  
  hist(local_aveg_context_of_CNVs$Ratio[local_aveg_context_of_CNVs$CNV == "Del"],
       main=paste("Fraction of chromosome altered", "Del", tumour, sep="\n"))
  abline(v=quantile(local_aveg_context_of_CNVs$Ratio[local_aveg_context_of_CNVs$CNV == "Del"], 0.25),
         col="red")
  
  first_quantile_cutoff <- rbind(first_quantile_cutoff,
                                 c(tumour, "Amp", "0.25", quantile(local_aveg_context_of_CNVs$Ratio[local_aveg_context_of_CNVs$CNV == "Amp"], 0.25)),
                                 c(tumour, "Del", "0.25", quantile(local_aveg_context_of_CNVs$Ratio[local_aveg_context_of_CNVs$CNV == "Del"], 0.25)),
                                 c(tumour, "Amp", "0.10", quantile(local_aveg_context_of_CNVs$Ratio[local_aveg_context_of_CNVs$CNV == "Amp"], 0.10)),
                                 c(tumour, "Del", "0.10", quantile(local_aveg_context_of_CNVs$Ratio[local_aveg_context_of_CNVs$CNV == "Del"], 0.10)))
  
}
dev.off()

first_quantile_cutoff <- as.data.frame(first_quantile_cutoff)
colnames(first_quantile_cutoff) <- c("Tumour", "CNV", "Quantile", "Ratio_cutoff")
save(first_quantile_cutoff, file="first_quantile_cutoff.Rdata")

load("first_quantile_cutoff.Rdata")
aveg_context_of_CNVs <- aveg_context_of_CNVs[!is.na(aveg_context_of_CNVs$Genes),]

CNVs_above_fraction_0.25 <- list()
CNVs_above_fraction_0.10 <- list()
for(tumour in tumours){
  for(CNV_type in c("Amp", "Del")){
    local_aveg_context_of_CNVs <- aveg_context_of_CNVs[aveg_context_of_CNVs$Tumour == tumour,]
    local_aveg_context_of_CNVs <- local_aveg_context_of_CNVs[local_aveg_context_of_CNVs$CNV == CNV_type,]
    local_aveg_context_of_CNVs <- unique(local_aveg_context_of_CNVs)
    local_cutoff <- first_quantile_cutoff[intersect(which(first_quantile_cutoff$Tumour == tumour),
                                                    which(first_quantile_cutoff$CNV == CNV_type)),]
    cutoff_0.25 <- as.numeric(as.character(local_cutoff[local_cutoff$Quantile == "0.25", "Ratio_cutoff"]))
    cutoff_0.10 <- as.numeric(as.character(local_cutoff[local_cutoff$Quantile == "0.10", "Ratio_cutoff"]))
    genes_0.25 <- local_aveg_context_of_CNVs[local_aveg_context_of_CNVs$Ratio <= cutoff_0.25,"Genes"]
    genes_0.10 <- local_aveg_context_of_CNVs[local_aveg_context_of_CNVs$Ratio <= cutoff_0.10,"Genes"]
    
    CNVs_above_fraction_0.25[[tumour]][[CNV_type]] <- genes_0.25
    CNVs_above_fraction_0.10[[tumour]][[CNV_type]] <- genes_0.10
  }
}

save(CNVs_above_fraction_0.25, file="CNVs_above_fraction_0.25.Rdata")
save(CNVs_above_fraction_0.10, file="CNVs_above_fraction_0.10.Rdata")
