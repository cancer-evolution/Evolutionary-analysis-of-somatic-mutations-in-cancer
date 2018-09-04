##CNA classification by fraction of genome altered

library(Homo.sapiens)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(clinfun)

load("CNVs_curated.Rdata")
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

tumours <- c("LUAD", "LUSC", "BRCA", "PRAD", "LIHC", "COAD", "STAD")

CNV_chr_sample_df <- vector()
for(tumour in tumours){
  
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
      
      if(sum(is.na(local_CNVs_sample$Gene_age)) != nrow(local_CNVs_sample)){
        CNV_chr_sample <- aggregate(Gene ~ Chr+CNV, local_CNVs_sample, length)
        colnames(CNV_chr_sample)[3] <- "n_CNV"
        
        CNV_chr_sample$Total_genes <- n_genes_per_chr[match(CNV_chr_sample$Chr, n_genes_per_chr$Chromosome), "Gene"]
        
        CNV_chr_sample$Ratio <- CNV_chr_sample$n_CNV/CNV_chr_sample$Total_genes
        CNV_chr_sample$Patient <- sample
        
        CNV_chr_sample <- data.frame(local_CNVs_sample,
                                     CNV_chr_sample[match(paste(local_CNVs_sample$Chr, local_CNVs_sample$CNV),
                                                          paste(CNV_chr_sample$Chr, CNV_chr_sample$CNV)), 
                                                    c("n_CNV", "Total_genes", "Ratio")])
        
        CNV_chr_sample_mean <- aggregate(Ratio ~ CNV+Chr+Gene_age, CNV_chr_sample, mean)  ##Same value for all genes
        
        CNV_chr_sample_mean$Tumour <- tumour
        CNV_chr_sample_df <- rbind(CNV_chr_sample_df, CNV_chr_sample_mean) 
      }
    }  
  }
}

#save(CNV_chr_sample_df, file="CNV_chr_sample_df.Rdata")
load("CNV_chr_sample_df.Rdata")


CNV_chr_sample_df_mean_ratio <- aggregate(Ratio ~ CNV+Chr+Gene_age+Tumour, CNV_chr_sample_df, mean)

CNV_chr_sample_df_mean_ratio$Gene_age <- as.numeric(as.character(CNV_chr_sample_df_mean_ratio$Gene_age))

CNV_chr_sample_df_mean_ratio$Tumour <- factor(CNV_chr_sample_df_mean_ratio$Tumour, levels=tumours)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


pdf("Supp1_Fraction_chromosome_altered.pdf", width=9, height=23)
g_chr <- list()
for(chr in 1:22){
  local_chr <- paste("chr", chr, sep="")
  
  CNV_chr_sample_df_mean_ratio4 <- CNV_chr_sample_df_mean_ratio[CNV_chr_sample_df_mean_ratio$Chr == local_chr,]
  CNV_chr_sample_df_mean_ratio4$Gene_age <- factor(CNV_chr_sample_df_mean_ratio4$Gene_age, levels=1:16)
  CNV_chr_sample_df_mean_ratio4$Age <- ifelse(CNV_chr_sample_df_mean_ratio4$Gene_age %in% 1:3, "UC",
                                              ifelse(CNV_chr_sample_df_mean_ratio4$Gene_age %in% 4:9, "EM",
                                                     ifelse(CNV_chr_sample_df_mean_ratio4$Gene_age %in% 10:16, "MM", NA)))
  temp <- CNV_chr_sample_df_mean_ratio4
  temp$Gene_age <- as.numeric(as.character(temp$Gene_age))

  g_chr[[local_chr]] <- ggplot(CNV_chr_sample_df_mean_ratio4, aes(x=Gene_age, y =Ratio))+
    geom_smooth(data=temp, aes(x=Gene_age, y=Ratio,colour=Tumour), se=FALSE, size=1)+
    ylab(paste("Fraction copy-number aberrant",
               paste("Chr", chr), sep="\n"))+
    xlab("Phylostrata")+
    scale_x_continuous(breaks=1:16,limits=c(1,16))+
    facet_grid(.~CNV)+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.position="none",
          axis.text.y=element_text(size=7),
          axis.text.x=element_text(size=7))
}

do.call("grid.arrange", c(g_chr, ncol=2))

dev.off()



CNV_chr_sample_df_mean_ratio6 <- CNV_chr_sample_df_mean_ratio[CNV_chr_sample_df_mean_ratio$Chr == "chr6",]
CNV_chr_sample_df_mean_ratio6$Gene_age <- factor(CNV_chr_sample_df_mean_ratio6$Gene_age, levels=1:16)
CNV_chr_sample_df_mean_ratio6$Age <- ifelse(CNV_chr_sample_df_mean_ratio6$Gene_age %in% 1:3, "UC",
                                            ifelse(CNV_chr_sample_df_mean_ratio6$Gene_age %in% 4:9, "EM",
                                                   ifelse(CNV_chr_sample_df_mean_ratio6$Gene_age %in% 10:16, "MM", NA)))

CNV_chr_sample_df_mean_ratio6$Gene_age <- as.numeric(as.character(CNV_chr_sample_df_mean_ratio6$Gene_age))


pdf("Figure1_Fraction_Chr6_altered.pdf",
    height=2.5, width=4)
g <- ggplot(subset(CNV_chr_sample_df_mean_ratio6, CNV=="Amp"), aes(x=Gene_age, y =Ratio))+
  geom_smooth(aes(x=Gene_age, y=Ratio,colour=Tumour), se=FALSE, size=1)+
  ylab(paste("Fraction of chromosome altered",
             paste("Chr", 6), sep="\n"))+
  xlab("Phylostrata")+
  scale_x_continuous(breaks=1:15,limits=c(1,15))+
  #facet_grid(.~CNV)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.y=element_text(size=9),
        axis.text.x=element_text(size=9))
print(g)
dev.off()


groups <- temp$Gene_age
groups <- factor(groups,levels=1:15, ordered=TRUE) 
p_trend <- vector()
for(tumour in tumours){
  temp <- CNV_chr_sample_df_mean_ratio6[CNV_chr_sample_df_mean_ratio6$Tumour == tumour,]
  temp <- temp[temp$CNV == "Amp",]
  p_trend <- c(p_trend, jonckheere.test(g=groups, x=temp$Ratio, alternative="increasing")$p.value)
  
}

p_trend <- p.adjust(p_trend, method="BH")


