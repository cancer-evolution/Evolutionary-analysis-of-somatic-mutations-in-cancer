####Neighbourhood of regulators by somatic mutation class

library(ggplot2)
library(igraph)
library(reshape2)
library(gridExtra)
library(readr)
library(clinfun)
library(ggrepel)
library(ggsci)

source("helper_functions.R")

genes_phy <- read.csv("gene_phylostrata.txt")  ###File with gene name, entrez and phylostratum as column
genes_phy$Age <- ifelse(genes_phy$Phylostrata %in% 1:3, "UC",
                        ifelse(genes_phy$Phylostrata %in% 4:9, "EM",
                               ifelse(genes_phy$Phylostrata %in% 10:16, "MM", NA)))

UC_genes <- as.character(genes_phy[which(genes_phy$Phylostrata %in% 1:3),1])
EM_genes <- as.character(genes_phy[which(genes_phy$Phylostrata %in% 4:9),1])
MM_genes <- as.character(genes_phy[which(genes_phy$Phylostrata %in% 10:16),1])
tumours <- c("LUAD", "LUSC", "BRCA", "PRAD", "LIHC", "COAD", "STAD")

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
  miss_df <- rbind(miss_df, cbind(Tumour=tumour, Alt="Miss", Genes=missense))
  lof_df <- rbind(lof_df, cbind(Tumour=tumour, Alt="LoF", Genes=lof))
  
  all_miss <- c(all_miss, missense)
  all_lof <- c(all_lof, lof)
}



##Read in CNVs
##Using only focal ones
load("CNVs_curated2.Rdata")
load("patients_with_CNV_info.Rdata")

load("CNVs_above_fraction_0.25.Rdata")

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
all_amp_common <- names(table(all_amp))[table(all_amp)>=3]
all_del_common <- names(table(all_del))[table(all_del)>=3]
all_CNV_common_both <- all_amp_common[all_amp_common %in% all_del_common]
all_amp_common <- all_amp_common[!(all_amp_common %in% all_CNV_common_both)]
all_del_common <- all_del_common[!(all_del_common %in% all_CNV_common_both)]


all_miss_common <- names(table(all_miss))[table(all_miss)>=3]
all_lof_common <- names(table(all_lof))[table(all_lof)>=3]
all_mut_common_both <- all_miss_common[all_miss_common %in% all_lof_common]
all_miss_common <- all_miss_common[!(all_miss_common %in% all_mut_common_both)]
all_lof_common <- all_lof_common[!(all_lof_common %in% all_mut_common_both)]


####Percentage EM of downstream genes
load("percentage_UC_per_gene_directed_PC2.Rdata")

percentage_UC_per_gene_directed$Database <- "RN"
percentage_UC_per_gene_directed$Genes <- rownames(percentage_UC_per_gene_directed)
percentage_UC_per_gene_directed$Age <- percentage_UC_per_gene_directed$Gene_age

percentage_UC_EM_MM_alt <- add_mut_CNV_to_degree_df(percentage_UC_per_gene_directed[,c("Genes","Degree","Age","Database")], all_miss_common, all_lof_common,
                                                    all_amp_common, all_del_common)


##Classify targets into quantiles based on age

##Separate by %UC and %EM

#At least >2/3 UC to be UC, at least >2/3 EM to be EM, and to be mixed at least 10% of each

per_UC_EM_downstream <- percentage_UC_per_gene_directed[,c("Per_UC", "Per_EM", "Per_MM", "Gene_age", "Degree")]
colnames(per_UC_EM_downstream)[4] <- "Age"
per_UC_EM_downstream$Database <- "RN"
per_UC_EM_downstream$Genes <- rownames(per_UC_EM_downstream)

per_UC_EM_downstream$Alt <- percentage_UC_EM_MM_alt[match(per_UC_EM_downstream$Genes, percentage_UC_EM_MM_alt$Genes), "Alt"]

per_UC_EM_downstream$Alt_binary <- ifelse(per_UC_EM_downstream$Alt %in% c("Miss", "LoF"), "Point",
                                           ifelse(per_UC_EM_downstream$Alt %in% c("Amp", "Del"), "CNV", "WT"))

per_UC_EM_downstream$Regulator_class <- ifelse(per_UC_EM_downstream$Per_UC >= 200/3, "UC_downstream",
                                               ifelse(per_UC_EM_downstream$Per_EM >= 200/3, "EM_downstream", 
                                                      ifelse(per_UC_EM_downstream$Per_MM >= 200/3, "MM_downstream",
                                                        ifelse(per_UC_EM_downstream$Per_UC >= 10 & per_UC_EM_downstream$Per_EM >= 10 & per_UC_EM_downstream$Per_MM < 10, "Mixed_downstream", "None"))))



regulator_classification <- per_UC_EM_downstream[,c("Per_UC", "Per_EM", "Per_MM", "Regulator_class", "Alt_binary")]
regulator_classification$Regulator <- rownames(regulator_classification)

save(regulator_classification, file="regulator_classification.Rdata")

sum(regulator_classification$Regulator_class == "UC_downstream")
sum(regulator_classification$Regulator_class == "EM_downstream")
sum(regulator_classification$Regulator_class == "Mixed_downstream")
sum(regulator_classification$Regulator_class == "MM_downstream")
sum(regulator_classification$Regulator_class == "None")

regulator_classification_print <- regulator_classification[,c("Regulator", "Regulator_class")]
regulator_classification_print$Regulator_class[regulator_classification_print$Regulator_class == "UC_downstream"] <- "UC-t"
regulator_classification_print$Regulator_class[regulator_classification_print$Regulator_class == "EM_downstream"] <- "EM-t"
regulator_classification_print$Regulator_class[regulator_classification_print$Regulator_class == "Mixed_downstream"] <- "UC/EM-i"
regulator_classification_print$Regulator_class[regulator_classification_print$Regulator_class == "MM_downstream"] <- "MM-t"
regulator_classification_print$Regulator_class[regulator_classification_print$Regulator_class == "None"] <- NA

write.table(regulator_classification_print, file="Supp_table1_Regulator_classifications.txt",
            quote=FALSE, sep="\t", row.names=FALSE)



#Sort regulators
regulator_classification_UC <- regulator_classification[regulator_classification$Regulator_class == "UC_downstream",]
regulator_classification_EM <- regulator_classification[regulator_classification$Regulator_class == "EM_downstream",]
regulator_classification_MM <- regulator_classification[regulator_classification$Regulator_class == "MM_downstream",]
regulator_classification_Mixed <- regulator_classification[regulator_classification$Regulator_class == "Mixed_downstream",]

order_regulators_UC <- regulator_classification_UC$Regulator[order(regulator_classification_UC$Per_UC,
                                                                   regulator_classification_UC$Per_EM, decreasing=TRUE)]
order_regulators_EM <- regulator_classification_EM$Regulator[order(regulator_classification_EM$Per_UC,
                                                                   regulator_classification_EM$Per_EM, decreasing=TRUE)]
order_regulators_MM <- regulator_classification_MM$Regulator[order(regulator_classification_MM$Per_UC,
                                                                   regulator_classification_MM$Per_EM, decreasing=TRUE)]
order_regulators_Mixed <- regulator_classification_Mixed$Regulator[order(regulator_classification_Mixed$Per_UC,
                                                                         regulator_classification_Mixed$Per_EM, decreasing=TRUE)]

order_regulators_all <- rev(c(order_regulators_UC, order_regulators_Mixed, 
                              order_regulators_EM, order_regulators_MM))

regulator_classification2 <- rbind(regulator_classification_UC[match(order_regulators_UC, regulator_classification_UC$Regulator),],
                                   regulator_classification_Mixed[match(order_regulators_Mixed, regulator_classification_Mixed$Regulator),],
                                   regulator_classification_EM[match(order_regulators_EM, regulator_classification_EM$Regulator),],
                                   regulator_classification_MM[match(order_regulators_MM, regulator_classification_MM$Regulator),])

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

p.adjust(p_enrich_regulator_ages$value)

p_enrich_regulator_ages

colnames(per_each_age) <- c("Regulator_class", "Per_UC", "Per_EM", "Per_MM")




regulator_classification_melt <- melt(regulator_classification2)

colnames(regulator_classification_melt) <- c("Regulator_class", "Alt_binary", "Regulator", "Downstream_age", "Percentage")
regulator_classification_melt$Downstream_age <- factor(regulator_classification_melt$Downstream_age,
                                                       levels=c("Per_MM", "Per_EM", "Per_UC"))
regulator_classification_melt$Regulator <- factor(regulator_classification_melt$Regulator,
                                                  levels=order_regulators_all)
regulator_classification_melt$Regulator_age <- genes_phy[match(regulator_classification_melt$Regulator,
                                                               genes_phy$GeneID), "Age"]
regulator_classification_melt$Regulator_class <- factor(regulator_classification_melt$Regulator_class,
                                                        levels=c("UC_downstream", "Mixed_downstream", "EM_downstream", 
                                                                 "MM_downstream", "None"))

regulator_classification_melt$Alt_binary[regulator_classification_melt$Alt_binary == "WT"] <- NA
regulator_classification_melt$Alt_binary <- factor(regulator_classification_melt$Alt_binary,
                                                   levels=rev(c("CNV", "Point")))

regulator_classification_melt$Regulator <- factor(regulator_classification_melt$Regulator,
                                                  levels=rev(levels(regulator_classification_melt$Regulator)))

table(regulator_classification_melt[regulator_classification_melt$Alt_binary == "Point","Regulator_class"])
table(regulator_classification_melt[regulator_classification_melt$Alt_binary == "CNV","Regulator_class"])
table(regulator_classification_melt[is.na(regulator_classification_melt$Alt_binary),"Regulator_class"])


pdf("Figure3_Barplot_regulators.pdf",
    width=9, height=4)
g <- ggplot(subset(regulator_classification_melt, Regulator_class != "MM_downstream"), aes(y=Percentage, x=Regulator))+
  geom_bar(aes(fill=Downstream_age), stat='identity', position='stack', width=1)+
  #geom_area(aes(fill=Downstream_age), position = "stack", group=1)+
#  geom_point(aes(y=105, colour=Alt_binary), size=0.25)+
  geom_line(aes(y=-5, group=Regulator_class, colour=Regulator_class), size=3)+
  scale_fill_manual(values=c(Per_UC=gg_color_hue(3)[1],
                             Per_EM=gg_color_hue(3)[2],
                             Per_MM=gg_color_hue(3)[3]))+
  scale_colour_manual(values=c(CNV="red",
                               Point="blue",
                               UC_downstream="royalblue",
                               EM_downstream="grey",
                               MM_downstream="brown",
                               Mixed_downstream="gold"))+
  #ylim(-10,105)+
  xlab("")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.background = element_blank())
print(g)
dev.off()





regulator_classification_alt <- regulator_classification_melt
regulator_classification_alt <- regulator_classification_alt[regulator_classification_alt$Downstream_age == "Per_UC",]
regulator_classification_alt$Alt_binary <- as.character(regulator_classification_alt$Alt_binary)
regulator_classification_alt$Alt_binary[is.na(regulator_classification_alt$Alt_binary)] <- "WT"
regulator_classification_alt$Position <- nrow(regulator_classification_alt):1
regulator_classification_alt$N_point <- ifelse(regulator_classification_alt$Alt_binary == "Point", 1, 0)
regulator_classification_alt$N_CNV <- ifelse(regulator_classification_alt$Alt_binary == "CNV", 1, 0)

regulator_classification_alt_both <- regulator_classification_alt[regulator_classification_alt$Alt_binary %in% c("Point", "CNV"),]

regulator_classification_alt_both$Alt_binary <- factor(regulator_classification_alt_both$Alt_binary,
                                                       levels=c("Point", "CNV"))

regulator_classification_alt$Alt_binary[regulator_classification_alt$Alt_binary == "WT"] <- NA
regulator_classification_alt <- regulator_classification_alt[!is.na(regulator_classification_alt$Alt_binary),]

pdf("Figure3_Density_alterations.pdf",
    height=1.5, width=7)
g <- ggplot(subset(regulator_classification_alt, Regulator_class != "MM_downstream"), aes(x=Position))+
  geom_density(aes(y=..scaled..,colour=Alt_binary, fill=Alt_binary))+
  #coord_flip()+
  xlim(max(regulator_classification_alt_both$Position), 0)+
  scale_fill_manual(values=c(Point="grey48", CNV="grey87"))+
  scale_colour_manual(values=c(Point="grey48", CNV="grey87"))+
  #ylim(0.0002, 0.0003)+  
  #scale_x_reverse()+
  theme_bw()+
  theme(legend.position="none")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
print(g)
dev.off()





quantiles1 <- count_regulators_each_class(per_UC_EM_downstream)

quantiles1$Alt_binary <- factor(quantiles1$Alt_binary, levels=c("WT", "CNV", "Point"))

quantiles1$Regulator_class <- factor(quantiles1$Regulator_class,
                                     levels=c("EM_downstream", "Mixed_downstream", "UC_downstream"))

pdf("Figure3_Stacked_barplot_alterations.pdf",
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


temp <- subset(quantiles1, Regulator_class %in% c("UC_downstream", "Mixed_downstream","EM_downstream"))
temp <- temp[temp$Alt_binary != "WT",]


pdf("Supp3_Pie chart UC-EM-t.pdf",
    height=4, width=4)
g <- ggplot(subset(temp, Regulator_class=="Mixed_downstream"), aes(x=1, y=Number_of_regulators))+
  geom_bar(aes(fill=Alt_binary), stat='identity')+
  scale_fill_manual(values=c(Point="purple", CNV="orange"))+
  coord_polar("y")+
  ggtitle("UC/EM-t")+
  theme_void()
print(g)  
dev.off()

pdf("Supp3_Pie chart EM-t.pdf",
    height=4, width=4)
g <- ggplot(subset(temp, Regulator_class=="EM_downstream"), aes(x=1, y=Number_of_regulators))+
  geom_bar(aes(fill=Alt_binary), stat='identity')+
  scale_fill_manual(values=c(Point="purple", CNV="orange"))+
  coord_polar("y")+
  ggtitle("EM-t")+
  theme_void()
print(g)  
dev.off()

pdf("Supp3_Pie chart UC-t.pdf",
    height=4, width=4)
g <- ggplot(subset(temp, Regulator_class=="UC_downstream"), aes(x=1, y=Number_of_regulators))+
  geom_bar(aes(fill=Alt_binary), stat='identity')+
  scale_fill_manual(values=c(Point="purple", CNV="orange"))+
  coord_polar("y")+
  ggtitle("UC-t")+
  theme_void()
print(g)
dev.off()



median(per_UC_EM_downstream[per_UC_EM_downstream$Alt_binary == "Point","Per_UC"])
median(per_UC_EM_downstream[per_UC_EM_downstream$Alt_binary == "Point","Per_EM"], na.rm=TRUE)

median(per_UC_EM_downstream[per_UC_EM_downstream$Alt_binary == "CNV","Per_UC"])
median(per_UC_EM_downstream[per_UC_EM_downstream$Alt_binary == "CNV","Per_EM"], na.rm=TRUE)

median(per_UC_EM_downstream[per_UC_EM_downstream$Alt_binary == "WT","Per_UC"])
median(per_UC_EM_downstream[per_UC_EM_downstream$Alt_binary == "WT","Per_EM"], na.rm=TRUE)


wilcox.test(per_UC_EM_downstream[per_UC_EM_downstream$Alt_binary == "Point","Per_EM"],
            per_UC_EM_downstream[per_UC_EM_downstream$Alt_binary == "CNV","Per_EM"], alternative="less")

wilcox.test(per_UC_EM_downstream[per_UC_EM_downstream$Alt_binary == "Point","Per_UC"],
            per_UC_EM_downstream[per_UC_EM_downstream$Alt_binary == "CNV","Per_UC"], alternative="greater")


wilcox.test(per_UC_EM_downstream[per_UC_EM_downstream$Alt_binary == "Point","Per_EM"],
            per_UC_EM_downstream[per_UC_EM_downstream$Alt_binary == "WT","Per_EM"], alternative="less")

wilcox.test(per_UC_EM_downstream[per_UC_EM_downstream$Alt_binary == "Point","Per_UC"],
            per_UC_EM_downstream[per_UC_EM_downstream$Alt_binary == "WT","Per_UC"], alternative="greater")
