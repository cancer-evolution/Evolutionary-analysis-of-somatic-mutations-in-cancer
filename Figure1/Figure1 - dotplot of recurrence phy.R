##Create dotbplot summarizing recurrance of mut and CNAs

library(ggplot2)
library(reshape2)

load("enrichment_of_recurrent_mutations.Rdata")
load("enrichment_of_recurrent_CNVs.Rdata")

#Calcualte average across tumour types

summary_recurrent_miss <- aggregate(Fraction_genes_miss ~ Phylostrata, enrichment_of_recurrent_mutations, median)
summary_recurrent_lof <- aggregate(Fraction_genes_lof ~ Phylostrata, enrichment_of_recurrent_mutations, median)


summary_recurrent_amp <- aggregate(Fraction_genes_amp ~ Phylostrata, enrichment_of_recurrent_CNVs, median)
summary_recurrent_del <- aggregate(Fraction_genes_del ~ Phylostrata, enrichment_of_recurrent_CNVs, median)

summary_recurrent_miss$Rank_miss <- rank(-summary_recurrent_miss$Fraction_genes_miss)
summary_recurrent_lof$Rank_lof <- rank(-summary_recurrent_lof$Fraction_genes_lof)
summary_recurrent_amp$Rank_amp <- rank(-summary_recurrent_amp$Fraction_genes_amp)
summary_recurrent_del$Rank_del <- rank(-summary_recurrent_del$Fraction_genes_del)

rank_recurrent <- data.frame(Phylostrata=1:16,
                             Missense=summary_recurrent_miss$Rank_miss,
                             LoF=summary_recurrent_lof$Rank_lof,
                             Amp=summary_recurrent_amp$Rank_amp,
                             Del=summary_recurrent_del$Rank_del)

rank_recurrent$Phylostrata <- factor(rank_recurrent$Phylostrata, levels=1:16)
rank_recurrent_melt <- melt(rank_recurrent)

rank_recurrent_melt$Age <- ifelse(rank_recurrent_melt$Phylostrata %in% 1:3, "UC",
                                  ifelse(rank_recurrent_melt$Phylostrata %in% 4:9, "EM", "MM"))

rank_recurrent_melt$Age <- factor(rank_recurrent_melt$Age, levels=c("UC", "EM", "MM"))

rank_recurrent_melt$variable <- factor(rank_recurrent_melt$variable,
                                       levels=rev(c("Amp", "Del", "Missense", "LoF")))

pdf("Figure1_dotplot_recurrence.pdf",
    #height=2.5, width=5)  #with_legend
    height=1.5, width=5)
g <- ggplot(rank_recurrent_melt, aes(y=variable, x=Phylostrata))+
  geom_point(aes(size=-value, color=Age))+
  #geom_text(aes(label=value))+
  #scale_size_continuous(range = c(0.5,6))+
  scale_radius()+
  #scale_size_continuous(range = c(1, 5))+
  ylab("")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
print(g)
dev.off()
