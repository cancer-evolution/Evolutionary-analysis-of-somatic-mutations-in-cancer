##Create dotbplot summarizing recurrance of mut and CNAs

library(ggplot2)
library(reshape2)

load("enrichment_of_recurrent_mutations.Rdata")
load("enrichment_of_recurrent_CNVs.Rdata")

#Calcualte average across tumour types

summary_recurrent_mut <- aggregate(Rank_mut1 ~ Phylostrata, enrichment_of_recurrent_mutations, median)

summary_recurrent_amp <- aggregate(Rank_amp1 ~ Phylostrata, enrichment_of_recurrent_CNVs, median)
summary_recurrent_del <- aggregate(Rank_del1 ~ Phylostrata, enrichment_of_recurrent_CNVs, median)

rank_recurrent <- data.frame(Phylostrata=1:16,
                             Mutations=summary_recurrent_mut[match(1:16, summary_recurrent_mut$Phylostrata), "Rank_mut1"],
                             Amp=summary_recurrent_amp[match(1:16, summary_recurrent_amp$Phylostrata), "Rank_amp1"],
                             Del=summary_recurrent_del[match(1:16, summary_recurrent_del$Phylostrata), "Rank_del1"])

rank_recurrent$Phylostrata <- factor(rank_recurrent$Phylostrata, levels=1:16)
rank_recurrent_melt <- melt(rank_recurrent)

rank_recurrent_melt$Age <- ifelse(rank_recurrent_melt$Phylostrata %in% 1:3, "UC",
                                  ifelse(rank_recurrent_melt$Phylostrata %in% 4:9, "EM", "MM"))

rank_recurrent_melt$Age <- factor(rank_recurrent_melt$Age, levels=c("UC", "EM", "MM"))

rank_recurrent_melt$variable <- factor(rank_recurrent_melt$variable,
                                       levels=rev(c("Amp", "Del", "Mutations")))

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

