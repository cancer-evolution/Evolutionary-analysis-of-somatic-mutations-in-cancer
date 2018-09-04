library(ggplot2)
library(reshape2)
library(ggsci)

load("regulator_classification.Rdata")

cancer_census <- read.delim("Census_allTue Dec 12 03_33_11 2017.csv", sep=",")  #File from Cancer Census database

cancer_census$Regulator_class <- regulator_classification[match(cancer_census$Gene.Symbol, regulator_classification$Regulator), "Regulator_class"]


n_census <- length(unique(cancer_census$Gene.Symbol))
n_census_reg <- sum(unique(cancer_census$Gene.Symbol) %in% regulator_classification$Regulator)
n_reg <- nrow(regulator_classification)

sum(unique(cancer_census$Gene.Symbol) %in% regulator_classification$Regulator)/length(unique(cancer_census$Gene.Symbol))

n_census_class <- table(cancer_census$Regulator_class)

n_census_class <- data.frame(Class = names(n_census_class),
                             Number = as.vector(n_census_class))
n_census_class <- n_census_class[n_census_class$Class != "None",]
n_census_class <- n_census_class[n_census_class$Class != "MM_downstream",]

n_census_class$Total <- sum(n_census_class$Number)
n_census_class$Per <- (n_census_class$Number/n_census_class$Total)*100



n_census_class$Class <- factor(n_census_class$Class,
                               levels=c("UC_downstream", "Mixed_downstream", "EM_downstream"))

pdf("Figure5_Pie_chart_census.pdf",
    height=4, width=5)
g <- ggplot(n_census_class, aes(x="", y=Number, fill=Class))+
  geom_bar(width = 1, stat = "identity")+
  scale_fill_jco()+
  coord_polar("y", start=0)+
  theme(axis.text.x=element_blank())+
  theme_void()
print(g)
dev.off()

78/sum(n_census_class$Number)
90/sum(n_census_class$Number)
24/sum(n_census_class$Number)


#Null
n_null <- aggregate(Regulator ~ Regulator_class, regulator_classification, length)
n_null <- n_null[n_null$Regulator_class %in% c("UC_downstream", "Mixed_downstream", "EM_downstream"),]

n_null$Class <- factor(n_null$Regulator_class,
                               levels=c("UC_downstream", "Mixed_downstream", "EM_downstream"))

pdf("Figure5_Pie_chart_null.pdf",
    height=4, width=5)
g <- ggplot(n_null, aes(x="", y=Regulator, fill=Class))+
  geom_bar(width = 1, stat = "identity")+
  scale_fill_jco()+
  coord_polar("y", start=0)+
  theme(axis.text.x=element_blank())+
  theme_void()
print(g)
dev.off()


total_reg <- sum(n_null$Regulator)

n_null$Per <- (n_null$Regulator/total_reg)*100


n_total_class <- table(regulator_classification$Regulator_class)
n_total_class <- data.frame(Class = names(n_total_class),
                            Number = as.vector(n_total_class))
n_total_class <- n_total_class[n_total_class$Class != "None",]

n_regulators <- sum(n_total_class$Number)

n_census <- sum(n_census_class$Number)


fisher.test(cbind(c(n_census_class[n_census_class$Class == "EM_downstream", "Number"],
                    n_total_class[n_total_class$Class == "EM_downstream", "Number"]),
                  c(n_census, n_regulators)), alternative="greater")$p.value

fisher.test(cbind(c(n_census_class[n_census_class$Class == "UC_downstream", "Number"],
                    n_total_class[n_total_class$Class == "UC_downstream", "Number"]),
                  c(n_census, n_regulators)), alternative="greater")$p.value

fisher.test(cbind(c(n_census_class[n_census_class$Class == "Mixed_downstream", "Number"],
                    n_total_class[n_total_class$Class == "Mixed_downstream", "Number"]),
                  c(n_census, n_regulators)), alternative="greater")$p.value

fisher.test(cbind(c(n_census_class[n_census_class$Class == "MM_downstream", "Number"],
                    n_total_class[n_total_class$Class == "MM_downstream", "Number"]),
                  c(n_census, n_regulators)), alternative="greater")$p.value


genes_phy <- read.csv("gene_phylostrata.txt")  ###File with gene name, entrez and phylostratum as column
genes_phy$Age <- ifelse(genes_phy$Phylostrata %in% 1:3, "UC",
                        ifelse(genes_phy$Phylostrata %in% 4:9, "EM",
                               ifelse(genes_phy$Phylostrata %in% 10:16, "MM", NA)))
UC_genes <- as.character(genes_phy[genes_phy$Age == "UC", "GeneID"])
EM_genes <- as.character(genes_phy[genes_phy$Age == "EM", "GeneID"])
MM_genes <- as.character(genes_phy[genes_phy$Age == "MM", "GeneID"])


n_census_UC <- sum(cancer_census$Gene.Symbol %in% UC_genes)
n_census_EM <- sum(cancer_census$Gene.Symbol %in% EM_genes)
n_census_MM <- sum(cancer_census$Gene.Symbol %in% MM_genes)
n_census <- n_census_UC+n_census_EM+n_census_MM

n_total_UC <- length(UC_genes)
n_total_EM <- length(EM_genes)
n_total_MM <- length(MM_genes)
n_total <- n_total_UC+n_total_EM+n_total_MM

n_census_UC/n_census
n_census_EM/n_census
n_census_MM/n_census

fisher.test(cbind(c(n_census_UC, n_total_UC),
                  c(n_census, n_total)), alternative="greater")

fisher.test(cbind(c(n_census_EM, n_total_EM),
                  c(n_census, n_total)), alternative="greater")

fisher.test(cbind(c(n_census_MM, n_total_MM),
                  c(n_census, n_total)), alternative="greater")
