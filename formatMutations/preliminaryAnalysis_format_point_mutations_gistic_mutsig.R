###Filtering of variants. Select variants from MutSig

library(stringr)
library(readr)

tumours <- c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "ESCA", "GBM", 
             "HNSC", "KICH", "KIRC", "KIRP", "LGG",  "LIHC", "LUAD", "LUSC", 
             "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", 
             "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM")

mutsig2CV <- list()

for(tumour in tumours){
  file2CV <- paste("gdac.broadinstitute.org_",
                   tumour, "-TP.MutSigNozzleReport2CV.Level_4.2016012800.0.0/",
                   "sig_genes.txt", sep="")
  if(tumour == "SKCM"){
    file2CV <- paste("gdac.broadinstitute.org_",
                     tumour, "-TM.MutSigNozzleReport2CV.Level_4.2016012800.0.0/",
                     "sig_genes.txt", sep="")
  }
  
  if(file.exists(file2CV)){
    mutsig2CV[[tumour]] <- read.delim(file2CV)
  }else{
    print("MutSig2CV")
    print(tumour)
  }
}

#Select significant genes
mutsig2CV_sig <- list()

for(tumour in tumours){
  if(!is.null(nrow(mutsig2CV[[tumour]]))){
    mutsig2CV[[tumour]]$q <- as.numeric(as.character(gsub("<","", mutsig2CV[[tumour]]$q)))
    mutsig2CV_sig[[tumour]] <- subset(mutsig2CV[[tumour]], q < 0.05)
    colnames(mutsig2CV_sig[[tumour]])[2] <- "Hugo_Symbol"
  }
}

save(mutsig2CV_sig, file="mutsig2CV_sig.Rdata")

