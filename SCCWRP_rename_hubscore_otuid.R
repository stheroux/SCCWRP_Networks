setwd("/Users/Mel/Desktop/SCCWRP")
library(dplyr)
library(tidyverse)
library(phyloseq)

hubscores = read.csv("network_hubscores.csv", row.names = 1)
hubscores$TaxaName = NA
renamed = data.frame()
asv16S <- read.csv("16S_2022.csv", stringsAsFactors=FALSE, row.names = 1)
asv18S <- read.csv("18S_2022.csv", stringsAsFactors=FALSE, row.names = 1)
asvrbcL <- read.csv("rbcL_2022.csv", stringsAsFactors=FALSE, row.names = 1)

primers = c("asv16S","asv18S","asvrbcL")

for(primer in primers){
  #Format abundance table and create taxonomy table as inputs for phyloseq
  otumat <- get(primer)
  #print("got primer")
  taxmat <- otumat %>% select("ConsensusLineage")
  #print("taxmat select lineage")
  levels <- c("Domain","Phylum","Class","Order","Family","Genus","Species")
  taxmat <- taxmat %>% separate(ConsensusLineage, levels, sep = ";")
  #print("separated lineage")
  taxmat <- as.matrix(taxmat) 
  
  if(primer == "asv16S"){
    taxmat_16S <- taxmat
  }
  
  if(primer == "asv18S"){
    taxmat_18S <- taxmat
  }
  
  if(primer == "asvrbcL"){
    taxmat_rbcl <- taxmat
  }
}

for(i in 1:nrow(hubscores)){
  otu_id <- hubscores$OTU[i]
  primerset <- hubscores$Primer[i]
  taxalevel <- hubscores$Level[i]
  
  if(primerset == "asv16S"){
    hubscores$TaxaName[i] <- taxmat_16S[otu_id,taxalevel]

  }
  
  if(primerset == "asv18S"){
    hubscores$TaxaName[i] <- taxmat_18S[otu_id,taxalevel]
  }
  
  if(primerset == "asvrbcL"){
    hubscores$TaxaName[i] <- taxmat_rbcl[otu_id,taxalevel]
  }

    renamed <- rbind(renamed, hubscores[i,])
    
}

write.csv(renamed, "/Users/Mel/Desktop/renamed_network_hubscores.csv")


