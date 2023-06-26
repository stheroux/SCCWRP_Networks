setwd("/project/noujdine_61/mgosborn/SCCWRP")
library(dplyr)
library(xlsx)
library(tidyverse)
library(phyloseq)
library(SpiecEasi)
library(igraph)


## This script runs through network analysis for SCCWRP water stream samples. 

###################
###
### Loading metadata and abundance tables. 
### Preparing phyloseq object.
###
###################

#Load metadata
metadata <- read.xlsx("env8_forMelissa.xlsx", sheetIndex = 1, row.names = 1)

#Load abundance table
asv16S <- read.csv("16S_2022.csv", stringsAsFactors=FALSE, row.names = 1)#, check.names = FALSE)
asv18S <- read.csv("18S_2022.csv", stringsAsFactors=FALSE, row.names = 1)
asvrbcL <- read.csv("rbcL_2022.csv", stringsAsFactors=FALSE, row.names = 1)
#asvCO1 <- read.csv("CO1_blast_tax_table.csv", stringsAsFactors=FALSE, row.names = 1)


## Initialize dataframes that will be filled in later
pos_and_neg_edges = data.frame()
hubscores = data.frame()
numsamplesandtaxa = data.frame()

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
  taxmat <- as.matrix(taxmat) ### NOTE: Some ASVs/OTUs appear to be misclassified (i.e., classified as Cyanobacteria across all taxonomic levels, Genus classification incorrectly listed in Species column. Consider filtering out incorrect classifications)
  otumat <- subset(otumat, select = -c(ConsensusLineage))
  #remove samples that aren't present in metadata
  otumat <- otumat[, (colnames(otumat) %in% row.names(metadata) == TRUE)] 
  
  #Create phyloseq object
  OTU = otu_table(otumat, taxa_are_rows = TRUE)
  TAX = tax_table(taxmat)
  physeq = phyloseq(OTU, TAX)
  
  #Add metadata info to physeq object
  physeq@sam_data = sample_data(metadata)
  
  ###################
  ###
  ### Process phyloseq object (i.e., filtering samples and taxa)
  ###
  ###################
  
  #Remove "replicate" samples (i.e., 'eDNA', 'MB', and 'LB')
  physeq = subset_samples(physeq, Replicate != c("eDNA", "MB", "LB"))
  #Remove singletons and doubletons
  physeq <- prune_taxa(taxa_sums(physeq) > 2, physeq) 
  #Remove samples with <2000 reads
  physeq <- prune_samples(sample_sums(physeq) > 2000, physeq)
  
  ### NOTE: Did not normalize or apply alternate filter on abundance table. Consider adding here. 
  
  #Only keep assigned taxa
  physeq = subset_taxa(physeq, Domain!="Unassigned")
  
  ###################
  ###
  ### Create separate phyloseq objects for difference reference statuses. 
  ###
  ###################
  
  levels <- c("Family","Genus","Species")
  
  for (level in levels){
    #Conglomerate to particular taxonomic level (e.g., Order, Family... etc.)
    physeq_glom = tax_glom(physeq, level)
    #Get rid of ASVs/OTUs that are absent from all samples
    physeq_glom = prune_taxa(taxa_sums(physeq_glom) > 0, physeq_glom)
    ### NOTE: if there are unclassified or ambiguous taxa at a particular level (such as Species), consider adding a filter here to remove

    if (level == "Species"){
      # Only include taxa >0.1% relative abundance
      #print("Num taxa: ")
      #print(ntaxa(physeq_glom))
      rel_glom <- transform_sample_counts(physeq_glom, function(x) x / sum(x) )
      rel_glom = filter_taxa(rel_glom, function(x) sum(x) > 0.1, TRUE)
      physeq_glom = prune_taxa(taxa_names(physeq_glom) %in% taxa_names(rel_glom), physeq_glom)
      #print("Num taxa: ")
      #print(ntaxa(physeq_glom))
    }
    #print("double check, num taxa: ")
    #print(ntaxa(physeq_glom))

    
    #Create separate phyloseq object for each reference status (i.e., Reference, Intermediate, and Stressed)
    ref_statuses = c("physeq_Reference", "physeq_Intermediate", "physeq_Stressed")
    physeq_glom <- subset_samples(physeq_glom, RefStatus != "NA")
    physeq_Reference <- subset_samples(physeq_glom, RefStatus == "Reference")
    physeq_Intermediate <- subset_samples(physeq_glom, RefStatus == "Intermediate")
    physeq_Stressed <- subset_samples(physeq_glom, RefStatus == "Stressed")
    
    ###################
    ###
    ### Record number of samples and number of taxa 
    ###
    ###################
    
    ## primer, taxonomic level, num taxa, num all samples, num reference samples, num intermediate samples, num stressed samples
    samtax = c(primer, level, ntaxa(physeq_glom), nsamples(physeq_glom), nsamples(physeq_Reference), nsamples(physeq_Intermediate), nsamples(physeq_Stressed))
    numsamplesandtaxa <- rbind(samtax, numsamplesandtaxa)

    ###################
    ###
    ### SpiecEasi
    ###
    ###################

    for(ref_status in ref_statuses){
      temp = get(ref_status)
      pargs <- list(rep.num=50, ncores=36)
      se.mb <- spiec.easi(temp, method='mb', nlambda = 100, sel.criterion='bstars', pulsar.select=TRUE, pulsar.params=pargs)
      print(paste("Primer: ", primer, ", Reference Status:", ref_status, "& Taxonomic Level: ", level, sep = " "))
      print(getStability(se.mb))
      optbeta <- as.matrix(symBeta(getOptBeta(se.mb)))
      edge_cols <-  ifelse(optbeta>0, 'rgb(0.3,0.8,0,', 'rgb(1,0.25,1,')[upper.tri(optbeta) & optbeta!=0]

      optbeta <- symBeta(getOptBeta(se.mb))
      edge_weights <- abs((Matrix::summary(t(optbeta))[,3]))
      edge_weights <- edge_weights/max(edge_weights)
      ig2.mb <- adj2igraph(getRefit(se.mb),  rmEmptyNodes=TRUE,
                           vertex.attr=list(name=taxa_names(temp)),
                           edge.attr=list(color=edge_cols, curved = 0, weight = edge_weights))

      TotalNodes = gorder(ig2.mb)
      TotalEdges = gsize(ig2.mb)
      PositiveEdges = sum(optbeta@x>0)
      NegativeEdges = sum(optbeta@x<0)
      PropPosNeg = PositiveEdges/NegativeEdges
      PosTotal = PositiveEdges/TotalEdges
      NegTotal = NegativeEdges/TotalEdges
      AvgPathLength = average.path.length(ig2.mb) #Shortest paths between vertices
      Cluster = cluster_louvain(ig2.mb) #Louvain method is an algorithm to detect communities in large networks. It maximizes a modularity score for each community, where the modularity quantifies the quality of an assignment of nodes to communities.This means evaluating how much more densely connected the nodes within a community are, compared to how connected they would be in a random network.
      Modularity = modularity(Cluster) #Strength of division of network into modules
      DegreeDist = degree(ig2.mb)
      AvgDegree = mean(DegreeDist)
      ClusterCoeff = transitivity(ig2.mb) #Transitivity ("clustering coefficient") measures the probability that the adjacent vertices of a vertex are connected.

      #degree heterogeneity
      temp2 = 0
      for(i in 1:length(table(DegreeDist))){
        temp = sum((1-(table(DegreeDist)[i]/sum(table(DegreeDist))))^2)
        temp2 = temp2+temp
      }
      HetIndex = sqrt((1/TotalNodes)*temp2)
      MaxHetInt = sqrt(1 - (3/TotalNodes) + ((TotalNodes+2)/(TotalNodes^3)))
      HetMeasure = HetIndex/MaxHetInt

      temp <- c(primer, level, ref_status, "All Samples", TotalNodes, TotalEdges, PositiveEdges, NegativeEdges, PropPosNeg, PosTotal, NegTotal, AvgPathLength, Modularity, AvgDegree, HetMeasure, ClusterCoeff)
      pos_and_neg_edges <- rbind(pos_and_neg_edges, temp)
      colnames(pos_and_neg_edges) <- c("Primer","Taxonomic Level", "Reference Status", "Num Subsampled", "Total Nodes","Total Edges", "Positive Edges", "Negative Edges", "Pos/Neg","Pos/Total","Neg/Total", "Avg Path Length", "Modularity", "Avg Degree", "Heterogeneity", "Clustering Coefficient")

      ###HUB SCORES
      temp <- as.data.frame(head(sort(hub_score(ig2.mb)$vector, decreasing = TRUE), 10))
      colnames(temp)[1] <-  "Hub Score"
      temp <- tibble::rownames_to_column(temp, "OTU")
      temp$ReferenceStatus = ref_status
      temp$Level = level
      temp$Primer = primer
      hubscores <- rbind(hubscores, temp)

      for(i in 1:100){
        #randomly select 50 samples (100 times)
        temp = get(ref_status)
        temp <- prune_samples(sample(sample_names(temp), 50), temp)
        pargs <- list(rep.num=50, ncores=36)
        se.mb <- spiec.easi(temp, method='mb', nlambda = 100,sel.criterion='bstars', pulsar.select=TRUE,pulsar.params=pargs)
        #print(Sys.time())
        print(paste("Primer: ", primer, ", Reference Status: ", ref_status, ", Taxonomic Level: ", level, ", Iteration: ", i, " of 100", sep = ""))
        print(getStability(se.mb))
        optbeta <- as.matrix(symBeta(getOptBeta(se.mb)))
        edge_cols <-  ifelse(optbeta>0, 'rgb(0.3,0.8,0,', 'rgb(1,0.25,1,')[upper.tri(optbeta) & optbeta!=0]

        optbeta <- symBeta(getOptBeta(se.mb))
        edge_weights <- abs((Matrix::summary(t(optbeta))[,3]))
        edge_weights <- edge_weights/max(edge_weights)
        ig2.mb <- adj2igraph(getRefit(se.mb),  rmEmptyNodes=TRUE,
                             vertex.attr=list(name=taxa_names(temp)),
                             edge.attr=list(color=edge_cols, curved = 0, weight = edge_weights))

        TotalNodes = gorder(ig2.mb)
        TotalEdges = gsize(ig2.mb)
        PositiveEdges = sum(optbeta@x>0)
        NegativeEdges = sum(optbeta@x<0)
        PropPosNeg = PositiveEdges/NegativeEdges
        PosTotal = PositiveEdges/TotalEdges
        NegTotal = NegativeEdges/TotalEdges
        AvgPathLength = average.path.length(ig2.mb) #Shortest paths between vertices
        Cluster = cluster_louvain(ig2.mb) #Louvain method is an algorithm to detect communities in large networks. It maximizes a modularity score for each community, where the modularity quantifies the quality of an assignment of nodes to communities.This means evaluating how much more densely connected the nodes within a community are, compared to how connected they would be in a random network.
        Modularity = modularity(Cluster) #Strength of division of network into modules
        DegreeDist = degree(ig2.mb)
        AvgDegree = mean(DegreeDist)
        ClusterCoeff = transitivity(ig2.mb) #Transitivity ("clustering coefficient") measures the probability that the adjacent vertices of a vertex are connected.

        #degree heterogeneity
        temp2 = 0
        for(i in 1:length(table(DegreeDist))){
          temp = sum((1-(table(DegreeDist)[i]/sum(table(DegreeDist))))^2)
          temp2 = temp2+temp
        }
        HetIndex = sqrt((1/TotalNodes)*temp2)
        MaxHetInt = sqrt(1 - (3/TotalNodes) + ((TotalNodes+2)/(TotalNodes^3)))
        HetMeasure = HetIndex/MaxHetInt

        temp <- c(primer, level, ref_status, "50", TotalNodes, TotalEdges, PositiveEdges, NegativeEdges, PropPosNeg, PosTotal, NegTotal, AvgPathLength, Modularity, AvgDegree, HetMeasure, ClusterCoeff)
        pos_and_neg_edges <- rbind(pos_and_neg_edges, temp)
        colnames(pos_and_neg_edges) <- c("Primer","Taxonomic Level", "Reference Status", "Num Subsampled", "Total Nodes","Total Edges", "Positive Edges", "Negative Edges", "Pos/Neg","Pos/Total","Neg/Total", "Avg Path Length", "Modularity", "Avg Degree", "Heterogeneity", "Clustering Coefficient")

      }

      ###################
      ###
      ### Network graphs from SpiecEasi analysis
      ###
      ###################

      df <- igraph::as_data_frame(ig2.mb, 'both')
      check <- as.data.frame(taxmat[,"Order"])
      colnames(check) <- "Order"
      check$id <- row.names(check)

      df$vertices <- df$vertices %>%
        left_join(check, c('name'='id'))

      updated_g <- graph_from_data_frame(df$edges,
                                         directed = F,
                                         vertices = df$vertices)
      updated_g$layout <- layout_with_dh

      coul <- scales::viridis_pal(option = "turbo")(length(unique(tax_table(get(ref_status))[,4])))

      edge_function <- paste(edge_cols,edge_weights,')',sep = "")
      edge_col_val <- c()
      for(i in 1:length(edge_function)){
        edge_col_val[i] <- eval(parse(text=edge_function[i]))
      }

      plot(updated_g, edge.color = edge_col_val, edge.curved=.2, vertex.size = ((hub_score(ig2.mb)$vector)*10)+1, vertex.label = NA, vertex.color = coul[as.numeric(as.factor(vertex_attr(updated_g, "Order")))],vertex.frame.color = "black")

      legend("left", title = "Order", legend=levels(as.factor(V(updated_g)$Order))  , col = coul , bty = "n", pch=20 , pt.cex = 1, cex = 1 , horiz = FALSE, inset = c(-0.3, 0))
      sizeCut<- c(0.2,0.4,0.6,0.8,1.0)
      sizeCutScale <- sizeCut*10+1
      a <- legend('right',title = "Hub Score", legend=unique(sizeCut), pt.cex= sizeCutScale, inset = c(-0, 0), bty = "n", y.intersp=1.1)
      x <- (a$text$x + a$rect$left) / 2
      y <- a$text$y
      symbols(x,y,circles=sizeCutScale/200,inches=FALSE,add=TRUE,bg='black')

      if(ref_status == "physeq_Reference"){
        network_physeq_Reference <- recordPlot()
      }

      if(ref_status == "physeq_Intermediate"){
        network_physeq_Intermediate <- recordPlot()
      }

      if(ref_status == "physeq_Stressed"){
        network_physeq_Stressed <- recordPlot()
      }


    } # for(ref_status in ref_statuses)
  } # for (level in levels)
}


# write.csv(pos_and_neg_edges, "/Users/Mel/Desktop/SCCWRP/network_topology.csv")
# write.csv(hubscores, "/Users/Mel/Desktop/SCCWRP/network_hubscores.csv")     
colnames(numsamplesandtaxa) <- c("Primer","Taxa Level","Num Taxa","Num All Samples", "Num Reference Samples","Num Intermediate Samples","Num Stressed Samples")
write.csv(numsamplesandtaxa, "/project/noujdine_61/mgosborn/SCCWRP/num_samples_vs_taxa.csv")
write.csv(pos_and_neg_edges, "/project/noujdine_61/mgosborn/SCCWRP/network_topology_species.csv")
write.csv(hubscores, "/project/noujdine_61/mgosborn/SCCWRP/network_hubscores_species.csv")     


