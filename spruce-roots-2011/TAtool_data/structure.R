##################################################
## Project: Flakaliden spruce 2011
## Script purpose: Generate fullData for TAtool
## Date: 2000417
## Author: Alonso Serrano alonso.serrano@slu.se
##################################################
##################################################
library(RCurl)
library(dplyr)
#structure
options(stringsAsFactors = FALSE)

projectDir <- "/mnt/picea/projects/spruce/vhurry/14_SpruceRoots_Project/" 

descriptionFile <- "data/intro.html"
description <-  readChar(descriptionFile, file.info(descriptionFile)$size)

#start fullData
fullData <- list()

# main store general information about the experiment

fullData$main$description <- description
fullData$species <- "Picea abies"
fullData$speciesNickname <- "pabies"
fullData$genePattern <- "MA_\\d+g\\d+|G\\w+|C\\w+"

#is next is "" no external source will be called
fullData$externalWeb <-  "Congenie"
fullData$externalURL <-  "http://congenie.org/?genelist=enable&_term="
fullData$enrOptions <- c("go", "kegg", "pfam", "mapman")

# conditions is a quick way to control how subdataset the experiment has, if only one, name it control
fullData$conditions <- c("ND+NE", "ND", "NE")
fullData$profileColors[['ND']] <- c("#009444", "#009444", "#009444")
fullData$profileColors[['NE']] <- c("#BE9230", "#BE9230", "#BE9230")
fullData$profileColors[['ND+NE']] <- c("#009444", "#BE9230", "#00A5FF")
#names(fullData$profileColors) <- fullData$conditions

#prepare metadata
meta <- read.table(file.path(projectDir, "data/samples.tsv"), sep='\t', header=T, stringsAsFactors = F)
#meta <- as.data.frame(lapply(meta, as.character))
meta <- meta[-grep("W34-C4", meta$Sample),]
meta <- meta[,c("Week", "Treatment")]

meta$Treatment <- ifelse(meta$Treatment == 'C', "ND", "NE")
# in the preprocessing, stablish the levels in binary conditions

# meta will store specific metadata for each dataset
# TODO simplify this
fullData$meta[['ND']] <- meta[meta$Treatment=="ND",]
fullData$meta[['NE']] <- meta[meta$Treatment=="NE",]
fullData$meta[['ND+NE']] <- meta

# expData store the expression data for each dataset
# TOOO, simplify in combined datasets
load(file.path(projectDir, "DE/controlVsdData.RData"))
load(file.path(projectDir, "DE/fertilisedVsdData.RData"))
load(file.path(projectDir, "DE/combinedVsdData.RData"))

fullData$expData[['ND']] <- controlData
fullData$expData[['NE']] <- fertilisedData
fullData$expData[['ND+NE']] <- combinedData

load(file.path(projectDir, "spruce-networks/control/cluster/controlEnr.RData"))
load(file.path(projectDir, "spruce-networks/fertilised/cluster/fertilisedEnr.RData"))
load(file.path(projectDir, "spruce-networks/combined/cluster/combinedEnr.RData"))

# enr will store the enrichment results, better to have it precalculated for faster access
fullData$enr <- list() # initializing list prevents "more elements supplied than there are to replace" error
fullData$enr[['ND+NE']] <- combinedEnr
fullData$enr[['ND']] <- controlEnr
fullData$enr[['NE']] <- fertilisedEnr

names(fullData$enr[['ND+NE']]) <- as.numeric(gsub("Cluster", "", names(fullData$enr[['ND+NE']])))
names(fullData$enr[['ND']]) <- as.numeric(gsub("Cluster", "", names(fullData$enr[['ND']])))
names(fullData$enr[['NE']]) <- as.numeric(gsub("Cluster", "", names(fullData$enr[['NE']])))

for (x in names(fullData$enr$ND) ){
  fullData$enr[['ND']][[x]]<- fullData$enr[['ND']][[x]][c("go","mapman", "kegg", "pfam")]
}

for (x in names(fullData$enr$NE) ){
  fullData$enr[['NE']][[x]]<- fullData$enr[['NE']][[x]][c("go","mapman", "kegg", "pfam")]
}

for (x in names(fullData$enr[['ND+NE']] ) ){
  fullData$enr[['ND+NE']] [[x]]<- fullData$enr[['ND+NE']] [[x]][c("go","mapman", "kegg", "pfam")]
}

#fullData$clusters[['ND+NE']]$cluster <- as.numeric(gsub("Cluster", "", fullData$clusters[['ND+NE']]$cluster))



load(file.path(projectDir, "DE/allDEresults.RData"))

fullData$de <- list()
fullData$de[['ND+NE']]$description<- "With each timepoint, NE vs. ND"
fullData$de[['ND+NE']]$results <- combined.res

names(fullData$de$`ND+NE`$results) <- gsub("F", "NE", names(fullData$de$`ND+NE`$results))
names(fullData$de$`ND+NE`$results) <- gsub("C", "ND", names(fullData$de$`ND+NE`$results))
fullData$de[['ND']]$description <- "Nutrient-deficientonly, each time point vs. previous time point"
fullData$de[['ND']]$results <- control.res
names(fullData$de$`ND`$results) <- gsub("C", "ND", names(fullData$de$`ND`$results))

fullData$de[['NE']]$description <- "Nutrient-enriched only, each time point vs. previous time point"
fullData$de[['NE']]$results <- fertilised.res
names(fullData$de$`NE`$results) <- gsub("F", "NE", names(fullData$de$`NE`$results))



#TODO match and reorder rows in expData
#match(rownames(fullData$expData$combined), fullData$meta$combined$Sample)

#background

fullData$background[['ND']] <- read.table(paste0(projectDir, "DE/background.txt"), sep='\t', header=T, stringsAsFactors = F)[,1]
fullData$background[['NE']] <- read.table(paste0(projectDir, "DE/background.txt"), sep='\t', header=T, stringsAsFactors = F)[,1]
fullData$background[['ND+NE']] <- read.table(paste0(projectDir, "DE/background.txt"), sep='\t', header=T, stringsAsFactors = F)[,1]
#cluster
#2 columns, gene and cluster

fullData$clusters[['ND']] <- read.table(paste0(projectDir, "spruce-networks/control/cluster/InfomapClusters.tsv"), sep='\t', header=T, stringsAsFactors = F)
fullData$clusters[['NE']] <- read.table(paste0(projectDir, "spruce-networks/fertilised/cluster/InfomapClusters.tsv"), sep='\t', header=T, stringsAsFactors = F)
fullData$clusters[['ND+NE']] <- read.table(paste0(projectDir, "spruce-networks/combined/cluster/InfomapClusters.tsv"), sep='\t', header=T, stringsAsFactors = F)

fullData$clusters[['ND']]$cluster <- as.numeric(gsub("Cluster", "", fullData$clusters[['ND']]$cluster))
fullData$clusters[['NE']]$cluster <- as.numeric(gsub("Cluster", "", fullData$clusters[['NE']]$cluster))
fullData$clusters[['ND+NE']]$cluster <- as.numeric(gsub("Cluster", "", fullData$clusters[['ND+NE']]$cluster))


# sample points, time, stage, % humidity, these will be represented in the x axis
fullData$samplePoints[['ND']] <- fullData$meta[['ND']]$Week
fullData$samplePoints[['NE']] <- fullData$meta[['NE']]$Week
fullData$samplePoints[['ND+NE']] <- fullData$meta[['ND+NE']]$Week

# variables, tissues, mutants, control vs fertilised. Each one of this will have a line in the plots
fullData$variables[['ND']] <- fullData$meta[['ND']]$Treatment
fullData$variables[['NE']] <- fullData$meta[['NE']]$Treatment
fullData$variables[['ND+NE']] <- fullData$meta[['ND+NE']]$Treatment

#process images
fullData$networks[['ND+NE']]$name <- "ND+NE network"
fullData$networks[['ND+NE']]$description <- "ND+NE network generated with Seidr, clustering with Infomap, visualization with Cytoscape"
fullData$networks[['ND+NE']]$image <- sprintf('<img height="720" width="1280" src="data:image/png;base64,%s">', 
                                           base64Encode(readBin(("data/spruce_combined_network.png"), "raw", file.info(("data/spruce_combined_network.png"))[1, "size"]), "txt"))
fullData$networks[['ND+NE']]$stats <-read.table(file.path(projectDir, "spruce-networks/combined/aggregated", "seidrStatsReduced.tsv"), header=T, sep='\t')
rownames(fullData$networks[['ND+NE']]$stats) <- fullData$networks[['ND+NE']]$stats[,1]
fullData$networks[['ND+NE']]$stats <- fullData$networks[['ND+NE']]$stats[,-1]
fullData$networks[['ND+NE']]$edgeList <- read.table(file.path(projectDir, "spruce-networks/combined/aggregated", "combinedEdgeList.tsv"), header=T, sep='\t')


fullData$networks[['ND']]$name <- "ND network"
fullData$networks[['ND']]$description <- "ND network generated with Seidr, clustering with Infomap, visualization with Cytoscape"
fullData$networks[['ND']]$image <- sprintf('<img height="720" width="1280" src="data:image/png;base64,%s">', 
                                           base64Encode(readBin(("data/spruce_control_network.png"), "raw", file.info(("data/spruce_control_network.png"))[1, "size"]), "txt"))
fullData$networks[['ND']]$stats <-read.table(file.path(projectDir, "spruce-networks/control/aggregated", "seidrStatsReduced.tsv"), header=T, sep='\t')
rownames(fullData$networks[['ND']]$stats) <- fullData$networks[['ND']]$stats[,1]
fullData$networks[['ND']]$stats <- fullData$networks[['ND']]$stats[,-1]
fullData$networks[['ND']]$edgeList <- read.table(file.path(projectDir, "spruce-networks/control/aggregated", "controlEdgeList.tsv"), header=T, sep='\t')


fullData$networks[['NE']]$name <- "NE network"
fullData$networks[['NE']]$description <- "NE network generated with Seidr, clustering with Infomap, visualization with Cytoscape"
fullData$networks[['NE']]$image <- sprintf('<img height="720" width="1280" src="data:image/png;base64,%s">', 
                                           base64Encode(readBin(("data/spruce_fertilized_network.png"), "raw", file.info(("data/spruce_fertilized_network.png"))[1, "size"]), "txt"))
fullData$networks[['NE']]$stats <-read.table(file.path(projectDir, "spruce-networks/fertilised/aggregated", "seidrStatsReduced.tsv"), header=T, sep='\t')
rownames(fullData$networks[['NE']]$stats) <- fullData$networks[['NE']]$stats[,1]
fullData$networks[['NE']]$stats <- fullData$networks[['NE']]$stats[,-1]
fullData$networks[['NE']]$edgeList <- read.table(file.path(projectDir, "spruce-networks/fertilised/aggregated", "fertilisedEdgeList.tsv"), header=T, sep='\t')


# temp <- read.table(paste0(projectDir, "spruce-networks/control/aggregated/stats.txt"), header=T, sep='\t')
# temp2 <- temp %>% select(Target, PageRank_target, Betweenness_target, Strength_target, Eigenvector_target, Katz_target)
# temp3 <- temp2 %>% distinct()
# rownames(temp3) <- temp3$Target
# temp3 <- temp3[ , !(names(temp3) %in% "Target")]
# fullData$networks[['control']]$stats <- temp3


# Save data
save(fullData, file= "data/fullData.RData")
