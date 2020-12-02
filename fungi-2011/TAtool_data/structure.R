##################################################
## Project: Flakaliden fungi 2011
## Script purpose: Generate fullData for TAtool
## Date: 2000417
## Author: Alonso Serrano alonso.serrano@slu.se
##################################################

library(RCurl)
library(dplyr)

#structure
options(stringsAsFactors = FALSE)

projectDir <- "/mnt/picea/projects/spruce/vhurry/fungi-flakaliden-2011/" 
deFolder <- file.path(projectDir,"DE")
dataFolder <- file.path(projectDir,"data")

descriptionFile <- "data/intro.html"
description <-  readChar(descriptionFile, file.info(descriptionFile)$size)
#description <- str_replace_all(description, "(\\W)", "\\\\\\1")

#start fullData
fullData <- list()
fullData$species <- "2011Fungi"
fullData$speciesNickname <- "ko"
fullData$genePattern <- "K\\d+"

#is next is "" no external source will be called
fullData$externalWeb <-  ""
fullData$externalURL <-  ""
fullData$enrOptions <- c("ko_pathway")
# main store general information about the experiment

fullData$main$enrDate <- as.Date("19-11-22") # enrichment date of generation


# conditions is a quick way to control how subdataset the experiment has, if only one, name it control
fullData$conditions <- c("ND+NE", "ND", "NE")
fullData$profileColors[['ND']] <- c("#009444", "#009444", "#009444")
fullData$profileColors[['NE']] <- c("#BE9230", "#BE9230", "#BE9230")
fullData$profileColors[['ND+NE']] <- c("#009444", "#BE9230", "#00A5FF")


#prepare metadata
meta <- read.table(file.path(dataFolder, "meta.tsv"), header = TRUE, sep = '\t', comment.char = "", stringsAsFactors = TRUE)
#meta$Treatment <- ifelse(meta$Condition == 'control', 'C', 'F')
meta$Treatment <- ifelse(meta$Condition == 'control', "ND", "NE")
meta$Week <- meta$Sampling.date..week..
#meta$group <- paste0(meta$Week, gsub("T","F",meta$Plot))
#meta <- meta[order(meta$group), ]
meta$Treatment <- ifelse(meta$Condition == 'control', "ND", "NE")
# meta$group <- substr(meta$group, 1, 4)
# meta$group <- as.factor(meta$group)
# meta$Treatment <- as.factor(meta$Treatment)
# meta$Condition <- ifelse(meta$Condition == 'control', 'control', 'fertilised')

# in the preprocessing, stablish the levels in binary conditions

# meta will store specific metadata for each dataset

fullData$meta[['ND']] <- meta[meta$Treatment=="ND",][c("Treatment", "Week")]
fullData$meta[['NE']] <- meta[meta$Treatment=="NE",][c("Treatment", "Week")]
fullData$meta[['ND+NE']] <- meta[c("Treatment", "Week")]

# expData store the expression data for each dataset
# TOOO, simplify in combined datasets
load(file.path(deFolder, "controlVsdData.RData"))
load(file.path(deFolder, "fertilisedVsdData.RData"))
load(file.path(deFolder, "combinedVsdData.RData"))

fullData$expData[['ND+NE']] <- combinedData
fullData$expData[['ND']] <- controlData
fullData$expData[['NE']] <- fertilisedData


load(file.path(projectDir, "networks/control/cluster/controlEnr.RData"))
load(file.path(projectDir, "networks/fertilised/cluster/fertilisedEnr.RData"))
load(file.path(projectDir, "networks/combined/cluster/combinedEnr.RData"))

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


load(file.path(deFolder, "allDEresults.RData"))

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

fullData$background[['ND']] <- read.table(file.path(deFolder, "background.txt"), sep='\t', header=T, stringsAsFactors = F)[,1]
fullData$background[['NE']] <- read.table(file.path(deFolder, "background.txt"), sep='\t', header=T, stringsAsFactors = F)[,1]
fullData$background[['ND+NE']] <- read.table(file.path(deFolder, "background.txt"), sep='\t', header=T, stringsAsFactors = F)[,1]
#cluster
#2 columns, gene and cluster

fullData$clusters[['ND']] <- read.table(file.path(projectDir, "networks/control/cluster/InfomapClusters.tsv"), sep='\t', header=T, stringsAsFactors = F)
fullData$clusters[['NE']] <- read.table(file.path(projectDir, "networks/fertilised/cluster/InfomapClusters.tsv"), sep='\t', header=T, stringsAsFactors = F)
fullData$clusters[['ND+NE']] <- read.table(file.path(projectDir, "networks/combined/cluster/InfomapClusters.tsv"), sep='\t', header=T, stringsAsFactors = F)

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
                                              base64Encode(readBin(("data/fungi_combined_network.png"), "raw", file.info(("data/fungi_combined_network.png"))[1, "size"]), "txt"))
# fullData$networks[['ND+NE']]$stats <-read.table(file.path(projectDir, "spruce-networks/combined/aggregated", "seidrStatsReduced.tsv"), header=T, sep='\t')
# rownames(fullData$networks[['ND+NE']]$stats) <- fullData$networks[['ND+NE']]$stats[,1]
# fullData$networks[['ND+NE']]$stats <- fullData$networks[['ND+NE']]$stats[,-1]
fullData$networks[['ND+NE']]$edgeList <- read.table(file.path(projectDir, "networks/combined/cluster", "edgeList.txt"), header=T, sep='\t')


fullData$networks[['ND']]$name <- "ND network"
fullData$networks[['ND']]$description <- "ND network generated with Seidr, clustering with Infomap, visualization with Cytoscape"
fullData$networks[['ND']]$image <- sprintf('<img height="720" width="1280" src="data:image/png;base64,%s">', 
                                           base64Encode(readBin(("data/fungi_control_network.png"), "raw", file.info(("data/fungi_control_network.png"))[1, "size"]), "txt"))
# fullData$networks[['ND']]$stats <-read.table(file.path(projectDir, "spruce-networks/control/aggregated", "seidrStatsReduced.tsv"), header=T, sep='\t')
# rownames(fullData$networks[['ND']]$stats) <- fullData$networks[['ND']]$stats[,1]
# fullData$networks[['ND']]$stats <- fullData$networks[['ND']]$stats[,-1]
fullData$networks[['ND']]$edgeList <- read.table(file.path(projectDir, "networks/control/cluster", "edgeList.txt"), header=T, sep='\t')


fullData$networks[['NE']]$name <- "NE network"
fullData$networks[['NE']]$description <- "NE network generated with Seidr, clustering with Infomap, visualization with Cytoscape"
fullData$networks[['NE']]$image <- sprintf('<img height="720" width="1280" src="data:image/png;base64,%s">', 
                                           base64Encode(readBin(("data/fungi_fertilized_network.png"), "raw", file.info(("data/fungi_fertilized_network.png"))[1, "size"]), "txt"))
#fullData$networks[['NE']]$stats <-read.table(file.path(projectDir, "spruce-networks/fertilised/aggregated", "seidrStatsReduced.tsv"), header=T, sep='\t')
#rownames(fullData$networks[['NE']]$stats) <- fullData$networks[['NE']]$stats[,1]
#fullData$networks[['NE']]$stats <- fullData$networks[['NE']]$stats[,-1]
fullData$networks[['NE']]$edgeList <- read.table(file.path(projectDir, "networks/fertilised/cluster", "edgeList.txt"), header=T, sep='\t')


# Save data
save(fullData, file= file.path("data/fullData.RData"))


# custom data 1

load(file.path(dataFolder, "megaData.RData"))
customData1 <- list()
customData1$taxa <- c("superkingdom","kingdom","phylum","class","order","family","genus","species")

customData1$megaData <- megaData[c(1:120,124,125)]
# object.size(customData1)
#as.factor(customData1$megaData[113:122])

# customData1$megaData$superkingdom <- as.factor(customData1$megaData$species)
# customData1$megaData$kingdom <- as.factor(customData1$megaData$kingdom)
# customData1$megaData$phylum <- as.factor(customData1$megaData$phylum)
# customData1$megaData$class <- as.factor(customData1$megaData$class)
# customData1$megaData$order <- as.factor(customData1$megaData$order)
# customData1$megaData$family <- as.factor(customData1$megaData$family)
# customData1$megaData$genus <- as.factor(customData1$megaData$genus)
# customData1$megaData$species <- as.factor(customData1$megaData$species)
# customData1$megaData$Trophic.Mode <- as.factor(customData1$megaData$Trophic.Mode)
# customData1$megaData$Guild <- as.factor(customData1$megaData$Guild)
# object.size(customData1)

customData1$profileColors[['ND']] <- c("#009444", "#009444", "#009444")
customData1$profileColors[['NE']] <- c("#BE9230", "#BE9230", "#BE9230")
customData1$profileColors[['ND+NE']] <- c("#009444", "#BE9230", "#00A5FF")

customData1$rawData[['ND']] <- customData1$megaData[grep("W\\d+C",colnames(customData1$megaData))]
customData1$rawData[['NE']] <- megaData[grep("W\\d+F",colnames(megaData))]
customData1$rawData[['ND+NE']] <- megaData[grep("W\\d+[CF]",colnames(megaData))]

customData1$meta[['ND']] <- data.frame(Treatment=rep("ND", length(colnames(customData1$rawData[['ND']]))),
                                            Week = substr(colnames(customData1$rawData[['ND']]),2,3))

customData1$meta[['NE']] <- data.frame(Treatment=rep("NE", length(colnames(customData1$rawData[['NE']]))),
                                               Week = substr(colnames(customData1$rawData[['NE']]),2,3))

customData1$meta[['ND+NE']] <- data.frame(Treatment=substr(colnames(customData1$rawData[['ND+NE']]),4,4),
                                             Week = substr(colnames(customData1$rawData[['ND+NE']]),2,3))
customData1$meta[['ND+NE']]$Treatment <- ifelse(customData1$meta[['ND+NE']]$Treatment == "C", "ND", "NE")

save(fullData, customData1, file= file.path("data/fullData.RData"))