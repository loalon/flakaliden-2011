##################################################
## Project: Flakaliden fungi 2011
## Script purpose: Generate fullData for TAtool
## Date: 2000417
## Author: Alonso Serrano alonso.serrano@slu.se
##################################################

library(here)
library(RCurl)
library(dplyr)

#structure
options(stringsAsFactors = FALSE)

projectDir <- "/mnt/picea/projects/spruce/vhurry/fungi-flakaliden-2011/" 
deFolder <- file.path(projectDir,"DE")
dataFolder <- file.path(projectDir,"data")

descriptionFile <- here("data/intro.html")
description <-  readChar(descriptionFile, file.info(descriptionFile)$size)
#description <- str_replace_all(description, "(\\W)", "\\\\\\1")

#start fullData
fullData <- list()
speciesMap <- read.table(here("config/speciesMap.tsv"),sep='\t', stringsAsFactors = F, header=T, quote = "", allowEscapes = T )
# main store general information about the experiment
fullData$main$description <- description
fullData$species <- "2011Fungi"
speciesRow <- speciesMap[speciesMap$species == fullData$species,]
fullData$main$enrDate <- as.Date("19-11-22") # enrichment date of generation

fullData$speciesNickname <- speciesRow$nickname
fullData$genePattern <- speciesRow$genePattern
fullData$enrOptions <- unlist(strsplit(speciesRow$enrichment,','))
#is next is "" no external source will be called


# conditions is a quick way to control how subdataset the experiment has, if only one, name it control
fullData$conditions <- c("control", "fertilised", "combined")
fullData$profileColors[['control']] <- c("#009444", "#009444", "#009444")
fullData$profileColors[['fertilised']] <- c("#BE9230", "#BE9230", "#BE9230")
fullData$profileColors[['combined']] <- c("#00A5FF", "#00A5FF", "#00A5FF")
#names(fullData$profileColors) <- fullData$conditions

#prepare metadata
meta <- read.table(file.path(dataFolder, "meta.tsv"), header = TRUE, sep = '\t', comment.char = "", stringsAsFactors = TRUE)
meta$Treatment <- ifelse(meta$Condition == 'control', 'C', 'F')
meta$Week <- meta$Sampling.date..week..
meta$group <- paste0(meta$Week, gsub("T","F",meta$Plot))
meta <- meta[order(meta$group), ]
meta$group <- substr(meta$group, 1, 4)
meta$group <- as.factor(meta$group)
meta$Treatment <- as.factor(meta$Treatment)
meta$Condition <- ifelse(meta$Condition == 'control', 'control', 'fertilised')

# in the preprocessing, stablish the levels in binary conditions

# meta will store specific metadata for each dataset
# TODO simplify this
fullData$meta[['control']] <- meta[meta$Condition =="control",][c("Condition", "Week")]
fullData$meta[['fertilised']] <- meta[meta$Condition =="fertilised",][c("Condition", "Week")]
fullData$meta[['combined']] <- meta[c("Condition", "Week")]

# expData store the expression data for each dataset
# TOOO, simplify in combined datasets
load(file.path(deFolder, "controlVsdData.RData"))
load(file.path(deFolder, "fertilisedVsdData.RData"))
load(file.path(deFolder, "combinedVsdData.RData"))

fullData$expData[['control']] <- controlData
fullData$expData[['fertilised']] <- fertilisedData
fullData$expData[['combined']] <- combinedData

load(file.path(projectDir, "networks/control/cluster/controlEnr.RData"))
load(file.path(projectDir, "networks/fertilised/cluster/fertilisedEnr.RData"))
load(file.path(projectDir, "networks/combined/cluster/combinedEnr.RData"))

# enr will store the enrichment results, better to have it precalculated for faster access
fullData$enr <- list() # initializing list prevents "more elements supplied than there are to replace" error
fullData$enr[['control']] <- controlEnr
fullData$enr[['fertilised']] <- fertilisedEnr
fullData$enr[['combined']] <- combinedEnr

load(file.path(deFolder, "allDEresults.RData"))

fullData$de <- list()
fullData$de[['control']]$description <- "Control only, each time point vs. previous time point"
fullData$de[['control']]$results <- control.res
fullData$de[['fertilised']]$description <- "fertilised only, each time point vs. previous time point"
fullData$de[['fertilised']]$results <- fertilised.res
fullData$de[['combined']]$description<- "Within each timepoint, Fertilised vs. Control"
fullData$de[['combined']]$results <- combined.res

#TODO match and reorder rows in expData
#match(rownames(fullData$expData$combined), fullData$meta$combined$Sample)

#background

fullData$background[['control']] <- read.table(file.path(deFolder, "background.txt"), sep='\t', header=T, stringsAsFactors = F)[,1]
fullData$background[['fertilised']] <- read.table(file.path(deFolder, "background.txt"), sep='\t', header=T, stringsAsFactors = F)[,1]
fullData$background[['combined']] <- read.table(file.path(deFolder, "background.txt"), sep='\t', header=T, stringsAsFactors = F)[,1]
#cluster
#2 columns, gene and cluster

fullData$clusters[['control']] <- read.table(file.path(projectDir, "networks/control/cluster/InfomapClusters.tsv"), sep='\t', header=T, stringsAsFactors = F)
fullData$clusters[['fertilised']] <- read.table(file.path(projectDir, "networks/fertilised/cluster/InfomapClusters.tsv"), sep='\t', header=T, stringsAsFactors = F)
fullData$clusters[['combined']] <- read.table(file.path(projectDir, "networks/combined/cluster/InfomapClusters.tsv"), sep='\t', header=T, stringsAsFactors = F)

# sample points, time, stage, % humidity, these will be represented in the x axis
fullData$samplePoints[['control']] <- fullData$meta$control$Week
fullData$samplePoints[['fertilised']] <- fullData$meta$fertilised$Week
fullData$samplePoints[['combined']] <- fullData$meta$combined$Week

# variables, tissues, mutants, control vs fertilised. Each one of this will have a line in the plots
fullData$variables[['control']] <- fullData$meta$control$Condition
fullData$variables[['fertilised']] <- fullData$meta$fertilised$Condition
fullData$variables[['combined']] <- fullData$meta$combined$Condition

#process images
#txt <- base64Encode(readBin(here("toydata/control.png"), "raw", file.info(here("toydata/control.png"))[1, "size"]), "txt")

fullData$networks[['control']]$name <- "Control network"
fullData$networks[['control']]$description <- "Control network generated with Seidr, visualization with Infomap"
#fullData$networks[['control']]$image <- sprintf('<img height="720" width="1280" src="data:image/png;base64,%s">', txt)

# temp <- read.table(paste0(projectDir, "spruce-networks/control/aggregated/stats.txt"), header=T, sep='\t')
# temp2 <- temp %>% select(Target, PageRank_target, Betweenness_target, Strength_target, Eigenvector_target, Katz_target)
# temp3 <- temp2 %>% distinct()
# rownames(temp3) <- temp3$Target
# temp3 <- temp3[ , !(names(temp3) %in% "Target")]
# fullData$networks[['control']]$stats <- temp3
# fullData$networks[['control']]$edgeList <- read.table(file.path(projectDir, "spruce-networks/control/aggregated", "controlEdgeList.tsv"), header=T, sep='\t')

# Save data
save(fullData, file= file.path("/mnt/picea/home/bastian/shiny-apps/flakaliden-2011-fungi", "data/fullData.RData"))


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

customData1$profileColors[['control']] <- c("#009444", "#009444", "#009444")
customData1$profileColors[['fertilised']] <- c("#BE9230", "#BE9230", "#BE9230")
customData1$profileColors[['combined']] <- c("#009444", "#BE9230", "#00A5FF")


# customData1$rawData[['control']] <- megaData[grep("W\\d+C",colnames(megaData))]
# customData1$rawData[['fertilised']] <- megaData[grep("W\\d+F",colnames(megaData))]
# customData1$rawData[['combined']] <- megaData[grep("W\\d+[CF]",colnames(megaData))]

customData1$meta[['control']] <- data.frame(Treatment=rep("control", length(colnames(customData1$rawData[['control']]))),
                                            Week = substr(colnames(customData1$rawData[['control']]),2,3))

customData1$meta[['fertilised']] <- data.frame(Treatment=rep("fertilised", length(colnames(customData1$rawData[['fertilised']]))),
                                               Week = substr(colnames(customData1$rawData[['fertilised']]),2,3))

customData1$meta[['combined']] <- data.frame(Treatment=substr(colnames(customData1$rawData[['combined']]),4,4),
                                             Week = substr(colnames(customData1$rawData[['combined']]),2,3))
customData1$meta[['combined']]$Treatment <- ifelse(customData1$meta[['combined']]$Treatment == "C", "control", "fertilised")

save(fullData, customData1, file= file.path("/mnt/picea/home/bastian/shiny-apps/flakaliden-2011-fungi", "data/fullData.RData"))