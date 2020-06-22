##################################################
## Project: Flakaliden fungi 2011
## Script purpose: Obtain enrichment tests for fungi KO
## Date: 2000117
## Author: Alonso Serrano alonso.serrano@slu.se
##################################################

library(here)
library(RCurl)
library(dplyr)
source(here("UPSCb-common/src/R/gopher.R"))
source("~/scripts/koPathwayDiseaesCleanup/koProcessing.R")

#structure
options(stringsAsFactors = FALSE)

projectDir <- "/mnt/picea/projects/spruce/vhurry/fungi-flakaliden-2011/" 
deFolder <- file.path(projectDir,"DE")
dataFolder <- file.path(projectDir,"data")
controlClusters <- read.table(file.path(projectDir, "networks/control/cluster/InfomapClusters.tsv"), sep='\t', header=T, stringsAsFactors = F)
fertilisedClusters <- read.table(file.path(projectDir, "networks/fertilised/cluster/InfomapClusters.tsv"), sep='\t', header=T, stringsAsFactors = F)
combinedClusters <- read.table(file.path(projectDir, "networks/combined/cluster/InfomapClusters.tsv"), sep='\t', header=T, stringsAsFactors = F)

controlBackground <- controlClusters$gene
fertilisedBackground <- fertilisedClusters$gene
combinedBackground <- combinedClusters$gene
#load("/mnt/picea/projects/spruce/vhurry/spruce-flakaliden-root-2011/spruce-networks/control/cluster/controlEnr.RData")

controlEnr <- lapply(unique(controlClusters$cluster), function(x){
  genes <- controlClusters[controlClusters$cluster == x,]$gene
  
  enr <- gopher(genes, controlBackground, url="ko", task=c("ko_pathway"))
  if(!is.null(enr$ko_pathway)){ 
    enr$ko_pathway <- koPathwayDiseaseCleanup(koTranslate(enr$ko_pathway))
  }
  enr
})
names(controlEnr) <- unique(controlClusters$cluster)

fertilisedEnr<- lapply(unique(fertilisedClusters$cluster), function(x){
  genes <- fertilisedClusters[fertilisedClusters$cluster == x,]$gene
  
  enr <- gopher(genes, fertilisedBackground, url="ko", task=c("ko_pathway"))
  if(!is.null(enr$ko_pathway)){ 
    enr$ko_pathway <- koPathwayDiseaseCleanup(koTranslate(enr$ko_pathway))
  }
  enr
})
names(fertilisedEnr) <- unique(fertilisedClusters$cluster)

combinedEnr <- lapply(unique(combinedClusters$cluster), function(x){
  genes <- combinedClusters[combinedClusters$cluster == x,]$gene
  
  enr <- gopher(genes, combinedBackground, url="ko", task=c("ko_pathway"))
  if(!is.null(enr$ko_pathway)){ 
    enr$ko_pathway <- koPathwayDiseaseCleanup(koTranslate(enr$ko_pathway))
  }
  enr
})
names(combinedEnr) <- unique(combinedClusters$cluster)

save(controlEnr, file = file.path(projectDir, "networks/control/cluster/controlEnr.RData"))
save(fertilisedEnr, file = file.path(projectDir, "networks/fertilised/cluster/fertilisedEnr.RData"))
save(combinedEnr, file = file.path(projectDir, "networks/combined/cluster/combinedEnr.RData"))