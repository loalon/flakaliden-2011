#library(RLinuxModules)
library(data.table)
library(here)
library(tidyverse)
library(stringr)
library(reshape2)
module("load bioinfo-tools seidr-devel")
module("load bioinfo-tools InfoMap")
library(dplyr)
source(here("Rtoolbox/src/utilsDE.r"))
source(here("/Rtoolbox/src/plot3Dvector.R"))
source(here("/Rtoolbox/src/plotVectorPCA.R"))
source(here("Rtoolbox/src/plotEnrichedTreemap.R"))
source(here("Rtoolbox/src/plotUpAndDown.R"))
source(here("UPSCb-common/src/R/gopher.R"))
source(here("Rtoolbox/src/getEigengenes.R"))

source("~/Git/Rtoolbox/src/infomapTools.R")
source("~/Git/Rtoolbox/src/plotEigenGene.R")
source("~/Git/Rtoolbox/src/plotEnrichedTreemap.R")
source('~/Git/UPSCb-common/src/R/gopher.R')


runInfomap <- function(seidrFile, edgeIndexFile, treeFile, markovTime=1, folder) {
  
  projectDir <- "/mnt/picea/projects/spruce/vhurry/fungi-flakaliden-2011/fungi-species/piloderma_croceum/networks/control/results/aggregated"
  resolveFile <- file.path(projectDir, "seidrResolve.txt")
  edgeIndexFile <- file.path(projectDir, "indexEdgeList.txt")
  seidrFile <- file.path(projectDir, "aggregatedThreshold.sf")
  treeFile <- file.path(projectDir, "indexEdgeList.tree")
  markovTime = 1
  module("load bioinfo-tools seidr-devel")
  module("load bioinfo-tools InfoMap")
  
  system(paste("Infomap", edgeIndexFile,"-z --markov-time", markovTime, "\\."))
  infomapRes <- system(paste("seidr resolve -s", seidrFile, treeFile), intern=TRUE)
  infomapTable <-data.frame(do.call(rbind, strsplit(infomapRes, "\t")))
  infomapTable <- prepareData(infomapTable)
  infomapTable$gene <- gsub("\r","",infomapTable$gene)
  infomapTable$Level1 <- infomapTable$P1
  infomapTable$Level2 <- paste0(infomapTable$Level1,":",infomapTable$P2)
  infomapTable$Level2 <- ifelse(infomapTable$Level2 %like% "NA", NA, infomapTable$Level2)
  
  if("P3" %in% colnames(infomapTable)) {
    infomapTable$Level3 <- paste0(infomapTable$Level1,":",infomapTable$P2, ":",infomapTable$P3)
    infomapTable$Level3 <- ifelse(infomapTable$Level3 %like% "NA", NA, infomapTable$Level3)
  }
  return(infomapTable)
}


projectDir <- "/mnt/picea/projects/spruce/vhurry/fungi-flakaliden-2011/fungi-species/piloderma_croceum/networks/control/results/aggregated"
resolveFile <- file.path(projectDir, "seidrResolve.txt")
edgeFile <- file.path(projectDir, "edgeIndexFile.txt")
