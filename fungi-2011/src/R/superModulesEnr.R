# Simon's request for grouping modules
# #200120
# 1, 3, 4, 8, and 14 grouped together, I mean
# And then 2, 5, 6, 7, 9 and 11 together

source('~/Git/UPSCb-common/src/R/gopher.R')
source("~/Git/Rtoolbox/src/plotEnrichedTreemap.R")

projectFolder <- "/mnt/picea/projects/spruce/vhurry/fungi-flakaliden-2011"
networkDir <- file.path(projectFolder,"networks")
deDir <- file.path(projectFolder,"DE")

background <- read.table(file.path(deDir, "background.txt"), stringsAsFactors = F, header=F)[[1]]

infomapCombined <- read.table(file.path(networkDir, "combined/cluster/InfomapClusters.tsv"), 
                              header=T, stringsAsFactors = F)

superA <- c("Cluster1", "Cluster3", "Cluster4", "Cluster8", "Cluster14")
superB <- c("Cluster2", "Cluster5", "Cluster6", "Cluster7", "Cluster9", "Cluster11")

moduleA <- infomapCombined[infomapCombined$cluster %in% superA,]$gene
moduleB <- infomapCombined[infomapCombined$cluster %in% superB,]$gene

enrA <- gopher(moduleA, background = infomapCombined$gene, task = list('ko_pathway'), url='ko', alpha=0.05)
enrB <- gopher(moduleB, background = infomapCombined$gene, task = list('ko_pathway'), url='ko', alpha=0.05)

plotEnrichedTreemap(enrA, enrichment = 'ko_pathway', namespace = 'none',
                          clusterColor=clusterTreemapColors[1],
                          clusterText=clusterTreemapText[1], 
                          title = paste("Clusters: 1,3,4,8,14 - "," KO enrichment"))

plotEnrichedTreemap(enrB, enrichment = 'ko_pathway', namespace = 'none',
                    clusterColor=clusterTreemapColors[3],
                    clusterText=clusterTreemapText[3], 
                    title = paste("Clusters: 2,5,6,7,9,11 - "," KO enrichment"))
