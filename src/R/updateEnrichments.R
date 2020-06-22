##################################################
## Project: Flakaliden enrichent updater
## Script purpose: Obtain/update enrichment tests for fungi or spruce
## Date: 2000529
## Author: Alonso Serrano alonso.serrano@slu.se
##################################################

library(here)
library(RCurl)
library(dplyr)
source(here("UPSCb-common/src/R/gopher.R"))
source("~/scripts/koPathwayDiseaesCleanup/koProcessing.R")

#structure
options(stringsAsFactors = FALSE)

enrichClusters <- function(infomapFile, url, task, geneCol="gene", clusterCol="cluster", background="",
                           url2="", task2=""){
  
  clusters <- read.table(infomapFile, sep='\t', header=T, stringsAsFactors = F)
  if(background==""){
    background <- clusters[[geneCol]]
    
  }
  print(length(background))
  enrichedClusters <- lapply(unique(clusters[[clusterCol]]), function(x){
    genes <- clusters[clusters[[clusterCol]] == x,][[geneCol]]
    
    enr <- gopher(genes, background, url=url, task=task)
    # if(!is.null(enr$ko_pathway)) 
    #   enr$ko_pathway <- koPathwayDiseaseCleanup(koTranslate(enr$ko_pathway))
    # if(!is.null(enr$ko)) 
    #   enr$ko <- koTranslate(enr$ko)
    # if(!is.null(enr$kegg)) 
    #   enr$kegg <- koTranslate(enr$kegg)  
    
    if(url2!="" && task2!="") {
      enr2 <- gopher(genes, background, url=url2, task=task2)
      # if(!is.null(enr2$ko_pathway)) 
      #   enr2$ko_pathway <- koPathwayDiseaseCleanup(koTranslate(enr2$ko_pathway))
      # if(!is.null(enr2$ko)) 
      #   enr2$ko <- koTranslate(enr2$ko)
      # if(!is.null(enr2$kegg)) 
      #   enr2$kegg <- koTranslate(enr2$kegg)
      
      enr <- c(enr, enr2)
    }
   
    
    enr
  })
  names(enrichedClusters) <- unique(clusters[[clusterCol]])
  

  return(enrichedClusters)
}

#cenococcum update

infomapFile <- "/mnt/picea/projects/spruce/vhurry/fungi-flakaliden-2011/fungi-species/cenococcum_geophilum/networks/control/results/aggregated/infomapClusters.tsv" 
controlEnr <- enrichClusters(infomapFile, 
                              url="fungi2011", 
                              task=c("cog","ko","ko_pathway","kog"),
                              url2="fungi2011_cenococcum_geophilum",
                              task2="go")
save(controlEnr, file="/mnt/picea/projects/spruce/vhurry/fungi-flakaliden-2011/fungi-species/cenococcum_geophilum/networks/control/results/aggregated/enrichedClusters.RData")


#cortinarius update

infomapFile <- "/mnt/picea/projects/spruce/vhurry/fungi-flakaliden-2011/fungi-species/cortinarius_glaucopus/networks/control/results/aggregated/infomapClusters.tsv" 
controlEnr <- enrichClusters(infomapFile, 
                             url="fungi2011", 
                             task=c("cog","ko","ko_pathway","kog"),
                             url2="fungi2011_cortinarius_glaucopus",
                             task2="go")
save(controlEnr, file="/mnt/picea/projects/spruce/vhurry/fungi-flakaliden-2011/fungi-species/cortinarius_glaucopus/networks/control/results/aggregated/enrichedClusters.RData")

# spruce control update
infomapFile <- "/mnt/picea/projects/spruce/vhurry/spruce-flakaliden-root-2011/spruce-networks/control/cluster/InfomapClusters.tsv" 
controlEnr <- enrichClusters(infomapFile, 
                             url="pabies", 
                             task=c("go","kegg","pfam","mapman")
                             )
save(controlEnr, file="/mnt/picea/projects/spruce/vhurry/spruce-flakaliden-root-2011/spruce-networks/control/cluster/enrichedClusters.RData")

#spruce fertilised update

infomapFile <- "/mnt/picea/projects/spruce/vhurry/spruce-flakaliden-root-2011/spruce-networks/fertilised/cluster/InfomapClusters.tsv" 
fertilisedEnr <- enrichClusters(infomapFile, 
                             url="pabies", 
                             task=c("go","kegg","pfam","mapman")
)
save(fertilisedEnr, file="/mnt/picea/projects/spruce/vhurry/spruce-flakaliden-root-2011/spruce-networks/fertilised/cluster/enrichedClusters.RData")

#spruce combined update

infomapFile <- "/mnt/picea/projects/spruce/vhurry/spruce-flakaliden-root-2011/spruce-networks/combined/cluster/InfomapClusters.tsv" 
combinedEnr <- enrichClusters(infomapFile, 
                             url="pabies", 
                             task=c("go","kegg","pfam","mapman")
)
save(combinedEnr, file="/mnt/picea/projects/spruce/vhurry/spruce-flakaliden-root-2011/spruce-networks/combined/cluster/enrichedClusters.RData")
