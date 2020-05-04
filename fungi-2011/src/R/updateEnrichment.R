library(here)
source("~/Git/UPSCb-common/src/R/gopher.R")
source("~/Git/Rtoolbox/src/getEigengenes.R")

currentDate <- Sys.Date()

conditions <- c("control", "fertilised", "combined")
fungiDir <- "/mnt/picea/projects/spruce/vhurry/fungi-flakaliden-2011"

fungiDE <- file.path(fungiDir, "DE")
fungiNetworks <- file.path(fungiDir, "networks")

load(file.path(fungiDE, "controlVsdData.RData"))
load(file.path(fungiDE, "fertilisedVsdData.RData"))
load(file.path(fungiDE, "combinedVsdData.RData"))

fungiData <- list(controlData, fertilisedData, combinedData)
names(fungiData) <- conditions
# gopher parameters in case we need to change them
background <- read.table(file.path(fungiDE, "background.txt"), sep='\t', header=F, stringsAsFactors = F)[,1] 
url <- 'ko'
alpha <- 0.05
task <- list("ko_pathway")
# main loop it will create a list of lists
# spruceEnr
#   control
#     go -> tibble
#     mapman -> tibble
#     ...
#   fertilised
#     go -> tibble
#     ...
#   ...

fungiEnr <- lapply(conditions, function(cond){
  print(cond)
  clusterTable <- read.table(file.path(fungiNetworks, cond,  paste0(cond,"InfomapClusters.tsv")), sep='\t', header=T, stringsAsFactors = F)
  clusters <- unique(clusterTable$cluster)
  #print(length(clusters))
  enr <- lapply(clusters, function(clus){
    print(clus)
    genes <- clusterTable[clusterTable$cluster==clus,]$gene
    gopher(genes, task=task, background=background, url=url, alpha=alpha)
  })
  names(enr) <- clusters
  enr
})

names(fungiEnr) <- conditions

#save data with generated name
save(fungiEnr, file=file.path(fungiNetworks, paste0("fungiEnr_", as.character(currentDate), ".RData")))


fungiSplitEnr <- lapply(conditions, function(cond){
  print(cond)
  clusterTable <- read.table(file.path(fungiNetworks, cond,  paste0(cond,"InfomapClusters.tsv")), sep='\t', header=T, stringsAsFactors = F)
  clusters <- unique(clusterTable$cluster)
  #print(length(clusters))
  enr <- lapply(clusters, function(clus){
    print(clus)
    genes <- clusterTable[clusterTable$cluster==clus,]$gene
    genes <- getEigenGenes(fungiData[[cond]], genes)
    prof <- lapply(c("positive", "negative"), function(pr){
      res <- NULL
      if(length(genes[[pr]])>0){
        print(length(genes[[pr]]))
        res <- gopher(genes[[pr]], task=task, background=background, url=url, alpha=alpha)
      } 
     res
    })
    names(prof) <- c("positive", "negative")
    prof
  })
  names(enr) <- clusters
  enr
})

names(fungiSplitEnr) <- conditions

save(fungiSplitEnr, file=file.path(fungiNetworks, paste0("fungiSplitEnr_", as.character(currentDate), ".RData")))
