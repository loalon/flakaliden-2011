---
  title: "Flakaliden spruce roots + fungi cluster enrichment analysis"
author: "Alonso Serrano"
date: "`r format(Sys.time(), '%d %B %Y')`"
output: 
  html_document: 
  number_sections: yes
toc: yes
---
  
  <style>
  * {
    box-sizing: border-box;
  }

/* Create two equal columns that floats next to each other */
  .column {
    float: left;
    width: 50%;
    padding: 10px;
    
  }

/* Clear floats after the columns */
  .row:after {
    content: "";
    display: table;
    clear: both;
  }
</style>
  
  ```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

# Introduction
This file contains the script and results of the data from spruce roots and fungi from the 2011 Flakaliden dataset. 

## Prerequisites
The variance stabilization transformed dataset, either as TSV or RData obtained from de differential expression analysis step.
Infomap clustering results from Seidr network
Pre enriched data using gofer

## Setup
```{r generalSetup, message=FALSE, warning=FALSE}
library(data.table)
library(dplyr)
library(here)

# preload packages for depencies
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(tidyr)
library(tibble)
library(jsonlite)
library(treemap)
library(KEGGREST)
library(reshape2)

source("~/Git/Rtoolbox/src/plotEigenGene.R")
source("~/Git/Rtoolbox/src/getEigengenes.R")
source("~/Git/Rtoolbox/src/plotEnrichedTreemap.R")
source("~/Git/UPSCb-common/src/R/gopher.R")
source("~/scripts/koPathwayDiseaesCleanup/koProcessing.R")
#source("corAuxFunctions.R")

spruceDir <- "/mnt/picea/projects/spruce/vhurry/14_SpruceRoots_Project"
fungiDir <- "/mnt/picea/projects/spruce/vhurry/fungi-flakaliden-2011"
cortinariusDir <- "/mnt/picea/projects/spruce/vhurry/fungi-flakaliden-2011/fungi-species/cortinarius_glaucopus"
cenococcumDir <- "/mnt/picea/projects/spruce/vhurry/fungi-flakaliden-2011/fungi-species/cenococcum_geophilum"

spruceDE <- file.path(spruceDir,"DE")
fungiDE <- file.path(fungiDir,"DE")
cortiDE <- file.path(spruceDir,"DE")
cenoDE <- file.path(fungiDir,"DE")

checkFile <- function(fileName) {
  if (!file.exists(fileName)) stop(paste("File", fileName, "doesn't exist"))
}
```

## Check files

### VSD data file
The files are tab separated, rows contain sample meta information (week and treatment),
columns contain gene or KEGG orthologs names. The values are expression data after variance stabilization transformation

### Cluster files
They must contain two columns gene and cluster and be tab separated<br />
  ex.<br />
  gene	cluster<br />
  MA_5320g0010	Cluster1<br />
  MA_32763g0010	Cluster1<br />
  MA_66378g0010	Cluster1<br />
  
  
  ```{r loadRdata}
load(file.path(spruceDE, "controlVsdData.RData"))
controlSpruceData <- controlData
load(file.path(spruceDE, "fertilisedVsdData.RData"))
fertilisedSpruceData <- fertilisedData

load(file.path(fungiDE, "controlVsdData.RData"))
controlFungiData <- controlData
load(file.path(fungiDE, "fertilisedVsdData.RData"))
fertilisedFungiData <- fertilisedData
```

```{r}
# load spruce data files
controlSpruceClusterFile <- file.path(spruceDir,"spruce-networks/control/cluster","InfomapClusters.tsv")
checkFile(controlSpruceClusterFile)

# load fungi data files
controlFungiClusterFile <- file.path(fungiDir,"networks/control/cluster","InfomapClusters.tsv")
checkFile(controlFungiClusterFile)
```

## Load files
```{r loadFiles}
#name correction, spruce rownames, have an extra dash
rownames(controlSpruceData) <- gsub("-", "", rownames(controlSpruceData))
rownames(fertilisedSpruceData) <- gsub("-", "", rownames(fertilisedSpruceData))

controlSpruceClusters <- read.table(controlSpruceClusterFile, stringsAsFactors=F, header=T, sep='\t', comment.char="")
controlFungiClusters <- read.table(controlFungiClusterFile, stringsAsFactors=F, header=T, sep='\t', comment.char="")
```

## Create metadata based on the information on sample names
```{r}
controlSpruceGroup <- substr(rownames(controlSpruceData),4,4)
controlSpruceTime <- as.integer(substr(rownames(controlSpruceData),2,3))

controlFungiGroup <- substr(rownames(controlFungiData),4,4)
controlFungiTime <- as.integer(substr(rownames(controlFungiData),2,3))

fertilisedSpruceGroup <- substr(rownames(fertilisedSpruceData),4,4)
fertilisedSpruceTime <- as.integer(substr(rownames(fertilisedSpruceData),2,3))

fertilisedFungiGroup <- substr(rownames(fertilisedFungiData),4,4)
fertilisedFungiTime <- as.integer(substr(rownames(fertilisedFungiData),2,3))
```

# Control profile matching
```{r controlMatching, echo=TRUE, message=FALSE}
control.res <- profileSimilarity(controlSpruceData, controlFungiData, controlSpruceClusters, controlFungiClusters, controlSpruceGroup, controlFungiGroup, controlSpruceTime, controlFungiTime, 
                                 "Spruce-control", "Fungi-control")

control.res.sort <- control.res[order(control.res$score, decreasing = T),]
```

## Best correlated
```{r controlHead}
head(control.res.sort)
```
## Best reverse correlated
```{r controlTail}
tail(control.res.sort)
```
## Save results to file
```{r controlSave}
# control.res.sort$`Spruce-control`<- paste0("spruce",control.res.sort$`Spruce-control`)
# control.res.sort$`Fungi-control`<- paste0("fungi",control.res.sort$`Fungi-control`)
write.table(control.res.sort, file=here("results",paste0("controlProfileMatching_",specie,"_",Sys.Date(),".tsv")), sep='\t', quote=F, row.names=F)
```

```{r}
spruceControl <- cluster2Network(controlSpruceData, controlSpruceClusters, "spruce")
fungiControl <- cluster2Network(controlFungiData, controlFungiClusters, "fungi")
write.table(rbind(spruceControl, fungiControl), 
            file=here("results",paste0("controlClusterNetworkMeta",specie,"_",Sys.Date(),".tsv")),
            sep='\t', row.names = F, quote=F)

# spruceFertilised <- cluster2Network(fertilisedSpruceData, fertilisedSpruceClusters, "spruce")
# fungiFertilised <- cluster2Network(fertilisedFungiData, fertilisedFungiClusters, "fungi")
# 
# write.table(rbind(spruceFertilised, fungiFertilised),
#             file=here("results",paste0("fertilisedClusterNetworkMeta",specie,"_",Sys.Date(),".tsv"))
#             , sep='\t', row.names = F, quote=F)
```

## Obtain profile degradation matrix
```{r}
load(here("env", paste0(specie, "_controlDEG.RData")))
```


```{r controlDEF, eval=F}
controlDEG <- getProfileDegradationMatrix(control.res.sort, 
                                          conditionSpruceClusters=controlSpruceClusters, 
                                          conditionFungiClusters=controlFungiClusters,
                                          data1Spruce=controlSpruceData, 
                                          data1Fungi=controlFungiData, 
                                          data2Spruce=fertilisedSpruceData, 
                                          data2Fungi=fertilisedFungiData,
                                          group1Spruce=controlSpruceGroup, 
                                          group1Fungi=controlFungiGroup,
                                          group2Spruce=fertilisedSpruceGroup , 
                                          group2Fungi=fertilisedFungiGroup ,
                                          time1Spruce=controlSpruceTime , 
                                          time1Fungi=controlFungiTime,
                                          time2Spruce=fertilisedSpruceTime, 
                                          time2Fungi=fertilisedFungiTime ,
                                          conditionNames = c("control", "fertilised")
                                          
)


```

```{r, eval=F}
save(controlDEG, file=here("env", paste0(specie, "_controlDEG.RData")))
```

```{r}
write.table(controlDEG, file=here("results", paste0("controlDEGfull_",specie,"_",Sys.Date(),".tsv")), quote = F, sep='\t', row.names = F)
```

## Reduced violin plot
```{r controlViolin, warning=FALSE}
controlDEG.df <- controlDEG[controlDEG$control_score>0,]
controlDEG.df <- controlDEG.df[order(controlDEG.df$degradation, decreasing = F),]
controlDEG.df <- melt(controlDEG.df[,c("control_score","fertilised_score")])

p <- ggplot(controlDEG.df, aes(x=variable, y=value  )) + 
  geom_violin() +
  stat_summary(fun.y=mean, geom="point", shape=23, size=2) + 
  stat_summary(fun.y=median, geom="point", size=2, color="red") +
  geom_boxplot(width=0.1)

print(p)
```

```{r eval=FALSE, include=FALSE}
png(here("results", paste0("controlDEGplot_",specie,"_",Sys.Date(),".png")),
    units="px", width=1920, height=1080, res=300)
print(p)
dev.off()

write.table(controlDEG.df, file=here("results", paste0("controlDEG_",specie,"_",Sys.Date(),"_df_sort.tsv")), quote = F, sep='\t', row.names = F)

```

## Full violin plot
```{r controlViolin2, warning=FALSE}
controlDEG.df <- controlDEG
controlDEG.df <- controlDEG.df[order(controlDEG.df$degradation, decreasing = F),]
controlDEG.df <- melt(controlDEG.df[,c("control_score","fertilised_score")])

p <- ggplot(controlDEG.df, aes(x=variable, y=value  )) + 
  geom_violin() +
  stat_summary(fun.y=mean, geom="point", shape=23, size=2) + 
  stat_summary(fun.y=median, geom="point", size=2, color="red") +
  geom_boxplot(width=0.1)

print(p)
```

```{r eval=FALSE, include=FALSE}
# png(here("results", paste0("controlDEGplot_",specie,"_",Sys.Date(),"_full.png")),
#units="px", width=1920, height=1080, res=300)
print(p)
# dev.off()
```

Average degradation
```{r}
controlAvg <- mean(controlDEG$degradation)
print(controlAvg)
```

## Significant degradations
Get the maximum difference

```{r}
res.pos <- controlDEG %>% 
  #filter((control_score > 0) & (fertilised_score > 0)) #%>% 
  #filter((degradation<1) & (degradation > -1))
  filter(degradation>0)

res.pos <- res.pos[order(res.pos$degradation, decreasing = T),]

res.neg <- controlDEG %>% 
  #filter((control_score < 0) & (fertilised_score < 0)) #%>% 
  #filter((degradation<1) & (degradation > -1))
  filter(degradation<0)

res.neg <- res.neg[order(res.neg$degradation, decreasing = T),]

write.table(rbind(res.pos, res.neg), file=here("results", paste0("maxDifference_",specie,"_",Sys.Date(),".tsv")), quote = F, sep='\t', row.names = F)
```


Get the top 20 most important changes(10 degradations, 10 recovers)
```{r}
mostDegraded <- rbind(head(res.pos, n=10),
                      tail(res.neg, n=10))
mostDegraded
```

Enrich & plot results


```{r, eval=T, include=F}
load(file=here("env", paste0(specie, "_enrRes.RData")))
```

```{r superLoop, eval=F}
enrRes <- lapply(1:20, function(x) {
  spruceClus <- mostDegraded[x,]$`Spruce-control`
  fungiClus <- mostDegraded[x,]$`Fungi-control`
  
  res <- list()
  spruceMulti <- enrichThis(controlSpruceData,  fertilisedSpruceData,
                            controlSpruceClusters, spruceClus,
                            controlSpruceClusters$gene, url="pabies", type=c("go","pfam","mapman","kegg")
  )
  # only work with the top 5 species
  fungiMulti <- enrichThis(controlFungiData,  fertilisedFungiData,
                           controlFungiClusters, fungiClus,
                           controlFungiClusters$gene, url="ko", type="ko_pathway"
  )
  
  # fungiMultiGeneric <- enrichThis(controlFungiData,  fertilisedFungiData,
  #          controlFungiClusters, fungiClus,
  #          controlFungiClusters$gene, url="fungi2011", type=c("cog","ko","ko_pathway","kog") 
  #          )
  
  res$spruce <- spruceMulti
  res$fungi <- fungiMulti
  res
})
```

```{r, eval=F,include=F}
save(enrRes, file=here("env", paste0(specie, "_enrRes.RData")))
```





## Enrich & plot
```{r superLoop2, fig.hold='hold', message=FALSE, warning=FALSE, out.width="50%", results='asis'}
None <- lapply(1:20, function(x) {
  spruceClus <- mostDegraded[x,]$`Spruce-control`
  fungiClus <- mostDegraded[x,]$`Fungi-control`
  cat("\n")
  cat("<hr>")
  cat(paste0("<H3>Spruce cluster: ", spruceClus, ". Fungi cluster ",fungiClus, ". Degradation: ", round(mostDegraded[x,]$degradation, 4),"</H3>"))
  cat("\n\n")
  spruceMulti <- plotThis(controlSpruceData,  fertilisedSpruceData,
                          controlSpruceClusters,              controlSpruceGroup, 
                          controlSpruceTime, 
                          fertilisedSpruceGroup, 
                          fertilisedSpruceTime, 
                          spruceClus, 
                          
                          title1 = paste("Spruce control cluster", spruceClus, "with control data"),
                          title2 = paste("Spruce control cluster", spruceClus, "with fertilised data")
                          
  )
  
  fungiMulti <- plotThis(controlFungiData,  fertilisedFungiData,
                         controlFungiClusters, 
                         controlFungiGroup, 
                         controlFungiTime, 
                         fertilisedFungiGroup, 
                         fertilisedFungiTime, 
                         fungiClus, 
                         
                         title1 = paste("Fungi control cluster", fungiClus, "with control data"),
                         title2 = paste("Fungi control cluster", fungiClus, "with fertilised data")
  )
  
  cat(paste("<H4>Profiles plot</H4>"))
  cat("\n\n")
  cat(paste("<center>Spruce - fungi control profiles correlation",round(mostDegraded[x,]$control_score,4), "</center>"))
  cat("\n\n")
  par(mar = c(4, 4, 0.1, 0.1))
  print(spruceMulti$plot1)
  print(fungiMulti$plot1)
  cat(paste("<div class='row'>",
            "<div class='column'><center>Spruce correlation", round(spruceMulti$cor1, 4), "</center></div>",
            "<div class='column'><center>Fungi correlation", round(fungiMulti$cor1, 4), "</center></div></div>"
  ))
  
  print(spruceMulti$plot2)
  print(fungiMulti$plot2)
  par(mfrow=c(1,1))
  cat("\n\n")
  cat("\n\n")
  cat(paste("<center>Spruce - fungi fertilised profiles correlation",round(mostDegraded[x,]$fertilised_score,4),"</center>"))
  cat("\n\n")
  cat("\n\n")
  
  cat(paste("<H4>Spruce treemaps</H4>"))
  cat("\n\n")
  par(mar = c(4, 4, 0.1, 0.1))
  None <- sapply(c("go","pfam","mapman","kegg"), function(y){
    if(!is.null(enrRes[[x]]$spruce[[y]])) {
      if(y=="kegg"){
        enrRes[[x]]$spruce[[y]] <- koTranslate(enrRes[[x]]$spruce[[y]])
      }
      plotEnrichedTreemap(enrRes[[x]]$spruce, enrichment = y, namespace = "none", title=paste(y,"enrichment"))
    }
  })
  par(mfrow=c(1,1))
  cat("\n\n")
  cat("<hr>")
  
  cat(paste("<H4>Fungi treemaps</H4>"))
  cat("\n\n")
  par(mar = c(4, 4, 0.1, 0.1))
  
  None <- sapply(c("ko_pathway"), function(y){
    if(!is.null(enrRes[[x]]$fungi[[y]])) {
      if(y=="ko_pathway"){
        enrRes[[x]]$fungi[[y]] <- koPathwayDiseaseCleanup(koTranslate(enrRes[[x]]$fungi[[y]]))
      }
      plotEnrichedTreemap(enrRes[[x]]$fungi, enrichment = y, namespace = "none",title=paste(y,"enrichment"))
    }
    
  })
  par(mfrow=c(1,1))
  cat("\n\n")
  cat("<hr>")
  cat("<hr>")
})

```
```{r splitClusters}

controlFungiClusters

splitClusters <- lapply(unique(controlFungiClusters$cluster), function(cluster){
  print(cluster)
})

splitClusters <- lapply(clusterNames, function(cluster){
  theCluster <- InfomapClusters[InfomapClusters$cluster==cluster,]$gene
  cat(paste0("<H3>", cluster, "","</H3>"))
  cat("\n\n")
  res <- getEigenGenes(theData, theCluster)
  pos <-length(res$positive)
  neg <- length(res$negative)
  posRatio <- round((pos/(pos+neg)) * 100, digits = 2)
  
  ratio <- paste0(posRatio,"-",100-posRatio)
  cat(paste("Genes in main profile", pos, "\n"))
  cat(paste("Genes in reverse profile", neg, "\n"))
  cat(paste("Main-reverse ratio", ratio, "\n"))
  res
})
names(splitClusters) <- clusterNames
```

```{r, eval=F, include=F}
save.image(file.path(here("env"),"env.RData"))
```

# Session information
```{r}
sessionInfo()
```
