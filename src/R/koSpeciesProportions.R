library(dplyr)
library(tidyr)
library(ggplot2)

allFile <- "/mnt/picea/projects/spruce/vhurry/fungi-flakaliden-2011/data/annotation_results.emapper.annotations"
taxonomyFile <- "/mnt/picea/projects/spruce/vhurry/fungi-flakaliden-2011/data/gene_taxonomy_update.tsv"

data <- read.table(allFile, stringsAsFactors=F, header=F, sep='\t', comment.char="", quote = "")
taxo <- read.table(taxonomyFile, header = TRUE, sep = '\t', comment.char = "", stringsAsFactors = FALSE)

# 1. query_name
# 2. seed eggNOG ortholog
# 3. seed ortholog evalue
# 4. seed ortholog score
# 5. Predicted taxonomic group
# 6. Predicted protein name
# 7. Gene Ontology terms 
# 8. EC number
# 9. KEGG_ko
# 10. KEGG_Pathway
# 11. KEGG_Module
# 12. KEGG_Reaction
# 13. KEGG_rclass
# 14. BRITE
# 15. KEGG_TC
# 16. CAZy 
# 17. BiGG Reaction
# 18. tax_scope: eggNOG taxonomic level used for annotation
# 19. eggNOG OGs 
# 20. bestOG (deprecated, use smallest from eggnog OGs)
# 21. COG Functional Category
# 22. eggNOG free text description
# 

tempData <- data.frame(gene=data[[1]], ko=data[[9]])
tempData$ko <- gsub("ko:","", tempData$ko)
tempData <- tempData[tempData$ko!="",]

tempData <- left_join(tempData, taxo, by="gene")

#ko2species <- data.frame(ko=tempData$ko, species=tempData$species)
#temp <- separate_rows(ko2species[1:7,], ko, species, sep=',', convert = TRUE)

write.table(tempData, file="/mnt/picea/projects/spruce/vhurry/fungi-flakaliden-2011/data/ko2Species.tsv",
            sep='\t', row.names = F, col.names = F, quote=F)

#Run python3 script
#python ko2species.py 

mapFile <- "/mnt/picea/projects/spruce/vhurry/fungi-flakaliden-2011/data/ko2speciesMap.tsv"

mapData <- read.table(mapFile, stringsAsFactors=F, header=F, sep='\t', comment.char="", quote = "")

temp <- separate_rows(mapData, V2, sep=',', convert = TRUE)

kos <- unique(temp$V1)
theList <- lapply(kos, function(ko){
  res <- temp[temp$V1==ko,]$V2
  res
})

names(theList) <- kos

load("splitClusters.RData")

#using sort(table(megaData$species), decreasing = T)[1:50][!grepl("Unclassified",names(sort(table(megaData$species), decreasing = T)[1:50]))]
# we get the top 40 species with more representation
# in gene number without unclassified

# Piloderma croceum           Hyaloscypha variabilis 
# 38665                            13388 
# Hyaloscypha bicolor             Cenococcum geophilum 
# 13180                            12902 
# Thelephora ganbajun            Cortinarius glaucopus 
# 9030                             8540 
# Mycena galopus          Rhodocollybia butyracea 
# 6802                             4725 
# Hebeloma cylindrosporum              Hygrophorus russula 
# 3681                             2039 
# Acephala macrosclerotiorum                 Hydnum rufescens 
# 1892                             1429 
# Mycena amicta        Archaeorhizomyces finlayi 
# 1351                             1021 
# Trichophaea hybrida               Mycena epipterygia 
# 932                              652 
# Chalara longipes Fibularhizoctonia sp. CBS 109695 
# 458                              444 
# Pezoloma ericae       Archaeorhizomyces borealis 
# 429                              394 

top20Species<-c("Piloderma croceum","Hyaloscypha variabilis",
                "Hyaloscypha bicolor","Cenococcum geophilum",
                "Thelephora ganbajun","Cortinarius glaucopus",
                "Mycena galopus","Rhodocollybia butyracea",
                "Hebeloma cylindrosporum","Hygrophorus russula",
                "Acephala macrosclerotiorum","Hydnum rufescens",
                "Mycena amicta","Archaeorhizomyces finlayi",
                "Trichophaea hybrida","Mycena epipterygia",
                "Chalara longipes","Fibularhizoctonia sp. CBS 109695",
                "Pezoloma ericae","Archaeorhizomyces borealis",
                "Tulasnella calospora","Gymnopus earleae",
                "Glonium stellatum","Russula ochroleuca",
                "Oidiodendron maius","Wilcoxina mikolae",
                "Septobasidium sp.","Galerina marginata",
                "Venturia inaequalis","Xenasmatella vaga")



# allClusters <- c(rbind(unique(names(splitClusters)), paste0(unique(names(splitClusters)),"i")))

# 
# test <- lapply(unique(names(splitClustersEnr)), function(cluster){
#   
#   lapply(c('positive','negative'), function(profile)  {
#     y <- splitClustersEnr[[cluster]][[profile]]
#     if(!is.null(y$ko2species)) {
#       if(profile=="positive") {
#         y$ko2species$cluster <- rep(cluster, nrow(y$ko2species))
#       }else{
#         y$ko2species$cluster <- rep(paste0(cluster,"i"), nrow(y$ko2species))
#       }
#       #  
#     }
#     y$ko2species
#     
#   })  
# })


newDF <- data.frame(cluster=character(), species=character(), n=double(), stringsAsFactors=FALSE)

for(cluster in unique(names(splitClusters)) ) {
  print(cluster)
  
  if(!is.null(splitClusters[[cluster]]$positive)) {
    listPos <- theList[splitClusters[[cluster]]$positive]
    
    posSpecies <- list()
    
    for (ln in listPos){
      tablePos <- table(ln)
      for(item in 1:length(tablePos)) {
        specie <- names(tablePos)[item]
        if(!is.null(posSpecies[[specie]])) {
          posSpecies[[specie]] <- posSpecies[[specie]] + tablePos[item]
        } else {
          posSpecies[[specie]] <- 0
        }

      }
    }
    
    posSums <- 0 
    for(i in posSpecies){ posSums <- posSums + i }

    for(item in 1:length(posSpecies)) {
      newLine <- data.frame(cluster=cluster, species=names(posSpecies)[item], n=round((posSpecies[[item]]*100)/posSums),
                            stringsAsFactors = F)
      newDF <- rbind(newDF, newLine)

    }
        
    #posSums <- 0 
    #for(i in listPos){ posSums <- posSums + (length(i))}
    
    # for (ln in listPos){
    #   tablePos <- table(ln)  
    #   for(item in 1:length(tablePos)) {
    #     newLine <- data.frame(cluster=cluster, species=names(tablePos)[item], n=round((tablePos[item]*100)/posSums),
    #                           stringsAsFactors = F)
    #     newDF <- rbind(newDF, newLine)
    #     
    #   }
    # }
  }
  
  if(!is.null(splitClusters[[cluster]]$negative)) {
    listNeg <- theList[splitClusters[[cluster]]$negative]
    
    negSpecies <- list()
    
    for (ln in listNeg){
      tableNeg <- table(ln)
      for(item in 1:length(tableNeg)) {
        specie <- names(tableNeg)[item]
        if(!is.null(negSpecies[[specie]])) {
          negSpecies[[specie]] <- negSpecies[[specie]] + tableNeg[item]
        } else {
          negSpecies[[specie]] <- 0
        }
        
      }
    }
    
    negSums <- 0 
    for(i in negSpecies){ negSums <- negSums + i }
    
    for(item in 1:length(negSpecies)) {
      newLine <- data.frame(cluster=paste0(cluster,"i"), species=names(negSpecies)[item], n=round((negSpecies[[item]]*100)/negSums),
                            stringsAsFactors = F)
      newDF <- rbind(newDF, newLine)
      
    }
    

  }
  
  # if(!is.null(splitClusters[[cluster]]$negative)){
  #   listNeg <- theList[splitClusters[[cluster]]$negative]
  #   
  #   negSums <- 0 
  #   for(i in listNeg){ negSums <- negSums + (length(i))}
  # 
  #   for (ln in listNeg){
  #     tableNeg <- table(ln)  
  #     for(item in 1:length(tableNeg)) {
  #       newLine <- data.frame(cluster=paste0(cluster,"i"), species=names(tableNeg)[item], n=round((tableNeg[item]*100)/negSums),
  #                             stringsAsFactors = F)
  #       newDF <- rbind(newDF, newLine)
  #     }
  #   }
  # }
  
}

head(newDF)

defDF <- newDF[newDF$species %in% top20Species,]

newDF[newDF$species == "Piloderma croceum" & newDF$cluster == "Cluster3",]
defDF[defDF$species == "Piloderma croceum" & defDF$cluster == "Cluster3",]

defDF$cluster <-  factor(defDF$cluster, 
       levels = rev(c("Cluster1", "Cluster1i", "Cluster2", "Cluster2i", "Cluster3", "Cluster3i", "Cluster4", "Cluster4i",
                  "Cluster5", "Cluster5i", "Cluster6", "Cluster6i", "Cluster7", "Cluster7i", "Cluster8", "Cluster8i",  
                  "Cluster9", "Cluster9i","Cluster10" , "Cluster10i", "Cluster11" , "Cluster11i" ,"Cluster12" , "Cluster12i", "Cluster13" , "Cluster13i", "Cluster14" , "Cluster14i",
                "Cluster15",  "Cluster15i", "Cluster16" , "Cluster16i" ,"Cluster17" , "Cluster17i" ,"Cluster18",  "Cluster19" , "Cluster19i",
                "Cluster20" , "Cluster20i", "Cluster21" , "Cluster21i", "Cluster22" , "Cluster22i" ,"Cluster23" , "Cluster23i", "Cluster24" , "Cluster24i" ,"Cluster25", 
                "Cluster25i", "Cluster26","Cluster26i")))
ggplot(defDF, aes(x=species, 
                  y=as.factor(cluster), 
                  group=n
                  )) +
  geom_point(aes(color=n, size=n)) +
  #scale_size_manual(values=c(0, 1, 2, 3, 5, 4, 8, 7))+
  #scale_color_manual(values=c('#999999','#E69F00', '#56B4E9'))+
  expand_limits(x=0) +
  labs(x="Species", y="Clusters", size="% species-transcripts", color="% species-transcripts") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

