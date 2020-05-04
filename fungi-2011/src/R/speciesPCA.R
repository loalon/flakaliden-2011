# SPECIES PCA
# 
# The purpose of this is script is to obtain a PCA based on the species
library(dplyr)
library(DESeq2)
source("~/Git/Rtoolbox/src/plot2Dpca.R")
library(plotly)
library(scales)
library(htmlwidgets)


projectFolder <- "/mnt/picea/projects/spruce/vhurry/fungi-flakaliden-2011"
dataFolder<- file.path(projectFolder, "data")
resultsFolder<- file.path(projectFolder, "reports")

load(file.path(dataFolder,"megaData.RData"))


speciesData <- megaData[,c(1:112, 120)]
colnames(speciesData)
#remove genes
speciesData <-speciesData[,-1]
speciesData <- speciesData %>%
  select(species, everything())

#if control
# speciesDataName <- speciesData$species
# speciesData <- speciesData[,grepl("C", colnames(speciesData))]
# speciesData$species <- speciesDataName
# 
# #if fert
# speciesDataName <- speciesData$species
# speciesData <- speciesData[,grepl("F", colnames(speciesData))]
# speciesData$species <- speciesDataName

bySpecies <- speciesData %>%
  group_by(species)  %>%
  summarise_if(is.numeric, sum, na.rm = TRUE)

speciesDF <- data.frame(bySpecies)
rownames(speciesDF) <- speciesDF$species

speciesDF <- speciesDF[,-1]

#subsetting

speciesMeta <- data.frame(time = substr(colnames(speciesDF), 2,3),
                          treatment = substr(colnames(speciesDF), 4,4))
#if control
#speciesMeta <- speciesMeta[grepl("C", speciesMeta$treatment),]
#
#if fert
#speciesMeta <- speciesMeta[grepl("F", speciesMeta$treatment),]

#export a summary
#write.table(rowSums(speciesDF), file="speciesPCAsummary.tsv", sep='\t', quote = F, col.names = F)

#trophic mode and species within
# TvsS <- data.frame(trophic=megaData$Trophic.Mode, species=megaData$species, stringsAsFactors = F)
# 
# TvsS.res <- sapply(unique(TvsS$trophic),function(x){
#   print(x)
#   y <- TvsS[TvsS$trophic==x, ]$species
#   length(unique(y))
# })



SvsT <- data.frame(species=megaData$species, trophic=megaData$Trophic.Mode, stringsAsFactors = F)
SvsT.res <- SvsT %>% group_by(species) %>% filter(row_number(species) == 1) %>% arrange(species)

speciesDF <-speciesDF[!(row.names(speciesDF) %in% c("Unclassified.Ascomycota", "Unclassified.Basidiomycota")), ]
#TvsS.res <-TvsS.res[!(TvsS.res$species %in% c("Unclassified.Ascomycota", "Unclassified.Basidiomycota")), ]
SvsT.res <-SvsT.res[!(SvsT.res$species %in% c("Unclassified.Ascomycota", "Unclassified.Basidiomycota")), ]

SvsTrait <- data.frame(species=megaData$species, trait=megaData$Trait, stringsAsFactors = F)
SvsTrait.res <- SvsTrait %>% group_by(species) %>% filter(row_number(species) == 1) %>% arrange(species)

speciesDF <-speciesDF[!(row.names(speciesDF) %in% c("Unclassified.Ascomycota", "Unclassified.Basidiomycota")), ]
#TvsS.res <-TvsS.res[!(TvsS.res$species %in% c("Unclassified.Ascomycota", "Unclassified.Basidiomycota")), ]
SvsT.res <-SvsT.res[!(SvsT.res$species %in% c("Unclassified.Ascomycota", "Unclassified.Basidiomycota")), ]
SvsTrait.res <-SvsTrait.res[!(SvsTrait.res$species %in% c("Unclassified.Ascomycota", "Unclassified.Basidiomycota")), ]


trophicColors=c("#DB3832", "#E26D38",
                "#3F2313", "#764693",
                "#FEDD4F", "#00FF00",
                "#459E55","#264294", "#000000")
names(trophicColors) <- c("Pathotroph","Pathotroph-Saprotroph",
                          "Pathotroph-Saprotroph-Symbiotroph","Pathotroph-Symbiotroph",
                          "Saprotroph","Saprotroph-Pathotroph-Symbiotroph",
                          "Saprotroph-Symbiotroph","Symbiotroph", "-")

traitColors=c("#D3D3D3", "#000000", "#E26D38","#3F2313", "#764693",
              "#FEDD4F", "#00FF00", "#459E55","#264294",  "#800080")

names(traitColors) <- c("NULL", "-","White Rot", "Soft Rot", "Brown Rot; White Rot", "Brown Rot",
                        "Poisonous",  "Hypogeous", "Brown Rot-White Rot", "Blue-Staining")



ddsSpecies <- DESeqDataSetFromMatrix((speciesDF), colData=speciesMeta, design = ~1)

vsd.QA <- varianceStabilizingTransformation(ddsSpecies, blind=TRUE)

assay(vsd.QA) <- assay(vsd.QA) - min(assay(vsd.QA))

#vsd.QA.t <- (assay(vsd.QA))
pca <- prcomp(assay(vsd.QA))

speciesSize <- rescale(rowSums(speciesDF), c(1, 20))

rowSums(speciesDF)
# floor (log10 (abs (30))) + 1


df3D <- data.frame(PC1= pca$x[,1],
                   PC2 = pca$x[,2],
                   PC3 = pca$x[,3],
                   group = rownames(pca$x),
                   color= SvsT.res$species,
                   size= speciesSize)
percents <- round(summary(pca)$importance[2,]*100)



p <- plot_ly(df3D, x = ~PC1, y = ~PC2, 
             mode = 'markers', size = ~size, color=~color,
             sizes = c(10, 100), colors= trophicColors,text = ~group, marker = list(opacity = 0.625, sizemode = 'diameter')) %>%
  layout(title = title,
         xaxis = list(title = paste("PC1 (",percents[1],"%)",sep="")),
         yaxis = list(title = paste("PC2 (",percents[2],"%)",sep=""))
         
  )
#htmlwidgets::saveWidget(p, "speciesControlPC1vsPC2_noUnclassified_Basidio_or_Asco.html")
#htmlwidgets::saveWidget(p, "speciesFertilisedPC1vsPC2_noUnclassified_Basidio_or_Asco.html")

p <- plot_ly(df3D, x = ~PC2, y = ~PC3, 
             mode = 'markers', size = ~size,color=~color,
             sizes = c(10, 100), colors= trophicColors,text = ~group, marker = list(opacity = 0.625, sizemode = 'diameter')) %>%
  layout(title = title,
         xaxis = list(title = paste("PC2 (",percents[2],"%)",sep="")),
         yaxis = list(title = paste("PC3 (",percents[3],"%)",sep=""))
         
  )
#htmlwidgets::saveWidget(p, "speciesControlPC2vsPC3_noUnclassified_Basidio_or_Asco.html")
htmlwidgets::saveWidget(p, "speciesFertilisedPC2vsPC3_noUnclassified_Basidio_or_Asco.html")

p <- plot_ly(df3D, x = ~PC1, y = ~PC3, z = ~PC2, 
              mode = 'markers', size = ~size,color=~color,
             sizes = c(10, 100), colors= trophicColors, text = ~group, marker = list(opacity = 0.625, sizemode = 'diameter')) %>%
  add_markers() %>%
  layout(title = title,
         scene=list(xaxis = list(title = paste("PC1 (",percents[1],"%)",sep="")),
                    zaxis = list(title = paste("PC2 (",percents[2],"%)",sep="")),
                    yaxis = list(title = paste("PC3 (",percents[3],"%)",sep=""))
         )
  )
#htmlwidgets::saveWidget(p, "speciesControl3DPCA_noUnclassified_Basidio_or_Asco.html")
htmlwidgets::saveWidget(p, "speciesFertilised3DPCA_noUnclassified_Basidio_or_Asco.html")


###############splitting species
###############

speciesControlDF <- speciesDF[,grepl("C", colnames(speciesDF))]

avgControl <- lapply(unique(substr(colnames(speciesControlDF),1,4)) , function(x){
  ind<-grepl(x, colnames(speciesControlDF))
  
  res <- apply(speciesControlDF[,ind], 1, mean)
  res
})
names(avgControl) <- unique(substr(colnames(speciesControlDF),1,3))

avgControlDF <- as.data.frame(avgControl)
rownames(avgControlDF) <- paste0("control_", rownames(avgControlDF))


speciesFertilisedDF <- speciesDF[,grepl("F", colnames(speciesDF))]

avgFertilised <- lapply(unique(substr(colnames(speciesFertilisedDF),1,4)) , function(x){
  ind<-grepl(x, colnames(speciesFertilisedDF))
  
  res <- apply(speciesFertilisedDF[,ind], 1, mean)
  res
})
names(avgFertilised) <- unique(substr(colnames(speciesFertilisedDF),1,3))

avgFertilisedDF <- as.data.frame(avgFertilised)
rownames(avgFertilisedDF) <- paste0("fertilised_", rownames(avgFertilisedDF))

avgCombinedDF <- rbind(avgControlDF, avgFertilisedDF)

avgCombinedDF <- ceiling(avgCombinedDF)
ddsCombinedSpecies <- DESeqDataSetFromMatrix(avgCombinedDF, colData=data.frame(time=substr(colnames(avgCombinedDF), 2, 3)) , design = ~1)

vsd.QA <- varianceStabilizingTransformation(ddsCombinedSpecies, blind=TRUE)

assay(vsd.QA) <- assay(vsd.QA) - min(assay(vsd.QA))

#vsd.QA.t <- (assay(vsd.QA))
pca <- prcomp(assay(vsd.QA))

#floor (log10 (abs (rowSums(assay(ddsCombinedSpecies))))) + 1
# speciesSize <- floor (log10 (abs (rowSums(assay(ddsCombinedSpecies))))) + 1
# speciesSize <- ifelse(speciesSize ==-Inf, 1, speciesSize)

#speciesSize <- rescale(rowSums(avgCombinedDF), c(1, 20))
speciesSize <- rowSums(avgCombinedDF)
# speciesSize <- ifelse(speciesSize==1, 10,speciesSize)
# speciesSize <- ifelse(speciesSize==2, 10,speciesSize)
# speciesSize <- ifelse(speciesSize==3, 10,speciesSize)
# speciesSize <- ifelse(speciesSize==4, 20,speciesSize)
# speciesSize <- ifelse(speciesSize==5, 20,speciesSize)
# speciesSize <- ifelse(speciesSize==6, 50,speciesSize)
# speciesSize <- ifelse(speciesSize==7, 100,speciesSize)


df3D <- data.frame(PC1= pca$x[,1],
                   PC2 = pca$x[,2],
                   PC3 = pca$x[,3],
                   group = rownames(pca$x),
                   color= rep(SvsT.res$trophic,2),
                   size= speciesSize,
                   line=substr(rownames(pca$x),1,1))
percents <- round(summary(pca)$importance[2,]*100)

#original setting

pc1vs2 <- plot_ly(df3D, x = ~PC1, y = ~PC2, 
             mode = 'markers', 
             size = ~size,
             color=~color,
             sizes = c(2,  100), 
             colors= trophicColors,
             text = ~group, marker = list(opacity = 0.5, sizemode = 'diameter')) %>%
  layout(title = title,
         xaxis = list(title = paste("PC1 (",percents[1],"%)",sep="")),
         yaxis = list(title = paste("PC2 (",percents[2],"%)",sep="")),
         showlegend = FALSE
         
  )

htmlwidgets::saveWidget(pc1vs2, 
                        file.path(resultsFolder, paste0("speciesPCA_PC1vsPC2_",format(Sys.time(),"%Y%m%d"),".html")))

##PC2 vs PC3
pc2vs3 <- plot_ly(df3D, x = ~PC2, y = ~PC3, 
             mode = 'markers', 
             size = ~size,
             color=~color,
             sizes = c(2,100),
             colors= trophicColors,
             text = ~group, marker = list(opacity = 0.5, sizemode = 'diameter')) %>%
  layout(title = title,
         xaxis = list(title = paste("PC2 (",percents[2],"%)",sep="")),
         yaxis = list(title = paste("PC3 (",percents[3],"%)",sep="")),
         showlegend = FALSE
  )

htmlwidgets::saveWidget(pc2vs3, 
                        file.path(resultsFolder, paste0("speciesPCA_PC2vsPC3_",format(Sys.time(),"%Y%m%d"),".html")))

## species interactive explorer
pcaInteractive <- plot_ly(df3D, x = ~PC2, y = ~PC3, 
             mode = 'markers', 
             size = ~size,
             color=~group,
             sizes = c(1,100),
             
             #colors= trophicColors,
             text = ~group, marker = list(opacity = 0.5, sizemode = 'diameter')) %>%
  layout(title = title,
         xaxis = list(title = paste("PC2 (",percents[2],"%)",sep="")),
         yaxis = list(title = paste("PC3 (",percents[3],"%)",sep="")),
         showlegend = TRUE
  )

#htmlwidgets::saveWidget(pcaInteractive, file.path(resultsFolder, "speciesPCA_interactive.html"), selfcontained = T)
htmlwidgets::saveWidget(pcaInteractive, 
                        file.path(resultsFolder, paste0("speciesPCA_interactive_",format(Sys.time(),"%Y%m%d"),".html")))

##3D version
pca3D <- plot_ly(df3D, x = ~PC1, y = ~PC3, z=~PC2,
                  mode = 'markers', 
                  size = ~size,
                  color=~color,
                  sizes = c(2,100),
                  colors= trophicColors,
                  text = ~group, marker = list(opacity = 0.5, sizemode = 'diameter')) %>%
  layout(title = title, scene=list(
         xaxis = list(title = paste("PC1 (",percents[1],"%)",sep="")),
         yaxis = list(title = paste("PC3 (",percents[3],"%)",sep="")),
         zaxis = list(title = paste("PC2 (",percents[2],"%)",sep=""))),
         showlegend = FALSE
  )
#htmlwidgets::saveWidget(pca3D, file.path(resultsFolder, "pca3D.html"), selfcontained = T)
htmlwidgets::saveWidget(pca3D, 
                        file.path(resultsFolder, paste0("pca3D_",format(Sys.time(),"%Y%m%d"),".html")))

