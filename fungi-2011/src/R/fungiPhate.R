##################################################
## Project: Flakaliden 2011
## Script purpose: Dimensin reduction with PHATE
## Date: 200114
## Author: Alonso Serrano
##################################################

# Phate script

install.packages("phateR")

library("phateR")
library(ggplot2)
library(dplyr)
library(DESeq2)
library(scales)

setwd("C:/projects/phate/")


#obta
load("C:/projects/phate/vsd.Group.QA.RData")
load("C:/projects/phate/megaData.RData")
load("C:/projects/phate/combined_res_filter2_200221.RData")

# execute for faster loading

source("Rtoolbox/src/plotEigenGene.R")
source("Rtoolbox/src/plotEnrichedTreemap.R")
source("UPSCb-common/src/R/gopher.R")

data <- vsd.Group.matrix
controlData <- data[,grepl("C", colnames(data))]
fertilisedData <- data[,grepl("F", colnames(data))]

data.combined.phate <- phate(data, n.jobs = -2)
# execute when done
#save(data.combined.phate, "data_combined_phate.RData")
data.control.phate <- phate(controlData, n.jobs = -2)
data.fertilised.phate <- phate(fertilisedData, n.jobs = -2)

### AESTHETICS

colorRampPalette( c( "#BC9241", "#949598","#039453" ) )( 15 )

species <- megaData$species

speciesColors <- rep("#000000", length(unique(species)))
names(speciesColors) <- unique(species)

speciesColors["Piloderma croceum"]<-  "#60B266"
speciesColors["Cortinarius glaucopus"]<- "#FEF65B"
speciesColors["Cenococcum geophilum"]<-  "#009AA9"

speciesColors["Archaeorhizomyces borealis"]<-  "#9BAEC6"
speciesColors["Archaeorhizomyces finlayi"]<-  "#009AA9"


#phate is calculated we can transpose the vsd matrix for plotting purposes
vsd.Group.matrix <- t(vsd.Group.matrix)

# General phate plot
g <- ggplot(data.combined.phate, aes(x=PHATE1, y=PHATE2)) +
  geom_point(size=0.01) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 

png(paste0("phateCombinedGeneral.png"), 
    units="px", width=1920, height=1080, res=300)
print(g)
dev.off()

#########################
#DE results
##########################
lapply(names(combined.res.filter2), function(x){
  query <- combined.res.filter2[[x]]
  
  
  upGenes <- intersect(rownames(query[query$log2FoldChange>0,]), rownames(data.combined.phate$embedding))
  downGenes <- intersect(rownames(query[query$log2FoldChange<0,]), rownames(data.combined.phate$embedding))
  #query <- c("Cortinarius","Cenococcum","Piloderma","Thelephora","Hygrophorus","Hyaloscypha")
  
  geneColors <- rep("#000000", length(rownames(data.combined.phate$embedding)))
  names(geneColors) <- rownames(data.combined.phate$embedding)
  geneColors[upGenes]<- "#FF0000"
  geneColors[downGenes]<- "#0000FF"
  
  g <- ggplot(data.combined.phate, aes(x=PHATE1, y=PHATE2, color=rownames(data.combined.phate$embedding))) +
    geom_point(size=0.01) +
    
    scale_color_manual(values=geneColors) +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()) +
    theme(legend.position = "none")
  
  png(paste0("phate_upanddown_",x,".png"), 
      units="px", width=1920, height=1080, res=300)
  print(g)
  dev.off()
})

#########################
#treatment
##########################

treatmentColors <- c("#009444","#BE9230")
names(treatmentColors) <- c("C", "F")

g<- ggplot(data.combined.phate, aes(x=PHATE1, y=PHATE2, color=megaData$treatment)) +
  geom_point(size=0.01) +
  scale_color_manual(values=treatmentColors) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(legend.position = "none")

png(paste0("phate_treatmentAproximation.png"), 
    units="px", width=1920, height=1080, res=300)
print(g)
dev.off()

#########################
#time
##########################

timeColors <- c("#e6e9bd","#b4d7bc","#7dc6bb","#1cb7ba","#00a7b9","#008b9a",
                "#849fba","#a181a6","#c76b97","#ee2e82","#8e5766","#8a7b73",
                "#869b7f","#42885e","#1aa64a","#75b443","#9bbe3b","#e1d51d",
                "#ffe200")
names(timeColors) <- c("19", "20", "21", "22", "23", "24", "25", "26", "28", "31", "32", "34" ,"35", "36", "37" ,"38", "39", "40", "41")

g<- ggplot(data.combined.phate, aes(x=PHATE1, y=PHATE2, color=megaData$time)) +
  geom_point(size=0.01) +
  scale_color_manual(values=timeColors) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(legend.position = "none")

png(paste0("phate_timeAproximation.png"), 
    units="px", width=1920, height=1080, res=300)
print(g)
dev.off()




#########################
#trophic
##########################

#lets check that the order of megadata is the same one as in the vst data

all(megaData$gene == rownames(data))

# all matches

trophicColors=c("#DB3832", "#E26D38","#3F2313", "#764693","#FEDD4F", "#00FF00","#459E55","#264294", "#000000")
names(trophicColors) <- c("Pathotroph","Pathotroph-Saprotroph","Pathotroph-Saprotroph-Symbiotroph","Pathotroph-Symbiotroph",
                          "Saprotroph","Saprotroph-Pathotroph-Symbiotroph","Saprotroph-Symbiotroph","Symbiotroph", "-")

g<- ggplot(data.combined.phate, aes(x=PHATE1, y=PHATE2, color=megaData$Trophic.Mode)) +
  geom_point(size=0.01) +
  scale_color_manual(values=trophicColors) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(legend.position = "none")

png(paste0("phate_trophicMode.png"), 
    units="px", width=1920, height=1080, res=300)
print(g)
dev.off()

#########################
# genus
##########################
# using taxa
Taxo <- megaData$genus

genusColors <- rep("#000000", length(unique(Taxo)))
names(genusColors) <- unique(Taxo)
genusColors["Cortinarius"]<- "#A4D171"
genusColors["Cenococcum"]<-  "#009AA9"
genusColors["Piloderma"]<-  "#60B266"
genusColors["Thelephora"]<-  "#C9BBDC"
genusColors["Hygrophorus"]<-  "#A084BD"
genusColors["Hyaloscypha"]<-  "#9BAEC6"

queryGenus <- c("Cortinarius","Cenococcum","Piloderma","Thelephora","Hygrophorus","Hyaloscypha")

lapply(queryGenus, function(query){
  
  #reset genusTrasnperency to all 0
  genusTransparency <- vector(mode = "integer", length(unique(Taxo)))
  names(genusTransparency) <- unique(Taxo)
  genusTransparency[query] <- 1
  
  g <- ggplot(data.combined.phate, aes(x=PHATE1, y=PHATE2, color=Taxo, alpha= Taxo)) +
    geom_point(size=0.01) +
    scale_alpha_manual(values = genusTransparency) +
    scale_color_manual(values=genusColors) +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()) +
    theme(legend.position = "none")
  
  png(paste0("phateCombined_genus_",query,".png"), 
      units="px", width=1920, height=1080, res=300)
  print(g)
  dev.off()
})

#########################
# species
##########################
species <- megaData$species

speciesColors <- rep("#000000", length(unique(species)))
names(speciesColors) <- unique(species)
speciesColors["Hyaloscypha bicolor"]<- "#A4D171"
speciesColors["Hyaloscypha variabilis"]<-  "#009AA9"
speciesColors["Piloderma croceum"]<-  "#60B266"
speciesColors["Unclassified.Piloderma"]<-  "#C9BBDC"
speciesColors["Fibularhizoctonia sp. CBS 109695"]<-  "#A084BD"
speciesColors["Unclassified.Fibularhizoctonia"]<-  "#9BAEC6"

speciesColors["Cortinarius glaucopus"]<- "#A4D171"
speciesColors["Hebeloma cylindrosporum"]<-  "#009AA9"
speciesColors["Unclassified.Cortinarius"]<-  "#60B266"
speciesColors["Hygrophorus russula"]<-  "#C9BBDC"
speciesColors["Unclassified.Hygrophorus"]<-  "#A084BD"
speciesColors["Mycena galopus"]<-  "#9BAEC6"

speciesColors["Unclassified.Mycena"]<- "#A4D171"
speciesColors["Rhodocollybia butyracea"]<-  "#009AA9"
speciesColors["Unclassified.Lactarius"]<-  "#60B266"
speciesColors["Hygrophorus russula"]<-  "#C9BBDC"
speciesColors["Thelephora ganbajun"]<-  "#A084BD"
speciesColors["Unclassified.Thelephora"]<-  "#9BAEC6"

speciesColors["Acephala macrosclerotiorum"]<- "#A4D171"
speciesColors["Hydnum rufescens"]<-  "#009AA9"
speciesColors["Mycena alexandri"]<-  "#60B266"
speciesColors["Mycena amicta"]<-  "#C9BBDC"
speciesColors["Mycena epipterygia"]<-  "#A084BD"
speciesColors["Septobasidium sp."]<-  "#9BAEC6"
speciesColors["Unclassified.Septobasidium"]<-  "#9BAEC6"


speciesColors["Archaeorhizomyces borealis"]<-  "#9BAEC6"
speciesColors["Archaeorhizomyces finlayi"]<-  "#009AA9"



querySpecies <- c("Hyaloscypha bicolor","Hyaloscypha variabilis","Piloderma croceum","Unclassified.Piloderma",
                  "Fibularhizoctonia sp. CBS 109695","Unclassified.Fibularhizoctonia","Cortinarius glaucopus",
                  "Hebeloma cylindrosporum","Unclassified.Cortinarius","Hygrophorus russula","Unclassified.Hygrophorus",
                  "Mycena galopus","Unclassified.Mycena","Rhodocollybia butyracea","Unclassified.Lactarius", "Hygrophorus russula",
                  "Thelephora ganbajun","Unclassified.Thelephora", "Acephala macrosclerotiorum","Hydnum rufescens","Mycena alexandri","Mycena amicta",
                  "Mycena epipterygia","Septobasidium sp.","Unclassified.Septobasidium")

querySpecies <- c("Archaeorhizomyces borealis", "Archaeorhizomyces finlayi")

querySpecies <- c("Hyaloscypha variabilis", "Piloderma croceum")

lapply(querySpecies, function(query){
  #query <- 'Hyaloscypha bicolor'
  print(query)
  speciesTransparency <- vector(mode = "integer", length(unique(species)))
  names(speciesTransparency) <- unique(species)
  
  speciesTransparency[query] <- 1
  
  g <- ggplot(data.combined.phate, aes(x=PHATE1, y=PHATE2, color=species, alpha= species)) +
    geom_point(size=0.01) +
    scale_alpha_manual(values = speciesTransparency) +
    scale_color_manual(values=speciesColors) +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()) +
    theme(legend.position = "none")
  
  png(paste0("phateCombined_species_",query,".png"), 
      units="px", width=1920, height=1080, res=300)
  print(g)
  dev.off()
  
})



tempPCA <- prcomp(t(vsd.Group.matrix))

tempDF <- data.frame(PC1= tempPCA$x[,1], PC2 = tempPCA$x[,2], PC3 = tempPCA$x[,3], group = rownames(tempPCA$x),
                     color=megaData$treatment )
ceno.percents <- round(summary(ceno.pca)$importance[2,]*100)

p <- plot_ly(cenoDF, x = ~PC1, y = ~PC2, 
             mode = 'markers', color=~color, colors=c("#009444","#BE9230"),
             text = ~group, marker = list(opacity = 0.5, sizemode = 'diameter')) %>%
  
  #add_markers() %>%
  layout(title = title,
         xaxis = list(title = paste("PC1 (",ceno.percents[1],"%)",sep="")),
         yaxis = list(title = paste("PC2 (",ceno.percents[2],"%)",sep=""))
         
  )




#########################
# species specific arms
##########################

# Cortinarius glaucopus AND Unclassified.Cortinarius AND 34
# and
# Cortinarius glaucopus AND Unclassified.Cortinarius AND 36


species <- megaData$species

speciesColors <- rep("#000000", length(unique(species)))
names(speciesColors) <- unique(species)

speciesColors["Piloderma croceum"]<-  "#60B266"
speciesColors["Cortinarius glaucopus"]<- "#FEF65B"
speciesColors["Cenococcum geophilum"]<-  "#009AA9"

#querySpecies <- c("Piloderma croceum", "Cortinarius glaucopus", "Cenococcum geophilum")


querySpecies <- "Cortinarius glaucopus"
queryTime <- "36" #max genes in that time
#queryURL <- "fungi2011_piloderma_croceum"
queryURL <- "fungi2011_cortinarius_glaucopus"

getSpeciesTimePlot <- function(querySpecies, queryTime, queryURL){
  require(dplyr)

  if (queryTime=="") {
    speciesTime <- megaData %>% filter(species==querySpecies)
  } else {
    speciesTime <- megaData %>% filter(time==queryTime, species==querySpecies)  
  }
  
  
  speciesFullGenes <- megaData %>% filter(species==querySpecies) %>% select(gene)
  queryGenes <- speciesTime$gene
  
  genesTransparency <- vector(mode = "integer", length(megaData$gene))
  names(genesTransparency) <- unique(megaData$gene)
  
  genesTransparency[queryGenes] <- 1
  
  
  g <-ggplot(data.combined.phate, aes(x=PHATE1, y=PHATE2,
                                       color=species,
                                       alpha= megaData$gene
                                       )) +
    geom_point(size=0.01) +
    scale_alpha_manual(values = genesTransparency) +
    scale_color_manual(values=speciesColors) +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()) +
    theme(legend.position = "none")



  png(paste0("phateCombined_species_",querySpecies,"_time_",queryTime, ".png"),
      units="px", width=1920, height=1080, res=300)
  print(g)
  dev.off()
 
  enr <- gopher(queryGenes, url=queryURL, task = "go")
  enr2 <- gopher(queryGenes, background= speciesFullGenes$gene, url="fungi2011", task = c("kog", "cog","ko", "ko_pathway") )

  enrTotal <- c(enr, enr2)

  lapply(c("go","kog", "cog","ko", "ko_pathway" ), function(x){
    if (!is.null(enrTotal[[x]]) ) {
      print(x)
      png(paste0("phateCombined_species_",querySpecies,"_time_",queryTime, "_",x,"_treemap.png"),
          units="px", width=1920, height=1080, res=300)
      plotEnrichedTreemap(enrTotal, enrichment = x, namespace = "none", clusterColor = speciesColors[querySpecies])
      dev.off()

    }


  })

  enr3 <- gopher(queryGenes, url="fungi2011", task = c("kog", "cog","ko", "ko_pathway") )

  lapply(c("kog", "cog","ko", "ko_pathway" ), function(x){
    if (!is.null(enr3[[x]]) ) {
      print(x)
      png(paste0("phateCombined_species_",querySpecies,"_time_",queryTime, "_",x,"_vs_AllGenes_treemap.png"),
          units="px", width=1920, height=1080, res=300)
      plotEnrichedTreemap(enr3, enrichment = x, namespace = "none", clusterColor = speciesColors[querySpecies])
      dev.off()

    }
  })

  
  png(paste0("phateCombined_species_",querySpecies,"_time_",queryTime, "_eigengene.png"),
      units="px", width=1920, height=1080, res=300)
  
  
  p <- plotEigengene(vsd.Group.matrix, 
                            queryGenes, 
                            substr(rownames(vsd.Group.matrix),4,4), 
                            as.integer(substr(rownames(vsd.Group.matrix),2,3)), 
                            timeUnits = "Time",
                            inverse = F, title = "", noGrid = T, colors = c("#009444", "#BE9230"),
                            noLegend=T, multiline = F, legendTitle="",
                            plotBothProfiles=F) 
  print(p)
    
  dev.off()

}


getSpeciesTimePlot("Cortinarius glaucopus", "36", "fungi2011_cortinarius_glaucopus")
getSpeciesTimePlot("Cortinarius glaucopus", "34", "fungi2011_cortinarius_glaucopus")

getSpeciesTimePlot("Piloderma croceum", "36", "fungi2011_piloderma_croceum")
getSpeciesTimePlot("Piloderma croceum", "39", "fungi2011_piloderma_croceum")
getSpeciesTimePlot("Piloderma croceum", "20", "fungi2011_piloderma_croceum")
getSpeciesTimePlot("Piloderma croceum", "22", "fungi2011_piloderma_croceum")
getSpeciesTimePlot("Piloderma croceum", "26", "fungi2011_piloderma_croceum")
getSpeciesTimePlot("Piloderma croceum", "28", "fungi2011_piloderma_croceum")

#isolate cenococcum
getSpeciesTimePlot("Cenococcum geophilum", "", "fungi2011_cenococcum_geophilum")



############# DE log2foldchage plot

# Gold - BC9241
# Grey - 949598
# Green - 039453

colorRampPalette( c( "#BC9241", "#949598","#039453" ) )( 15 )


treatmentDF <- data.frame(lfc= dds_treatment[,2], stringsAsFactors = F)

treatmentDF$scale <- rescale(treatmentDF$lfc)

g <-ggplot(data.combined.phate, aes(x=PHATE1, y=PHATE2,color=treatmentDF$lfc)) +
  geom_point(size=0.01) +

  # scale_colour_gradient2(low = "#039453", mid = "#949598",
  #                        high = "#BC9241", midpoint = 0, space = "Lab",
  #                        na.value = "black", guide = "colourbar", aesthetics = "colour")  +
  scale_colour_gradientn(colours=c("#BC9241", "#B6924D",
                                   "grey", 

                                   "#17945C", "#039453"),
                         values=c(1, 0.68,  0.55, 0.42, 0)
                         ) +
  
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(legend.position = "none")


##########################################################################################
## species and treatment
#########################################################################################

getSpeciesTreatmentPlot <- function(querySpecies, queryTreatment, queryURL, lfc= 0.5, padj=0.99){
  require(dplyr)
  
  if (queryTreatment == 'C')
    speciesTreat <- megaData %>% filter(treatment_lfc< (-1*lfc), species==querySpecies)
  else if (queryTreatment == 'F')
    speciesTreat <- megaData %>% filter(treatment_lfc> lfc, species==querySpecies)
  else if (queryTreatment == 'None')
    speciesTreat <- megaData %>% filter(treatment_lfc< lfc & treatment_lfc> (-1*lfc), species==querySpecies)  
  
  speciesFullGenes <- megaData %>% filter(species==querySpecies) %>% select(gene)
  
  queryGenes <- speciesTreat$gene
  print(querySpecies)
  print(length(queryGenes))
  
  genesTransparency <- vector(mode = "integer", length(megaData$gene))
  names(genesTransparency) <- unique(megaData$gene)
  
  genesTransparency[queryGenes] <- 1
  
  
  g <-ggplot(data.combined.phate, aes(x=PHATE1, y=PHATE2,
                                      color=species,
                                      alpha= megaData$gene
  )) +
    geom_point(size=0.01) +
    scale_alpha_manual(values = genesTransparency) +
    scale_color_manual(values=speciesColors) +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()) +
    theme(legend.position = "none")
  
  

  png(paste0("phateCombined_species_",querySpecies,"_treatment_",queryTreatment, ".png"),
      units="px", width=1920, height=1080, res=300)
  print(g)
  dev.off()

  # enr <- gopher(queryGenes, url=queryURL, task = "go")
  # enr2 <- gopher(queryGenes, background= speciesFullGenes$gene, url="fungi2011", task = c("kog", "cog","ko", "ko_pathway") )
  # 
  # enrTotal <- c(enr, enr2)
  # 
  # lapply(c("go","kog", "cog","ko", "ko_pathway" ), function(x){
  #   if (!is.null(enrTotal[[x]]) ) {
  #     print(x)
  #     png(paste0("phateCombined_species_",querySpecies,"_treatment_",queryTreatment, "_",x,"_treemap.png"),
  #         units="px", width=1920, height=1080, res=300)
  #     plotEnrichedTreemap(enrTotal, enrichment = x, namespace = "none", clusterColor = speciesColors[querySpecies])
  #     dev.off()
  # 
  #   }
  # 
  # 
  #  })
  # 
  enr3 <- gopher(queryGenes, url="fungi2011", task = c('cog', 'ko', 'ko_pathway','kog') )
  
  print(enr3)

  lapply(c("kog", "cog","ko", "ko_pathway" ), function(x){
    if (!is.null(enr3[[x]]) ) {
      print(x)
      png(paste0("phateCombined_species_",querySpecies,"_treatment_",queryTreatment, "_",x,"_vs_AllGenes_treemap.png"),
          units="px", width=1920, height=1080, res=300)
      plotEnrichedTreemap(enr3, enrichment = x, namespace = "none", clusterColor = speciesColors[querySpecies])
      dev.off()

    }
  })
  
  
  png(paste0("phateCombined_species_",querySpecies,"_treatment_",queryTreatment, "_eigengene.png"),
      units="px", width=1920, height=1080, res=300)
  
  
  p <- plotEigengene(vsd.Group.matrix, 
                     queryGenes, 
                     substr(rownames(vsd.Group.matrix),4,4), 
                     as.integer(substr(rownames(vsd.Group.matrix),2,3)), 
                     timeUnits = "Time",
                     inverse = F, title = "", noGrid = T, colors = c("#009444", "#BE9230"),
                     noLegend=T, multiline = F, legendTitle="",
                     plotBothProfiles=F) 
  print(p)
  
  dev.off()
}

getSpeciesTreatmentPlot("Piloderma croceum", "C", "fungi2011_piloderma_croceum")
getSpeciesTreatmentPlot("Piloderma croceum", "F", "fungi2011_piloderma_croceum")

getSpeciesTreatmentPlot("Cortinarius glaucopus", "C", "fungi2011_cortinarius_glaucopus")
getSpeciesTreatmentPlot("Cortinarius glaucopus", "F", "fungi2011_cortinarius_glaucopus")
getSpeciesTreatmentPlot("Cortinarius glaucopus", "None", "fungi2011_cortinarius_glaucopus")

getSpeciesTreatmentPlot("Cenococcum geophilum", "C", "fungi2011_cenococcum_geophilum")
getSpeciesTreatmentPlot("Cenococcum geophilum", "F", "fungi2011_cenococcum_geophilum")
getSpeciesTreatmentPlot("Cenococcum geophilum", "None", "fungi2011_cenococcum_geophilum")

getSpeciesTreatmentPlot("Archaeorhizomyces borealis", "C", "fungi2011_cenococcum_geophilum")
getSpeciesTreatmentPlot("Archaeorhizomyces borealis", "F", "fungi2011_cenococcum_geophilum")

getSpeciesTreatmentPlot("Archaeorhizomyces finlayi", "C", "fungi2011_cenococcum_geophilum")
getSpeciesTreatmentPlot("Archaeorhizomyces finlayi", "F", "fungi2011_cenococcum_geophilum")



#################################################
##SPECIES SPECIFIC TREATMENT DE
#################################################

getSpeciesDEPlot <- function(querySpecies, queryTreatment, queryURL, lfc=0.5){
  # if (queryTreatment == 'C')
  #   speciesTreat <- megaData %>% filter(treatment_lfc< (-1*lfc), species==querySpecies)
  # else if (queryTreatment == 'F')
  #   speciesTreat <- megaData %>% filter(treatment_lfc> lfc, species==querySpecies)
  
  speciesFullGenes <- megaData %>% filter(species==querySpecies) %>% select(gene)
  
  queryGenes <- speciesFullGenes$gene
  print(querySpecies)
  print(length(queryGenes))
  
  genesTransparency <- vector(mode = "integer", length(megaData$gene))
  names(genesTransparency) <- unique(megaData$gene)
  
  genesTransparency[queryGenes] <- 1
  
  colorRampPalette( c( "#BC9241", "#949598","#039453" ) )( 15 )
  
  
  treatmentDF <- data.frame(lfc= dds_treatment[,2], stringsAsFactors = F)
  
  treatmentDF$scale <- rescale(treatmentDF$lfc)
  
  png(paste0("phateCombined_species_",querySpecies,"_DEtreatment.png"),
      units="px", width=1920, height=1080, res=300)
  
  g <-ggplot(data.combined.phate, aes(x=PHATE1, y=PHATE2,color=treatmentDF$lfc, alpha=megaData$gene)) +
    geom_point(size=0.01) +
    
    # scale_colour_gradient2(low = "#039453", mid = "#949598",
    #                        high = "#BC9241", midpoint = 0, space = "Lab",
    #                        na.value = "black", guide = "colourbar", aesthetics = "colour")  +
    scale_alpha_manual(values = genesTransparency) +
    scale_colour_gradientn(colours=c("#BC9241", "#B6924D",
                                     "grey", 
                                     
                                     "#17945C", "#039453"),
                           values=c(1, 0.68,  0.55, 0.42, 0)
    ) +
    
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()) +
    theme(legend.position = "none")
  
  print(g)
  dev.off()

}

getSpeciesDEPlot("Piloderma croceum", "C", "fungi2011_piloderma_croceum")
getSpeciesDEPlot("Cortinarius glaucopus", "C", "fungi2011_piloderma_croceum")
getSpeciesDEPlot("Cenococcum geophilum", "C", "fungi2011_piloderma_croceum")
getSpeciesDEPlot("Archaeorhizomyces borealis", "C", "fungi2011_piloderma_croceum")
getSpeciesDEPlot("Archaeorhizomyces finlayi", "C", "fungi2011_piloderma_croceum")

getSpeciesTreatmentPlot("Archaeorhizomyces borealis", "F", "fungi2011_cenococcum_geophilum")

getSpeciesTreatmentPlot("Archaeorhizomyces finlayi", "C", "fungi2011_cenococcum_geophilum")

dim(megaData2 %>% filter(species=="Piloderma croceum") %>% select(gene))
dim(megaData2 %>% filter(species=="Piloderma croceum", treatment_lfc>0.5) %>% select(gene))
dim(megaData2 %>% filter(species=="Piloderma croceum", treatment_lfc<(-1*0.5)) %>% select(gene))

dim(megaData2 %>% filter(species=="Cortinarius glaucopus") %>% select(gene))
dim(megaData2 %>% filter(species=="Cortinarius glaucopus", treatment_lfc>0.5) %>% select(gene))
dim(megaData2 %>% filter(species=="Cortinarius glaucopus", treatment_lfc<(-1*0.5)) %>% select(gene))

dim(megaData2 %>% filter(species=="Cenococcum geophilum") %>% select(gene))
dim(megaData2 %>% filter(species=="Cenococcum geophilum", treatment_lfc>0.5) %>% select(gene))
dim(megaData2 %>% filter(species=="Cenococcum geophilum", treatment_lfc<(-1*0.5)) %>% select(gene))

dim(megaData %>% filter(species=="Archaeorhizomyces borealis") %>% select(gene))
dim(megaData %>% filter(species=="Archaeorhizomyces finlayi") %>% select(gene))


#################################################
##SPECIES SPECIFIC CLUSTER DE
#################################################

getSpeciesModulesPlot <- function(querySpecies, clusters, col_vector){
  library(RColorBrewer)
  #speciesFullGenes <- megaData %>% filter(species==querySpecies) %>% select(gene)
  
  #queryGenes <- speciesFullGenes$gene
  
  genesTransparency <- vector(mode = "integer", length(megaData$gene))
  names(genesTransparency) <- unique(megaData$gene)
  
  genesTransparency[clusters$gene] <- 1
  
  #colorRampPalette( c( "#BC9241", "#949598","#039453" ) )( 15 )
  
  megaDataTemp <- left_join(megaData, clusters, by="gene")


  
  megaDataTemp$cluster <- ifelse(is.na(megaDataTemp$cluster), "noCluster", megaDataTemp$cluster)
  
  moduleColors <- rep("white",length(unique(megaDataTemp$cluster)))
  
  names(moduleColors) <- unique(megaDataTemp$cluster)
  

  
  topClusters <- names(sort(table(clusters$cluster), decreasing = T)[1:max(20,length(unique(clusters$cluster)))])
  #moduleColors[topClusters] <- rainbow(length(topClusters))
  
  # n <- length(topClusters)*3
  # qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  # col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  # pie(rep(1,n), col=sample(col_vector, n))
  
  # moduleColors[topClusters] <- col_vector[1:length(topClusters)]
  moduleColors[topClusters] <- col_vector[1:length(topClusters)]
  
  g <-ggplot(data.combined.phate, aes(x=PHATE1, y=PHATE2,color=megaDataTemp$cluster, alpha=megaData$gene)) +
    geom_point(size=0.01) +
    

  #   g <- ggplot(df, aes(x = x, y = y, group = group))
  # 
  # g <- g + geom_line(aes(colour = group))
  # 
  # g <- g + geom_point(aes(colour = group), alpha = 0.8)
    # scale_colour_gradient2(low = "#039453", mid = "#949598",
    #                        high = "#BC9241", midpoint = 0, space = "Lab",
    #                        na.value = "black", guide = "colourbar", aesthetics = "colour")  +
    scale_alpha_manual(values = genesTransparency) +
    scale_colour_manual(values = moduleColors) +
    guides(alpha=FALSE)+
    #scale_colour_gradientn(colours=rainbow(30)) +
    # scale_colour_gradientn(colours=c("#BC9241", "#B6924D",
    #                                  "grey", 
    #                                  
    #                                  "#17945C", "#039453"),
    #                        values=c(1, 0.68,  0.55, 0.42, 0)
    # ) +
    
    theme(axis.line = element_line(colour = "black"),
          legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()) 
  
  png(paste0("phateCombined_species_",querySpecies,"_top20ClustersLevel1.png"),
   units="px", width=1920, height=1080, res=300)
  print(g)
  dev.off()
  
}

cenococcumInfomapClusters <- read.delim("C:/projects/phate/cenococcum_geophilum/infomapClustersLevel1.tsv", stringsAsFactors=FALSE)
cortinariusInfomapClusters <- read.delim("C:/projects/phate/cortinarius_glaucopus/infomapClustersLevel1.tsv", stringsAsFactors=FALSE)
hbicolorInfomapClusters <- read.delim("C:/projects/phate/hyaloscypha_bicolor/infomapClustersLevel1.tsv", stringsAsFactors=FALSE)
hvariabilisInfomapClusters <- read.delim("C:/projects/phate/hyaloscypha_variabilis/infomapClustersLevel1.tsv", stringsAsFactors=FALSE)
pilodermaInfomapClusters <- read.delim("C:/projects/phate/piloderma_croceum/infomapClustersLevel1.tsv", stringsAsFactors=FALSE)

getSpeciesModulesPlot("cenococcum_geophilum", cenococcumInfomapClusters, newPalette)

getSpeciesModulesPlot("cortinarius_glaucopus", cortinariusInfomapClusters, newPalette)

getSpeciesModulesPlot("hyaloscypha_bicolor", hbicolorInfomapClusters, newPalette)

getSpeciesModulesPlot("hyaloscypha_variabilis", hvariabilisInfomapClusters, newPalette)

getSpeciesModulesPlot("piloderma_croceum", pilodermaInfomapClusters, newPalette)

###########################
# SAVE SESSION INFO
# All phate calculations used 3.6.8
writeLines(capture.output(sessionInfo()), paste0("sessionInfo_",Sys.Date(),".txt"))
save.image(paste0("session_",Sys.Date(),".RData"))

