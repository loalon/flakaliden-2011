library(here)
library(ggplot2)
library(reshape2)
library(DESeq2)
library(dplyr)

data <- read.table(here("env/tempSpecies.txt"), sep='\t')
data$species <- rownames(data)

melted <-  melt(data, id = c("species"))

ggplot(melted, aes(x = species, y = variable)) + 
  geom_point(aes(fill=variable, size = value), alpha = 0.5) +
  #scale_fill_manual(values = c("#00AFBB", "#E7B800", "#FC4E07", "#86BBD8")) +
  scale_size(range = c(0.5, 12)) +  # Adjust the range of points size
  theme(#legend.key=element_blank(), 
        axis.text.x = element_text(colour = "black", size =8, face = "bold", angle = 90, vjust = 0.3, hjust = 1), 
        axis.text.y = element_text(colour = "black", face = "bold", size = 8), 
        legend.text = element_text(size = 8, face ="bold", colour ="black"), 
        legend.title = element_text(size = 8, face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
        legend.position = "right")


projectFolder <- "/mnt/picea/projects/spruce/vhurry/fungi-flakaliden-2011"
networkDir <- file.path(projectFolder,"networks")
deDir <- file.path(projectFolder,"DE")
dataDir <- file.path(projectFolder,"data")

load(file.path(deDir, "ddsTreatment_ko.RData"))
combined.res <- results(ddsFungiTreatment)

clusters <- read.table(file.path(networkDir, "combined/cluster/InfomapClusters.tsv"), sep='\t', header=T)

res <- sapply(unique(clusters$cluster), function(x){
  kos <- clusters[clusters$cluster==x,]$gene
  
  daMean <- mean(combined.res[rownames(combined.res) %in% kos,]$log2FoldChange)
  
  res <- "ND"
  if (daMean > 0 ) {
    res <- "NE"
  }
  res
  
})

kos <- read.table(file.path(dataDir, "koTaxo.tsv"), sep='\t', header=T)

cluster_hemis <- left_join(clusters, data.frame(cluster=names(res), hemisphere=res), by="cluster" )

totalData <- left_join(cluster_hemis, kos, by=c("gene" = "ko"))

finalTable <- melt(table(totalData$species, totalData$hemisphere))

soi <- c("Piloderma croceum",
         "Thelephora ganbajun",
         "Cortinarius glaucopus",
         "Unclassified.Piloderma",
         "Cenococcum geophilum",
         "Hyaloscypha bicolor",
         "Unclassified.Mycena",
         "Hyaloscypha variabilis",
         "Mycena galopus",
         "Hebeloma cylindrosporum",
         "Unclassified.Cortinarius",
         "Hygrophorus russula",
         "Unclassified.Hyaloscypha"
)

reFinalTable <- finalTable[finalTable$Var1 %in% soi,]
colours = c("#009444", "#BE9230")

ggplot(reFinalTable, aes(x = Var1, y = Var2)) + 
  geom_point(aes(fill=Var2, size = value), alpha = 1, shape = 21) +
  scale_fill_manual(values = colours, guide = FALSE) +
  #scale_size_continuous(limits = c(1, 10000), range = c(0.1,17)) + 
  scale_size(range = c(0.5, 20), breaks = c(1,10,100,500,1000,2000,5000,10000)) +  # Adjust the range of points size
  theme(#legend.key=element_blank(), 
    axis.text.x = element_text(colour = "black", size =8, face = "bold", angle = 90, vjust = 0.3, hjust = 1), 
    axis.text.y = element_text(colour = "black", face = "bold", size = 8), 
    legend.text = element_text(size = 8, face ="bold", colour ="black"), 
    legend.title = element_text(size = 8, face = "bold"), 
    panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
    legend.position = "right")
