##################################################
## Project: Flakaliden 2011
## Script purpose: Dimensin reduction with PHATE
## and funGUILD analysis
## Date: 200114
## Author: Alonso Serrano
##################################################

install.packages("phateR")

library("phateR")
library(ggplot2)
library(viridis)

setwd("C:/projects/phate/")

#data <- read.delim('fullCountsMatrix.tsv', header = T, row.names = 1, stringsAsFactors = F, sep='\t')
#save(data, file="data.RData")

load("C:/projects/phate/guildVstData_200121.RData")

data <- combinedGuildVstData
controlData <- controlGuildVstData
fertilisedData <- fertilisedGuildVstData

#data.phate <- phate(data)
data.control.phate <- phate(controlData)
data.fertilised.phate <- phate(fertilisedData)
data.combined.phate <- phate(data)

data.phate <- data.fertilised.phate

summary(data.phate)

group <- substr(colnames(data),2,4)
week <- substr(colnames(data),2,3)
treatment <- substr(colnames(data),4,4)

palette(rainbow(10))
plot(data.combined.phate, col = data.phate$branches)

modes <- c("Pathotroph","Pathotroph-Saprotroph",
           "Pathotroph-Saprotroph-Symbiotroph","Pathotroph-Symbiotroph",
           "Saprotroph","Saprotroph-Pathotroph-Symbiotroph",
           "Saprotroph-Symbiotroph","Symbiotroph")
modeColors=c("#DB3832", "#E26D38",
             "#3F2313", "#764693",
             "#FEDD4F", "#00FF00",
             "#459E55","#264294")
names(modeColors) <- modes

# General phate plot
ggplot(data.combined.phate, aes(x=PHATE1, y=PHATE2, color=data.phate$branches)) +
  geom_point(size=0.01) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 

# using trophic mode
Trophic.Mode <- localguilds$Trophic.Mode
ggplot(data.phate, aes(x=PHATE1, y=PHATE2, color=Trophic.Mode)) +
  geom_point(size=0.01) +
  scale_color_manual(values=modeColors) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 

# using trophic mode
Guild <- localguilds$Guild
ggplot(data.phate, aes(x=PHATE1, y=PHATE2, color=Guild)) +
  geom_point(size=0.01) +
  #scale_color_manual(values=modeColors) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
 theme(legend.position = "none")


library(stringr)
taxos <- str_split_fixed(localguilds$taxonomy, ";", n=Inf)
colnames(taxos) <- c("superkingdom",	"kingdom",	"phylum",	"class",	"order",	"family",	"genus",	"species")

fullguilds <- cbind(localguilds, taxos)

# using taxa
Taxo <- fullguilds$genus

modeColors <- rep("#000000", length(unique(Taxo)))
names(modeColors) <- unique(Taxo)
modeColors["Cortinarius"]<- "#A4D171"
modeColors["Cenococcum"]<-  "#009AA9"
modeColors["Piloderma"]<-  "#60B266"
modeColors["Thelephora"]<-  "#C9BBDC"
modeColors["Hygrophorus"]<-  "#A084BD"
modeColors["Hyaloscypha"]<-  "#9BAEC6"

#reset modeTrasnperency to all 0
modeTransparency <- vector(mode = "integer", length(unique(Taxo)))
names(modeTransparency) <- unique(Taxo)

#look 4 specific
query <- 'Hyaloscypha'
#query <- c("Cortinarius","Cenococcum","Piloderma","Thelephora","Hygrophorus","Hyaloscypha")
modeTransparency[query] <- 1

g <- ggplot(data.phate, aes(x=PHATE1, y=PHATE2, color=Taxo, alpha= Taxo)) +
  geom_point(size=0.01) +
  scale_alpha_manual(values = modeTransparency) +
  scale_color_manual(values=modeColors) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(legend.position = "none")

png(paste0("phate_fertilised_",query,"_genus",".png"), 
    units="px", width=1920, height=1080, res=300)
g
dev.off()

#####################
####################
#Paint by species

species <- fullguilds$species

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

querySpecies <- c("Hyaloscypha bicolor","Hyaloscypha variabilis","Piloderma croceum","Unclassified.Piloderma",
                  "Fibularhizoctonia sp. CBS 109695","Unclassified.Fibularhizoctonia","Cortinarius glaucopus",
                  "Hebeloma cylindrosporum","Unclassified.Cortinarius","Hygrophorus russula","Unclassified.Hygrophorus",
                  "Mycena galopus","Unclassified.Mycena","Rhodocollybia butyracea","Unclassified.Lactarius", "Hygrophorus russula",
                  "Thelephora ganbajun","Unclassified.Thelephora", "Acephala macrosclerotiorum","Hydnum rufescens","Mycena alexandri","Mycena amicta",
                  "Mycena epipterygia","Septobasidium sp.","Unclassified.Septobasidium")



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
  
  png(paste0("phate_",query,"_species",".png"), 
      units="px", width=1920, height=1080, res=300)
  print(g)
  dev.off()
  
})

#reset modeTrasnperency to all 0
speciesTransparency <- vector(mode = "integer", length(unique(species)))
names(speciesTransparency) <- unique(species)

#look 4 specific
query <- 'Hyaloscypha bicolor'
#query <- c("Cortinarius","Cenococcum","Piloderma","Thelephora","Hygrophorus","Hyaloscypha")
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

png(paste0("phate_",query,"_species",".png"), 
    units="px", width=1920, height=1080, res=300)
g
dev.off()


#######################
########################

#look 4 specific
query <- combined.res.filter2$`19Fvs19C`

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
#time
##########################

timeColors <- c("#e6e9bd","#b4d7bc","#7dc6bb","#1cb7ba","#00a7b9","#008b9a",
                "#849fba","#a181a6","#c76b97","#ee2e82","#8e5766","#8a7b73",
                "#869b7f","#42885e","#1aa64a","#75b443","#9bbe3b","#e1d51d",
                "#ffe200")
names(timeColors) <- c("19", "20", "21", "22", "23", "24", "25", "26", "28", "31", "32", "34" ,"35", "36", "37" ,"38", "39", "40", "41")



parData <-data.frame(gene = rownames(data.combined.phate$embedding), stringsAsFactors = F) 

parData <- left_join(parData, megaData, by="gene")

g<- ggplot(data.combined.phate, aes(x=PHATE1, y=PHATE2, color=parData$time)) +
  geom_point(size=0.01) +
  scale_color_manual(values=timeColors) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(legend.position = "none")

png(paste0("phateGuild_timeAproximation.png"), 
    units="px", width=1920, height=1080, res=300)
print(g)
dev.off()


#######################
########################













######################
# fullguilds[fullguilds$family == 'Cortinariaceae',]$genus
# 
# 
# 
# 
# 
# 
# ggplot(data.phate) +
#   geom_point(aes(PHATE1, PHATE2, size=0.01, color=localguilds$Trophic.Mode))
# # +
# #   labs(color=group) #+
# #scale_color_viridis(option="B")
# 
# ggplot(data.phate) +
#   geom_point(aes(PHATE1, PHATE2, size=5, color=substr(colnames(data),4,4))) +
#   labs(color=substr(colnames(data),4,4)) #+
