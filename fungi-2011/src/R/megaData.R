##################################################
## Project: Flakaliden 2011
## Script purpose: Generate megadata for fungi 2011
## data. Megadata facilitates the downstream analysis
## by adding the information into a single RData file.
## 
## Date: 200420
## Author: Alonso Serrano
## Mail: alonso.serrano@slu.se
##       bioinformatics@loalon.com
##################################################
##################################################

library(dplyr)

# load data
projectFolder <- "/mnt/picea/projects/spruce/vhurry/fungi-flakaliden-2011"
deFolder <- file.path(projectFolder,"DE")
dataFolder <- file.path(projectFolder,"data")

rawFile <- file.path(dataFolder, "kingdom.Fungi.2011.raw.tsv")
metaFile <- file.path(dataFolder, "meta.tsv")

data <- read.table(rawFile, header = T, sep='\t', comment.char = "", quote="", row.names = 1, stringsAsFactors = F)
colnames(data)<-gsub("X", "", colnames(data))

# load metadata and cleanup
meta <- read.table(metaFile, header = TRUE, sep = '\t', comment.char = "", stringsAsFactors = F)

meta$Treatment <- ifelse(meta$Condition == 'control', 'C', 'F')
meta$Week <- meta$Sampling.date..week..
meta$Time <- substr(meta$Sampling.date..week.., 2,3)
meta$group <- paste0(meta$Week, gsub("T","F",meta$Plot))
meta <- meta[order(meta$group), ]
meta$group <- substr(meta$group, 1, 4)
meta$group <- as.factor(meta$group)
meta$Treatment <- as.factor(meta$Treatment)

# remove non useful samples
data <- data[,colnames(data) %in% meta$Sample.ID]

# Order them by group
data <- data[, match(meta$Sample.ID, colnames(data))]

# change the names of the samples to group
colnames(data) <- paste0(meta$Sampling.date..week..,meta$Plot)

#control data
controlData <- data[,grepl('C',(colnames(data)))]

# fertilised data
fertilisedData <- data[,grepl('F',(colnames(data)))]

# create the megadata
megaData <- data
megaData$gene <- rownames(megaData)
megaData <- megaData %>%  select(gene, everything())

taxonomyFile <- file.path(dataFolder, "gene_taxonomy_update.tsv")
taxo <- read.table(taxonomyFile, header = TRUE, sep = '\t', comment.char = "", stringsAsFactors = FALSE)
megaData <- left_join(megaData, taxo, by="gene")

# add guild information
guildFile <- file.path(dataFolder, "taxo4funguild.guilds.txt")
guilds <- read.table(guildFile, header = TRUE, sep = '\t', comment.char = "", stringsAsFactors = F)
megaData <- left_join(megaData, guilds, by="gene")

# add DE information
load(file.path(deFolder, "resFungiTreatment.RData"))
treatmentDF <- data.frame(gene = rownames(resFungiTreatment), 
                          treatment_lfc= resFungiTreatment[,"log2FoldChange"], 
                          treatment_padj = resFungiTreatment[,"padj"],
                          stringsAsFactors = F)

megaData <- left_join(megaData, treatmentDF, by="gene")

# save a TSV file and a Rdata file for future uses
write.table(megaData, file=file.path(dataFolder, "megaData2.tsv"), quote=F, row.names=F, sep='\t')
save(megaData, file=file.path(dataFolder, "megaData2.RData"))

writeLines(capture.output(sessionInfo()), file.path(dataFolder,
                                                    paste0("megaData_sessionInfo_",Sys.Date(), ".txt"))
)
