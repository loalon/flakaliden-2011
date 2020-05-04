#generate megadata
library("DESeq2")
# library(vegan)
# library(ape)
library(dplyr)

#
projectFolder <- "/mnt/picea/projects/spruce/vhurry/fungi-flakaliden-2011"
deFolder <- file.path(projectFolder,"DE")
dataFolder <- file.path(projectFolder,"data")

#guildFile <- file.path(dataFolder, "fungi.guilds.tsv")
rawFile <- file.path(dataFolder, "kingdom.Fungi.2011.raw.tsv")
metaFile <- file.path(dataFolder, "meta.tsv")
taxonomyFile <- file.path(dataFolder, "gene_taxonomy_arc.tsv")

data <- read.table(rawFile, header = T, sep='\t', comment.char = "", quote="", row.names = 1, stringsAsFactors = F)
colnames(data)<-gsub("X", "", colnames(data))

meta <- read.table(metaFile, header = TRUE, sep = '\t', comment.char = "", stringsAsFactors = F)

meta$Treatment <- ifelse(meta$Condition == 'control', 'C', 'F')
meta$Week <- meta$Sampling.date..week..
meta$Time <- substr(meta$Sampling.date..week.., 2,3)
meta$group <- paste0(meta$Week, gsub("T","F",meta$Plot))
meta <- meta[order(meta$group), ]
meta$group <- substr(meta$group, 1, 4)
meta$group <- as.factor(meta$group)
meta$Treatment <- as.factor(meta$Treatment)

#' #remove non usefull samples
data <- data[,colnames(data) %in% meta$Sample.ID]

#' #Order them by group
data <- data[, match(meta$Sample.ID, colnames(data))]

#' change the names of the samples to group
colnames(data) <- paste0(meta$Sampling.date..week..,meta$Plot)

#control data
controlData <- data[,grepl('C',(colnames(data)))]

# fertilised data
fertilisedData <- data[,grepl('F',(colnames(data)))]

#data <- data[(rowSums(data) > 0),]
#guilds <- read.table(guildFile, header = TRUE, sep = '\t', comment.char = "", stringsAsFactors = F)
#
####
#####CREATING MEGADATA
megaData <- data
megaData$gene <- rownames(megaData)
megaData <- megaData %>%  select(gene, everything())

taxo <- read.table(taxonomyFile, header = TRUE, sep = '\t', comment.char = "", stringsAsFactors = FALSE)
#data for funGuild
#taxo_df <- data.frame(gene = taxo$gene, taxonomy =apply(taxo[,2:8],1,paste,  collapse = ";" ) )
#write.table(taxo_df, file=file.path(dataFolder, "taxo4funguild.tsv"), quote=F, row.names=F, sep='\t')

megaData <- left_join(megaData, taxo, by="gene")
