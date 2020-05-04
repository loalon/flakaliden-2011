# script for proportion plotting of fungi
# Alonso, Simon and Andreas
# 200122
library(here)
library(stringr)
library(reshape2)
library(matrixStats)
library(ggplot2)
library(scales)
library(dplyr)

# thanks to andreas for most of the code
# 
ggformat <- theme_classic()+
  theme(axis.text = element_text(face = "bold", size = 12),
        axis.text.x = element_text(vjust = 0.5, angle = 45),
        axis.line = element_line(size = 1, linetype = "solid"),
        axis.title = element_text(size = 13, face = "bold"),
        #axis.ticks = element_line(size = 1, linetype = "solid"),
        legend.text = element_text(size=13, face = "bold"),
        legend.title = element_text(size = 15, face = "bold"),
        plot.title = element_text(size = 19, face = "bold"),
        strip.text = element_text(size = 17, face = "bold"))


ggformat_pca <- theme_classic()+
  theme(axis.text = element_text(face = "bold", size = 10),
        axis.text.x = element_text(vjust = 1),
        axis.line = element_line(size = 1, linetype = "solid"),
        axis.title = element_text(size = 10, face = "bold"),
        axis.ticks = element_line(size = 1, linetype = "solid"),
        legend.text = element_text(size=13, face = "bold"),
        legend.title = element_text(size = 15, face = "bold"),
        plot.title = element_text(size = 19, face = "bold"),
        strip.text = element_text(size = 17, face = "bold"))

##########################

projectFolder <- "/mnt/picea/projects/spruce/vhurry/fungi-flakaliden-2011"
deFolder <- file.path(projectFolder,"DE")
dataFolder <- file.path(projectFolder,"data")
#guildFile <- file.path(dataFolder, "fungi.guilds.tsv")
rawFile <- file.path(dataFolder, "kingdom.Fungi.2011.raw.tsv")
metaFile <- file.path(dataFolder, "meta.tsv")
taxonomyFile <- file.path(dataFolder, "gene_taxonomy.tsv")

data <- read.table(rawFile, header = T, sep='\t', comment.char = "", quote="", row.names = 1, stringsAsFactors = F)
colnames(data)<-gsub("X", "", colnames(data))

meta <- read.table(metaFile, header = TRUE, sep = '\t', comment.char = "", stringsAsFactors = F)

meta$Treatment <- ifelse(meta$Condition == 'control', 'C', 'F')
meta$Week <- meta$Sampling.date..week..
meta$Time <- substr(meta$Sampling.date..week.., 2,3)
meta$group <- paste0(meta$Week, gsub("T","F",meta$Plot))
meta <- meta[order(meta$group), ]
meta$group <- substr(meta$group, 1, 4)
# meta$group <- as.factor(meta$group)
# meta$Treatment <- as.factor(meta$Treatment)

#' #remove non usefull samples
data <- data[,colnames(data) %in% meta$Sample.ID]

#' #Order them by group
data <- data[, match(meta$Sample.ID, colnames(data))]

#' change the names of the samples to group
colnames(data) <- paste0(meta$Sampling.date..week..,meta$Plot)

# #control data
# controlData <- data[,grepl('C',(colnames(data)))]
# 
# # fertilised data
# fertilisedData <- data[,grepl('F',(colnames(data)))]

#data <- data[(rowSums(data) > 0),]
#guilds <- read.table(guildFile, header = TRUE, sep = '\t', comment.char = "", stringsAsFactors = F, row.names=1)
taxo <- read.table(taxonomyFile, header = TRUE, sep = '\t', comment.char = "", stringsAsFactors = FALSE, row.names=1)
taxo$gene <- row.names(taxo)


data.proportion <- prop.table(as.matrix(data), margin = 2)


melt_filt2 <- function (mat) {
  melt1 <- melt(mat, varnames = c("gene", "group"), value.name = "Count")
  melt2 <- melt1[melt1$Count>0,]
  return(melt2)
}
data.melt <- melt_filt2(data.proportion)
data.melt$gene <- as.character(data.melt$gene)
data.melt$group <- as.character(data.melt$group)

#data.melt2 <- left_join(data.melt, meta, by = "group")
data.melt2 <- left_join(data.melt, taxo, by = "gene")



data.family <- aggregate(data.melt2$Count, list(paste(data.melt2$group, data.melt2$family, sep = "-")), sum)
colnames(data.family) <- c("group_family", "Count")

data.family[,c(3,4)] <- str_split_fixed(data.family$group_family, "-", 2)
data.family <- data.family[,-1]

data.family.matrix <- acast(data.family, V4~V3, value.var = "Count")
data.family.matrix[is.na(data.family.matrix)] <- 0
#Extract top 12 families for plotting
data.family.matrix.top12 <- rownames(data.family.matrix)[order(rowMedians(data.family.matrix), decreasing = TRUE)][1:12]

write.table(data.family.matrix, file="data_fungi_family.tsv", sep='\t', quote=F)
