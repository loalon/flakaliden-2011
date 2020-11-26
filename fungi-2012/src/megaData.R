##################################################
## Project: Flakaliden 2011
## Script purpose: Generate megadata for fungi 2012
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
projectFolder <- "/mnt/picea/projects/spruce/vhurry/fungi-flakaliden-2012"
deFolder <- file.path(projectFolder,"DE")
dataFolder <- file.path(projectFolder,"data")

rawFile <- file.path(dataFolder, "Assembly_2012.raw.tsv")
metaFile <- file.path(dataFolder, "meta.tsv")

data <- read.table(rawFile, header = T, sep='\t', comment.char = "", quote="", row.names = 1, stringsAsFactors = F)

# load metadata and cleanup
meta <-read.table(metaFile, header = T, sep=';', comment.char = "", stringsAsFactors = F)

meta$date[meta$date == 'Early_June'] <- '23'
meta$date[meta$date == 'Late_June'] <- '25'
meta$date[meta$date == 'August'] <- '32'
meta$date[meta$date == 'October'] <- '41'

meta$treatment <- ifelse(meta$treatment=="25_year", "Fertilised", meta$treatment)
meta <- meta[meta$treatment != "5_year",]
meta$treatment <- factor(meta$treatment, levels = c("Control", "Fertilised") )

meta$group <- paste0("W",meta$date, meta$plot.1, "P", meta$plot)
meta$group <- gsub("N","F",meta$group)

meta$group <- as.factor(meta$group)
meta$treatment <- as.factor(meta$treatment)
meta$date <- as.factor(meta$date)

# remove non useful samples
data <- data[,colnames(data) %in% meta$SciLifeID]

# Order them by group
data <- data[, match(meta$SciLifeID, colnames(data))]

# change the names of the samples to group
colnames(data) <- meta$group

# control data
controlData <- data[,grepl('C',(colnames(data)))]

# fertilised data
fertilisedData <- data[,grepl('F',(colnames(data)))]

# create the megadata
megaData <- data
megaData$gene <- rownames(megaData)
megaData <- megaData %>%  select(gene, everything())

taxonomyFile <- file.path(dataFolder, "gene_taxonomy.tsv")
taxo <- read.table(taxonomyFile, header = TRUE, sep = '\t', comment.char = "", stringsAsFactors = FALSE)
megaData <- left_join(megaData, taxo, by="gene")

# save a TSV file and a Rdata file for future uses
write.table(megaData, file=file.path(dataFolder, "megaData.tsv"), quote=F, row.names=F, sep='\t')
save(megaData, file=file.path(dataFolder, "megaData.RData"))

writeLines(capture.output(sessionInfo()), file.path(dataFolder,
                                                    paste0("megaData_sessionInfo_",Sys.Date(), ".txt"))
)

