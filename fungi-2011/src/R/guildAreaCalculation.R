library(here)
library(stringr)
library(reshape2)
library(matrixStats)
library(ggplot2)
library(scales)
library(dplyr)

projectFolder <- "/mnt/picea/projects/spruce/vhurry/fungi-flakaliden-2011"
guildFolder <- file.path(projectFolder,"funguild")
otufile <- file.path(guildFolder,"juliasOTUs.guilds.txt")
taxos <- read.table(file.path(guildFolder, "Supp_table_15.csv"), header=T, sep=',', stringsAsFactors = F )
taxos$OTU_ID <- rownames(taxos)

guilds <- read.table(otufile, header=T, sep='\t', stringsAsFactors = F)

data <- guilds[,2:109]
rownames(data) <- guilds[,1]

meta <- guilds[c(1, 110:119)]
colnames(meta)[1] <- "OTU_ID"

data.proportion <- prop.table(as.matrix(data), margin = 2)

taxo <- taxos[,c("OTU_ID", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")]

melt_filt2 <- function (mat) {
  melt1 <- melt(mat, varnames = c("OTU_ID", "group"), value.name = "Count")
  melt2 <- melt1[melt1$Count>0,]
  return(melt2)
}
data.melt <- melt_filt2(data.proportion)
data.melt$OTU_ID <- as.character(data.melt$OTU_ID)
data.melt$group <- as.character(data.melt$group)

data.melt2 <- left_join(data.melt, meta, by = "OTU_ID")
data.melt2 <- left_join(data.melt2, taxo, by = "OTU_ID")

## familes
## 
data.family <- aggregate(data.melt2$Count, list(paste(data.melt2$group, data.melt2$Family, sep = "-")), sum)
colnames(data.family) <- c("group_family", "Count")

data.family[,c(3,4)] <- str_split_fixed(data.family$group_family, "-", 2)
data.family <- data.family[,-1]

data.family.matrix <- acast(data.family, V4~V3, value.var = "Count")
data.family.matrix[is.na(data.family.matrix)] <- 0
#Extract top 12 families for plotting
data.family.matrix.top12 <- rownames(data.family.matrix)[order(rowMedians(data.family.matrix), decreasing = TRUE)][1:12]

write.table(data.family.matrix, file="julias_OTU_family.tsv", sep='\t', quote=F, col.names=NA)


## trophic.modes
## 
data.trophic <- aggregate(data.melt2$Count, list(paste(data.melt2$group, data.melt2$Trophic.Mode, sep = "-")), sum)
colnames(data.trophic) <- c("group_trophic", "Count")

data.trophic[,c(3,4)] <- str_split_fixed(data.trophic$group_trophic, "-", 2)
data.trophic <- data.trophic[,-1]

data.trophic.matrix <- acast(data.trophic, V4~V3, value.var = "Count")
data.trophic.matrix[is.na(data.trophic.matrix)] <- 0
#Extract top 12 families for plotting
data.trophic.matrix.top12 <- rownames(data.trophic.matrix)[order(rowMedians(data.trophic.matrix), decreasing = TRUE)][1:12]

write.table(data.trophic.matrix, file="julias_OTU_trophic.tsv", sep='\t', quote=F, col.names=NA)

