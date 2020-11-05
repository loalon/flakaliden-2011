### ==============================
## load the necessary libraries
### ==============================
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(vsn))
suppressPackageStartupMessages(library(scatterplot3d))
suppressPackageStartupMessages(library(LSD))
supressPackageStartupMessages(library("ggplot2", lib.loc="/usr/local/lib/R/site-library"))
library("gplots", lib.loc="/usr/local/lib/R/site-library")

source('http://bioconductor.org/biocLite.R')
library(arrayQualityMetrics)
source("~/UPSCb/src/R/plot.multidensity.R")

### ==============================
## set the working dir
### ==============================
setwd("/mnt/picea/projects/spruce/14_SpruceRoots_Project/")

### ==============================
## read the samples details
### ==============================
samples <- read.csv("~/UPSCb/projects/spruce-roots/doc/sample.csv")

### ==============================
## read the HTSeq files in a matrix
## names are set according to the sample.csv!
### ==============================
res <- mclapply(dir("HTSeq",pattern="*.txt",full.names=TRUE),function(fil){
  read.delim(fil,header=FALSE,stringsAsFactors=FALSE)
},mc.cores=9)
names(res) <- gsub(".*CXX_|_index[0-9]+_trimmomatic_sortmerna_STAR\\.txt","",dir("HTSeq",pattern="*.txt"))
names(res) <- samples$Name[match(names(res),samples$ID)]

### ==============================
## get the count table 
### ==============================
addInfo <- c("no_feature","ambiguous","too_low_aQual","not_aligned","alignment_not_unique")
sel <- match(addInfo,res[[1]][,1])
count.table <- do.call(cbind,lapply(res,"[",3))[-sel,]
colnames(count.table) <- names(res)
rownames(count.table) <- res[[1]][,1][-sel]

### ==============================
## according to the previous graphs two samples are not good quality:
## W25-F3 and W34-C4
## so in the following these are removed from the count table
### ==============================
count.table <- count.table[,order(colnames(count.table))]
count.table <- count.table[,!colnames(count.table) %in% c("W25-F3","W34-C4")]

### ==============================
## For visualization, the data is
## submitted to a variance stabilization
## transformation 
### ==============================
conditions <- sub("-[C,F][1-6]m?","",colnames(count.table))
dds <- DESeqDataSetFromMatrix(countData = count.table,
                              colData = data.frame(condition=conditions),
                              design = ~ condition)
colData(dds)$condition <- factor(colData(dds)$condition,
                                 levels=unique(conditions))
vsd <- varianceStabilizingTransformation(dds, blind=TRUE) ## blind=TRUE is good for testing, if I want to use it for analysis
## then it is better to use blind=FORCE (in this case the replicates are not ignored)

## hcluster
vsd.matrix <- assay(vsd) ## extracts the information to a simple matrix
vsd.matrix <- t(vsd.matrix)
dist.matrix <- dist(vsd.matrix)
plot(hclust(dist.matrix),method="average")

## heatmap
dist.matrix <- as.matrix(dist.matrix)
heatmap.2(dist.matrix)

##
vsd.matrix <- t(vsd.matrix)
eSet <- ExpressionSet(assayData=vsd.matrix)
arrayQualityMetrics(eSet)