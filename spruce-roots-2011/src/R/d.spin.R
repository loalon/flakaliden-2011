#' ---
#' title: "Spruce roots project - Differential expression analysis"
#' author: "Alonso Serrano"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---

#' # Introduction
#' 

#' # Setup
#' 
#' Load libraries
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(tximport))
suppressPackageStartupMessages(library(BiocParallel))
suppressPackageStartupMessages(library(plotly))
suppressPackageStartupMessages(library(dplyr))
source("~/Rtoolbox/utilsDE.r")


#' Set the project folder

projectFolder <- "/mnt/picea/projects/spruce/vhurry/14_SpruceRoots_Project"

#' Read metadata
#' ```{r cache=TRUE}
meta <- read.table(file.path(projectFolder, "samples.tsv"), header = TRUE, sep = '\t', stringsAsFactors = TRUE)
#' ```


#' Add the interaction group to the metadata
meta$group <- paste0(meta$Week, meta$Treatment)
meta$group <- factor(gsub("W", "", meta$group))
head(meta, n=3)

#' Get the Salmon files
files <- dir(paste0(projectFolder, "/Salmon"),'quant.sf', recursive = TRUE, full.names = TRUE)
head(files)

#' Extract Sample ID from filenames and check if all files are present
mvec <- basename(dirname((files)))
if (! all(meta$Name %in% mvec)) {
  stop("Metadata doesn't match samples")
} else {
  print("All files present")
}

#' Get the proper paths for the each Salmon quant.sf file
files <- paste0(projectFolder, "/Salmon/", meta$Name, "/quant.sf")

#' Import transcript id to gene id table
#' ```{r cache=TRUE}
tx2gene <- read.table(paste0(projectFolder,"/tx2gene.tsv"), header = TRUE, sep = '\t', stringsAsFactors = FALSE)
txi <- tximport(files, type="salmon", tx2gene=tx2gene)


#' Import data as a DESeq object and place the sample names to the columns

ddsMatrix <- DESeqDataSetFromTximport(txi, colData=meta, ~group)
colnames(ddsMatrix) <- meta$Sample

#' Some genes (5887 )have 0 expression accross all samples
#' we can remove them
ddsMatrix <- ddsMatrix[(rowSums(counts(ddsMatrix)) > 0),]

#' # QA
#' ## PCA
#' Use vst to plot a PCA
vsd.QA <- vst(ddsMatrix, blind=TRUE)
#' ```
#' 
geneNumber <- length(rownames(assay(vsd.QA))) #66360
plotPCA(vsd.QA, intgroup = c("Treatment", "Week"), ntop=geneNumber)

#' The PCA shows and outlier in C:W34 
#' To get which replicate is the outlier we use the pairs function
pairs(~W34.C2+W34.C3+W34.C4,
      data=data.frame(t(t(assay(vsd.QA))[c('W34-C2','W34-C3','W34-C4'),])), 
      lower.panel = panel.smooth, 
      upper.panel = panel.cor,
      main="Replicate scatterplot",
      gap=0, row1attop=FALSE)

#' Sample W34.C4 has a low correlation with the other replicates, we remove it from the dds dataset and from the metadata
#' 
meta <- meta[-grep("W34-C4", meta$Sample),]
ddsMatrix <- ddsMatrix[,-grep("W34-C4", colnames(ddsMatrix))]

#' And we apply another cleanup, 3 more genes now have 0 expression
ddsMatrix <- ddsMatrix[(rowSums(counts(ddsMatrix)) > 0),]

#' We use vst again and then plot a PCA
#' ```{r cache=TRUE}
vsd.QA <- vst(ddsMatrix, blind=TRUE)
#' ```
geneNumber=length(rownames(assay(vsd.QA)))
plotPCA(vsd.QA, intgroup = c("Treatment", "Week"), ntop=geneNumber) 

#' The PCA seems correct
#' 
#' ## PCA for treatment group
plotPCA(vsd.QA, intgroup = c("Treatment"), ntop=geneNumber) 

#' ## PCA for time group
plotPCA(vsd.QA, intgroup = c("Week"), ntop=geneNumber) 

#' ## PCA with letters instead of numbers for increase clarity
vsd.QA.t <- t(assay(vsd.QA))
pca <- prcomp(vsd.QA.t)
percents <- round(summary(pca)$importance[2,]*100)
pchs <- LETTERS[as.integer(factor(substr(rownames(vsd.QA.t), 2, 3)))]
cols <- c("#fc8d62", "#66c2a5")[as.integer(factor(substr(rownames(vsd.QA.t),5,5)))]
PCAlegend <- data.frame(pchs, substr(rownames(vsd.QA.t), 2, 3))

par(oma = c(4, 1, 1, 4))
plot(pca$x[,1], pca$x[,2],
     xlab=paste("Comp. 1 (",percents[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percents[2],"%)",sep=""),
     col=cols,pch=pchs, cex = 1.0, main="")
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("topright", legend=unique(PCAlegend$pchs:meta$Week), cex = 0.8, xpd=TRUE)
legend("bottomright", legend=unique(meta$Treatment), fill=unique(cols),cex = 0.8, xpd=TRUE)

#' ## Scree plot, to see how PC proportions are distributed
screePlot(pca, 10)

#' The first 5 components explain 54 % of the variance
summary(pca)$importance[2,][1:5] *100
sum(summary(pca)$importance[2,][1:5]) * 100

#' ## PCA bi-plot comparition
#' As seen before, we can focus on the 5 components and plot a comparative bi-plot
par(cex=1.0, cex.axis=0.8, cex.main=0.8)
pairs(pca$x[,1:5], col=cols, main="Principal components analysis bi-plot\nPCs 1-5", pch=pchs)


#' ## 3D PCA
#' We represent the first 3 compoments as a 3D interactive object
#' 

colors <- c("#F8766D", "#00BFC4", "#56B4E9")

treatment <- substr(rownames(vsd.QA.t), 5, 5)
treatment[treatment == 'C'] <- 'Control'
treatment[treatment == 'F'] <- 'Fertilised'

df3D <- data.frame(treatment)
df3D$time <- substr(rownames(vsd.QA.t),2,3)
df3D$colors <- colors[as.numeric(df3D$treatment)]
df3D$PC1 <- pca$x[,1]
df3D$PC2 <- pca$x[,2]
df3D$PC3 <- pca$x[,3] 


plot_ly(df3D, x = ~PC1, y = ~PC3, z = ~PC2, 
        color = ~treatment, colors = c("#F8766D", "#00BFC4"),
        text = ~paste("Week:", time, '<br>Treatment:', treatment) ) %>% 
  add_markers() %>%
  layout(scene = list(xaxis = list(title = paste("PC1 (",percents[1],"%)",sep="")),
                      yaxis = list(title = paste("PC3 (",percents[3],"%)",sep="")),
                      zaxis = list(title = paste("PC2 (",percents[2],"%)",sep=""))))

#' # Differential Expression
#' The dataset is ready for DESeq2
#' We check the design of the DE object 
design(ddsMatrix)
#' We are interested on future comparitions like DE genes in each timepoints or between timepoints
#' that is why we created the interaction group
#' 
#' ## Run DESeq
#' ```{r cache=TRUE}
dds <- DESeq(ddsMatrix)
#' ```

#'
#' We check the result names
resultsNames(dds)

#' In order to make all the comparitions automated, we prepare some vector with the elements that will be used
#' in the result extraction
controlVector <- as.character(unique(meta$group[grep("C", meta$group)]))
fertilisedVector <- as.character(unique(meta$group[grep("F", meta$group)]))

