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
source("~/spruceRoots/0.common/utilsDE.r")
#register(MulticoreParam())

#' ```{r echo=FALSE}
#' load("spruceDE2.RData")
#' ```

#' Set the project folder

projectFolder <- "/mnt/picea/projects/spruce/vhurry/14_SpruceRoots_Project"

#' Read metadata
meta <- read.table(file.path(projectFolder, "samples.tsv"), header = TRUE, sep = '\t', stringsAsFactors = TRUE)

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
vsd.QA <- vst(ddsMatrix, blind=TRUE)
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
#' ```{r eval=FALSE}
#' dds <- DESeq(ddsMatrix)
#' ```


#' We check the result names
resultsNames(dds)

#' In order to make all the comparitions automated, we prepare some vector with the elements that will be used
#' in the result extraction
controlVector <- as.character(unique(meta$group[grep("C", meta$group)]))
fertilisedVector <- as.character(unique(meta$group[grep("F", meta$group)]))



#' ## Result extraction
#' We want to make 3 different comparations, each one composed of several outputs from DESeq2
#' In order to do that, we will extract all the results for each comparation type, store them
#' as a list and then we can process them.
#' 
#' ### Control week vs week
#' For this comparition we take consider only the control condition and we compare each pair of weeks
#' to determine which genes are up or down regaluated from one week to the previous week
#' 
#' ```{r eval=FALSE}
#' control.res <- mapply(getRes, controlVector[-1], 
#' controlVector[-length(controlVector)], 
#' MoreArgs = list(localDDS=dds, group="group"))
#' ```
#' 
#' ```{r eval=FALSE}
#' ### Fertilised week vs week
#' For this comparition we take consider only the fertilised condition and we compare each pair of weeks
#' to determine which genes are up or down regaluated from one week to the previous week
#' 
#' fertilised.res <- mapply(getRes, fertilisedVector[-1], 
#' fertilisedVector[-length(fertilisedVector)], 
#' MoreArgs = list(localDDS=dds, group="group"))
#' ```
#' 
#' ### Timepoint Fertilised vs Control
#' For this comparition we take consider only a week each time and then compare DE genes due to fertilisation
#' 
#' ```{r eval=FALSE}
#' combined.res <- mapply(getRes, fertilisedVector, controlVector, 
#' MoreArgs = list(localDDS=dds, group="group"))
#' ```
#' 
#' ```{r eval=FALSE}
#' names(control.res) <- c("20Cvs19C", "21Cvs20C", "22Cvs21C", "23Cvs22C",
#'                      "24Cvs23C", "25Cvs24C", "26Cvs25C", "28Cvs26C",
#'                      "31Cvs28C", "32Cvs31C", "34Cvs32C", "35Cvs34C",
#'                      "36Cvs35C", "37Cvs36C", "38Cvs37C", "39Cvs38C",
#'                      "40Cvs39C", "41Cvs40C")

#' names(fertilised.res) <- c("20Fvs19F", "21Fvs20F", "22Fvs21F", "23Fvs22F",
#'                         "24Fvs23F", "25Fvs24F", "26Fvs25F", "28Fvs26F",
#'                         "31Fvs28F", "32Fvs31F", "34Fvs32F", "35Fvs34F",
#'                         "36Fvs35F", "37Fvs36F", "38Fvs37F", "39Fvs38F",
#'                         "40Fvs39F", "41Fvs40F")
#' 
#' names(combined.res) <- c("19Fvs19C", "20Fvs20C", "21Fvs21C", "22Fvs22C",
#'                       "23Fvs23C", "24Fvs24C", "25Fvs25C", "26Fvs26C",
#'                       "28Fvs28C", "31Fvs31C", "32Fvs32C", "34Fvs34C",
#'                       "35Fvs35C", "36Fvs36C", "37Fvs37C", "38Fvs38C",
#'                       "39Fvs39C", "40Fvs40C", "41Fvs41C")
#' ```
#' 
#' ## Filtering
#' The results contain all the genes, we need to filter them to get only those genes
#' that are stadistical significant and that have a strong condtition effect
#' For that we discard all the genes that has a padj higher than 0.01 and
#' a log2FoldChange lower than 0.5
#' 
control.res.filter <- lapply(control.res, filterDE)
fertilised.res.filter <- lapply(fertilised.res, filterDE)
combined.res.filter <- lapply(combined.res, filterDE)

#' It is a good idea to do a check and see if the filtering was correct
#' ```{r echo=FALSE}
#' lapply(control.res.filter, function(x) print(all(abs(x$log2FoldChange)>0.5 & x$padj<0.01)) )
#' lapply(fertilised.res.filter, function(x) print(all(abs(x$log2FoldChange)>0.5 & x$padj<0.01)) )
#' lapply(combined.res.filter, function(x) print(all(abs(x$log2FoldChange)>0.5 & x$padj<0.01)) )
#' ```

#' ## Plot DE genes
#' We can plot the results as positive and negative bars to reflect the amount of genes up or down
#' regulated in each comparition 
#' For that we create a couple more vector to present the results
mixVector <- c(rbind(controlVector, fertilisedVector))
#' will order the results by week, and then treatment. 19C, 19F, 20C, ...
concatVector <- c(cbind(controlVector, fertilisedVector))
#' will plot the results by treatment. 19C, 20C, ..., 41C, 19F, 20F, ...

#' Control week+1 vs week plot
plotUpAndDown(control.res.filter)
#' Fertilised week+1 vs week plot
plotUpAndDown(fertilised.res.filter)
#' Fertilised vs Control per week plot
plotUpAndDown(combined.res.filter)

#' ## Top genes
#' To check the top genes in each comparition we create a dot plot for the top 6 genes with two lines (read and blue)
#' that highligths the comparition pair that we want to check and store it as a PDF
#' 
#' ```{r eval=FALSE}
#plotTop(control.res.filter, dds, plotGroup=c('Treatment','Week'), infoVector=concatVector, name = "controlFilter.pdf")
#plotTop(fertilised.res.filter, dds, plotGroup=c('Treatment','Week'), infoVector=concatVector, name = "fertilisedFilter.pdf")
#plotTop(combined.res.filter, dds, plotGroup=c('Week','Treatment'), infoVector=mixVector, name = "combinedFilter.pdf")
#```
#' Result export for the combined network

vsdCombined <- vst(dds, blind=FALSE)
combinedData <- t(assay(vsdCombined))
combinedData <- combinedData[,which(colMads(as.matrix(combinedData)) > 0)]
#' For the network construction we need a headless files with the data and a file only with the gene names
write.table(combinedData, file ="spruceRoots2011Combined.tsv", sep='\t', row.names=FALSE, col.names=FALSE)
writeLines(colnames(combinedData), con="spruceRoots2011CombinedNames.txt")

ddsSplit <- dds
design(ddsSplit) <- ~Treatment
ddsSplit <- DESeq(ddsSplit)
vsdSplit <- vst(ddsSplit, blind=F)
vsdSplitData <- t(assay(vsdSplit))
controlData <- vsdSplitData[grep("C",rownames(vsdSplitData)),]
controlData <- controlData[,which(colMads(as.matrix(controlData)) > 0)]

fertilisedData <- vsdSplitData[grep("F",rownames(vsdSplitData)),]
fertilisedData <- fertilisedData[,which(colMads(as.matrix(fertilisedData)) > 0)]

write.table(controlData, file ="spruceRoots2011Control.tsv", sep='\t', row.names=FALSE, col.names=FALSE)
writeLines(colnames(controlData), con="spruceRoots2011ControlNames.txt")

write.table(fertilisedData, file ="spruceRoots2011Fertilised.tsv", sep='\t', row.names=FALSE, col.names=FALSE)
writeLines(colnames(fertilisedData), con="spruceRoots2011FertilisedNames.txt")


#' # Session information
sessionInfo()



#save(control.res, combined.res, fertilised.res, dds,ddsMatrix, control.res.filter, combined.res.filter, fertilised.res.filter, file="spruceDE2.RData")
writeLines(capture.output(sessionInfo()), paste("sessionInfo",Sys.Date() ,".txt", sep=""))

#' # References
