### ==============================
## load the necessary libraries
### ==============================
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(vsn))
suppressPackageStartupMessages(library(scatterplot3d))
suppressPackageStartupMessages(library(LSD))
source("~/UPSCb/src/R/plot.multidensity.R")
source("~/UPSCb/src/R/rmd.R")

### ==============================
## set the working directory
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
## get the last stat lines
### ==============================
count.stats <- do.call(cbind,lapply(res,"[",2))[sel,]
colnames(count.stats) <- names(res)
rownames(count.stats) <- addInfo
count.stats <- rbind(count.stats,aligned=colSums(count.table))
count.stats <- count.stats[rowSums(count.stats) > 0,]


### ==============================
## according to the previous graphs two samples are not good quality:
## W25-F3 and W34-C4
## also according to the first PCA two samples seems to be mixed up, because
## they don't cluster with the other samples in their group: 
## W37-F3 and W40-F1
## so in the following these are removed from the count table
### ==============================
count.table <- count.table[,order(colnames(count.table))]
count.table <- count.table[,!colnames(count.table) %in% c("W25-F3","W34-C4","W37-F3","W40-F1")]

### ==============================
## display the per-gene mean expression after removal of the "bad" samples
## i.e. the mean raw count of every 
## gene across samples is calculated
## and displayed on a log10 scale
### ==============================
pal=brewer.pal(8,"Dark2")
pal2=brewer.pal(12,"Paired")

plot(density(log10(rowMeans(count.table))),col=pal[1],
     main="mean raw counts distribution",
     xlab="mean raw counts (log10)")


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
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)

### ==============================
## plot the mean against the variance (sd)
## it should be flat if the VST worked properly
## it is not but it is likely due to the 
## variability of the gene expression
## This is not a concern for visualization
## but might be for normalization
### ==============================
meanSdPlot(assay(vsd)[rowSums(count.table)>0,], ylim = c(0,2.5))
title(main="Mean vs Variance - no cutoff")

sel2 <- rowMeans(assay(vsd)) > 2
meanSdPlot(assay(vsd)[rowMeans(assay(vsd)) > 2,])
title(main="Mean vs Variance - cutoff:2")

### ==============================
## Perform a Principal Component Analysis
### ==============================
pc <- prcomp(t(assay(vsd)))
percent <- round(summary(pc)$importance[2,]*100)
pchs <- LETTERS[as.integer(factor(conditions))]
cf <- substr(x=colnames(count.table),5,5)
cols <- c("darkred","darkgreen")[as.integer(factor(cf))]
PCAlegend <- data.frame(pchs,conditions)

### ==============================
## plot the first 3 dimensions
## First component: 25%
## Second component: 11%
## Third component: 5%
### ==============================
### ==============================
## weeks are represented by letters (A-S)
## red: Control, green: fertilized
### ==============================
scatterplot3d(pc$x[,1],
              pc$x[,2],
              pc$x[,3],
              xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
              ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
              zlab=paste("Comp. 3 (",percent[3],"%)",sep=""),
              color=cols,pch=pchs)
legend("topright",
       legend=unique(PCAlegend$pchs:PCAlegend$conditions),
       cex=0.4)

### ==============================
## Dimensions are plotted separately
### ==============================
par(mar=c(5.1,4.1,4.1,2.1))

### ==============================
## Component1-Component2
### ==============================
### ==============================
## weeks are represented by letters (A-S)
## red: Control, green: fertilized
### ==============================
plot(pc$x[,1],
     pc$x[,2],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     col=cols,pch=pchs,
     main="Principal Component Analysis")
legend("topright",
       legend=unique(PCAlegend$pchs:PCAlegend$conditions),
       cex=0.4)

### ==============================
## Component2-Component3
### ==============================
### ==============================
## weeks are represented by letters (A-S)
## red: Control, green: fertilized
### ==============================
plot(pc$x[,2],
     pc$x[,3],
     xlab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     ylab=paste("Comp. 3 (",percent[3],"%)",sep=""),
     col=cols,pch=pchs,
     main="Principal Component Analysis")
legend("topright",
       legend=unique(PCAlegend$pchs:PCAlegend$conditions),
       cex=0.4)

### ==============================
## Component1-Component3
### ==============================
### ==============================
## weeks are represented by letters (A-S)
## red: Control, green: fertilized
### ==============================
plot(pc$x[,1],
     pc$x[,3],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 3 (",percent[3],"%)",sep=""),
     col=cols,pch=pchs,
     main="Principal Component Analysis")
legend("topright",
       legend=unique(PCAlegend$pchs:PCAlegend$conditions),
       cex=0.4)


###================================
## SELECTING WHICH GENES SHOULD BE CONSIDERED FOR FURTHER ANALSIS
###===============================

### ==============================
## Set cutoff for lowl-expressed genes
### ==============================
df <- as.data.frame(assay(vsd))
colnames(df) <- names(count.table)

### ==============================
## the vst adds a min value for non-expressed 
## gene. The following line of code is just used to 
## substract this value from all expression data.
## then all genes that are 0 can be considered 
## as non-expressed, which will be used for the
## analysis below 
## introduce rmd first !!! (from ("~/UPSCb/src/R/rmd.R"))
### ============================
rmd <- function(x){
  mean(abs(x-mean(x)))/mean(x)}

df <- df - min(df)
ncols <- ncol(df)
df$expressed <- rowSums(df) > 0
df$mean.expression <- rowMeans(df[,1:ncols])
df$sample.expressed <- rowSums(df[,1:ncols]>0)
df$rmd <- apply(df[,1:ncols],1,rmd)

### ==============================
## Table of genes that are expressed 
## in all samples (TRUE) or not 
## expressed in all samples (FALSE)
## Here: 59098 out of 70736 genes are expressed
## in all samples
### ==============================

table(df$expressed)
sum(table(df$expressed))

### ==============================
## Table of genes that are expressed 
## in the given number of samples 
## (in none=0, in 1 sample=1 etc)
### ==============================

table(df$sample.expressed)
barplot(table(df$sample.expressed)[2:109])
title(ylab = "number of genes", xlab = "number of samples (2:109)", 
      font.lab = 1)

### ==============================
## Plot of mean expression accross 
## all samples (Y-axis) against
## number of samples in which
## these genes are expressed (x-axis)
### ==============================

par(mar=c(4.1,4.1,1.1,0.1))
boxplot(split(df$mean.expression,df$sample.expressed))
title(ylab = "Mean expression", 
      xlab = "Number of samples with expressed genes", font.lab = 1) 

par(mar=c(4.1,4.1,1.1,0.1))
boxplot(split(df$mean.expression,df$sample.expressed)[0:109])
title(ylab = "Mean expression", 
      xlab = "Number of samples with expressed genes (0:109)", font.lab = 1)

### ==============================
## Plot of mean variance accross 
## all samples (Y-axis) against
## number of samples in which
## these genes are expressed (x-axis)
### ==============================
boxplot(split(df$rmd,df$sample.expressed))
title(ylab = "Mean variance", 
      xlab = "Number of samples with expressed genes", font.lab = 1)

### ==============================
## Plot the mean against the variance 
## after variance stabilization transformation
## and with mean expression values
## ordere from small to large
## X-axis indicates number of genes
## add two lines into the plot indicating
## 1. the first expressed gene 
## (in at least 1 sample) and
## 2. the first gene 
## expressed in all samples (this is the
## gene out of the 70736 all-sample expressed 
## ones with lowest mean expression
## accross all samples)
## Note that on the right side of this line
## there are still genes that are not 
## expressed in all samples (x-axis ordered)
## by mean expression, not by number of samples
## in which each gene is expressed
## We won't work with the ca. 35000 genes on
## on the left side of this line. But we will
## consider the 35000 genes on the right side
## of this line (! not necessarily expressed
## in all samples! )
### ==============================

meanSdPlot(assay(vsd))
nrow(assay(vsd))
abline(v=sum(!df$expressed))
abline(v=match("109",df[order(df$mean.expression),"sample.expressed"]))

### ==============================
## Get expression information for the
## gene that is indicated by the 
## right-most vertical line in the aforementioned
## plot. The mean expression value of
## this gene (1.67407) will be our cut-off
## for considering genes that are (non-)expressed
## the rmd value here is Relative Mean Difference 
## (variance) and has a commons scale accross
## all genes
### ==============================

df[order(df$mean.expression),][match("109",df[order(df$mean.expression),
                                              "sample.expressed"]),]

### ==============================
## Now we set a threshold at that position
## and plot once again the mean expression 
## accross all samples (Y-axis) against
## number of samples in which
## these genes are expressed (x-axis) 
### ==============================

pos <- match("109",df[order(df$mean.expression),"sample.expressed"])
df <- df[order(df$mean.expression),]
dfs <- df[pos:nrow(df),]
boxplot(split(dfs$mean.expression,dfs$sample.expressed))
title(ylab = "Mean expression", 
      xlab = "Number of samples with expressed genes", font.lab = 1)

### ==============================
## This plot below indicates
## relative mean difference (variance) 
## accross all samples against
## number of samples in which
## these genes are expressed (x-axis)
## Compared to what we had before 
## (when considering all genes), the
## variance has decreased and the 
## quality of the considered data
## therefore increased
### ==============================
boxplot(split(dfs$rmd,dfs$sample.expressed))
title(ylab = "Relative Mean Difference", 
      xlab = "Number of samples with expressed genes", font.lab = 1)

### ==============================
## Perform a Principal Component Analysis for the genes that are expressed above
## the level of mean expression
### ==============================
vsd109 <- vsd[order(df$mean.expression)]
pc109 <- prcomp(t(assay(vsd109)))
percent109 <- round(summary(pc2)$importance[2,]*100)

meanSdPlot(assay(vsd109))


scatterplot3d(pc2$x[,1],
              pc2$x[,2],
              pc2$x[,3],
              xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
              ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
              zlab=paste("Comp. 3 (",percent[3],"%)",sep=""),
              color=cols,pch=pchs)


## subset the lower expression


