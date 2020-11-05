### ==============================
## load the necessary libraries
### ==============================
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(vsn))
suppressPackageStartupMessages(library(scatterplot3d))
suppressPackageStartupMessages(library(LSD))
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
names(res) <- gsub(".*CXX_|_index[0-9]+_trimmomatic_sortmerna_STAR\\.txt","",
                   dir("HTSeq",pattern="*.txt"))
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
## plot the stats
## stats are sorted by week 
### ==============================
pal=brewer.pal(6,"Dark2")[1:nrow(count.stats)]
par(mar=c(5.1,4.1,4.1,2.1))
weeks <- sub("-[C,F][1-6]m?","",names(res))
sapply(sort(unique(weeks)),function(week,count.stats){
  dat <- count.stats[,grep(week,colnames(count.stats))]
  barplot(as.matrix(dat),col=pal,beside=TRUE,las=2,
          main=paste("week",week,"read proportion"),
          ylim=c(1,1e7))
  legend("top",fill=pal,legend=rownames(dat),bty="n",cex=0.5,horiz=TRUE)
},count.stats)

## 16.5 % of the genes are not expressed
## from a total of 70,736 genes
### ==============================
sel <- rowSums(count.table) == 0
sum(sel) / nrow(count.table)
length(sel)

### ==============================
## display the per-gene mean expression
## i.e. the mean raw count of every 
## gene across samples is calculated
## and displayed on a log10 scale
### ==============================
plot(density(log10(rowMeans(count.table))),col=pal[1],
     main="mean raw counts distribution",
     xlab="mean raw counts (log10)")

### ==============================
## The same is done for the individual
## samples colored by sample type
## separated into weeks
### ==============================

pal=brewer.pal(8,"Dark2")
pal2=brewer.pal(12,"Paired")
sapply(sort(unique(sub("-[C,F][1-6]m?","",names(res)))),function(week,count.table){
  dat <- count.table[,grep(week,colnames(count.table))]
  plot.multidensity(log10(dat),col=pal,
                    legend.x="topright",legend.cex=0.5,
                    main=paste("week",week,"sample raw counts distribution"),
                    xlab="per gene raw counts (log10)",lwd=3)
},count.table)

### ==============================
## according to the previous graphs two samples are not good quality:
## W25-F3 and W34-C4
## so in the following these are removed from the count table
### ==============================
count.table <- count.table[,order(colnames(count.table))]
count.table <- count.table[,!colnames(count.table) %in% c("W25-F3","W34-C4")]
sizeF <- estimateSizeFactorsForMatrix(count.table)
min(sizeF)
max(sizeF)

### ==============================
## display the per-gene mean expression after removal of the "bad" samples
## i.e. the mean raw count of every 
## gene across samples is calculated
## and displayed on a log10 scale
### ==============================
plot(density(log10(rowMeans(count.table))),col=pal[1],
     main="mean raw counts distribution",
     xlab="mean raw counts (log10)")

### ==============================
## For visualization, the data is
## submitted to a variance stabilisation
## transformation (conditions are not taken into account)
### ==============================
conditions <- sub("-[C,F][1-6]m?","",colnames(count.table))

dds <- DESeqDataSetFromMatrix(countData = count.table,
                              colData = data.frame(condition=conditions),
                              design = ~ condition)
colData(dds)$condition <- factor(colData(dds)$condition,
                                 levels=unique(conditions))
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)

### ==============================
## dds designed by combined week&treatment
### ==============================
treatment <- substr(x=colnames(count.table),1,5)
dds2 <- DESeqDataSetFromMatrix(countData = count.table,
                              colData = data.frame(condition=treatment),
                              design = ~ condition)
colData(dds2)$condition <- factor(colData(dds2)$condition,
                                 levels=unique(treatment))

### replicates are collapsed - Only for visualisation!!!!

ddsColl <- collapseReplicates(dds2, groupby = dds2$condition)
vsdColl <- varianceStabilizingTransformation(ddsColl, blind=TRUE)

### ==============================
## plot the mean against the variance (sd)
## it should be flat if the VST worked properly
## it is not but it is likely due to the 
## variability of the gene expression
## This is not a concern for visualization
## but might be for normalization
### ==============================
par(mfrow=c(1,2))
meanSdPlot(assay(vsd)[rowSums(count.table)>0,], ylim = c(0,2.5))
title(main="Mean vs Variance - vsd")

meanSdPlot(assay(vsdColl)[rowSums(count.table)>0,], ylim = c(0,2.5))
title(main="Mean vs Variance - vsdColl")

### ==============================
## Perform a Principal Component Analysis
### ==============================
pc <- prcomp(t(assay(vsd)))
percent <- round(summary(pc)$importance[2,]*100)
pchs <- LETTERS[as.integer(factor(conditions))]
cf <- substr(x=colnames(count.table),5,5)
cols <- c("darkred","darkgreen")[as.integer(factor(cf))]
PCAlegend <- data.frame(pchs,conditions)

pc2 <- prcomp(t(assay(vsdColl)))
percent <- round(summary(pc2)$importance[2,]*100)
tr <- substr(x=unique(treatment),1,3)
pchs2 <- LETTERS[as.integer(factor(tr))]
cf2 <- substr(x=unique(treatment),5,5)
cols2 <- c("darkred","darkgreen")[as.integer(factor(cf2))]
PCAlegend2 <- data.frame(pchs2,tr)

### ==============================
## Component1-Component2
### ==============================
### ==============================
## weeks are represented by letters (A-S)
## red: Control, green: fertilized
### ==============================
par(mar=c(5.1,4.1,4.1,2.1))
par(mfrow=c(1,1))
plot(pc$x[,1],
     pc$x[,2],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     col=cols,pch=pchs,
     main="Principal Component Analysis - vsd")
legend("topright",
       legend=unique(PCAlegend$pchs:PCAlegend$conditions),
       cex=0.4)

plot(pc2$x[,1],
     pc2$x[,2],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     col=cols2,pch=pchs2,
     main="Principal Component Analysis - vsdColl")
legend("topright",
       legend=unique(PCAlegend2$pchs2:PCAlegend2$tr),
       cex=0.4)

plot(pc$x[,1],
     pc$x[,3],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 3 (",percent[3],"%)",sep=""),
     col=cols,pch=pchs,
     main="Principal Component Analysis - vsd")
legend("topright",
       legend=unique(PCAlegend$pchs:PCAlegend$conditions),
       cex=0.4)

plot(pc2$x[,1],
     pc2$x[,3],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 3 (",percent[3],"%)",sep=""),
     col=cols2,pch=pchs2,
     main="Principal Component Analysis - vsdColl")
legend("topright",
       legend=unique(PCAlegend2$pchs2:PCAlegend2$tr),
       cex=0.4)

plot(pc$x[,2],
     pc$x[,3],
     xlab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     ylab=paste("Comp. 3 (",percent[3],"%)",sep=""),
     col=cols,pch=pchs,
     main="Principal Component Analysis - vsd")
legend("topright",
       legend=unique(PCAlegend$pchs:PCAlegend$conditions),
       cex=0.4)

plot(pc2$x[,2],
     pc2$x[,3],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     col=cols2,pch=pchs2,
     main="Principal Component Analysis - vsdColl")
legend("topright",
       legend=unique(PCAlegend2$pchs2:PCAlegend2$tr),
       cex=0.4)

### ==============================
## From the PCA it looks like two samples have been switched
## in W37-F and W40-F
## comparing the correlation between the samples to see which one is the outlier
### ==============================
## cor=0,96 ##
heatscatter(assay(vsd)[,86],assay(vsd)[,87],
            xlab=colnames(count.table)[86],
            ylab=colnames(count.table)[87])

## cor=0,95 ##
heatscatter(assay(vsd)[,87],assay(vsd)[,88],
            xlab=colnames(count.table)[87],
            ylab=colnames(count.table)[88])

## cor=0,95 ##
heatscatter(assay(vsd)[,86],assay(vsd)[,88],
            xlab=colnames(count.table)[86],
            ylab=colnames(count.table)[88])

## cor=0,94 ##
heatscatter(assay(vsd)[,104],assay(vsd)[,105],
            xlab=colnames(count.table)[104],
            ylab=colnames(count.table)[105])

## cor=0,95 ##
heatscatter(assay(vsd)[,104],assay(vsd)[,106],
            xlab=colnames(count.table)[104],
            ylab=colnames(count.table)[106])

## cor=0,96 ##
heatscatter(assay(vsd)[,106],assay(vsd)[,105],
            xlab=colnames(count.table)[106],
            ylab=colnames(count.table)[105])

### ==============================
## W40-F1 (104) is outlier, W37-F3 (88) is outlier 
## W40-F1 has a higher correlation to W37-F1 and W37-F2
## W40-F1 and W37-F1 could be exchanged, but there is no real basis for doing that
## So we could remove them completely, but for now they will stay in the analysis
### ==============================
## cor=0,94 ##
heatscatter(assay(vsd)[,104],assay(vsd)[,87],
            xlab=colnames(count.table)[104],
            ylab=colnames(count.table)[87])

## cor=0,94 ##
heatscatter(assay(vsd)[,104],assay(vsd)[,86],
            xlab=colnames(count.table)[104],
            ylab=colnames(count.table)[86])

## cor=0,95 ##
heatscatter(assay(vsd)[,88],assay(vsd)[,106],
            xlab=colnames(count.table)[88],
            ylab=colnames(count.table)[106])

## cor=0,95 ##
heatscatter(assay(vsd)[,88],assay(vsd)[,105],
            xlab=colnames(count.table)[88],
            ylab=colnames(count.table)[105])


###================================
## SELECTING WHICH GENES SHOULD BE CONSIDERED FOR FURTHER ANALSIS
###===============================

### ==============================
## Set cutoff for lowly-expressed genes
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
  mean(abs(x-mean(x)))/mean(x)
}

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
barplot(table(df$sample.expressed)[3:110])
title(ylab = "number of genes", xlab = "number of samples (3-110)", font.lab = 1)

### ==============================
## Plot of mean expression accross 
## all samples (Y-axis) against
## number of samples in which
## these genes are expressed (x-axis)
### ==============================

par(mar=c(4.1,4.1,1.1,0.1))
boxplot(split(df$mean.expression,df$sample.expressed))
title(ylab = "Mean expression", xlab = "Number of samples with expressed genes", font.lab = 1) 

par(mar=c(4.1,4.1,1.1,0.1))
boxplot(split(df$mean.expression,df$sample.expressed)[0:110])
title(ylab = "Mean expression", xlab = "Number of samples with expressed genes (0:110)", 
      font.lab = 1)

### ==============================
## Plot of mean variance accross 
## all samples (Y-axis) against
## number of samples in which
## these genes are expressed (x-axis)
### ==============================
boxplot(split(df$rmd,df$sample.expressed))
title(ylab = "Mean variance", xlab = "Number of samples with expressed genes", font.lab = 1)

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
abline(v=match("111",df[order(df$mean.expression),"sample.expressed"]))

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

df[order(df$mean.expression),][match("111",df[order(df$mean.expression),"sample.expressed"]),]

### ==============================
## Now we set a threshold at that position
## and plot once again the mean expression 
## accross all samples (Y-axis) against
## number of samples in which
## these genes are expressed (x-axis) 
### ==============================

pos <- match("111",df[order(df$mean.expression),"sample.expressed"])
df <- df[order(df$mean.expression),]
dfs <- df[pos:nrow(df),]
boxplot(split(dfs$mean.expression,dfs$sample.expressed))
title(ylab = "Mean expression", xlab = "Number of samples with expressed genes", font.lab = 1)

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
title(ylab = "Relative Mean Difference", xlab = "Number of samples with expressed genes", 
      font.lab = 1)

## save.image(file="~/UPSCb/projects/spruce-roots/src/R/basicAnalysis_image1.RData")



