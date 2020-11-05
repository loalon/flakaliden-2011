### ============= 1. Start and set-up of the data =================
## 1.1 load the necessary libraries
### ==============================
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(vsn))
suppressPackageStartupMessages(library(scatterplot3d))
suppressPackageStartupMessages(library(LSD))
suppressPackageStartupMessages(library(knitr))
suppressPackageStartupMessages(library(gplots))
source("/mnt/picea/home/zsofia/UPSCb/src/R/plot.multidensity.R")
source("/mnt/picea/home/zsofia/UPSCb/src/R/rmd.R")
library(limma)

### ==============================
## 1.2 set the working directory
### ==============================
setwd("/mnt/picea/projects/spruce/14_SpruceRoots_Project/")

### ==============================
## 1.3 read the samples details
### ==============================
samples <- read.csv("/mnt/picea/home/zsofia/UPSCb/projects/spruce-roots/doc/sample.csv")

### ==============================
## 1.4 read the HTSeq files in a matrix
## names are set according to the sample.csv!
### ==============================
res <- mclapply(dir("HTSeq",pattern="*.txt",full.names=TRUE),function(fil){
  read.delim(fil,header=FALSE,stringsAsFactors=FALSE)
},mc.cores=9)
names(res) <- gsub(".*CXX_|_index[0-9]+_trimmomatic_sortmerna_STAR\\.txt","",dir("HTSeq",pattern="*.txt"))
names(res) <- samples$Name[match(names(res),samples$ID)]

### ==============================
## 1.5 get the count table 
### ==============================
addInfo <- c("no_feature","ambiguous","too_low_aQual","not_aligned","alignment_not_unique")
sel <- match(addInfo,res[[1]][,1])
count.table <- do.call(cbind,lapply(res,"[",3))[-sel,]
colnames(count.table) <- names(res)
rownames(count.table) <- res[[1]][,1][-sel]

## according to the previous analysis W34-C4 is not good quality and is removed
## count table is subset to control treatment only

count.table <- count.table[,order(colnames(count.table))]
count.table <- count.table[,!colnames(count.table) %in% c("W34-C4")]
count.table <- count.table[,grep("C",colnames(count.table))]

### ==============================
## 18.77 % of the genes are not expressed from a total of 70,736 genes
### ==============================
sel <- rowSums(count.table) == 0
sum(sel) / nrow(count.table)
length(sel)

### ==============================
## display the per-gene mean expression i.e. the mean raw count of every 
## gene across samples is calculated and displayed on a log10 scale
### ==============================
pal=brewer.pal(8,"Dark2")
plot(density(log10(rowMeans(count.table))),col=pal[1],
     main="mean raw counts distribution",
     xlab="mean raw counts (log10)")

### ==============================
## For visualization, the data is submitted to a variance stabilization
## transformation 
### ==============================

conditions <- sub("-[C][1-6]m?","",colnames(count.table))
dds <- DESeqDataSetFromMatrix(countData = count.table,
                              colData = data.frame(condition=conditions),
                              design = ~ condition)
colData(dds)$condition <- factor(colData(dds)$condition,
                                 levels=unique(conditions))
## vsd <- varianceStabilizingTransformation(dds, blind=TRUE) #blind=TRUE takes the replicates into account! #

rld <- rlogTransformation(dds,blind = TRUE) # more resistant to size factor differences

save.image(file="control_rld_image.RData")

### ==============================
## plot the mean against the variance (sd) - it should be flat if the VST worked properly
## it is not but it is likely due to the variability of the gene expression
## This is not a concern for visualization but might be for normalization
### ==============================
# meanSdPlot(assay(vsd)[rowSums(count.table)>0,], ylim = c(0,2.5))
# title(main="Mean vs Variance - no cutoff")

### =============================
## Cutoff 2.3 seems to be sufficient enough for getting an equal variance
###============================
# vsd2 <- vsd[rowMeans(assay(vsd)) > 2.3,]
# meanSdPlot(assay(vsd2))
# title(main="Mean vs Variance - cutoff:2.3")

### ==============================
## Perform a Principal Component Analysis (with cutoff)
### ==============================
pc <- prcomp(t(assay(rld)))
percent <- round(summary(pc)$importance[2,]*100)
cols <- rep(brewer.pal(7,"Dark2"),3)[as.integer(factor(conditions))]
pchs <- LETTERS[as.integer(factor(conditions))]


### ==============================
## plot the first 3 dimensions
## First component: 28%
## Second component: 13%
## Third component: 6%
### ==============================
### ==============================
## weeks are represented by letters (A-S)
### ==============================
scatterplot3d(pc$x[,1],
              pc$x[,2],
              pc$x[,3],
              xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
              ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
              zlab=paste("Comp. 3 (",percent[3],"%)",sep=""),
              color=cols,pch=pchs)
### ==============================
## Dimensions are plotted separately
### ==============================
par(mar=c(5.1,4.1,4.1,2.1))

### ==============================
## Component1-Component2
### ==============================
### ==============================
## weeks are represented by letters (A-S)
### ==============================
plot(pc$x[,1],
     pc$x[,2],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     col=cols,pch=pchs,
     main="Principal Component Analysis",sub="variance stabilized counts")

### ==============================
## Component1-Component3
### ==============================
### ==============================
## weeks are represented by letters (A-S)
### ==============================
plot(pc$x[,1],
     pc$x[,3],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 3 (",percent[3],"%)",sep=""),
     col=cols,pch=pchs,
     main="Principal Component Analysis",sub="variance stabilized counts")


### ==============================
## Component2-Component3
### ==============================
### ==============================
## weeks are represented by letters (A-S)
### ==============================
plot(pc$x[,2],
     pc$x[,3],
     xlab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     ylab=paste("Comp. 3 (",percent[3],"%)",sep=""),
     col=cols,pch=pchs,
     main="Principal Component Analysis",sub="variance stabilized counts")

###================================
## SELECTING WHICH GENES SHOULD BE CONSIDERED FOR FURTHER ANALSIS
###===============================

### ==============================
## Set cutoff for lowly-expressed genes
### ==============================
df <- as.data.frame(assay(vsd))
colnames(df) <- names(count.table)

### ==============================
## the vst adds a min value for non-expressed gene. The following line of code is just used to 
## substract this value from all expression data.
## then all genes that are 0 can be considered as non-expressed, which will be used for the
## analysis below introduce rmd first !!! (from ("~/UPSCb/src/R/rmd.R"))
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
## Table of genes that are expressed in all samples (TRUE) or not 
## expressed in all samples (FALSE)
## Here: 57462 out of 70736 genes are expressed in all samples
### ==============================

table(df$expressed)
sum(table(df$expressed))

### ==============================
## Table of genes that are expressed in the given number of samples 
## (in none=0, in 1 sample=1 etc)
### ==============================

table(df$sample.expressed)
barplot(table(df$sample.expressed)[0:56])
title(ylab = "number of genes", xlab = "number of samples (0-56)", font.lab = 1)

barplot(table(df$sample.expressed)[4:53])
title(ylab = "number of genes", xlab = "number of samples (4-53)", font.lab = 1)


### ==============================
## Plot of mean expression accross all samples (Y-axis) against
## number of samples in which these genes are expressed (x-axis)
### ==============================

par(mar=c(4.1,4.1,1.1,0.1))
boxplot(split(df$mean.expression,df$sample.expressed))
title(ylab = "Mean expression", xlab = "Number of samples with expressed genes", font.lab = 1) 

par(mar=c(4.1,4.1,1.1,0.1))
boxplot(split(df$mean.expression,df$sample.expressed)[0:55])
title(ylab = "Mean expression", xlab = "Number of samples with expressed genes (0:55)", font.lab = 1)

### ==============================
## Plot of mean variance accross all samples (Y-axis) against
## number of samples in which these genes are expressed (x-axis)
### ==============================
boxplot(split(df$rmd,df$sample.expressed))
title(ylab = "Mean variance", xlab = "Number of samples with expressed genes", font.lab = 1)

### ==============================
## Plot the mean against the variance after variance stabilization transformation
## and with mean expression values ordere from small to large
## X-axis indicates number of genes add two lines into the plot indicating
## 1. the first expressed gene (in at least 1 sample) and
## 2. the first gene expressed in all samples (this is the
## gene out of the 70736 all-sample expressed ones with lowest mean expression
## accross all samples)
## Note that on the right side of this line there are still genes that are not 
## expressed in all samples (x-axis ordered) by mean expression, not by number of samples
## in which each gene is expressed
## We won't work with the ca. 35000 genes on on the left side of this line. But we will
## consider the 35000 genes on the right side of this line (! not necessarily expressed
## in all samples! )
### ==============================

meanSdPlot(assay(vsd))
nrow(assay(vsd))
abline(v=sum(!df$expressed))
abline(v=match("56",df[order(df$mean.expression),"sample.expressed"]))

### ==============================
## Get expression information for the
## gene that is indicated by the 
## right-most vertical line in the aforementioned
## plot. The mean expression value of
## this gene (1.66885) will be our cut-off
## for considering genes that are (non-)expressed
## the rmd value here is Relative Mean Difference 
## (variance) and has a commons scale accross
## all genes
### ==============================

df[order(df$mean.expression),][match("56",df[order(df$mean.expression),"sample.expressed"]),]

### ==============================
## Now we set a threshold at that position
## and plot once again the mean expression 
## accross all samples (Y-axis) against
## number of samples in which
## these genes are expressed (x-axis) 
### ==============================

pos <- match("56",df[order(df$mean.expression),"sample.expressed"])
df <- df[order(df$mean.expression),]
dfs <- df[pos:nrow(df),]
boxplot(split(dfs$mean.expression,dfs$sample.expressed))
title(ylab = "Mean expression", xlab = "Number of samples with expressed genes", font.lab = 1)

### ==============================
## This plot below indicates relative mean difference (variance) 
## accross all samples against number of samples in which
## these genes are expressed (x-axis)
## Compared to what we had before (when considering all genes), the
## variance has decreased and the quality of the considered data
## therefore increased
### ==============================
boxplot(split(dfs$rmd,dfs$sample.expressed))
title(ylab = "Relative Mean Difference", xlab = "Number of samples with expressed genes", font.lab = 1)

### ==============================
## Heat map (sample to sample clustering)
### ==============================
hmcol <- colorRampPalette(brewer.pal(9,"GnBu"))(100)
dist.matrix <- dist(t(assay(vsd)))
mat <- as.matrix(dist.matrix)
rownames(mat) <- colnames(mat) <- with(colData(dds),
                                       paste(condition))
heatmap.2(mat,trace="none",col=rev(hmcol))

### ==============================
## rlog transformation
## heat map with rlog transformed data (rld looks better so that is used in wgcna)
### ==============================

dist.matrix <- dist(t(assay(rld)))
mat <- as.matrix(dist.matrix)
rownames(mat) <- colnames(mat) <- with(colData(dds),
                                       paste(condition))
heatmap.2(mat,trace="none",col=rev(hmcol))

### ==============================
## create rld data files for WGCNA analysis (done already)
### ==============================

## write.table(t(assay(vsd)),
           # "/mnt/picea/projects/spruce/14_SpruceRoots_Project/data/network/control_rld.tsv",sep="\t",
           # col.names = TRUE, row.names = FALSE, quote = FALSE)
## write.table(colnames(count.table),
           # "/mnt/picea/projects/spruce/14_SpruceRoots_Project/data/network/control_rld_sampeInfo.tsv",sep="\n",
           # col.names = FALSE, row.names = FALSE, quote = FALSE)

### ==============================
## filter the genes for at least 5 reads in 3 samples
### ==============================
plot(density(assay(vsd)))

library(genefilter)
kOverA(5,3)
filterfun(kOverA(5,3))
rldfilter_5_3 <- genefilter(assay(rld),filterfun(kOverA(5,3)))
rldfilter_5_3[1:6]

lines(density(assay(rld)[rldfilter_5_3,]), col="red")

> sum(rldfilter_5_3)
> rld_f <- rld[t(rldfilter_5_3),]
> dim(rld_f)
> dim(rld)







