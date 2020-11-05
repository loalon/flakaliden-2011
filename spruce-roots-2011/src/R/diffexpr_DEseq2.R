### ============= 1. Start and set-up of the data =================
## 1.1 load the necessary libraries
### ==============================
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(vsn))
suppressPackageStartupMessages(library(scatterplot3d))
suppressPackageStartupMessages(library(LSD))
suppressPackageStartupMessages(library(knitr))
source("~/UPSCb/src/R/plot.multidensity.R")
source("~/UPSCb/src/R/rmd.R")
library(limma)

### ==============================
## 1.2 set the working directory
### ==============================
setwd("/mnt/picea/projects/spruce/14_SpruceRoots_Project/")

### ==============================
## 1.3 read the samples details
### ==============================
samples <- read.csv("~/UPSCb/projects/spruce-roots/doc/sample.csv")

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

## according to the previous analysis two samples are not good quality:
## W25-F3 and W34-C4
## so in the following these are removed from the count table

count.table <- count.table[,order(colnames(count.table))]
count.table <- count.table[,!colnames(count.table) %in% c("W25-F3","W34-C4","W37-F3","W40-F1")]

### ===========2. Creat the frame of work ===================
## 2.1 Create filters
### ==============================
conditions <- substr(x=colnames(count.table),1,5)
weeks <- sub("-[C,F][1-6]m?","",colnames(count.table))
treatments <- substr(x=colnames(count.table),5,5)
sgroups <- c(rep("G1",3),rep("G2",4),rep("G3",5),rep("G4",2),rep("G5",5))[as.integer(factor(weeks))]

### ==============================
## 2.2 create the design matrix - this is the design for my biological question
### ==============================
df <- data.frame (
  condition=conditions,
  week=weeks,
  treatment=treatments,
  sgroup=sgroups)
head(df)

### ==============================
## 2.2 create dds object
### ==============================
dds <- DESeqDataSetFromMatrix(countData = count.table,
                              colData = df,
                              design = ~ condition)

## check the size factors

dds <- estimateSizeFactors(dds)
sizeFactors(dds)

## estimate the disperison (takes a long time!)

dds <- estimateDispersions(dds)
plotDispEsts(dds)

### ===========3. Negative binomial Wald Test===================
## run the test - negative binomial Wald test (takes a long time!)
### ==============================
dds <- nbinomWaldTest(dds)

### ==============================
## create ddsC (for comparison by condition)
### ==============================
ddsC <- dds

### ===========4. Methods and Functions===================
## 
### ==============================
setGeneric(name="VolcanoPlotMA",def=function(object,alpha=0.01){
  standardGeneric("VolcanoPlotMA")})

setMethod(f="VolcanoPlotMA",
          signature="DataFrame",
          definition=function(object,alpha=0.001){
            
            ## lib
            require(LSD)
            
            ## selectors
            sel <- ! is.na(object$padj)
            sel2 <- object$padj[sel]<=alpha
            
            ## plot
            heatscatter(object$log2FoldChange[sel],
                        -log10(object$padj[sel]),
                        main="Volcano",xlab="Log2 Fold Change",
                        ylab="- log(10) adj. p-value")
            
            ## legend
            legend("topleft",bty="n",paste("cutoff @",alpha),lty=2,col="gray")
            
            ## points
            points(object$log2FoldChange[sel][sel2],-log10(object$padj[sel][sel2]),col="lightblue",pch=19)
            points(object$log2FoldChange[sel][sel2],-log10(object$padj[sel][sel2]),col="dodgerblue3",pch=19,cex=0.5)
            
            ## circle the points for the dot plot
            abline(h=-log10(alpha),lty=2,col="gray")
          })

setMethod(f="plotMA",
          signature="DataFrame",
          definition=function(object,alpha=0.001){
            
            ## lib
            require(LSD)
            
            ## selectors
            sel <- ! is.na(object$padj)
            sel2 <- object$padj[sel]<=alpha
            
            ## graphic params
            orig.par <- par(no.readonly=TRUE)
            par(mfrow=c(2,1))
            
            ## plots
            kde2dplot(log10(object$baseMean[sel]),
                      object$log2FoldChange[sel],
                      grid=250,ncol=30,nlevels=10,
                      main="MA density estimation"
            )
            
            heatscatter(log10(object$baseMean[sel]),
                        object$log2FoldChange[sel],
                        add.contour=TRUE,main="MA")
            
            mtext(paste(sum(sel2),"sig. feats. @",alpha,"cutoff"),
                  side=1,line=2)
            
            points(log10(object$baseMean[sel][sel2]),
                   object$log2FoldChange[sel][sel2],
                   col="darkred",pch=19,cex=.5)
            
            par(orig.par,no.readonly=TRUE)
            invisible(TRUE)
          })
### ============5. Comparison of seasonal groups ==================
## According to the PCA plots I divided the data set into five seasonal groups:
## G1. W19-W21 (A-C on the plots)
## G2. W22-W25 (D-G on the plots)
## G3. W26-W34 (H-L on the plots)
## G4. W35-W36 (M-N on the plots)
## G5. W37-W41 (O-S on the plots)
## G4 is a group of "outliers", because they are distributed between G3 and G5
### ==============================
design(dds) <- ~sgroup
ddsG <- DESeq(dds)

### ==============================
## 5.1 comparing G1 to G2
### ==============================

res1 <- results(ddsG,contrast=c("sgroup","G1","G2"))
head(res1)

VolcanoPlotMA(res1,alpha=0.01)
sum(res1$padj <= 0.01,na.rm=TRUE)
plot(density(res1$log2FoldChange[res1$padj <= .01 & !is.na(res1$padj)]))
table(sign(res1$log2FoldChange[res1$padj <= .01 & !is.na(res1$padj)]))
table(sign(res1$log2FoldChange[res1$padj <= .01 & !is.na(res1$padj) & abs(res1$log2FoldChange)>=2]))

### ==============================
## 5.2 comparing G2 to G3
### ==============================

res2 <- results(ddsG,contrast=c("sgroup","G2","G3"))
head(res2)

VolcanoPlotMA(res2,alpha=0.01)
sum(res2$padj <= 0.01,na.rm=TRUE)
plot(density(res2$log2FoldChange[res2$padj <= .01 & !is.na(res2$padj)]))
table(sign(res2$log2FoldChange[res2$padj <= .01 & !is.na(res2$padj)]))
table(sign(res2$log2FoldChange[res2$padj <= .01 & !is.na(res2$padj) & abs(res2$log2FoldChange)>=2]))

### ==============================
## 5.3 comparing G3 to G5
### ==============================

res3 <- results(ddsG,contrast=c("sgroup","G3","G5"))
head(res3)

VolcanoPlotMA(res3,alpha=0.01)
sum(res3$padj <= 0.01,na.rm=TRUE)
plot(density(res3$log2FoldChange[res3$padj <= .01 & !is.na(res3$padj)]))
table(sign(res3$log2FoldChange[res3$padj <= .01 & !is.na(res3$padj)]))
table(sign(res3$log2FoldChange[res3$padj <= .01 & !is.na(res3$padj) & abs(res3$log2FoldChange)>=2]))

### =============6. Treatment comparison in different weeks=================
## 6.1 comparing W35 F and C
### ==============================
res4 <- results(ddsC,contrast=c("condition","W35-C","W35-F"))
head(res4)

VolcanoPlotMA(res4,alpha=0.01)
sum(res4$padj <= 0.01,na.rm=TRUE)
plot(density(res4$log2FoldChange[res4$padj <= .01 & !is.na(res4$padj)]))
table(sign(res4$log2FoldChange[res4$padj <= .01 & !is.na(res4$padj)]))
table(sign(res4$log2FoldChange[res4$padj <= .01 & !is.na(res4$padj) & abs(res4$log2FoldChange)>=2]))

### ==============================
## 6.2 comparing W36 F and C
### ==============================
res5 <- results(ddsC,contrast=c("condition","W36-C","W36-F"))
head(res5)

VolcanoPlotMA(res5,alpha=0.01)
sum(res5$padj <= 0.01,na.rm=TRUE)
plot(density(res5$log2FoldChange[res5$padj <= .01 & !is.na(res5$padj)]))
table(sign(res5$log2FoldChange[res5$padj <= .01 & !is.na(res5$padj)]))
table(sign(res5$log2FoldChange[res5$padj <= .01 & !is.na(res5$padj) & abs(res5$log2FoldChange)>=2]))

### ==============================
## 6.3 comparing W32 F and W34 F (after fertilisation if off)
### ==============================
res6 <- results(ddsC,contrast=c("condition","W32-F","W34-F"))
head(res6)

VolcanoPlotMA(res6,alpha=0.01)
sum(res6$padj <= 0.01,na.rm=TRUE)
plot(density(res6$log2FoldChange[res6$padj <= .01 & !is.na(res6$padj)]))
table(sign(res6$log2FoldChange[res6$padj <= .01 & !is.na(res6$padj)]))
table(sign(res6$log2FoldChange[res6$padj <= .01 & !is.na(res6$padj) & abs(res6$log2FoldChange)>=2]))



