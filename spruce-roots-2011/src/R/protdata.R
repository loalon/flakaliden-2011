### ==============================
## load the necessary libraries
### ==============================
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(vsn))
suppressPackageStartupMessages(library(scatterplot3d))
suppressPackageStartupMessages(library(LSD))
source("~/UPSCb/src/R/rmd.R")

### ==============================
## set the working directory
### ==============================
setwd("~/UPSCb/projects/spruce-roots/src/R/")

### ==============================
## read the samples details 
## data is free amino acid as umol/g DW (5 reps for each time and treatment)
### ==============================
samples <- read.csv("~/UPSCb/projects/spruce-roots/doc/prot_amino.csv", sep = ";")

### ==============================
## creating vectors
### ==============================
names <- samples[,1]
conditions <- sub("[C,F][1-6]?","",samples[,1])
cf <- substr(x=samples[,1],4,4)

### =============================
## scatterplot matrices
## ==============================
pairs(samples[2:18])

### ==============================
## Perform a Principal Component Analysis
## concetrations are standardised first
### ==============================
stdconc <- as.data.frame(scale(samples[2:18]))
pca <- prcomp(stdconc)
summary(pca)
percent <- round(summary(pca)$importance[2,]*100)
percent
pchs <- LETTERS[as.integer(factor(conditions))]
cols <- c("darkred","darkgreen")[as.integer(factor(cf))]
PCAlegend <- data.frame(pchs,conditions)

### ==============================
## plot the first 3 dimensions
## First component: 86%
## Second component: 7%
## Third component: 3%
### ==============================
### ==============================
## weeks are represented by letters (A-S)
## red: Control, green: fertilized
### ==============================
scatterplot3d(pca$x[,1],
              pca$x[,2],
              pca$x[,3],
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
plot(pca$x[,1],
     pca$x[,2],
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
plot(pca$x[,2],
     pca$x[,3],
     xlab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     ylab=paste("Comp. 3 (",percent[3],"%)",sep=""),
     col=cols,pch=pchs,
     main="Principal Component Analysis")
legend("topright",
       legend=unique(PCAlegend$pchs:PCAlegend$conditions),
       cex=0.4)
