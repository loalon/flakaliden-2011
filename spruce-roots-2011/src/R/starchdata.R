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
## data is starch as %/g DW (5 reps for each time and treatment)
### ==============================
samples <- read.csv("~/UPSCb/projects/spruce-roots/doc/starch.csv", sep = ";")

### ==============================
## creating vectors
### ==============================
names <- samples[,1]
conditions <- sub("[C,F][1-6]","",samples[,1])
cf <- substr(x=samples[,1],4,4)
reps <- substr(x=samples[,1],1,4)

### ==============================
## Perform a Principal Component Analysis
## concetrations are standardised first
### ==============================
stdconc <- as.data.frame(scale(samples[2]))
pca <- prcomp(stdconc)
summary(pca)
percent <- round(summary(pca)$importance[2,]*100)
percent
pchs <- LETTERS[as.integer(factor(conditions))]
cols <- c("darkred","darkgreen")[as.integer(factor(cf))]
PCAlegend <- data.frame(pchs,conditions)

### ==============================
## There is only one component
## First component: 100%
## weeks are represented by letters (A-S)
## red: Control, green: fertilized
### ==============================
plot(pca$x[,1],
              xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
              col=cols,pch=pchs)
legend("topright",
       legend=unique(PCAlegend$pchs:PCAlegend$conditions),
       cex=0.4)
