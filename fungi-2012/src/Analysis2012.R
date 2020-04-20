#
# this new version includes the changes John made in may, with fungi exclusive KO
library(here)
library(DESeq2)
library("vsn")
source("~/Git/Rtoolbox/src/utilsDE.r")
source("~/Git/Rtoolbox/src/plotUpAndDown.R")

projectDir <- "/mnt/picea/projects/spruce/vhurry/fungi-flakaliden-2012"
dataDir <- file.path(projectDir, "data")

###Import data
###
data <- read.delim(file.path(dataDir, "kos.raw.tsv"), header = T, sep='\t', stringsAsFactors = F, row.names = 1)
#data = read.delim("kos.raw_2012.tsv")

#there is a tow of 'unclassified data', it should be removed
data <- data[!rownames(data)=='Unclassified',]

#meta = read.csv2("Roots (1).csv")
meta <- read.delim(file.path(dataDir, "meta.tsv"), header = T, sep=';')

# order the levels for DE, both for treatment and date
print(levels(meta$treatment))
meta$treatment <- factor(meta$treatment, levels = c("Control","5_year","25_year" ))
meta$date <- factor(meta$date, levels = c("Early_June", "Late_June", "August", "October"))
meta$group <- as.factor(paste(meta$treatment, meta$date, sep="_"))
meta$group <- factor(meta$group, levels = c("Control_Early_June","Control_Late_June","Control_August","Control_October",
                                            "5_year_Early_June",  "5_year_Late_June", "5_year_August", "5_year_October",
                                            "25_year_Early_June","25_year_Late_June", "25_year_August", "25_year_October"))


#Filter out not used samples
data_rel = data[,which(colnames(data) %in% meta$SciLifeID)]
#OPTIONAL, clean rownames, to remove extra information and keep only Kxxxxx format for clarity
#rownames(data_rel) <- substr(rownames(data_rel), 1, 6)

#Design. We are going to look at the DE between treaments

#Check that samples are in cols and KO or genes in rows
dim(data_rel) #5553 KO in rows, 107 samples in cols

#import data as dds object
dds <- DESeqDataSetFromMatrix(data_rel, colData=meta, design= ~group)

#clean the data by looking for 0 expressed genes or KO
#82 KO have 0 expression
dds <- dds[(rowSums(counts(dds)) > 0),]

#QA
#normalize the data to do the QA, blind is true, no model awareness for QA
vsd.QA <- vst(dds, blind=TRUE)

#plotPCA by default only uses 500 genes, let's use all of them
geneNumber <- length(rownames(assay(vsd.QA)))

#now quick PCA plot
plotPCA(vsd.QA, intgroup = c("treatment"), ntop=geneNumber)
plotPCA(vsd.QA, intgroup = c("date"), ntop=geneNumber)
plotPCA(vsd.QA, intgroup = c("group"), ntop=geneNumber)
#there seems a more clear treatment effect than seasonal effect

#OPTIONAL - from PCA no obvious outliers, but we can check mean/sd distribution
meanSdPlot(assay(vsd.QA))

# PCA data extraction
pca <- prcomp(t(assay(vsd.QA)))
percents <- round(summary(pca)$importance[2,]*100)
summary(pca)

#OPTIONAL scree plot
screePlot <- function(pca, number=10) {
  proportionvariances <- ((pca$sdev^2) / (sum(pca$sdev^2)))*100
  barplot(proportionvariances[1:number], cex.names=1, 
          xlab=paste("Principal component (PC), 1-", number), 
          ylab="Proportion of variation (%)", main="Scree plot", ylim=c(0,100))
}
screePlot(pca, 10)

#OPTIONAL other lots, biplot, 3dPCA, etc...


##DE
#for starters, design is only the treatment
#design(dds) <- ~treatment

#execute dds
dds <- DESeq(dds)

#results from DE
resultsNames(dds)

controlVector <- as.character(unique(meta$group[grep("Control", meta$group)]))
fertilised5Vector <- as.character(unique(meta$group[grep("^5", meta$group)]))
fertilised25Vector <- as.character(unique(meta$group[grep("25", meta$group)]))



fert5_vs_control.res <- mapply(getRes, fertilised5Vector, 
                      controlVector, 
                      MoreArgs = list(localDDS=dds, group="group"))

fert25_vs_control.res <- mapply(getRes, fertilised25Vector, 
                               controlVector, 
                               MoreArgs = list(localDDS=dds, group="group"))

fert25_vs_fert5.res <- mapply(getRes, fertilised25Vector, 
                              fertilised5Vector, 
                               MoreArgs = list(localDDS=dds, group="group"))


names(fert5_vs_control.res) <- c("5_year_Early_June vs Control_Early_June", 
                                 "5_year_Late_June vs Control_Late_June",                        
                                 "5_year_August vs Control_August", 
                                 "5_year_October vs Control_October")

names(fert25_vs_control.res) <- c("25_year_Early_June vs Control_Early_June", 
                                 "25_year_Late_June vs Control_Late_June",                        
                                 "25_year_August vs Control_August", 
                                 "25_year_October vs Control_October")

names(fert25_vs_fert5.res) <- c("25_year_Early_June vs 5_year_Early_June", 
                                 "25_year_Late_June vs 5_year_Late_June",                        
                                 "25_year_August vs 5_year_August", 
                                 "25_year_October vs 5_year_October")


fert5_vs_control.res.filter <- lapply(fert5_vs_control.res, filterDE, p=0.05)
fert25_vs_control.res.filter <- lapply(fert25_vs_control.res, filterDE, p=0.05)
fert25_vs_fert5.res.filter <- lapply(fert25_vs_fert5.res, filterDE, p=0.05)

names(fert5_vs_control.res.filter) <- c("5T1 vs CT1", 
                                 "5T2 vs CT2",                        
                                 "5T3 vs CT3", 
                                 "5T4 vs CT4")

names(fert25_vs_control.res.filter) <- c("25T1 vs CT1", 
                                         "25T2 vs CT2",                        
                                         "25T3 vs CT3", 
                                         "25T4 vs CT4")

names(fert25_vs_fert5.res.filter) <- c("25T1 vs 5T1", 
                                       "25T2 vs 5T2",                        
                                       "25T3 vs 5T3", 
                                       "25T4 vs 5T4")

plotUpAndDown(fert5_vs_control.res.filter)
plotUpAndDown(fert25_vs_control.res.filter)
plotUpAndDown(fert25_vs_fert5.res.filter)

de4network(fert5_vs_control.res, here("results", "fert5_vs_control.tsv"))
de4network(fert25_vs_control.res, here("results", "fert25_vs_control.tsv"))
de4network(fert25_vs_fert5.res, here("results", "fert25_vs_fert5.tsv"))

##EXAMPLES

# source("~/Rtoolbox/plotEnrichedTreemap.R")
# source("~/Git/UPSCb/src/R/gopher.R")
# 
# up <- fert25_vs_control.res.filter$`25T3 vs CT3`[fert25_vs_control.res.filter$`25T3 vs CT3`$log2FoldChange>0,]
# down <- fert25_vs_control.res.filter$`25T3 vs CT3`[fert25_vs_control.res.filter$`25T3 vs CT3`$log2FoldChange<0,]
# up.enr <- gopher(rownames(up), task=list('ko_pathway'), url='ko')
# down.enr <- gopher(rownames(down), task=list('ko_pathway'), url='ko')
# plotEnrichedTreemap(up.enr, enrichment = 'ko_pathway', de='up')
# plotEnrichedTreemap(down.enr, enrichment = 'ko_pathway', de='down')



#read.table("~/fungi2012/fert5_vs_control.tsv", sep='\t', header = T)
## extract time inside each conditon??

#temp <- read.table("~/fungi2012/fert5_vs_control.tsv")
#result examples 5year vs. control

# res5 <- results(dds, contrast= c("treatment", "5_year", "Control"),
#                              filter = rowMedians(counts(dds)),
#                              parallel = TRUE)

#TASK, filter by lfc, padj, read the Schurch paper for filter values only with 3 replicates

#TASK plot the results

#OPTIONAL repeat DE with time and maybe with time and treatment

#If we are happy with the final results, we make a vst model aware and export that
vsd <- vst(dds, blind=FALSE)

#export the results in the proper format (samples in rows)
write.table(t(assay(vsd)), file="vstData.tsv", sep="\t", col.names = NA, quote = F)

#export for seidr
#combined design = group

#control


#####ADITIONAL
#####
#extract 25vs C genes, only down regulated
ko_25T1vsCT1 <- rownames(fert25_vs_control.res.filter$`25T1 vs CT1`[fert25_vs_control.res.filter$`25T1 vs CT1`$log2FoldChange<0,])
ko_25T2vsCT2 <- rownames(fert25_vs_control.res.filter$`25T2 vs CT2`[fert25_vs_control.res.filter$`25T2 vs CT2`$log2FoldChange<0,])
ko_25T3vsCT3 <- rownames(fert25_vs_control.res.filter$`25T3 vs CT3`[fert25_vs_control.res.filter$`25T3 vs CT3`$log2FoldChange<0,])
ko_25T4vsCT4 <- rownames(fert25_vs_control.res.filter$`25T4 vs CT4`[fert25_vs_control.res.filter$`25T4 vs CT4`$log2FoldChange<0,])

write.table(ko_25T1vsCT1, file=here("results", "ko_25T1vsCT1.tsv"), sep="\t", col.names = F, quote = F, row.names = F)
write.table(ko_25T2vsCT2, file=here("results", "ko_25T2vsCT2.tsv"), sep="\t", col.names = F, quote = F, row.names = F)
write.table(ko_25T3vsCT3, file=here("results", "ko_25T3vsCT3.tsv"), sep="\t", col.names = F, quote = F, row.names = F)
write.table(ko_25T4vsCT4, file=here("results", "ko_25T4vsCT4.tsv"), sep="\t", col.names = F, quote = F, row.names = F)
