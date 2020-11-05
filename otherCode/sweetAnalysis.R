library("devtools")
devtools::install_github("loalon/rBiobox")
library(rBiobox)
plotEigengene(toyData$expressionData, toyData$geneCluster, toyData$time, toyData$conditions)

sweetGenes <- c("MA_10064529g0010","MA_101378g0010","MA_10364691g0010",
"MA_10385255g0010","MA_10431590g0010","MA_10435683g0020",
"MA_10657g0010","MA_1083g0010","MA_111938g0010",
"MA_165733g0010","MA_17102g0010","MA_17777g0010",
"MA_178397g0010","MA_178981g0010","MA_18013g0010",
"MA_180850g0010","MA_181386g0010","MA_185837g0010",
"MA_188016g0010","MA_194437g0010","MA_199290g0010","MA_203675g0010",
"MA_2249g0010","MA_275250g0010","MA_317333g0010","MA_320782g0010",
"MA_322191g0010","MA_34530g0010","MA_36895g0010","MA_405662g0010",
"MA_420390g0010","MA_43131g0010","MA_44746g0010","MA_4590656g0010",
"MA_459631g0010","MA_481649g0010",
"MA_48670g0010","MA_55437g0010",
"MA_60853g0010","MA_61830g0010",
"MA_689129g0010","MA_69041g0010",
"MA_710181g0010","MA_72344g0010",
"MA_775888g0010","MA_776914g0010",
"MA_81780g0010","MA_8293037g0010",
"MA_8356859g0010","MA_854569g0010",
"MA_8809860g0010","MA_889646g0010",
"MA_91381g0010","MA_930638g0010",
"MA_93514g0020","MA_94599g0010")

load("/mnt/picea/projects/spruce/vhurry/14_SpruceRoots_Project/DE/combinedVsdData.RData")

commonGenes <- intersect(colnames(combinedData), sweetGenes)

plotEigengene(combinedData, 
              commonGenes,
              as.numeric(substr(rownames(combinedData), 2, 3)),
              substr(rownames(combinedData), 5, 5))


load("/mnt/picea/projects/spruce/vhurry/14_SpruceRoots_Project/DE/vsd_QA.RData")

rawData <- t(assay(vsd.QA))

nada <- plotEigengene(rawData, 
              sweetGenes,
              
              substr(rownames(rawData), 5, 5), as.numeric(substr(rownames(rawData), 2, 3)),multiline = T)

miniData <- rawData[,sweetGenes]

source("~/Git/Rtoolbox/src/plotEigenGene.R")


newTable <- read.table("sweetdata.tsv", header = T, sep='\t')
rownames(newTable) <- newTable$Gene
newTable <- newTable[,-1]
newTable <- t(newTable)

rBiobox::plotEigengene(newTable, 
              colnames(newTable),
              as.numeric(substr(rownames(newTable), 2, 3)),
              substr(rownames(newTable), 4, 4))

table(colMads(newTable) >0)


ndTable <- newTable[1:19,]
colMads[, which(colMads(ndTable) >0)]

ndTable <- ndTable[,which(colMads(as.matrix(ndTable)) > 0)]

plotEigengene(ndTable, 
              intersect(sweetGenes, colnames(ndTable)),
              
              substr(rownames(ndTable), 4, 4), 
              as.numeric(substr(rownames(ndTable), 2, 3)),multiline = T)

plotEigengene(ndTable, 
              intersect(sweetGenes, colnames(ndTable)),
              
              substr(rownames(ndTable), 4, 4), 
              as.numeric(substr(rownames(ndTable), 2, 3)),multiline = F)

neTable <- newTable[20:38,]
colMads[, which(colMads(neTable) >0)]

neTable <- neTable[,which(colMads(as.matrix(neTable)) > 0)]

plotEigengene(neTable, 
              intersect(sweetGenes, colnames(neTable)),
              
              substr(rownames(neTable), 4, 4), 
              as.numeric(substr(rownames(neTable), 2, 3)),multiline = T)

plotEigengene(neTable, 
              intersect(sweetGenes, colnames(neTable)),
              
              substr(rownames(neTable), 4, 4), 
              as.numeric(substr(rownames(neTable), 2, 3)),multiline = F, inverse = T)
