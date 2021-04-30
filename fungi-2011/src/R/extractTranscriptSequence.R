# Extract sequence information from specific clusters and generate fasta files

# Lets start with cenococcum geophilum
# 
# control clusters 14, 21, 28, and 4 are of high interest

library(ShortRead)

projectFolder <- "/mnt/picea/projects/spruce/vhurry/fungi-flakaliden-2011/transcriptSequences"

clusterFile <- read.table("/mnt/picea/projects/spruce/vhurry/fungi-flakaliden-2011/fungi-species/piloderma_croceum/networks/control/results/aggregated/infomapClusters.tsv",
                          sep='\t', header=T)

translationFile <- read.table(file.path(projectFolder, "final.reformat.gff"), sep='\t')


translationTable <- data.frame(transcriptID=translationFile$V1,
                               geneID=gsub("gene_id=", "", translationFile$V9))

#clustersOI <- c("Cluster14", "Cluster21", "Cluster4", "Cluster28")
clustersOI <- unique(clusterFile$cluster)

#fastaTranscripts <- readFasta(file.path(projectFolder, "final.fa"))

# remove irrelevant info from the names
s = readDNAStringSet(file.path(projectFolder, "final.fa"))
names(s) <- sub(" .*","",names(s))

# s[names(s) == "k117_951"]
# as.character(s[[1]])

# extract cluster gene lists
genesPerCluster <- sapply(clustersOI, function(x){
  clusterFile[clusterFile$cluster == x,]$gene  
})

# names(geneString) %in% translationTable$geneID

# get a reduced list
reducedFasta <-  s[match(translationTable$transcriptID, names(s))]

# all names match the translation table?
all(names(reducedFasta) == translationTable$transcriptID)

# if yes them we overwrite the names on the sequence data to 
names(reducedFasta) <- translationTable$geneID

clusterSequences <- sapply(genesPerCluster, function(x){
  reducedFasta[names(reducedFasta) %in% x]
})
# sanity check
#as.character(clusterSequences$Cluster64[names(clusterSequences$Cluster64) == "2011.gene297004"]) ==
#as.character(reducedFasta[names(reducedFasta) == "2011.gene297004"])

# get AA translation
clusterAA <- sapply(clusterSequences, function(x){
  translate(x)
})

#as.character(clusterAA$Cluster1[names(clusterAA$Cluster1) == "2011.gene103284"])

#as.character(translate(reducedFasta[names(reducedFasta) == "2011.gene103284"]))

sapply(names(clusterSequences), function(x){
  
  #fileConn <- file(paste0(x, ".fasta"))
  sink(paste0(x, ".fasta"))
  #print(length(clusterAA[[x]]))
  for (i in 1:length(clusterSequences[[x]])){
    currentName <- (names(clusterSequences[[x]])[i])
    cat(paste0(">", currentName))
    cat("\n")
    cat(as.character((clusterSequences[[x]][i])))
    cat("\n")
  }
  #writeFasta(clusterAA[[x]], paste0(x, ".fasta"))
  #close(fileConn)
  sink()
})

sapply(names(clusterAA), function(x){
  
  #fileConn <- file(paste0(x, ".fasta"))
  sink(paste0(x, "_AA.fasta"))
  #print(length(clusterAA[[x]]))
  for (i in 1:length(clusterAA[[x]])){
    currentName <- (names(clusterAA[[x]])[i])
    cat(paste0(">", currentName))
    cat("\n")
    cat(as.character((clusterAA[[x]][i])))
    cat("\n")
  }
  #writeFasta(clusterAA[[x]], paste0(x, ".fasta"))
  #close(fileConn)
  sink()
})

# get the full aa table
fullAA <- translate(reducedFasta)

#names(fullAA) <- translationTable[names(fullAA) %in% translationTable$transcriptID,]$geneID

dfAA <- as.data.frame(fullAA)
colnames(dfAA) <- "sequence"
dfAA$gene <- rownames(dfAA)
  
library(dplyr)
fullResult <- left_join(clusterFile, dfAA, by=c("gene"))

write.table(fullResult, 
            file="full_control_cluster_AAsequence.tsv", 
            sep='\t', 
            quote=FALSE,
            row.names = F
            )
    