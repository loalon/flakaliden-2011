setwd("~/Git/flakaliden-2011/others")

# sup table from Pereira et al., 2018


miniFS <- function(X, y, ntree=2000, rounds=10){
  
  require(randomForest)
  
  # TODO check input types
  
  X <- scale(X) # 
  y <- as.factor(y) # RF needs this as a factor
  
  # execute rounds and store them in a list
  # TODO parallelism
  rf_list <- lapply(seq(1, rounds), function(x) { 
    randomForest(x=X, 
                 y=y, 
                 ntree=ntree,
                 importance=T)
  })
  
  # Per RF round, get the ranking and sort it
  importance_list <- lapply(rf_list, function(x) { 
    fimp <- importance(x, type=1)
    res <- fimp[order(fimp, decreasing = T),]
    res[res>0]
  })

  # create a list of intersected candidates, it starts with all candidates
  # each intersection will decrease the output  
  finalGOI <- colnames(X)
  
  for (item in importance_list) {
    finalGOI <- intersect(finalGOI, names(item))
  }
  
  # final round of RF with intersected candidates
  final_rf <- randomForest(x=X[, finalGOI], 
                           y=y, 
                           ntree = ntree,
                           importance = T)
  
  final_imp <- importance(final_rf, type=1)
  final_imp <- final_imp[order(final_imp, decreasing = T),]
  
  return(final_imp)
}

# CENOCOCCUM
load("/mnt/picea/projects/spruce/vhurry/fungi-flakaliden-2011/fungi-species/cenococcum_geophilum/DE/combinedVsdData.RData")
cenoData <- combinedData

X_ceno <- scale(cenoData)
y_ceno <- as.factor(substr(rownames(X_ceno), 4, 4))

cenococcumRank <- miniFS(X_ceno, y_ceno)
write.table(cenococcumRank, file="cenococcum_rank.tsv", sep="\t", col.names = F, quote=F)

# PILODERMA
load("/mnt/picea/projects/spruce/vhurry/fungi-flakaliden-2011/fungi-species/piloderma_croceum/DE/combinedVsdData.RData")
piloData <- combinedData

X_pilo <- scale(piloData)
y_pilo <- as.factor(substr(rownames(X_pilo), 4, 4))

pilodermaRank <- miniFS(X_pilo, y_pilo)
write.table(pilodermaRank, file="piloderma_rank.tsv", sep="\t", col.names = F, quote=F)

#CORTINARIUS
load("/mnt/picea/projects/spruce/vhurry/fungi-flakaliden-2011/fungi-species/cortinarius_glaucopus/DE/combinedVsdData.RData")
cortiData <- combinedData

X_corti <- scale(cortiData)
y_corti <- as.factor(substr(rownames(X_corti), 4, 4))

cortinariusRank <- miniFS(X_corti, y_corti)
write.table(cortinariusRank, file="cortinarius_rank.tsv", sep="\t", col.names = F, quote=F)

#SPRUCE
load("/mnt/picea/projects/spruce/vhurry/spruce-flakaliden-root-2011/DE/combinedVsdData.RData")
spruceData <- combinedData

X_spruce <- scale(spruceData)
y_spruce <- as.factor(substr(rownames(X_spruce), 5, 5))

spruceRank <- miniFS(X_spruce, y_spruce)
write.table(spruceRank, file="spruce_rank.tsv", sep="\t", col.names = F, quote=F)
