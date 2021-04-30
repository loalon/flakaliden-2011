###CORRELATION ANALYSIS AUX FUNCTIONS
###
###
###
##### Auxiliary functions

getPlotData <-function(data1, genes1, group1, time1) {
  p1 <- plotEigengene(data1, genes1, group1, time1)
  p1.data <- p1$data %>%
    group_by(x) %>%
    summarise(n = n(), y = mean(y, na.rm = TRUE) )
  return(p1.data)
}

#TODO fix this using the minifunction
profileSimilarity <- function(data1, data2, clusterData1, clusterData2, 
                              group1, group2, time1, time2, name1, name2, minSize=20){
  
  # cluster1F -> cluster1S -> 0.34
  clusters1 <- unique(clusterData1$cluster)
  clusters2 <- unique(clusterData2$cluster)
  
  res <- data.frame(name1 = character(),
                    genes1 = character(),
                    name2 = character(), 
                    genes2 = character(),
                    score = double(), 
                    stringsAsFactors = FALSE) 
  
  sapply(clusters1, function(clus1){
    
    clus1Genes <- clusterData1[clusterData1$cluster==clus1,]$gene
    clus1Genes_split <- getEigenGenes(data1, clus1Genes)
    clus1i <- paste0(clus1, "i")
    lenClus1 <- length(clus1Genes_split$positive)
    lenClus1i <- length(clus1Genes_split$negative)
    
    p1 <- plotEigengene(data1, clus1Genes, group1, time1)
    p1.data <- p1$data %>%
      group_by(x) %>%
      summarise(n = n(), y = mean(y, na.rm = TRUE) )
    eigen1 <- p1.data$y
    eigen1i <- p1.data$y * -1
    
    sapply(clusters2, function(clus2){
      #print(paste("comparing",clus1, clus2))
      
      clus2Genes <- clusterData2[clusterData2$cluster==clus2,]$gene
      clus2Genes_split <- getEigenGenes(data2, clus2Genes)
      clus2i <- paste0(clus2, "i")
      lenClus2 <- length(clus2Genes_split$positive)
      lenClus2i <- length(clus2Genes_split$negative)
      
      p2 <- plotEigengene(data2,clus2Genes, group2, time2)
      p2.data <- p2$data %>%
        group_by(x) %>%
        summarise(n = n(), y = mean(y, na.rm = TRUE) )
      eigen2 <- p2.data$y
      eigen2i <- p2.data$y * -1
      
      #add comparitions
      #1vs2
      if((lenClus1 > minSize) && (lenClus2 > minSize))
        res <<- rbind(res, c(clus1, lenClus1, clus2, lenClus2, 
                             cor(eigen1, eigen2, method = c("spearman"))), stringsAsFactors = FALSE)
      #1ivs2
      if((lenClus1i > minSize) && (lenClus2 > minSize))
        res <<- rbind(res, c(clus1i, lenClus1i, clus2, lenClus2, 
                             cor(eigen1i, eigen2, method = c("spearman"))), stringsAsFactors = FALSE)
      #1vs2i
      if((lenClus1 > minSize) && (lenClus2i > minSize))
        res <<- rbind(res, c(clus1, lenClus1, clus2i, lenClus2i,
                             cor(eigen1, eigen2i, method = c("spearman"))), stringsAsFactors = FALSE)
      #1ivs2i
      if((lenClus1i > minSize) && (lenClus2i > minSize))
        res <<- rbind(res, c(clus1i, lenClus1i, clus2i, lenClus2i, 
                             cor(eigen1i, eigen2i, method = c("spearman"))), stringsAsFactors = FALSE)
    })
    
  })
  
  names(res) <- c(name1, paste0(name1,"-genes"),name2, paste0(name2,"-genes"), "score")
  res$score <- as.numeric(res$score)
  return(res)
}

# future uses
cluster2Network <- function(data, clusterData, species){
  
  temp <- lapply(unique(clusterData$cluster), function(localCluster) {
    spruceGenes <- getEigenGenes(data, clusterData[clusterData$cluster==localCluster,]$gene)
  })
  #spruceCluster <- "Cluster1"
  names(temp) <- paste0(species, unique(clusterData$cluster) )
  
  posTemp <- lapply(temp, function(x){
    x$positive
  })
  
  negTemp <- lapply(temp, function(x){
    x$negative
  })
  names(negTemp) <- paste0(names(negTemp),"i")
  
  posResult <- sapply(posTemp, function(x) {
    length(x)
  })
  negResult <- sapply(negTemp, function(x) {
    length(x)
  })
  
  posDF <- data.frame(clusters=names(posResult), geneNumber=posResult)
  negDF <- data.frame(clusters=names(negResult), geneNumber=negResult)
  res <- rbind(posDF, negDF)
  return (res)
}



getProfileDegradationMatrix <- function(pairedData, conditionSpruceClusters, conditionFungiClusters, 
                                        data1Spruce, data1Fungi,data2Spruce, data2Fungi,
                                        group1Spruce, group1Fungi,group2Spruce, group2Fungi,
                                        time1Spruce, time1Fungi,time2Spruce, time2Fungi,  
                                        
                                        conditionNames = c("control", "fertilised") ){
  
  condition.pos <- pairedData #[pairedData$score>0,]
  
  score1 <- vector(mode = "numeric", length = length(rownames(condition.pos)))
  score2 <- vector(mode = "numeric", length = length(rownames(condition.pos)))
  
  sapply(1:length(rownames(condition.pos)), function(l) {
    
    
    spruceCluster <- condition.pos[l,1]
    fungiCluster <- condition.pos[l,3]
    spruceReverse <- ifelse(endsWith(spruceCluster, "i"), TRUE, FALSE) 
    fungiReverse <- ifelse(endsWith(fungiCluster, "i"), TRUE, FALSE)
    
    spruceCluster <- ifelse(endsWith(spruceCluster, "i"), substr(spruceCluster, 1, nchar(spruceCluster)-1), spruceCluster) 
    fungiCluster <- ifelse(endsWith(fungiCluster, "i"), substr(fungiCluster, 1, nchar(fungiCluster)-1), fungiCluster) 
    
    c1SpruceClusterTotalGenes <- conditionSpruceClusters[conditionSpruceClusters$cluster==spruceCluster,]$gene
    c1FungiClusterTotalGenes <- conditionFungiClusters[conditionFungiClusters$cluster==fungiCluster,]$gene
    
    spruceGenes <- getEigenGenes(data1Spruce, c1SpruceClusterTotalGenes)
    c1SpruceClusterGenes <- if(spruceReverse) spruceGenes$negative else spruceGenes$positive
    
    fungiGenes <- getEigenGenes(data1Fungi, c1FungiClusterTotalGenes)
    c1FungiClusterGenes <- if(fungiReverse) fungiGenes$negative else fungiGenes$positive 
    
    # we extract only the genes that exist in both conditions and forget about uniques for the moment
    commonSpruceClusterGenes <- intersect(c1SpruceClusterGenes, colnames(data2Spruce))
    commonFungiClusterGenes <- intersect(c1FungiClusterGenes, colnames(data2Fungi))
    
    # TODO extract unique genes that exist in c1 but not in c2
    
    #initialise cor1 and cor2 for those cases that don't have genes in common
    cor1 <- 0
    cor2 <- 0
    if(length(commonSpruceClusterGenes)>0 && length(commonFungiClusterGenes)>0) {
      p1Spruce <- getPlotData(data1Spruce, commonSpruceClusterGenes, group1Spruce, time1Spruce)
      p1Fungi <- getPlotData(data1Fungi, commonFungiClusterGenes, group1Fungi, time1Fungi)
      
      p2Spruce <- getPlotData(data2Spruce, commonSpruceClusterGenes, group2Spruce, time2Spruce)
      p2Fungi <-getPlotData(data2Fungi, commonFungiClusterGenes, group2Fungi, time2Fungi)
      cor1 <- cor(p1Spruce$y, p1Fungi$y, method = c("spearman"))
      cor2 <- cor(p2Spruce$y, p2Fungi$y, method = c("spearman"))
    }
    
    score1[l] <<- cor1
    score2[l] <<- cor2
    
  }) #end of sapply
  
  res <- condition.pos[,1:4]
  res[[paste0(conditionNames[1],"_score")]] <- score1
  res[[paste0(conditionNames[2],"_score")]] <- score2
  res$degradation <- score1 - score2
  res.sort <- res[order(res$degradation, decreasing = T),]
  
  #write.table(fertilised2control.sort, file="fertilised2control.sort.tsv", quote = F, sep='\t', row.names = F)
  #
  res.sort
  
}

enrichThis <- function(data1, data2, clusterData1, clus1, background1, url ='pabies', type="go") {
  
  eigen <- getEigenGenes(data1, clusterData1[clusterData1$cluster==gsub("i","", clus1),]$gene)
  
  genes1 <- if(grepl("i", clus1, fixed = TRUE)){
    eigen$negative
  } else {
    eigen$positive
  }
  
  genes2 <- intersect(genes1, colnames(data2))
  genes3 <- intersect(genes1, genes2)
  
  if(length(genes3)>1) {
    enr1 <- gopher(genes3, background=background1, url=url, task = type)
  } else {
    enr1 <- NULL
  }
  return(enr1)
}


plotThis <- function(data1, data2, clusterData1, group1, time1, group2, time2, clus1, title1="", title2="") {
  
  results <- list()
  eigen <- getEigenGenes(data1, clusterData1[clusterData1$cluster==gsub("i","", clus1),]$gene)
  
  genes1 <- if(grepl("i", clus1, fixed = TRUE)){
    eigen$negative
  }else {
    eigen$positive
  }
  
  genes2 <- intersect(genes1, colnames(data2))
  
  p1 <- plotEigengene(data1, genes1, group1, time1, title=title1)
  p2 <- plotEigengene(data2, genes2, group2, time2, title=title2, colors="#BE9230")
  genes3 <- intersect(genes1, genes2)
  #data3 <- rbind(data1[,genes3], data2[,genes3])
  
  # p3 <- plotEigengene(data3, genes3, c(group1, group2),
  #                     c(time1,time2), title="", noLegend = F)
  
  if(length(genes3)>0) {
    plotData1 <- getPlotData(data1, genes3, group1, time1)
    plotData2 <- getPlotData(data2, genes3, group2, time2)
    cor1 <- cor(plotData1$y, plotData2$y, method = c("spearman"))
  } #else { 0 }
  
  results$plot1 <- p1
  results$plot2 <- p2
  results$cor1 <- cor1
  #results$plot3 <- p3
  
  return(results)
}
