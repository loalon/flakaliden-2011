library(dplyr)
library(here)
library(reshape2)
source("~/Git/Rtoolbox/src/plotEigenGene.R")
source("~/Git/Rtoolbox/src/getEigengenes.R")
source("~/Git/Rtoolbox/src/plotVectorPCA.R")

spruceDir <- "/mnt/picea/projects/spruce/vhurry/14_SpruceRoots_Project"
fungiDir <- "/mnt/picea/projects/spruce/vhurry/fungi-flakaliden-2011"

spruceDE <- file.path(spruceDir,"DE")
fungiDE <- file.path(fungiDir,"DE")

load(file.path(spruceDE, "allDEresults.RData"))
spruce.combined.res <- combined.res

load(file.path(fungiDE, "allDEresults.RData"))
fungi.combined.res <- combined.res

names(spruce.combined.res) <- paste0("spruce",names(spruce.combined.res))
names(fungi.combined.res) <- paste0("fungi",names(fungi.combined.res))

ifelse(is.na(spruce.combined.res$spruce19Fvs19C$padj), 1, spruce.combined.res$spruce19Fvs19C$padj)
spruce.combined.res$spruce19Fvs19C[spruce.combined.res$spruce19Fvs19C$padj<0.01,]

spruceDF <- lapply(spruce.combined.res, function(x){
  x$padj <- ifelse(is.na(x$padj), 1, x$padj)
  #x <- x[x$padj<0.01,]
  
  x[x$log2FoldChange>0,]$log2FoldChange
  })
fungiDF <- lapply(fungi.combined.res, function(x){
  x$padj <- ifelse(is.na(x$padj), 1, x$padj)
  #x <- x[x$padj<0.01,]
  
  x[x$log2FoldChange>0,]$log2FoldChange
  })
# 
# spruceDF2 <-  as.data.frame(spruceDF)
# fungiDF2 <-  as.data.frame(fungiDF) 
# 
# boxplot(spruceDF)
# boxplot(fungiDF)

mergedList <- list()
for(i in 1:19){
  x.name<- names(spruceDF)[i]
  y.name<- names(fungiDF)[i]
  
  mergedList[[x.name]] <- spruceDF[[i]]
  mergedList[[y.name]] <- fungiDF[[i]]
}
temp <- melt(mergedList)

myColors <- ifelse(startsWith(names(mergedList),"spruce"), "#00FF00", "#FF0000")
                              

# boxplot(mergedList,ylim=c(0,5), col=myColors )
# 
# ggplot(melt(mergedList[3:4]), aes(x=L1, y=value  )) + 
#   #geom_violin() +
#   stat_summary(fun.y=mean, geom="point", shape=23, size=2) + 
#   stat_summary(fun.y=median, geom="point", size=2, color="red") +
#   geom_boxplot(width=0.1) +
#   ylim(-5, 5)
# 
# ggplot(melt(fungiDF[1:3]), aes(x=L1, y=value  )) + 
#   #geom_violin() +
#   stat_summary(fun.y=mean, geom="point", shape=23, size=2) + 
#   stat_summary(fun.y=median, geom="point", size=2, color="red") +
#   geom_boxplot(width=0.1) +
#   ylim(-5, 5)
# 
# spruceDF3 <- melt(spruceDF2[,c(1:3)])
# fungiDF3 <- melt(fungiDF2[,c(1:3)])

# ggplot(spruceDF3, aes(x=variable, y=value  )) + 
#   #geom_violin() +
#   stat_summary(fun.y=mean, geom="point", shape=23, size=2) + 
#   stat_summary(fun.y=median, geom="point", size=2, color="red") +
#   geom_boxplot(width=0.1) +
#   ylim(-5, 5)
# 
# ggplot(fungiDF3, aes(x=variable, y=value  )) + 
#   #geom_violin() +
#   stat_summary(fun.y=mean, geom="point", shape=23, size=2) + 
#   stat_summary(fun.y=median, geom="point", size=2, color="red") +
#   geom_boxplot(width=0.1) +
#   ylim(-5, 5)


load("/mnt/picea/projects/spruce/vhurry/14_SpruceRoots_Project/DE/ddsTreatment.RData")
load("/mnt/picea/projects/spruce/vhurry/fungi-flakaliden-2011/DE/ddsTreatment.RData")

resFungiTreatment$log2FoldChange <- ifelse(is.na(resFungiTreatment$log2FoldChange), 0, resFungiTreatment$log2FoldChange)

resSpruceTreatment$log2FoldChange <- ifelse(is.na(resSpruceTreatment$log2FoldChange), 0, resSpruceTreatment$log2FoldChange)


treatmentList <- list(SpruceFvsC =resSpruceTreatment[resSpruceTreatment$log2FoldChange>0,]$log2FoldChange,
                      FungiFvsC=resFungiTreatment[resFungiTreatment$log2FoldChange>0,]$log2FoldChange)
treatmentList <- list(SpruceFvsC =resSpruceTreatment$log2FoldChange,
                      FungiFvsC=resFungiTreatment$log2FoldChange
                      )
boxplot(treatmentList,ylim=c(-5,5), col=c("#00FF00", "#FF0000"), ylab="log2FoldChange" )

ggplot(melt(treatmentList), aes(x=L1, y=value, col=L1  )) +
  #geom_violin() +
  #stat_summary(fun.y=mean, geom="point", shape=23, size=2) +

  geom_boxplot(width=0.1) +
  stat_summary(fun.y=median, geom="point", size=2, color="black") +
  ylim(-5, 5)
#   
