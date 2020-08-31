#taxo4funguild.tsv
library(dplyr)
load(file.path(deFolder, "ddsTreatment.RData"))
load(file.path(dataFolder, "megaData.RData"))
#taxo <- read.table(file.path(dataFolder, "taxo4funguild.tsv"), header=T, quote = "" ,sep='\t')

#COMBINED stuff

temp <- data.frame(species=megaData$species, 
                   mode=megaData$Trophic.Mode,
                   confidence=megaData$Confidence.Ranking,
                   guild=megaData$Guild)

temp <- temp %>% dplyr::filter(confidence %in% c("Highly Probable","Probable"))

deduped.data <- unique( temp[ , 1:2 ] )
finalCombined <- data.frame(table(deduped.data$mode))

guild.deduped.data <- unique( temp[ , c(1,2,4) ] )


guild.finalCombined <- table(guild.deduped.data[,2:3])

write.table(guild.finalCombined, "trophicMode_guild_Proportions_combined.tsv", sep='\t', row.names = T, quote=F, col.names = NA)

resDF <- dplyr::left_join(megaData, 
                 data.frame(resFungiTreatment, gene=rownames(data.frame(resFungiTreatment))),
                 by="gene")

#FERT stuff

resFert <- data.frame(species=resDF[resDF$FertilisedvsControl.log2FoldChange>0,]$species, 
                      mode=resDF[resDF$FertilisedvsControl.log2FoldChange>0,]$Trophic.Mode,
                      confidence=resDF[resDF$FertilisedvsControl.log2FoldChange>0,]$Confidence.Ranking)

resFert <- resFert %>% dplyr::filter(confidence %in% c("Highly Probable","Probable"))

deduped.fert <- unique( resFert[ , 1:2 ] )
finalFert <- data.frame(table(deduped.fert$mode))

#CONTROL
resControl <- data.frame(species=resDF[resDF$FertilisedvsControl.log2FoldChange<0,]$species, 
                      mode=resDF[resDF$FertilisedvsControl.log2FoldChange<0,]$Trophic.Mode,
                      confidence=resDF[resDF$FertilisedvsControl.log2FoldChange<0,]$Confidence.Ranking)

resControl <- resControl %>% dplyr::filter(confidence %in% c("Highly Probable","Probable"))
deduped.control <- unique( resControl[ , 1:2 ] )
finalControl <- data.frame(table(deduped.control$mode))

finalTable <- dplyr::left_join(dplyr::left_join(finalCombined, finalControl, by="Var1"), finalFert, by="Var1")
colnames(finalTable) <- c("Mode", "All", "ND_de", "NE_de")
finalTable

write.table(finalTable, "trophicModeProportions.tsv", sep='\t', row.names = F, quote=F)


# miniData  <- megaData[,2:112]
# 
# miniData$X19 <- floor(rowMeans(miniData[,grepl("19", colnames(miniData))]))
# miniData$X20 <- floor(rowMeans(miniData[,grepl("20", colnames(miniData))]))
# miniData$X21 <- floor(rowMeans(miniData[,grepl("21", colnames(miniData))]))
# miniData$X22 <- floor(rowMeans(miniData[,grepl("22", colnames(miniData))]))
# miniData$X23 <- floor(rowMeans(miniData[,grepl("23", colnames(miniData))]))
# miniData$X24 <- floor(rowMeans(miniData[,grepl("24", colnames(miniData))]))
# miniData$X25 <- floor(rowMeans(miniData[,grepl("25", colnames(miniData))]))
# miniData$X26 <- floor(rowMeans(miniData[,grepl("26", colnames(miniData))]))
# miniData$X28 <- floor(rowMeans(miniData[,grepl("28", colnames(miniData))]))
# miniData$X31 <- floor(rowMeans(miniData[,grepl("31", colnames(miniData))]))
# miniData$X32 <- floor(rowMeans(miniData[,grepl("32", colnames(miniData))]))
# miniData$X34 <- floor(rowMeans(miniData[,grepl("34", colnames(miniData))]))
# miniData$X35 <- floor(rowMeans(miniData[,grepl("35", colnames(miniData))]))
# miniData$X36 <- floor(rowMeans(miniData[,grepl("36", colnames(miniData))]))
# miniData$X37 <- floor(rowMeans(miniData[,grepl("37", colnames(miniData))]))
# miniData$X38 <- floor(rowMeans(miniData[,grepl("38", colnames(miniData))]))
# miniData$X39 <- floor(rowMeans(miniData[,grepl("39", colnames(miniData))]))
# miniData$X40 <- floor(rowMeans(miniData[,grepl("40", colnames(miniData))]))
# miniData$X41 <- floor(rowMeans(miniData[,grepl("41", colnames(miniData))]))
# 
# miniData$ND <- floor(rowMeans(miniData[,grepl("C", colnames(miniData))]))
# miniData$NE <- floor(rowMeans(miniData[,grepl("F", colnames(miniData))]))
# miniData$TOTAL <- rowMeans(miniData[,112:130])
# 
# maxiData <- cbind(megaData, miniData)
# 
# sum(maxiData[maxiData$Confidence.Ranking %in% c("Highly Probable","Probable"),]$NE)
# 
# sum(maxiData[maxiData$Confidence.Ranking %in% c("Highly Probable","Probable"),]$ND)
# sum(maxiData[maxiData$Confidence.Ranking %in% c("Highly Probable","Probable"),]$TOTAL)
