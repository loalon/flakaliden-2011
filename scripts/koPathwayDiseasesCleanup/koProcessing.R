

# diseases <- c("map05200","map05202","map05206","map05205","map05204","map05203","map05230","map05231","map05235","map05210",
#               "map05212","map05225","map05226","map05214","map05216","map05221","map05220","map05217","map05218","map05211",
#               "map05219","map05215","map05213","map05224","map05222","map05223","map05310","map05322","map05323","map05320",
#               "map05321","map05330","map05332","map05340","map05010","map05012","map05014","map05016","map05020","map05030",
#               "map05031","map05032","map05033","map05034","map05418","map05410","map05412","map05414","map05416","map04930",
#               "map04940","map04950","map04932","map04931","map04933","map04934","map05110","map05120","map05130","map05132",
#               "map05131","map05135","map05133","map05134","map05150","map05152","map05100","map05166","map05170","map05162",
#               "map05164","map05161","map05160","map05168","map05163","map05167","map05169","map05165","map05146","map05144",
#               "map05145","map05140","map05142","map05143","map01501","map01502","map01503","map01521","map01524","map01523",
#               "map01522","map03460")

koTranslate <- function(enrData, name="name") {
  require(KEGGREST)
  require(dplyr)
  
  if(!is.data.frame(enrData))
    stop("enrData is not a data.frame")
  
  #there is a limit of keggLink querys
  querySize <- 100
  
   if(length(rownames(querySize)) > querySize) { 
     
     enrData <- enrData[1:100,]
     # splits <- ceiling(nrow(enrData))
     # tempSplit <- split(enrData, sample(1:splits, nrow(enrData), replace=T))
     # 
     # lapply(tempSplit) <- function(x){
     #   x$name <- keggList(x$name)
     # }
     
   }
  enrData <- filter(enrData, def !='UNKNOWN')
  enrData[[name]] <- ifelse(enrData[[name]] != 'None', keggList(enrData[[name]]), 'None') 
  
  return(enrData)
}


koPathwayDiseaseCleanup <- function(enrData, name="name", id="id", diseasesFile="~/scripts/koPathwayDiseasesCleanup/diseases2ko.txt") {
  diseases <- read.table(diseasesFile, header=F)
  diseases <- as.character(diseases$V1)
  
  enrData[[name]][(enrData[[id]] %in% diseases)] <- "mis-annotated"  
  return(enrData)
}
  
# TEST  
# source('~/Git/UPSCb-common/src/R/gopher.R', echo=TRUE)
# enr <- gopher(c("K00001","K00002"),url='ko', task=list('ko_pathway'))
# 
# enr$ko_pathway <- koPathwayDiseaseCleanup(enr$ko_pathway)
# enr$ko_pathway <- koPathwayDiseaseCleanup(enr$ko_pathway)