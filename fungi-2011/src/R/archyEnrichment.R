##################################################
## Project: Flakaliden fungi 2011
## Script purpose: Obatin enrichment for the species Archaeorhizomycetaceae
## and generate treemaps for the enrichment results
## Date: 200403
## Author: Alonso Serrano alonso.serrano@slu.se
##################################################


library(here)
source(here("UPSCb-common/src/R/gopher.R"))
source(here("Rtoolbox/src/plotEnrichedTreemap.R"))

projectFolder <- "/mnt/picea/projects/spruce/vhurry/fungi-flakaliden-2011"
#deFolder <- file.path(projectFolder, "DE")
dataFolder <- file.path(projectFolder, "data")
load(file.path(dataFolder, "megaData.RData") )


megaData[megaData$genus == "Archaeorhizomyces",]
megaData[megaData$family == "Archaeorhizomycetaceae",]
  
a.finlayi <- megaData[megaData$species == "Archaeorhizomyces finlayi",]$gene
a.borealis <- megaData[megaData$species == "Archaeorhizomyces borealis",]$gene

enr.finlayi <- gopher(a.finlayi, url='fungi2011', task=list('cog', 'ko', 'ko_pathway','kog'))
enr.borealis <- gopher(a.borealis, url='fungi2011', task=list('cog', 'ko', 'ko_pathway','kog'))


lapply(c('cog', 'ko', 'ko_pathway','kog' ), function(x){
  if (!is.null(enr.finlayi[[x]]) ) {
    print(x)
    png(paste0("Archaeorhizomyces finlayi_",x,"_vs_AllGenes_treemap.png"),
        units="px", width=1920, height=1080, res=300)
    plotEnrichedTreemap(enr.finlayi, enrichment = x, namespace = "none")
    dev.off()
    
  }
})


lapply(c('cog', 'ko', 'ko_pathway','kog' ), function(x){
  if (!is.null(enr.borealis[[x]]) ) {
    print(x)
    png(paste0("Archaeorhizomyces borealis_",x,"_vs_AllGenes_treemap.png"),
        units="px", width=1920, height=1080, res=300)
    plotEnrichedTreemap(enr.borealis, enrichment = x, namespace = "none")
    dev.off()
    
  }
})
