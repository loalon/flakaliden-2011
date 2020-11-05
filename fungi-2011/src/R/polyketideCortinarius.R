
library("rBiobox")
source("~/Git/UPSCb-common/src/R/gopher.R")



spruceDir <- "/mnt/picea/projects/spruce/vhurry/14_SpruceRoots_Project"
fungiDir <- "/mnt/picea/projects/spruce/vhurry/fungi-flakaliden-2011"


fungiData <- file.path(fungiDir,"data")
fungiDE <- file.path(fungiDir,"DE")

load(file.path(fungiData, "megaData.RData"))

anno <- read.table(file.path(fungiData, "annotation_results.emapper.annotations"), header = F, sep='\t', quote = "", comment.char = "")

# GO:0016218 obsolete polyketide synthase activity
# GO:0034081 polyketide synthase complex
# GO:0034082 type II polyketide synthase complex
# GO:0034083 type III polyketide synthase complex
# GO:0030639 polyketide biosynthetic process

miniData <- megaData[megaData$species == "Cenococcum geophilum",]
mergedData <- dplyr::left_join(miniData, anno, by=c("gene"="V19") )

# COG3315 O-Methyltransferase involved in polyketide biosynthesis
# COG3319	Thioesterase domain of type I polyketide synthase or non-ribosomal peptide synthetase	Q
# COG3321	Acyl transferase domain in polyketide synthase (PKS) enzymes
# COG3559	Putative exporter of polyketide antibiotics	U
# 
polTerms <- c("COG3315", "COG3319", "COG3321", "COG3559")

res <- gopher("COG3321", task=("cog"), url="fungi2011", endpoint="term-to-gene")

load("/mnt/picea/projects/spruce/vhurry/fungi-flakaliden-2011/DE/transcriptsCombinedVsdData.RData")
load("/mnt/picea/projects/spruce/vhurry/fungi-flakaliden-2011/DE/ddsTreatment_transcripts.RData")

vsd.def <- varianceStabilizingTransformation(ddsFungiTreatment, blind = F)
vsd.mat <- t(assay(vsd.def))

plotMeanProfile(vsd.mat, 
                res$COG3321,
                substr(rownames(vsd.mat), 4, 4), 
                as.numeric(substr(rownames(vsd.mat), 2, 3)))

plotMeanProfile(vsd.mat, 
                res$COG3315,
                substr(rownames(vsd.mat), 4, 4), 
                as.numeric(substr(rownames(vsd.mat), 2, 3)))

full <- read.table("/mnt/picea/projects/spruce/vhurry/fungi-flakaliden-2011/data/kingdom.Fungi.2011.raw.tsv", header=T, sep='\t')
