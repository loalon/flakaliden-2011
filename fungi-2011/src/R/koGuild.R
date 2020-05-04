library(dplyr)

projectFolder <- "/mnt/picea/projects/spruce/vhurry/fungi-flakaliden-2011"
deFolder <- file.path(projectFolder,"DE")
dataFolder <- file.path(projectFolder,"data")
guildFile <- file.path(dataFolder, "fungi.guilds.tsv")
parsedFile <- file.path(dataFolder, "kos.parsed.tab")
taxoFile <- file.path(dataFolder, "gene_taxonomy.tsv")

guilds <- read.table(guildFile, header = TRUE, sep = '\t', comment.char = "", 
                     stringsAsFactors = F)
colnames(guilds)[1] <- "orf"

parsed <- read.table(parsedFile, header = TRUE, sep = '\t', comment.char = "", 
                     stringsAsFactors = F, quote = "")



koGuild <- left_join(parsed, guilds, by="orf")

koGuild_confidence <-filter(koGuild, Confidence.Ranking %in% c("Highly Probable","Probable") )

koGuild_res <- data.frame(Trophic.Mode = koGuild_confidence$Trophic.Mode,
                          Guild = koGuild_confidence$Guild,
                          KO = koGuild_confidence$ko)
write.table(koGuild_res, file="koGuild.tsv", sep='\t', row.names = F, quote = F)

# left join taxonomy
taxo <- read.table(taxoFile, header = TRUE, sep = '\t', comment.char = "", 
                   stringsAsFactors = F, quote = "")
colnames(taxo)[1] <- "orf"

koTaxo <- left_join(parsed, taxo, by="orf")

koTaxo_res <- data.frame(ko = koTaxo$ko,
                         superkingdom=koTaxo$superkingdom,
                         kingdom=koTaxo$kingdom,
                         phylum=koTaxo$phylum,
                         class=koTaxo$class,
                         order=koTaxo$order,
                         family=koTaxo$family,
                         genus=koTaxo$genus,
                         species=koTaxo$species
                         )
write.table(koTaxo_res, file="koTaxo.tsv", sep='\t', row.names = F, quote = F)

