
library(reshape2)
library(ggplot2)
source("~/Git/Rtoolbox/src/gopher.R")
projectFolder <- "/mnt/picea/projects/spruce/vhurry//spruce-flakaliden-root-2011"
deFolder <- file.path(projectFolder,"DE")
dataFolder <- file.path(projectFolder,"data")

combinedDataFolder <- file.path(projectFolder, "spruce-networks")

load(file.path(combinedDataFolder, "spruceEnr_2019-11-22.RData"))

infomapClusters <- read.table(file.path(projectFolder, "spruce-networks/combined/cluster/InfomapClusters.tsv"),
                              sep='\t', header=T)

# cluster 3, 5, 2, 12 on combined
# 
spruceEnr$combined$Cluster3$go
spruceEnr$combined$Cluster5$go
spruceEnr$combined$Cluster2$go
spruceEnr$combined$Cluster12$go

toi <- c("cell cycle", 
         "cell population proliferation",
         "signaling",
        "cytoskeleton",
        "regulation of cellular response to oxidative stress",
        "response to oxidative stress",
        "defense response signaling pathway, resistance gene-dependent",
        "LRR domain binding",
         "cell wall biogenesis",
         "plant-type cell wall biogenesis",
         "plant-type primary cell wall biogenesis",
         "plant-type secondary cell wall biogenesis",
         "plant-type cell wall organization or biogenesis")

group1 <- c("Cytoskeleton.microtubular network.Kinesin microtubule-based motor protein activities.Kinesin-14 motor protein",
"Cell wall.pectin.modification and degradation.pectate lyase",
"Protein modification.phosphorylation.TKL kinase superfamily.LRR-III kinase",
"Cell cycle.regulation.cyclin-dependent kinase complex.CDKB catalytic component",
"Cell cycle.regulation.cyclins.CYCA-type cyclin",
"Cell cycle.mitosis and meiosis.chromatin condensation.condensin I complex.CAP-H component",
"Cell cycle.regulation.cyclin-dependent kinase complex.CKS scaffolding component",
"Cytoskeleton.microtubular network.Kinesin microtubule-based motor protein activities.Kinesin-7 motor protein",
"Cell cycle.mitosis and meiosis.chromatin condensation.condensin I complex.CAP-G component",
"Cell wall.pectin.modification and degradation.pectin acetylesterase",
"Cell wall.hemicellulose.heteromannan.modification and degradation.endo-beta-1,4-mannanase",
"Cell cycle.cytokinesis.phragmoplast microtubule organization.EB1 microtubule plus-end-tracking protein",
"Cell wall.pectin.homogalacturonan.modification and degradation.pectin methylesterase",
"Cell cycle.interphase.DNA replication.elongation.DNA polymerase alpha complex.POLA3 primase component",
"Cell cycle.interphase.DNA replication.preinitiation.TOP2 DNA topoisomerase",
"Cell wall.hemicellulose.xyloglucan.modification and degradation.1,2-alpha-fucosidase",
"Cell cycle.mitosis and meiosis.chromatin condensation.condensin I complex.CAP-D2A component",
"Cell wall.pectin.modification and degradation.polygalacturonase activities.PGX1 polygalacturonase",
"Cytoskeleton.microtubular network.Kinesin microtubule-based motor protein activities.Kinesin-12 motor protein",
"Cell cycle.mitosis and meiosis.chromatin condensation.condensin II complex.CAP-D2B component",
"Cell cycle.cytokinesis.cell-plate formation.alpha-Aurora kinase",
"Cell cycle.mitosis and meiosis.chromatin condensation.condensin II complex.CAP-H2 component",
"Cytoskeleton.microtubular network.Kinesin microtubule-based motor protein activities.Kinesin-5 motor protein",
"Cell cycle.interphase.DNA replication.preinitiation.MCM replicative DNA helicase complex.MCM3 component",
"Cell cycle.mitosis and meiosis.meiotic recombination.meiotic exit.MS5/TDM1 meiotic exit regulator",
"Cell cycle.interphase.DNA replication.preinitiation.CDT1 helicase auxiliary factor",
"Cytoskeleton.microtubular network.Kinesin microtubule-based motor protein activities.Kinesin-10 motor protein",
"Cell cycle.interphase.DNA replication.elongation.DNA polymerase epsilon complex.POL2/POLE1 catalytic component",
"Cell cycle.mitosis and meiosis.sister chromatid separation.cohesin dissociation.ESP1 separase",
"Cytoskeleton.microtubular network.Kinesin microtubule-based motor protein activities.Kinesin-8 motor protein",
"Cell cycle.regulation.cyclins.CYCB-type cyclin",
"Cell cycle.cytokinesis.preprophase microtubule organization.Kinesin-14 microtubule-based motor protein",
"Cell cycle.cytokinesis.cell-plate formation.SNARE cell-plate vesicle fusion complex.SYP71 Qc-SNARE component",
"Cell cycle.mitosis and meiosis.chromatin condensation.condensin II complex.CAP-G2 component",
"Cell cycle.mitosis and meiosis.chromosome segregation.Kinesin-13 microtubule destabilizing motor protein",
"Cell wall.cell wall proteins.expansins.alpha-type expansin",
"Cytoskeleton.microfilament network.actin filament protein",
"Cell cycle.regulation.cyclins.CYCD-type cyclin",
"Cytoskeleton.microtubular network.Kinesin microtubule-based motor protein activities.Kinesin-13 motor protein",
"Cell wall.pectin.rhamnogalacturonan I.modification and degradation.beta-galactosidase",
"Cell wall.cell wall proteins.hydroxyproline-rich glycoproteins.extensins (EXTs).glycoproteins.LRR-domain extensin",
"Cell cycle.mitosis and meiosis.meiotic recombination.DNA strand exchange.RPA presynaptic filament assembly factor complex.RPA1 component",
"Cytoskeleton.cp-actin-dependent plastid movement.PMIR cp-actin stability co-factor",
"Cytoskeleton.microtubular network.microtubule dynamics.WDL microtubule-stabilizing factor",
"Cell wall.pectin.rhamnogalacturonan I.modification and degradation.alpha-L-arabinofuranosidase activities.bifunctional BXL-type alpha-L-arabinofuranosidase and beta-D-xylosidase",
"Cell wall.cell wall proteins.hydroxyproline-rich glycoproteins.extensins (EXTs).glycoproteins.LRR-domain extensin",
"Cell wall.hemicellulose.heteromannan.synthesis.mannan synthase activities.CSLD-type mannan synthase",
"Cell wall.hemicellulose.xylan.synthesis.galacturonosyltransferase",
"Cell wall.cellulose.synthesis.cellulose synthase complex (CSC).CSC components.CesA-type catalytic component",
"Cell wall.hemicellulose.heteromannan.modification and degradation.endo-beta-1,4-mannanase",
"Cell cycle.regulation.cyclin-dependent kinase inhibitor activities.SIM-type inhibitor",
"Cell wall.lignin.monolignol conjugation and polymerization.lignin peroxidase",
"Cell wall.hemicellulose.xylan.synthesis.glucuronosyltransferase activities.GUX-type glucuronosyltransferase",
"Cell wall.cell wall proteins.hydroxyproline-rich glycoproteins.arabinogalactan proteins (AGPs).glycoproteins.fasciclin-type arabinogalactan protein",
"Cytoskeleton.microfilament network.actin polymerisation.Arp2/3 actin polymerization initiation complex.ArpC2 component",
"Cell wall.cell wall proteins.hydroxyproline-rich glycoproteins.arabinogalactan proteins (AGPs).glycoproteins.xylogen-type arabinogalactan protein",
"Cell wall.hemicellulose.xylan.synthesis.galacturonosyltransferase",
"Cell wall.pectin.rhamnogalacturonan I.modification and degradation.beta-galactosidase",
"Cell wall.pectin.homogalacturonan.modification and degradation.pectin methylesterase",
"Cell wall.cell wall proteins.expansins.alpha-type expansin")

group2 <- c(
"Chromatin organisation.histones.H3-type histone",
"Cell cycle.mitosis and meiosis.chromosome segregation.centromere assembly and maintenance.WYRD inner centromere protein",
"Chromatin organisation.histone modifications.histone lysine methylation/demethylation.class I/Ez histone methyltransferase component",
"Chromatin organisation.histone modifications.histone lysine methylation/demethylation.PRC2 histone methylation complex.VRN/FIS/EMF core complexes.EZ-like catalytic component (CLF,SWN,MEA)",
"Chromatin organisation.histone modifications.histone phosphorylation.Aurora kinase",
"Chromatin organisation.histones.H1 linker histone",
"Chromatin organisation.histone chaperone activities.CAF1 histone chaperone complex.CAF1c/MSI component",
"Cell wall.lignin.monolignol conjugation and polymerization.lignin laccase",
"Cell wall.hemicellulose.heteromannan.synthesis.mannan synthesis accessory protein")

group3 <- c("RNA biosynthesis.transcriptional activation.AP2/ERF superfamily.DREB-type transcription factor",
"RNA biosynthesis.transcriptional activation.AP2/ERF superfamily.ERF-type transcription factor",
"RNA biosynthesis.transcriptional activation.TIFY transcription factor",
"RNA biosynthesis.transcriptional activation.WRKY transcription factor",
"RNA biosynthesis.transcriptional activation.AS2/LOB transcription factor",
"RNA biosynthesis.transcriptional activation.C2H2 zinc finger transcription factor",
"RNA biosynthesis.transcriptional activation.B3 superfamily.RAV/NGATHA transcription factor",
"RNA biosynthesis.transcriptional activation.MYB superfamily.MYB transcription factor",
"RNA biosynthesis.transcriptional activation.NAC transcription factor",
"RNA biosynthesis.transcriptional activation.GRF-GIF transcriptional complex.GRF transcription factor component")


group4 <- c("External stimuli response.biotic stress.pathogen effector.NLR effector receptor",
"Protein modification.phosphorylation.TKL kinase superfamily.L-lectin kinase",
"Protein modification.phosphorylation.STE kinase superfamily.MAP3K-MEKK kinase",
"Protein modification.phosphorylation.TKL kinase superfamily.LRK10-1-like kinase",
"Protein modification.phosphorylation.TKL kinase superfamily.RLCK-Os kinase",
"Protein modification.phosphorylation.TKL kinase superfamily.G-Lectin kinase families.SD-2 kinase",
"Protein modification.phosphorylation.TKL kinase superfamily.Crinkly-like kinase",
"Protein modification.phosphorylation.TKL kinase superfamily.LRR-XV kinase",
"Protein modification.phosphorylation.TKL kinase superfamily.LRR-X kinase families.LRR-Xa kinase",
"Phytohormones.jasmonic acid.synthesis.oxophytodienoate r")

group_protein_modification <- c(
  "Protein modification.phosphorylation.TKL kinase superfamily.LRR-III kinase",
  "Protein modification.phosphorylation.TKL kinase superfamily.L-lectin kinase",
  "Protein modification.phosphorylation.STE kinase superfamily.MAP3K-MEKK kinase",
  "Protein modification.phosphorylation.TKL kinase superfamily.LRK10-1-like kinase",
  "Protein modification.phosphorylation.TKL kinase superfamily.RLCK-Os kinase",
  "Protein modification.phosphorylation.TKL kinase superfamily.G-Lectin kinase families.SD-2 kinase",
  "Protein modification.phosphorylation.TKL kinase superfamily.Crinkly-like kinase",
  "Protein modification.phosphorylation.TKL kinase superfamily.LRR-XV kinase",
  "Protein modification.phosphorylation.TKL kinase superfamily.LRR-X kinase families.LRR-Xa kinase"
)

group_rna <- c(
  "RNA biosynthesis.transcriptional activation.AP2/ERF superfamily.DREB-type transcription factor",
  "RNA biosynthesis.transcriptional activation.AP2/ERF superfamily.ERF-type transcription factor",
  "RNA biosynthesis.transcriptional activation.TIFY transcription factor",
  "RNA biosynthesis.transcriptional activation.WRKY transcription factor",
  "RNA biosynthesis.transcriptional activation.AS2/LOB transcription factor",
  "RNA biosynthesis.transcriptional activation.C2H2 zinc finger transcription factor",
  "RNA biosynthesis.transcriptional activation.B3 superfamily.RAV/NGATHA transcription factor",
  "RNA biosynthesis.transcriptional activation.MYB superfamily.MYB transcription factor",
  "RNA biosynthesis.transcriptional activation.NAC transcription factor",
  "RNA biosynthesis.transcriptional activation.GRF-GIF transcriptional complex.GRF transcription factor component"
  )
terms <- as.factor(c(group1, group2, group3, group4))



newList <- spruceEnr$combined[c("Cluster3", "Cluster5", "Cluster2", "Cluster12")]

newList2 <- lapply(names(newList), function(x) {
  res <- newList[[x]]$mapman[newList[[x]]$mapman$name %in% terms,][c('name', 'nt', 'padj')]
  res$group <- res$name
  res$group <- ifelse(res$name %in% group1, "group1", 
                      ifelse(res$name %in% group2, "group2", 
                             ifelse(res$name %in% group3, "group3", "group4")))
  res
})

names(newList2) <- names(newList)

newMelt <- reshape2::melt(newList2,id=c("name","nt","padj"))
newMelt <- newMelt[newMelt$padj <0.01,]


newMelt$padj <- abs(log2(newMelt$padj))

#newMelt$name <- factor(newMelt$name, levels = levels(terms))
newMelt <- newMelt[order(newMelt$value),]

# Cell cycle - cell cycle
# Cell proliferation - cell population proliferation
# Cell wall biosynthesis - cell wall biogenesis?? plant-type cell wall biogenesis???
# plant-type primary cell wall biogenesis ?? plant-type secondary cell wall biogenesis??
# plant-type cell wall organization or biogenesis?? This one is more generic
# Cytoskeleton - cytoskeleton
# Cell wall remodeling is this part of biogenesis too? the definition says A process that results in the biosynthesis of constituent macromolecules, assembly, arrangement of constituent parts, or disassembly of a cellulose- and pectin-containing cell wall."
# Response to oxidative stress - regulation of cellular response to oxidative stress?
# Effector triggered immunity - defense response signaling pathway, resistance gene-dependent - synonym: "effector triggered immunity" 
# LRR - LRR domain binding??
# Signalling - signaling

library(ggplot2)

ggplot2::ggplot(newMelt, aes(x = L1, y = name)) + 
  #
  geom_point(aes(fill=padj, size = nt), alpha = 1, shape = 21) +
  #scale_fill_manual(values = colours, guide = FALSE) +
  #scale_size_continuous(limits = c(1, 10000), range = c(0.1,17)) + 
  scale_size(range = c(0.5, 10), breaks = c(1,10,100,500,1000,2000,5000,10000)) +  # Adjust the range of points size

  facet_grid(rows = vars(value), scales = "free_y", space = "free_y") + 
  
  theme(#legend.key=element_blank(), 
    #axis.text.x = element_text(colour = "black", size =8, face = "bold", angle = 90, vjust = 0.3, hjust = 1), 
    #axis.text.x=element_blank(),
    axis.text.y = element_text(colour = "black", face = "bold", size = 4), 
    #legend.text = element_text(size = 8, face ="bold", colour ="black"), 
    legend.title = element_text(size = 8, face = "bold"), 
    panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
    legend.position = "none") 
  















"Chromatin organisation.histones.H3-type histone",
"Chromatin organisation.histone modifications.histone lysine methylation/demethylation.class I/Ez histone methyltransferase component",
"Chromatin organisation.histone modifications.histone lysine methylation/demethylation.PRC2 histone methylation complex.VRN/FIS/EMF core complexes.EZ-like catalytic component (CLF,SWN,MEA)",
"Chromatin organisation.histone modifications.histone phosphorylation.Aurora kinase",
"Chromatin organisation.histones.H1 linker histone",
"Chromatin organisation.histone chaperone activities.CAF1 histone chaperone complex.CAF1c/MSI component",


# CELL CYCLE
group_cell_cycle <- c(
  "Cell cycle.regulation.cyclin-dependent kinase complex.CDKB catalytic component",
  "Cell cycle.regulation.cyclins.CYCA-type cyclin",
  "Cell cycle.mitosis and meiosis.chromatin condensation.condensin I complex.CAP-H component",
  "Cell cycle.regulation.cyclin-dependent kinase complex.CKS scaffolding component",
  "Cell cycle.mitosis and meiosis.chromatin condensation.condensin I complex.CAP-G component",
  "Cell cycle.interphase.DNA replication.elongation.DNA polymerase alpha complex.POLA3 primase component",
  "Cell cycle.interphase.DNA replication.preinitiation.TOP2 DNA topoisomerase",
  "Cell cycle.mitosis and meiosis.chromatin condensation.condensin II complex.CAP-D2B component",
  "Cell cycle.cytokinesis.cell-plate formation.alpha-Aurora kinase",
  "Cell cycle.mitosis and meiosis.chromatin condensation.condensin II complex.CAP-H2 component",
  "Cell cycle.cytokinesis.phragmoplast microtubule organization.EB1 microtubule plus-end-tracking protein",
  "Cell cycle.mitosis and meiosis.chromatin condensation.condensin I complex.CAP-D2A component",
  "Cell cycle.interphase.DNA replication.preinitiation.MCM replicative DNA helicase complex.MCM3 component",
  "Cell cycle.mitosis and meiosis.meiotic recombination.meiotic exit.MS5/TDM1 meiotic exit regulator",
  "Cell cycle.interphase.DNA replication.preinitiation.CDT1 helicase auxiliary factor",
  "Cell cycle.interphase.DNA replication.elongation.DNA polymerase epsilon complex.POL2/POLE1 catalytic component",
  "Cell cycle.mitosis and meiosis.sister chromatid separation.cohesin dissociation.ESP1 separase",
  "Cell cycle.regulation.cyclins.CYCB-type cyclin",
  "Cell cycle.cytokinesis.preprophase microtubule organization.Kinesin-14 microtubule-based motor protein",
  "Cell cycle.cytokinesis.cell-plate formation.SNARE cell-plate vesicle fusion complex.SYP71 Qc-SNARE component",
  "Cell cycle.mitosis and meiosis.chromatin condensation.condensin II complex.CAP-G2 component",
  "Cell cycle.mitosis and meiosis.chromosome segregation.Kinesin-13 microtubule destabilizing motor protein",
  "Cell cycle.mitosis and meiosis.chromosome segregation.centromere assembly and maintenance.WYRD inner centromere protein",
  "Cell cycle.regulation.cyclins.CYCD-type cyclin",
  "Cell cycle.mitosis and meiosis.meiotic recombination.DNA strand exchange.RPA presynaptic filament assembly factor complex.RPA1 component",
  "Cell cycle.regulation.cyclin-dependent kinase inhibitor activities.SIM-type inhibitor"
)

# CELL WALL
group_cell_wall <- c(
  "Cell wall.pectin.modification and degradation.pectate lyase",
  "Cell wall.pectin.modification and degradation.pectin acetylesterase",
  "Cell wall.hemicellulose.heteromannan.modification and degradation.endo-beta-1,4-mannanase",
  "Cell wall.pectin.homogalacturonan.modification and degradation.pectin methylesterase",
  "Cell wall.hemicellulose.xyloglucan.modification and degradation.1,2-alpha-fucosidase",
  "Cell wall.pectin.modification and degradation.polygalacturonase activities.PGX1 polygalacturonase",
  "Cell wall.cell wall proteins.expansins.alpha-type expansin",
  "Cell wall.pectin.rhamnogalacturonan I.modification and degradation.beta-galactosidase",
  "Cell wall.cell wall proteins.hydroxyproline-rich glycoproteins.extensins (EXTs).glycoproteins.LRR-domain extensin",
  "Cell wall.pectin.rhamnogalacturonan I.modification and degradation.alpha-L-arabinofuranosidase activities.bifunctional BXL-type alpha-L-arabinofuranosidase and beta-D-xylosidase",
  "Cell wall.cell wall proteins.hydroxyproline-rich glycoproteins.extensins (EXTs).glycoproteins.LRR-domain extensin",
  "Cell wall.hemicellulose.heteromannan.synthesis.mannan synthase activities.CSLD-type mannan synthase",
  "Cell wall.hemicellulose.xylan.synthesis.galacturonosyltransferase",
  "Cell wall.cellulose.synthesis.cellulose synthase complex (CSC).CSC components.CesA-type catalytic component",
  "Cell wall.hemicellulose.heteromannan.modification and degradation.endo-beta-1,4-mannanase",
  "Cell wall.lignin.monolignol conjugation and polymerization.lignin peroxidase",
  "Cell wall.hemicellulose.xylan.synthesis.glucuronosyltransferase activities.GUX-type glucuronosyltransferase",
  "Cell wall.cell wall proteins.hydroxyproline-rich glycoproteins.arabinogalactan proteins (AGPs).glycoproteins.fasciclin-type arabinogalactan protein",
  "Cell wall.cell wall proteins.hydroxyproline-rich glycoproteins.arabinogalactan proteins (AGPs).glycoproteins.xylogen-type arabinogalactan protein",
  "Cell wall.hemicellulose.xylan.synthesis.galacturonosyltransferase",
  "Cell wall.pectin.rhamnogalacturonan I.modification and degradation.beta-galactosidase",
  "Cell wall.pectin.homogalacturonan.modification and degradation.pectin methylesterase",
  "Cell wall.cell wall proteins.expansins.alpha-type expansin",
  "Cell wall.lignin.monolignol conjugation and polymerization.lignin laccase",
  "Cell wall.hemicellulose.heteromannan.synthesis.mannan synthesis accessory protein"
  
)

# CYTOSKELETON
group_cytoskeleton <- c(
  "Cytoskeleton.microtubular network.Kinesin microtubule-based motor protein activities.Kinesin-14 motor protein",
  "Cytoskeleton.microtubular network.Kinesin microtubule-based motor protein activities.Kinesin-7 motor protein",
  "Cytoskeleton.microtubular network.Kinesin microtubule-based motor protein activities.Kinesin-12 motor protein",
  "Cytoskeleton.microtubular network.Kinesin microtubule-based motor protein activities.Kinesin-10 motor protein",
  "Cytoskeleton.microtubular network.Kinesin microtubule-based motor protein activities.Kinesin-8 motor protein",
  "Cytoskeleton.microfilament network.actin filament protein",
  "Cytoskeleton.microtubular network.Kinesin microtubule-based motor protein activities.Kinesin-5 motor protein",
  "Cytoskeleton.microtubular network.Kinesin microtubule-based motor protein activities.Kinesin-13 motor protein",
  "Cytoskeleton.cp-actin-dependent plastid movement.PMIR cp-actin stability co-factor",
  "Cytoskeleton.microtubular network.microtubule dynamics.WDL microtubule-stabilizing factor",
  "Cytoskeleton.microfilament network.actin polymerisation.Arp2/3 actin polymerization initiation complex.ArpC2 component"
)

# RESPONSE TO OXIDATIVE STRESS
# PROTEIN MODIFICATION.PHOSPHORYLATION.TKL KINASE SUPERFAMILY
group_protein_modification <- c(
  "Protein modification.phosphorylation.TKL kinase superfamily.LRR-III kinase",
  "Protein modification.phosphorylation.TKL kinase superfamily.L-lectin kinase",
  "Protein modification.phosphorylation.STE kinase superfamily.MAP3K-MEKK kinase",
  "Protein modification.phosphorylation.TKL kinase superfamily.LRK10-1-like kinase",
  "Protein modification.phosphorylation.TKL kinase superfamily.RLCK-Os kinase",
  "Protein modification.phosphorylation.TKL kinase superfamily.G-Lectin kinase families.SD-2 kinase",
  "Protein modification.phosphorylation.TKL kinase superfamily.Crinkly-like kinase",
  "Protein modification.phosphorylation.TKL kinase superfamily.LRR-XV kinase",
  "Protein modification.phosphorylation.TKL kinase superfamily.LRR-X kinase families.LRR-Xa kinase"
)

# EXTERNAL STIMULI RESPONSE.BIOTIC STRESS. PATHOGEN EFFECTOR
group_external_stimuli <- c("External stimuli response.biotic stress.pathogen effector.NLR effector receptor")

# RNA BIOSYNTHESIS.TRANSCRIPTIONAL ACTIVATION
group_rna <- c(
  "RNA biosynthesis.transcriptional activation.AP2/ERF superfamily.DREB-type transcription factor",
  "RNA biosynthesis.transcriptional activation.AP2/ERF superfamily.ERF-type transcription factor",
  "RNA biosynthesis.transcriptional activation.TIFY transcription factor",
  "RNA biosynthesis.transcriptional activation.WRKY transcription factor",
  "RNA biosynthesis.transcriptional activation.AS2/LOB transcription factor",
  "RNA biosynthesis.transcriptional activation.C2H2 zinc finger transcription factor",
  "RNA biosynthesis.transcriptional activation.B3 superfamily.RAV/NGATHA transcription factor",
  "RNA biosynthesis.transcriptional activation.MYB superfamily.MYB transcription factor",
  "RNA biosynthesis.transcriptional activation.NAC transcription factor",
  "RNA biosynthesis.transcriptional activation.GRF-GIF transcriptional complex.GRF transcription factor component"
)

# GO extraction

clusterList <- spruceEnr$combined[c("Cluster3", "Cluster5", "Cluster2", "Cluster12")]
goTerms <- c("cell cycle",
             "cell wall organization or biogenesis",
             "structural constituent of cytoskeleton",
             "response to oxidative stress", 
             "phosphorylation",
             "cellular response to stress"
             
             )

newListGO <- lapply(names(clusterList), function(x) {
  res <- clusterList[[x]]$go[clusterList[[x]]$go$name %in% goTerms,][c('name', 'nt', 'padj')]
  res
})
names(newListGO) <- names(clusterList)

newMelt <- reshape2::melt(newListGO, id=c("name","nt","padj"))
newMelt$padj <- abs(log2(newMelt$padj))
newMelt$enr <- rep("go", nrow(newMelt))

ggplot2::ggplot(newMelt, aes(x = L1, y = name)) + 
  geom_point(aes(fill=nt, size = nt), alpha = 1, shape = 21) +
  #scale_fill_manual(values = colours, guide = FALSE) +
  #scale_size_continuous(limits = c(1, 10000), range = c(0.1,17)) + 
  scale_size(range = c(1, 20), breaks = c(1,10,100,500)) +  # Adjust the range of points size
  #facet_grid(rows = vars(value), scales = "free_y", space = "free_y") + 
  theme(#legend.key=element_blank(), 
    #axis.text.x = element_text(colour = "black", size =8, face = "bold", angle = 90, vjust = 0.3, hjust = 1), 
    #axis.text.x=element_blank(),
    axis.text.y = element_text(colour = "black", face = "bold", size = 14), 
    #legend.text = element_text(size = 8, face ="bold", colour ="black"), 
    legend.title = element_text(size = 8, face = "bold"), 
    panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
    legend.position = "right") 

# clustersOI <- c("Cluster3", "Cluster5", "Cluster2", "Cluster12")
# typeTerm <- "mapman"
# findTerm <- "RNA biosynthesis.transcriptional"

mergeTerms <- function(clustersOI, typeTerm, findTerm) {
  
  localList <- list()
  
  findTerm <- paste(findTerm, collapse='|')
  
  for (cluster in clustersOI) {
    print(cluster)
    
    # gene names in that cluster
    genesInCluster <- infomapClusters[infomapClusters$cluster==cluster,]$gene
    
    # indices with terms of relevance
    indices <- grep(findTerm, spruceEnr$combined[[cluster]][[typeTerm]]$name)
    
    # get those terms and print them to check
    localTerms <- spruceEnr$combined[[cluster]][[typeTerm]]$name[indices]
    print(localTerms)
    
    #get the ids
    localIDs <- spruceEnr$combined[[cluster]][[typeTerm]]$id[indices]
    
    if (length(localTerms) > 0) {
      # reverse gopher
      gopherTerms <- gopher(localIDs, task=typeTerm, url="pabies", endpoint = "term-to-gene")

      termGenes <- intersect(genesInCluster, unlist(gopherTerms))

      if(length(termGenes) > 0) {
        localList[[cluster]] <- termGenes
      }
    } 
  }
  
  return(localList)
}

# CELL CYCLE
cellcycle <- mergeTerms(clustersOI <- c("Cluster3", "Cluster5", "Cluster2", "Cluster12"),
                        typeTerm <- "mapman",
                        findTerm <- group_cell_cycle)

# CELL WALL
cellwall <- mergeTerms(clustersOI <- c("Cluster3", "Cluster5", "Cluster2", "Cluster12"),
                           typeTerm <- "mapman",
                           findTerm <- group_cell_wall)

# CYTOSKELETON
cytoskeleton <- mergeTerms(clustersOI <- c("Cluster3", "Cluster5", "Cluster2", "Cluster12"),
                           typeTerm <- "mapman",
                           findTerm <- group_cytoskeleton)

# RESPONSE TO OXIDATIVE STRESS
# PROTEIN MODIFICATION.PHOSPHORYLATION.TKL KINASE SUPERFAMILY
protein <- mergeTerms(clustersOI <- c("Cluster3", "Cluster5", "Cluster2", "Cluster12"),
                           typeTerm <- "mapman",
                           findTerm <- group_protein_modification)
# EXTERNAL STIMULI RESPONSE.BIOTIC STRESS. PATHOGEN EFFECTOR
external <- mergeTerms(clustersOI <- c("Cluster3", "Cluster5", "Cluster2", "Cluster12"),
                      typeTerm <- "mapman",
                      findTerm <- group_external_stimuli)
# RNA BIOSYNTHESIS.TRANSCRIPTIONAL ACTIVATION
rnabio <- mergeTerms(clustersOI <- c("Cluster3", "Cluster5", "Cluster2", "Cluster12"),
           typeTerm <- "mapman",
           findTerm <- group_rna)

# turn it into newMelt format
#  name nt       padj       L1 enr

newMelt2 <- newMelt

melting <- function(termDef, listOfTerms, dfToMelt) {
  for(cluster in names(listOfTerms)) {
    print(cluster)
    dfToMelt <- rbind(dfToMelt, c(termDef, 
                                  length(listOfTerms[[cluster]]), as.numeric(1), cluster, "mapman"))  
    print(dfToMelt)
  }
  return(dfToMelt)
}

newMelt2 <- melting("CELL CYCLE", cellcycle, newMelt2)
newMelt2 <- melting("CELL WALL", cellwall, newMelt2)
newMelt2 <- melting("CYTOSKELETON", cytoskeleton, newMelt2)
newMelt2 <- melting("PROTEIN MODIFICATION.PHOSPHORYLATION.TKL KINASE SUPERFAMILY", protein, newMelt2)
newMelt2 <- melting("EXTERNAL STIMULI RESPONSE.BIOTIC STRESS. PATHOGEN EFFECTOR", external, newMelt2)
newMelt2 <- melting("RNA BIOSYNTHESIS.TRANSCRIPTIONAL ACTIVATION", rnabio, newMelt2)


newMelt2$nt <- as.numeric(newMelt2$nt)


ggplot2::ggplot(newMelt2, aes(x = L1, y = name)) + 
  geom_point(aes(fill=nt, size = nt), alpha = 1, shape = 21) +
  #scale_fill_manual(values = colours, guide = FALSE) +
  #scale_size_continuous(limits = c(1, 10000), range = c(0.1,17)) + 
  scale_size(range = c(1, 20), breaks = c(1,10,100,500)) +  # Adjust the range of points size
  #facet_grid(rows = vars(value), scales = "free_y", space = "free_y") + 
  theme(#legend.key=element_blank(), 
    #axis.text.x = element_text(colour = "black", size =8, face = "bold", angle = 90, vjust = 0.3, hjust = 1), 
    #axis.text.x=element_blank(),
    axis.text.y = element_text(colour = "black", face = "bold", size = 10), 
    #legend.text = element_text(size = 8, face ="bold", colour ="black"), 
    legend.title = element_text(size = 8, face = "bold"), 
    panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
    legend.position = "right") 


###################################################

RNAlist <- unique(unlist(RNAlist))

print(RNAlist)

RNAtermlist <- sapply(clustersOI, function(cluster) {
  indices <- grep(findTerm, spruceEnr$combined[[cluster]][[typeTerm]]$name)
  spruceEnr$combined[[cluster]][[typeTerm]]$id[indices]
}, simplify = "array", USE.NAMES = F)
RNAtermlist <- unique(unlist(RNAtermlist))

gopher(RNAtermlist, task="mapman", url="pabies", endpoint = "term-to-gene")
