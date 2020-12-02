# Flakaliden 2011 analysis 
Source code for the *Metatranscriptomics captures the destabilised mycorrhizal economy in nutrient enriched boreal forests* 

## doc
Folder containing metadata and the ENA xml submission files 

## Spruce 2011  
Folder `spruce-roots-2011/pipeline`
* pre-processing scripts

Folder `spruce-roots-2011/src/R` 
* differentialExpression.Rmd  
  * Biological QA
  * Differential expression analysis
  * Data preparation for Seidr
	
* networkClustering.Rmd  
  * Networks clustering
	* Cluster enrichment
	
* TAtool_data
  * Scripts to format the data for the visualization app

## Spruce 2012  
Folder `spruce-roots-2012/src/R` 
* differentialExpression.Rmd  
  * Biological QA
  * Differential expression analysis

* differentialExpression_5year_included.Rmd  
  * Biological QA
  * Differential expression analysis
  * Includes 5 year fertilisation for future analysis

## Fungi 2011
* TAtool_data
  * Scripts to format the data for the visualization app
	
* differentialExpression_ko.Rmd  
  * Biological QA
  * Differential expression analysis
	* Data preparation for Seidr

* differentialExpression_transcripts.Rmd  
  * Biological QA
  * Differential expression analysis
	
* differentialExpression_fungi_species.Rmd  
  * Biological QA
  * Differential expression analysis
  * Data preparation for Seidr
	
* megaData.R
  * Creates a data frame with all the fungi information

* networkThreshold_fungi_ko.Rmd
	* KO hard threshold analysis
	
* networkClustering_fungi.Rmd
  * Networks clustering
	* Cluster enrichment
	
* prepareSpeciesNetwork.Rmd
  * Data preparation for Seidr

## Fungi 2012 
Folder `fungi-2012/src/R` 
* differentialExpression_ko.Rmd  
  * Biological QA
  * Differential expression analysis
* megaData.R
  * Creates a data frame with all the fungi information

## Common
Folder `src/R`
* correlationClusterAnalysis_*.Rmd
	* Correlation analysis based on module similarity
	
* bulkCorrelationAnalysis_*.Rmd
  * Correlation analysis results grouped and analysed
	
* differentialExpression_*.Rmd	
  * differential expression comparative between 2011 and 2012 datasets
	
* splitEnr*_Rmd
  * Split network modules analysis and enrichment
