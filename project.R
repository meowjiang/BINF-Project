source("https://bioconductor.org/biocLite.R"); ## try http:// if https:// URLs are not supported
biocLite("affy");
biocLite("limma");
biocLite("GEOquery");
biocLite("proxy")
library(affy);
library(limma);
library(GEOquery); 
library(proxy)

main <-function(diseaseIds){
  platforms <- basics.getPlatforms()
  topTables <- lapply(diseaseIds,disease.topGenes,platforms=platforms)
  topGenes <- lapply(topTables,disease.geneList)
  names(topGenes)<-lapply(topGenes,helper.getName)
  distanceMatrix <- outer(topGenes, topGenes,FUN = basics.tanimoto.vectorized)
  return (distanceMatrix)
}

goldStandard <-function(){
  #same diseases different datasets
  platforms <- basics.getPlatforms()
  diseaseIds <- c("GDS1321","GDS3472")
  topTables <- lapply(diseaseIds,disease.topGenes,platforms=platforms)
  topGenes <- lapply(topTables,disease.geneList)
  names(topGenes)<-lapply(topGenes,helper.getName)
  distanceMatrix <- outer(topGenes, topGenes,FUN = basics.tanimoto.vectorized)
  return (distanceMatrix)
}


basics.getPlatforms <-function(){
  platformNames <- c("GPL570","GPL96")
  platforms<-lapply(platformNames,helper.getPlatform)
  names(platforms)<-platformNames
  return (platforms)
}

basics.distance.vectorized <-Vectorize(basics.distance,c("topGenes1","topGenes2"))

basics.tanimoto.vectorized <-Vectorize(basics.tanimoto,c("topGenes1","topGenes2"))

basics.tanimoto <- function(topGenes1,topGenes2){
  t1 <- topGenes1
  t2 <- topGenes2
  
  t1Genes <-t1$genes
  t2Genes <-t2$genes
  
  intersect <- intersect(t1Genes,t2Genes)
  union <- unique(union(t1Genes,t2Genes))
  return (1-(length(intersect)/length(union)))
}


basics.distance <-function(topGenes1,topGenes2,platforms){
  t1 <-topGenes1
  t2 <-topGenes2
  
  t1Genes <-t1$genes
  t2Genes <-t2$genes
  
  #why do we have duplicate genes?
  #is it because the list is for all probles?
  t1PlatformGenes <-unique(Table(platforms[[t1$platform]])[,2])
  t2PlatformGenes <-unique(Table(platforms[[t2$platform]])[,2])
  
  common <- length(intersect(t1Genes,t2Genes))
  platformCommon <- length(intersect(t1PlatformGenes,t2PlatformGenes))
  bestIntersect <- min(length(t1Genes),length(t2Genes))
  bestmatch <- bestIntersect/platformCommon
  return(bestmatch-(common/platformCommon))

}

helper.getPlatform <- function(platformName){
  platform<- getGEO(platformName)
  return(platform)
}

helper.getName <- function(topGene){
  return(topGene$id)
}

disease.geneList <-function(toptable){
  return (list(genes=rownames(toptable$table),id=toptable$id,platform=toptable$platform))
}

disease.topGenes <- function(datasetID,platforms){
  geoSet <-getGEO(datasetID) 
  platformName <- Meta(geoSet)$platform  
  platform <- platforms[[platformName]]  
  geoData <- GDS2eSet(geoSet,do.log2 = TRUE,GPL=platform)
  geoMatrix <- as.matrix(geoData)
  mappings<- Table(platform)[,1:2]
  mappedMatrix <-   merge(geoMatrix,mappings,by.x="row.names",by.y="ID", all.x = TRUE)[,c(-1)] 
  aggregatedMatrix <-   aggregate(mappedMatrix[,-ncol(mappedMatrix)],by=list(GB_ACC=mappedMatrix$GB_ACC),FUN=mean)
  rownames(aggregatedMatrix)<-aggregatedMatrix[,1]
  aggregatedMatrix[,1]<-NULL
  design<-model.matrix(~ Columns(dataTable(geoSet))$disease.state)
  fittedMatrix <-   lmFit(aggregatedMatrix,design)
  bayesOut <-   eBayes(fittedMatrix) 
  topGenes <-   topTable(bayesOut,p.value=.01,number = nrow(aggregatedMatrix))
  print (datasetID)
  print ( head( topGenes))
  return(list(table=topGenes,id=datasetID,platform=platformName))
}

