source("https://bioconductor.org/biocLite.R"); ## try http:// if https:// URLs are not supported
biocLite("affy");
biocLite("limma");
biocLite("GEOquery");
biocLite("proxy")
library(affy);
library(limma);
library(GEOquery); 
library(proxy)
library(dplyr)





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

basics.tanimoto.vectorized <-Vectorize(basics.tanimoto,c("topGenes1","topGenes2"))

basics.tanimoto <- function(topGenes1,topGenes2){
  t1 <- topGenes1
  t2 <- topGenes2
  
  t1Positive <- t1$positive
  t1Negative <- t1$negative
  
  t2Positive <- t2$positive
  t2Negative <- t2$Negative
  
  positiveIntersect <- intersect(t1Positive,t2Positive)
  negativeIntersect <- intersect(t1Negative,t2Negative)
  
  union <- (union(t1Positive,t1Negative))
  union<-(union(union,t2Positive))
  union<-(unition(union.t2Negative))
            
  intersect = union(positveIntersert,negativeIntersect)
  return (1-(length(intersect)/length(union)))
}


helper.getPlatform <- function(platformName){
  platform<- getGEO(platformName)
  return(platform)
}

helper.getName <- function(topGene){
  return(topGene$id)
}

disease.geneList <-function(toptable){
  return (list(positive=toptable$positive$genes,negative=toptable$negative$genes,id=toptable$id,platform=toptable$platform))
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
  
  top_table_rownames <- topGenes
  top_table_rownames$gene = rownames(top_table_rownames)
  positive_genes = filter(top_table_rownames, t > 0)
  negative_genes = filter(top_table_rownames, t < 0)
  
  #save an output file of the topgenes 
  dir.create("data")
  outFile<- paste("data/",datasetID,"top_table.txt")
  write.table(topGenes, file = outFile, quote=F, sep="\t", row.names=T)

  return(list(positive=positive_genes,negative=negative_genese, id=datasetID,platform=platformName))
}

