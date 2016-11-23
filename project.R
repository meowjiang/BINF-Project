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


runMe <- function(){
  if(file.exists("data/error.txt")){
    file.remove("data/error.txt")
  }
  datasetIDS <- read.csv("ids",header = FALSE)
  IDList <-datasetIDS[,1]
  distanceMatrix <- main(levels(IDList))
  write.table(distanceMatrix,file="distanceMatrix.txt")
}

main <-function(diseaseIds){
  
  platforms <- basics.getPlatforms()
  topTables <- lapply(diseaseIds,disease.topGenes,platforms=platforms)
  topGenes <- lapply(topTables,disease.geneList)
  names(topGenes)<-lapply(topGenes,helper.getName)
  topGenes<-Filter(function(x)length(x$positive)>0||length(x$negative)>0,topGenes)
  distanceMatrix <- outer(topGenes, topGenes,FUN = basics.tanimoto.vectorized)
  return (distanceMatrix)
}

basics.getPlatforms <-function(){
  platformNames <- c("GPL570","GPL96")
  platforms<-lapply(platformNames,helper.getPlatform)
  names(platforms)<-platformNames
  return (platforms)
}

basics.tanimoto <- function(topGenes1,topGenes2){
  t1 <- topGenes1
  t2 <- topGenes2
  
  t1Positive <- t1$positive
  t1Negative <- t1$negative
  
  t2Positive <- t2$positive
  t2Negative <- t2$negative
  
  positiveIntersect <- intersect(t1Positive,t2Positive)
  negativeIntersect <- intersect(t1Negative,t2Negative)
  
  union <- (union(t1Positive,t1Negative))
  union<-(union(union,t2Positive))
  union<-(union(union,t2Negative))
  
  intersect = union(positiveIntersect,negativeIntersect)
  return (1-(length(intersect)/length(union)))
}

basics.tanimoto.vectorized <-Vectorize(basics.tanimoto,c("topGenes1","topGenes2"))

helper.getPlatform <- function(platformName){
  platform<- getGEO(platformName)
  return(platform)
}

helper.getName <- function(topGene){
  return(topGene$id)
}

helper.generateDesignMatrix <- function(geoSet,indices){
  controlPattern<-"benign nevi|high bone mineral density|uninfected|non-failing|control|normal|healthy|non-Alzheimer's|HIV-negative|non-obese"
  
  #No diseaseState, tissue instead
  if(unique(Meta(geoSet)$dataset_id)=="GDS4102")
    diseaseStates<-geoSet@dataTable@columns$tissue
  else
    diseaseStates<-geoSet@dataTable@columns$disease.state
  
  diseaseStates<-diseaseStates[indices]
  states<-unique(levels(diseaseStates))
  controlList <- grep(controlPattern,states,value=TRUE)
  diseaseList<-grep(controlPattern,states,value=TRUE,invert=TRUE)
  levels(diseaseStates) <- list(control=controlList, disease=diseaseList)
  design<-model.matrix(~ diseaseStates)
  
  return (design)
}

disease.geneList <-function(toptable){
  return (list(positive=toptable$positive$gene,negative=toptable$negative$gene,id=toptable$id,platform=toptable$platform))
}

helper.getHorribleIndices <-function(geoSet){
  
  id = unique(Meta(geoSet)$dataset_id)
  uglies<-list(GDS1615=list("normal","Crohn's disease"),
               GDS4358=list("control","HIV"),
               GDS1956=list("normal","amyotophic lateral sclerosis"),
               GDS4882=list("normal","hepatocellular carcinoma"),
               GDS5403=list("notmal","osteoarthritis"),
               GDS3874=list("healthy","type 1 diabetes"))
  
  #No diseaseState, tissue instead
  if(unique(Meta(geoSet)$dataset_id)=="GDS4102")
    diseaseStates<-geoSet@dataTable@columns$tissue
  else
    diseaseStates<-geoSet@dataTable@columns$disease.state
  
  if (id %in% names(uglies)){
    
    states <- uglies[[id]]
    indices <- which(diseaseStates %in% states  )
    
  }
  else{
    indices<- 1:length(diseaseStates)
  }
  return (indices)
  
}

disease.topGenes <- function(datasetID,platforms){
  dir.create("data")
  outFile<- paste("data/",datasetID,"top_table.txt",sep="")
  errorFile<- paste("data/","error.txt",sep="")
  
  geoSet <-getGEO(datasetID,getGPL = False) 
  platformName <- Meta(geoSet)$platform  
  platform <- platforms[[platformName]]  
  if(!file.exists(outFile)){
    geoData <- GDS2eSet(geoSet,do.log2 = TRUE,GPL=platform,getGPL = FALSE)
    geoMatrix <- as.matrix(geoData)
    indices<-helper.getHorribleIndices(geoSet)
    geoMatrix<-geoMatrix[,indices]
    mappings<- Table(platform)[,1:2]
    mappedMatrix <-   merge(geoMatrix,mappings,by.x="row.names",by.y="ID",all.x=TRUE)[,c(-1)] 
    aggregatedMatrix <-   aggregate(mappedMatrix[,-ncol(mappedMatrix)],by=list(GB_ACC=mappedMatrix$GB_ACC),FUN=mean)
    rownames(aggregatedMatrix)<-aggregatedMatrix[,1]

    #get rid of probes we couldn't match to a gene
    aggregatedMatrix<-filter(aggregatedMatrix,GB_ACC!="")
    
    aggregatedMatrix[,1]<-NULL
    design<-helper.generateDesignMatrix(geoSet,indices)
    fittedMatrix <-   lmFit(aggregatedMatrix,design)
    bayesOut <-   eBayes(fittedMatrix) 
    topGenes <-   topTable(bayesOut,p.value=.01,number = nrow(aggregatedMatrix))
    if(length(topGenes)==0)
      write.dcf(paste("uh oh, no top genes for ",datasetID),file = errorFile,append = TRUE)
    
    write.table(topGenes, file = outFile,quote=F, sep="\t", row.names=T,col.names = NA,)
  }
  else{
    #this is a bad idea, 3 is arbitrary
    if(file.info(outFile)$size<3){
      topGenes<-data.frame()
    }
    else{
      topGenes<- read.table(outFile)
      
    }
  }
  
  top_table_rownames <- topGenes
  top_table_rownames$gene = rownames(top_table_rownames)
  positive_genes = filter(top_table_rownames, t > 0)
  negative_genes = filter(top_table_rownames, t < 0)
  
  return(list(positive=positive_genes,negative=negative_genes, id=datasetID,platform=platformName))
}

