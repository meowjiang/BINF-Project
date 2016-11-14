source("https://bioconductor.org/biocLite.R"); ## try http:// if https:// URLs are not supported
biocLite("affy");
biocLite("limma");
biocLite("GEOquery");
biocLite("proxy")
library(affy);
library(limma);
library(GEOquery); 
library(proxy)


datasetIDS <- read.csv("ids",header = FALSE)
IDList <-datasetIDS[,1]


stats <-function(ids){
  fname="datasetStats.txt"
  fname2="datasetStats-2.txt"
  file.remove(fname)
  file.remove(fname2)
  fileConn1<-file(fname,open="a")
  fileConn2<-file(fname2,open="a")
  lapply(ids, outputStats,conn1=fileConn1,conn2=fileConn2)
  close(fileConn1)
  close(fileConn2)
}

outputStats <-function(id,conn1,conn2){
  geo<-getGEO(id)
  print(paste(id,":"))
  writeLines(paste(id,":"),con=conn1)
  writeLines(paste(id,":"),con=conn2)
  

  
  
  if(!(any(grep("disease.state",names(geo@dataTable@columns))))){
    print("FAIL: No disease.state" )
    writeLines("FAIL: No disease.state",con=conn1);
    writeLines("FAIL: No disease.state",con=conn2);
    writeLines("\n",con=conn1)
    writeLines("\n",con=conn2)
    
    
    return()
  }
  
  ds=geo@dataTable@columns$disease.state
  if(length(unique(ds))!=2){
    print("FAIL: Disease state non-binary")
    writeLines("FAIL: Disease state non-binary",con = conn1)
  }
  
  #writeLines(paste("Disease States:",toString(levels(unique(ds)))),con=conn)
  writeLines("Disease States:",con=conn1)
  writeLines(toString(names(summary(ds))),con=conn1)
  writeLines(toString(summary(ds)),con=conn1)


  #known list of patterns for controls
  controlPattern<-"benign nevi|high bone mineral density|uninfected|non-failing|control|normal|healthy|non-Alzheimer's|HIV-negative|non-obese"
  
  
  #If there was no control based on the known names we have write out an error
  if(!any(grep(controlPattern,ds))){
    writeLines("FAIL:No Control",con=conn1)
    writeLines("FAIL:No Control",con=conn2)
  }
  
  
  
  #Now write file for transformed disease states
  states<-unique(levels(ds))
  controlList <- grep(controlPattern,states,value=TRUE)
  diseaseList<-grep(controlPattern,states,value=TRUE,invert=TRUE)
  levels(ds) <- list(control=controlList, disease=diseaseList)
  
  
  if(length(unique(ds))!=2){
    print("FAIL: Disease state non-binary")
    writeLines("FAIL: Disease state non-binary",con = conn2)
  }
  writeLines("Transformed Disease States:",con=conn2)
  writeLines(toString(names(summary(ds))),con=conn2)
  writeLines(toString(summary(ds)),con=conn2)
  

  
  
  writeLines("\n",con=conn1)
  writeLines("\n",con=conn2)
}


