source("https://bioconductor.org/biocLite.R"); ## try http:// if https:// URLs are not supported
biocLite("affy");
biocLite("limma");
biocLite("GEOquery");
biocLite("proxy")
library(affy);
library(limma);
library(GEOquery); 
library(proxy)


GDS4513 <-getGEO("GDS4513")
GDS4491 <-getGEO("GDS4491")
GDS1615 <-getGEO("GDS1615")
datasetID<-"GDS1615"
geoSet <-getGEO(datasetID) 

design<-model.matrix(~ Columns(dataTable(geoSet))$disease.state)

# we need to be able to see if it's diseased or normal

#run through all of them
#print out the ones that don't have a control 
#print out the ones that aren't binary
datasetIDS <- read.csv("ids",header = FALSE)
IDList <-datasetIDS[,1]

design_matrix <- cbind(design_matrix[,1], pmax(design_matrix[,2], design_matrix[,3]))

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
  #This check doesn't make sense given that control is not the only word for a control
  
  #if(!any(grepl("control",ds))){
  #  print("FAIL:No Control")
  #  writeLines("FAIL:No Control",con=conn)
  #}
  #print(paste("Disease States:",toString(levels(unique(ds)))))
  
  #writeLines(paste("Disease States:",toString(levels(unique(ds)))),con=conn)
  writeLines("Disease States:",con=conn1)
  writeLines(toString(names(summary(ds))),con=conn1)
  writeLines(toString(summary(ds)),con=conn1)
  
  #GDS4513 no-relapse vs relapse?
  #GDS2205 non-failing vs failing
  #GDS5098 Seminoma, yolk sac tumor?
  #GDS2362 presymptomatic, experimentally acquired, symptomatic, naturally acquired, uninfected
  #GDS3353 high bone mineral density, low bone mineral density
  #GDS1287 high BMD, low BMD
  #GDS1523 clear cell adenocarcinoma, serous adenocarcinoma
  #GDS4336 cancer death: 0, cancer death: 1, cancer death: na
  #GDS3326 dedifferentiated liposarcoma, fibrosarcoma, leiomyosarcoma, lipoma, malignant fibrous histiocytoma, malignant peripheral nerve sheath tumor, myxofibrosarcoma, myxoid liposarcoma, synovial sarcoma, well-differentiated liposarcoma
  #GDS4337 non-diabetic, T2D
  #GDS4456 no recurrence/DOD, recurrence/DOD

  #known list of patterns for controls
  controlPattern<-"control|normal|healthy|non-Alzheimer's|HIV-negative|non-obese"
  
  
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
  
  #Can we take out some factors and associated values completely
  
  
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

#junk to figure out merging columns
#GDS1615 <-getGEO("GDS1615")
#originalDS<-GDS1615@dataTable@columns$disease.state
#ds  <-GDS1615@dataTable@columns$disease.state
#dsstates <-unique(levels(ds))
#this will set the levels
#we just need two lists one for things that will be control the other for things that are not control 
#control vs not control should be based on a list not one word
#controlPattern<-"control|normal|healthy"
#controlList <- grep(controlPattern,states,value=TRUE)
#diseaseList<-grep(controlPattern,states,value=TRUE,invert=TRUE)

#design<-model.matrix(~ Columns(dataTable(geoSet))$disease.state)


#levels(ds) <- list(control=controlList, disease=diseaseList)
#d1<-model.matrix(~originalDS)
#d2<-model.matrix(~ds)
# do it for a dataset that is already perfect
#see if it changes and whether it has just intercept and something 
#rerun stats and see which ones fail to be transformed by the new stuff
#add rest of things to control regex
#copy stuff to hozukimaru

