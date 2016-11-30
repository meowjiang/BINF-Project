# Detect bias in expression values between platforms

setwd("~/Desktop/Rachel's/BINFG4006/project")
all_diseases <- read.csv('diseases_1130.csv', header = T)

# Filter for diseases on GPL570 and GPL96 for now
library(tidyverse)

#diseases <- filter(all_diseases, Platform == 'GPL570' | Platform == 'GPL96')
diseases <- filter(all_diseases, Platform != '')

dataset_list <- unique(diseases$Accession.Number)

library(affy);
library(limma);
library(GEOquery);

main <- function(datasetID) {
  geoSet <-getGEO(datasetID) 
  platformName <- Meta(geoSet)$platform
  geoData <- GDS2eSet(geoSet, do.log2 = TRUE)
  geoMatrix <- as.matrix(geoData)
  return(geoMatrix)
}

# This aggregates all data
get_expression_values <- function(datasetID) {
  geoSet <-getGEO(datasetID) 
  platformName <- Meta(geoSet)$platform  
  geoData <- GDS2eSet(geoSet, do.log2 = TRUE)
  geoMatrix <- as.matrix(geoData)
  
  vector_platform <- data.frame(exp = as.vector(geoMatrix), platform = platformName)
  return(vector_platform)
}

# This returns the mean of each dataset
get_expression_values_mean <- function(datasetID) {
  geoSet <-getGEO(datasetID) 
  platformName <- Meta(geoSet)$platform  
  geoData <- GDS2eSet(geoSet, do.log2 = TRUE)
  geoMatrix <- as.matrix(geoData)
  mean_values <- rowMeans(geoMatrix)
  
  labeled_means <- data.frame(exp = as.vector(mean_values), platform = platformName)
  return(labeled_means)
}

# Put all expression values, labeled by platform, in a giant data frame
giant_exp_platform <- data.frame()

# Be careful about which function is used to process datasets
for (datasetID in dataset_list) {
  giant_exp_platform <- rbind(giant_exp_platform, get_expression_values_mean(datasetID))
}

write.table(giant_exp_platform, file = "~/Desktop/Rachel's/BINFG4006/project/giant_exp_platform_mean_all.txt", quote=F, sep="\t", row.names=F)
#giant_exp_platform <- read.table(file = "~/Desktop/Rachel's/BINFG4006/project/giant_exp_platform", sep="\t", header = T)

#dataset_list[1:10] %>%
#    lapply(foo) %>%
#    rbind(giant_exp_platform) -> giant_exp_platform

fit <- glm(exp ~ as.factor(platform), data = giant_exp_platform)

# Write output to file
sink('glm_mean_all.txt')
summary(fit)
sink()