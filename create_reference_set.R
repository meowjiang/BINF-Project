# Create reference set of the ATC codes (possibly with different character cut-offs) for each disease
dir <- "/Users/rachelhz/BINF-Project"
library(dplyr)

############# Variables to specify ##################
# How many ATC characters to keep (choose from: 1, 3, 4, 5, 7)
atc_cutoff <- 4
# The similarity cutoff for ATC codes shared
sim_cutoff <- 2

################# Load in data ######################
modify_indications <- function(indications, atc_cutoff) {
  # Only want disease-modifying indications
  indications<- filter(indications, category == 'DM')
  indications <- as.data.frame(indications)
  
  # Edit DOIDs to only the number
  indications$doid <- apply(indications, 1, function(x) sub('DOID:', '', x['doid_id']))
  
  # Edit ATCs to desired number of characters
  indications$atc_cut <- apply(indications, 1, function(x) substring(x['ATC.Code'], 1, atc_cutoff))
  
  return(indications)
}

############# Create and write reference set ##################
# Get the # of shared indications for DOIDs a, b
atc_intersection <- function(indications, a, b) {
  codes_a <- filter(indications, doid == a)$atc_cut
  codes_b <- filter(indications, doid == b)$atc_cut
  return(length(intersect(codes_a, codes_b)))
}

# Given the # of shared indications, assign 1 if num_shared >= cutoff, 0 o/w
disease_assoc_cutoff <- function(num_shared, cutoff) {
  if (num_shared >= cutoff) {
    return(1)
  }
  return(0)
}

create_refset <- function(indications, sim_cutoff, dir) {
  # Get all pair combinations for DOID
  all_pairs <- as.data.frame(t(combn(unique(indications$doid), 2)))
  colnames(all_pairs) <- c('node_a', 'node_b')
  
  # For each DOID pair (row), get the # of shared indications and indicate 0/1 for disease association
  all_pairs$shared <- apply(all_pairs, 1, function(x) 
    atc_intersection(indications, as.character(x['node_a']), as.character(x['node_b'])))
  all_pairs$ref <- apply(all_pairs, 1, function(x) disease_assoc_cutoff(as.numeric(x['shared']), sim_cutoff))
  
  # Write reference set to file with ATC code length and similarity cutoff
  file_name <- paste(paste('refset', as.character(atc_cutoff), as.character(sim_cutoff), sep = '_'), '.txt', sep = '')
  write.table(all_pairs, file = paste(dir, file_name, sep = '/'), quote=F, sep="\t", row.names=F)
}

################## Run me ########################
indications <- read.csv("~/Desktop/Rachel's/BINFG4006/project/catalog_RZ_indications.csv", header = T)
indications <- modify_indications(indications, atc_cutoff)
create_refset(indications, sim_cutoff, dir)

################## Other ########################
# How to read in a refset
foo <- read.table(file = paste(dir, 'reference_sets', 'refset_5_2.txt', sep = '/'), sep="\t", header = T, 
                  colClasses = c('character', 'character', 'integer', 'integer'))
