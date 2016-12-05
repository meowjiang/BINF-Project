# Evaluate the MST predictions against a reference set
library(dplyr)
dir <- "/Users/rachelhz/BINF-Project"

################# Specify reference set to use ######################
# How many ATC characters to keep (choose from: 1, 3, 4, 5, 7)
atc_cutoff <- 5
# The similarity cutoff for ATC codes shared
sim_cutoff <- 2

################# Run me for one refset ######################
main <- function(dir, atc_cutoff, sim_cutoff) {
  all_data <- load_combine_data(dir, atc_cutoff, sim_cutoff)
  ref_vs_pred <- select(all_data, ref, pred)
  prec <- precision(ref_vs_pred)
  rec <- recall(ref_vs_pred)
  fp_rate <- fpr(ref_vs_pred)
  f_1_score <- f_1(ref_vs_pred)
  return(c(prec, rec, fp_rate, f_1_score))
}

main(dir, atc_cutoff, sim_cutoff)

################# Run me for several refsets ######################
multi_refsets <- list(c(7,1), c(7,2), c(5,1), c(5,2), c(4,1), c(4,2))
# Results are: precision, recall, FPR, F1
for (x in multi_refsets) {
  print(x)
  print(main(dir, atc_cutoff = x[1], sim_cutoff = x[2]))
}

################# Read in and combine ref and preds ######################
match_nodes <- function(preds, a, b) {
  match_row <- filter(preds, g1 == a, g2 == b)
  if (nrow(match_row) > 0) {
    return(match_row[['pred']])
  } else {
    match_row <- filter(preds, g2 == a, g1 == b)
    if (nrow(match_row) > 0) {
      return(match_row[['pred']])
    }
    return('NA')
  }
}

load_combine_data <- function(dir, atc_cutoff, sim_cutoff) {
  refset_file <- paste(paste('refset', as.character(atc_cutoff), as.character(sim_cutoff), sep = '_'), '.txt', sep = '')
  
  refset <- read.table(file = paste(dir, 'reference_sets', refset_file, sep = '/'), sep="\t", header = T, 
                       colClasses = c('character', 'character', 'integer', 'integer'))
  
  preds <- read.table(file = paste(dir, 'all', 'predictions.txt', sep = '/'), header = T, 
                      colClasses = c('character', 'character', 'integer'))
  
  all_data <- refset[c('node_a', 'node_b', 'ref')]
  all_data$pred <- apply(all_data, 1, function(x) match_nodes(preds, x['node_a'], x['node_b']))
  # Remove NAs (disease pairs that didn't match)
  all_data <- filter(all_data, pred != 'NA')
  
  return(all_data)
}

################# Evaluation Metrics ######################
true_positive <- function(ref_vs_pred) {
  tp <- filter(ref_vs_pred, ref == 1, pred == 1)
  return(nrow(tp))
}

false_positive <- function(ref_vs_pred) {
  fp <- filter(ref_vs_pred, ref == 0, pred == 1)
  return(nrow(fp))
}

false_negative <- function(ref_vs_pred) {
  fn <- filter(ref_vs_pred, ref == 1, pred == 0)
  return(nrow(fn))
}

true_negative <- function(ref_vs_pred) {
  tn <- filter(ref_vs_pred, ref == 0, pred == 0)
  return(nrow(tn))
}

precision <- function(ref_vs_pred) {
  tp <- true_positive(ref_vs_pred)
  fp <- false_positive(ref_vs_pred)
  return(tp / (tp + fp))
}

recall <- function(ref_vs_pred) {
  tp <- true_positive(ref_vs_pred)
  fn <- false_negative(ref_vs_pred)
  return(tp / (tp + fn))
}

f_1 <- function(ref_vs_pred) {
  prec <- precision(ref_vs_pred)
  rec <- recall(ref_vs_pred)
  return((2*prec*rec)/(prec + rec))
}

fpr <- function(ref_vs_pred) {
  fp <- false_positive(ref_vs_pred)
  tn <- true_negative(ref_vs_pred)
  return(fp / (fp + tn))
}
