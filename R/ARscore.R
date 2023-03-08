####################################################################
## Upload packages
####################################################################

library(tidyverse)
library(fitdistrplus)
library(limma)
library(data.table)
library(janitor)

####################################################################
## Common functions
####################################################################

#' Left Right function to paste two strings
#'
#' @param left string
#' @param right string
#'
#' @return concatenated string
#' @export
#'
`%+%` <- function(left,right){
  return(paste0(left,right))
}

#' Left Right function to evaluate NOT in
#'
#' @param left vector
#' @param right vector
#'
#' @return boolean vector
#' @export
#'
`%nin%` <- function(left,right) {
  return(!(left %in% right))
}

#' function to evaluate NOT is.na
#'
#' @param target vector
#'
#' @return binary vector
#' @export
#'
nna <- function(target) {return(!is.na(target))}

####################################################################
## Functions for Score Calculation
####################################################################

#' Function to calculate aggregate reactivity scores
#'
#' Distributions are created by randomly selecting peptides
#' that were designed to viruses that are not the targets of
#' significant antibody responses
#' 
#' v 0.2.1 update - changed linear models to splines for interpolating
#' intermediate values of peptides
#'
#' @param norm_log dataframe containing intermediate reactivity metrics
#' @param all_peptide_fcs dataframe containing all peptide fold changes 
#' for an individual antibody profile
#' @param positives dataframe containing viruses determined to be 
#' significantly targeted by antibodies and therefore excluded from
#' random selection
#'
#' @return dataframe containing aggregate reactivity scores and data on 
#' random distributions
#' 
#'
#' @import tidyverse fitdistrplus
calc_scores <- function(norm_log, all_peptide_fcs, positives) {
  print("running ARscore algorithm")
  
  representations <- c(
    30, 42, 60, 85, 120, 170, 240,
    340, 480, 679, 960, 1358, 1960,
    2715, 3920, 5431, 7840
    #, 2715, 3920, 5431, 7840 not needed for monkey virus, toxome
    #, 15680 only needed for arboscan
  )
  
  all_peptide_fcs <- all_peptide_fcs %>% mutate(fc = log2(value))
  
  all_peptide_fcs <- all_peptide_fcs %>% 
    filter(taxon_genus %nin% positives$taxon_genus)
  
  
  distributions <- data.frame(matrix(nrow = length(representations), ncol = 1000))
  row.names(distributions) <- representations
  
  #calculate random distributions
  set.seed(120)
  for(i in c(1:length(representations))){
    n = row.names(distributions)[i]
    for(j in c(1:1000)){
      distributions[i, j] = mean(sample(all_peptide_fcs$fc, as.numeric(n), replace = F))
      if(distributions[i, j] == 0) {distributions[i,j] = .001} 
    }
    distributions[i, 1000] = mean(all_peptide_fcs$fc)
    #print("calculating distributions " %+% i %+% " of " %+% length(representations))
  }
  
  #model random distributions as a gamma distribution
  dist_info <- distributions[,c(1, 2, 3)]
  names(dist_info) <- c("total_peps", "shape", "rate")
  for(i in c(1:length(representations))){
    
    dist_info[i, 1] <- row.names(dist_info)[i] %>% as.numeric()
    
    gammafit  <-  fitdist(distributions[i,] %>% t() %>% as.numeric(), "gamma")
    
    dist_info[i, 2] <- gammafit[["estimate"]][["shape"]]
    dist_info[i, 3] <- gammafit[["estimate"]][["rate"]]
  }
  
  dist_info <- dist_info %>% mutate(mean = shape / rate) %>%
    mutate(variance = shape / rate^2)
  
  # gamma distribution parameters vary with n
  #shape_lm <- lm(log10(shape) ~ I(log10(total_peps)), data = dist_info)
  shape_spline <- smooth.spline(x = log10(dist_info$total_peps), y = log10(dist_info$shape), spar = .5)
  #rate_lm <- lm(log10(rate) ~ I(log10(total_peps)), data = dist_info)
  rate_spline <- smooth.spline(x = log10(dist_info$total_peps), y = log10(dist_info$rate), spar = .5)
  #mean_lm <- lm(log10(mean) ~ I(log10(total_peps)), data = dist_info)
  mean_spline <- smooth.spline(x = log10(dist_info$total_peps), y = log10(dist_info$mean), spar = .5)
  
  # add distribution info to norm_log
  norm_log <- norm_log %>%
    mutate(shape = 10^(predict(shape_spline, log10(total_peps))$y),
           rate = 10^(predict(rate_spline, log10(total_peps))$y)) %>%
    mutate(mean = shape / rate) %>%
    mutate(variance = shape / rate^2)
  
  # calculate scores and p values
  norm_log <- norm_log %>%
    mutate(vir_score = zscoreGamma(score_norm, shape = shape / sqrt(total_peps), rate = rate / sqrt(total_peps)),
           virus_fc = score_norm / mean) %>%
    mutate(vir_score = ifelse(vir_score < -10, -10, vir_score)) %>%
    mutate(p_val = -log10((1-pnorm(abs(zscoreGamma(score_norm, shape = shape, rate = rate)))))) %>%
    mutate(p_val = ifelse(p_val>15, 15, p_val))
  
  return(norm_log)
}


#' Shell function to iterate calculations of aggregate reactivity scores
#'
#' @param norm_log_1 dataframe containing intermediate reactivity metrics
#' @param all_peptide_fcs_1 dataframe containing all peptide fold changes 
#' for an individual antibody profile
#' @param max_iterations integer limiting the number of itererations to run
#' calc_scores. defaults to 10
#' @param p_cutoff numeric -log10 p value cutoff to determine significant  
#' targets of antibody reactivity. defaults to 4
#' @param score_cutoff numeric aggregate reactivity score cutoff to determine
#' signifiant targets of antibody reactivity. defaults to .72 the empircally 
#' determined threshold for VRC VirScan data. Recommended to adjust
#'
#' @return returns aggregate reactivity scores as well as a list denoting the
#' significant targets of reactivity detected during each iteration
#' 
#' @import tidyverse limma
iterative_scores <- function(norm_log_1, all_peptide_fcs_1, max_iterations = 10,
                              p_cutoff = -log10(.001), score_cutoff = .72) {
  iterations = 0
  positives_1 = norm_log_1 %>% filter(F)
  positives_output = list()
  library(limma)
  
  while(iterations < max_iterations) {
    print("iteration " %+% (iterations+1))
    scores <- calc_scores(norm_log = norm_log_1, all_peptide_fcs = all_peptide_fcs_1, positives = positives_1)
    
    #update variables
    positives_1 <- scores %>% filter(p_val > p_cutoff) %>% filter(vir_score > score_cutoff)
    positives_output[[iterations+1]] <- positives_1 
    
    #add break if no new positives
    if(iterations > 0) {
      if((positives_1$taxon_species %>% length()) <= (positives_output[[iterations]]$taxon_species %>% length())) {
        iterations <- iterations + 10
      }
    }
    
    iterations <- iterations + 1
  }
  
  outputs <- list()
  outputs[[1]] <- scores
  outputs[[2]] <- positives_output
  return(outputs) 
}

#' Calculate Aggregate Reactivity Scores
#'
#' @param hfc a dataframe containing hits fold change data from a PhIP-Seq
#' experiment. Used to exclude peptides that are reactive in beads only controls
#' @param fc a dataframe containing fold change data from a PhIP-Seq experiment.
#' Used to generate intermediate reactivity metrics, null distributions, and 
#' ARscores
#' @param set_max_iterations integer limiting the number of itererations to run
#' calc_scores. defaults to 10
#' @param set_p_cutoff numeric -log10 p value cutoff to determine significant  
#' targets of antibody reactivity. defaults to 4
#' @param set_score_cutoff numeric aggregate reactivity score cutoff to 
#' determine signifiant targets of antibody reactivity. defaults to .72 the 
#' empircally determined threshold for VRC VirScan data. Recommended to adjust
#' @param required_number_of_peptides integer that represents the minimum number
#' of peptides attributed to an antigen/pathogen required to calculate an 
#' ARscore
#'
#' @return a dataframe containing the aggregate reactivity scores for a PhIP run
#' @export
#'
#' @import tidyverse fitdistrplus limma
ARscore_algorithm <- function(hfc, fc, set_max_iterations = 10,
                               set_p_cutoff = -log10(.001), set_score_cutoff = .72,
                               required_number_of_peptides = 50) {
  
  ## peptides that have hits in beads
  bad_beads <- hfc %>% dplyr::select(pep_aa,contains('BEADS')) %>% pivot_longer(-pep_aa) %>%
    filter(value > 1) %>% dplyr::select(pep_aa) %>% pull
  
  ### fold change ###
  d1 <- fc
  
  ## make longform dataframe out of fc, floor fc at 1
  d2 <- d1 %>% filter(pep_aa %nin% bad_beads) %>% 
    dplyr::select(-contains('BEADS')) %>%
    #common pulldowns. Has to be edited in some screens
    pivot_longer(cols = c(contains('20A20G'), contains('20S'),
                          contains('ugIg'), contains('DPI')),names_to = 'sample_id') %>%
    mutate(value = ifelse(value<1,1,value))
  
  ids <- d2 %>% dplyr::select(sample_id) %>% unique()
  
  # need to group by taxon_species, sample_id and count number of reactivities. Then add back sample annotations 
  d4 <- d2 %>% group_by(taxon_species, taxon_genus, sample_id) %>% summarise(score = sum(log2(value))) %>% left_join(ids)
  d4_2 <- d2 %>% group_by(taxon_species, taxon_genus, sample_id) %>% summarise(total_peps = n()) %>%
    left_join(ids) %>% left_join(d4) %>% mutate(score = ifelse(is.na(score), 0, score))
  d4_2 <- d4_2 %>% mutate(score_norm = score / total_peps) %>% filter(total_peps >= required_number_of_peptides)
  
  library(tidyverse)
  library(fitdistrplus)
  library(limma)
  
  algorithm_output <- list()
  
  for(R in c(1:length(ids$sample_id))) {
    print(R)
    algorithm_output[[R]] <- iterative_scores(norm_log_1 = d4_2 %>% filter(sample_id == ids$sample_id[R]), 
                                              all_peptide_fcs_1 = d2 %>% filter(sample_id == ids$sample_id[R]),
                                              max_iterations = set_max_iterations,
                                              p_cutoff = set_p_cutoff, 
                                              score_cutoff = set_score_cutoff)
    
  }
  
  #algorithm output has the final output as well as a structure that tracks positive 
  scores <- algorithm_output[[1]][[1]]
  if(length(algorithm_output)>1){
  for(i in c(2:length(algorithm_output))){
    scores <- rbind(scores, algorithm_output[[i]][[1]])
  }
  }
  
  return(scores)
}

