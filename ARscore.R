####################################################################
## Upload packages
####################################################################

library(tidyverse)
library(fitdistrplus)
library(data.table)
library(janitor)

####################################################################
## Common functions
####################################################################

`%+%` <- function(left,right){
  return(paste0(left,right))
}

`%nin%` <- function(left,right) {
  return(!(left %in% right))
}

nna <- function(target) {return(!is.na(target))}

####################################################################
## Functions for Score Calculation
####################################################################

calc_scores <- function(norm_log, all_peptide_fcs, positives) {
  print("running ARscore algorithm")
  
  representations <- c(
    30, 60, 120, 240, 480, 960, 1960, 3920, 7840
    #, 3920, 7840 not needed for monkey virus, toxome
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
  
  # gamma distribution parameters change linearly with n
  shape_lm <- lm(log10(shape) ~ I(log10(total_peps)), data = dist_info)
  rate_lm <- lm(log10(rate) ~ I(log10(total_peps)), data = dist_info)
  mean_lm <- lm(log10(mean) ~ I(log10(total_peps)), data = dist_info)
  
  # add distribution info to norm_log
  norm_log <- norm_log %>%
    mutate(shape = 10^(log10(total_peps)*shape_lm$coefficients[2]+shape_lm$coefficients[1]),
           rate = 10^(log10(total_peps)*rate_lm$coefficients[2]+rate_lm$coefficients[1])) %>%
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
    pivot_longer(cols = c(contains('20A20G'), contains('20S'), contains('DPI')),names_to = 'sample_id') %>%
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
  for(i in c(2:length(algorithm_output))){
    scores <- rbind(scores, algorithm_output[[i]][[1]])
  }
  
  return(scores)
}

