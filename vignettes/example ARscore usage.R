
library(tidyverse)
library(ARscore)

#################
## data
#################
### Upload hits fold change ###
virscan_hfc <- ARscore::example_virscan_hfc
toxscan_hfc <- ARscore::example_toxscan_hfc

### Upload fold change ###
virscan_fc <- ARscore::example_virscan_fc
toxscan_fc <- ARscore::example_toxscan_fc

################
## calculate ARscores
################

VARscores <- ARscore_algorithm(hfc = virscan_hfc, fc = virscan_fc)
TARscores <- ARscore_algorithm(hfc = toxscan_hfc, fc = toxscan_fc)

#############
## make wide (if many samples were run)
#############

p_wide <- VARscores %>% ungroup() %>% dplyr::select(taxon_genus, taxon_species, total_peps, p_val, sample_id) %>% 
  pivot_wider(names_from = sample_id, values_from = p_val)
VARscore_wide <- VARscores %>% ungroup() %>% dplyr::select(taxon_genus, taxon_species, total_peps, vir_score, sample_id) %>% 
  pivot_wider(names_from = sample_id, values_from = vir_score)

