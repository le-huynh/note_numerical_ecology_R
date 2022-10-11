#'---
#' output:
#'     html_document:
#'       keep_md: TRUE
#'---

#+ message=FALSE
library(here)
library(tidyverse)
library(skimr)

#' The file Doubs.RData contains the following objects:  
load(here("data/Doubs.RData"))

#' ## `spe`: species (community) data frame (fish abundances)  
tibble(spe)
skim_without_charts(spe)

#' Minimum and maximum of abundance values in the whole data set 
range(spe)

#' Minimum and maximum value for each species 
apply(spe, 2, range)

spe %>% summarise(across(dplyr::everything(), range))

#' Number of absences 
sum(spe == 0)

#' Proportion of zeros in the community data set 
sum(spe == 0) / (nrow(spe) * ncol(spe))

#' At how many sites does each species occur? 
#' Calculate the relative occurrences of the species 
#' (proportions of the number of sites).

spe %>% pivot_longer(cols = everything(),
		     names_to = "species",
		     values_to = "value") %>%
	group_by(species) %>%
	nest() %>%
	mutate(site_presence = map_dbl(.x = data,
				       .f = function(x) {x %>% filter(value != 0) %>% nrow()}),
	       prob_site_presence = map_dbl(.x = site_presence,
	       			            .f = function(x) {x * 100 / nrow(spe)}))

#' how many species are present at each site (species richness)?
spe %>% rowwise() %>%
	transmute(presence_species = sum(c_across() > 0))


#' ## `env`: environmental data frame  
tibble(env)
skim_without_charts(env)


#' ## `spa`: spatial data frame – cartesian coordinates  
tibble(spa)
skim_without_charts(spa)


#' ## `fishtraits`: functional traits of fish species  
tibble(fishtraits)
skim_without_charts(fishtraits)


#' ## `latlong`: spatial data frame – latitude and longitude  
tibble(latlong)
skim_without_charts(latlong)


