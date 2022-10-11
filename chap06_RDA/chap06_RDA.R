#'---
#' output:
#'     html_document:
#'       keep_md: TRUE
#'---

#+ message=FALSE
library(here)
library(tidyverse)
library(ade4)
library(adegraphics)
library(adespatial)
library(vegan)
library(vegan3d)
library(MASS)
library(ellipse)
library(FactoMineR)
library(rrcov)

#' Load the data and functions  
load(here("data/Doubs.RData"))

source(here("functions/hcoplot.R"))
source(here("functions/triplot.rda.R"))
source(here("functions/plot.lda.R"))
source(here("functions/polyvars.R"))
source(here("functions/screestick.R"))

#' ## Prepare the data
#' Remove empty site 8
spe <- spe %>% slice(-8)
env <- env %>% slice(-8)
spa <- spa %>% slice(-8)

#' Set aside variable 'dfs' (distance from the source) for later use 
dfs <- env %>% dplyr::select(dfs)

#' Remove 'dfs' variable from `env` data frame
env2 <- env %>% dplyr::select(!dfs)

#' Recode the slope variable (slo) into a factor (qualitative) 
#' variable to show how these are handled in the ordinations 
slo2 <- rep(".very_steep", nrow(env))
slo2[env$slo <= quantile(env$slo)[4]] <- ".steep" 
slo2[env$slo <= quantile(env$slo)[3]] <- ".moderate" 
slo2[env$slo <= quantile(env$slo)[2]] <- ".low" 
slo2 <- factor(slo2, levels = c(".low", ".moderate", ".steep", ".very_steep"))
table(slo2)

#' Create `env3` data frame with slope as a qualitative variable 
env3 <- env2 
env3$slo <- slo2

#' Create two subsets of explanatory variables
#' Subset 1: Physiography (upstream-downstream gradient)
envtopo <- env2[, c(1 : 3)]
names(envtopo)

#' Subset 2: Water quality
envchem <- env2[, c(4 : 10)]
names(envchem)

#' Hellinger-transform the species data set 
spe.hel <- decostand(spe, "hellinger")

#' ## RDA using `vegan`
(spe.rda <- rda(spe.hel ~ ., env3))
summary(spe.rda) # Scaling 2 (default)

#' Canonical coefficients from the rda object
coef(spe.rda)

#' Unadjusted R^2 retrieved from the rda object 
(R2 <- RsquareAdj(spe.rda)$r.squared) 
#' Adjusted R^2 retrieved from the rda object
(R2adj <- RsquareAdj(spe.rda)$adj.r.squared)


#' ## plot RDA
#' Scaling 1 
plot(spe.rda,
     scaling = 1,
     display = c("sp", "lc", "cn"),
     main = "Triplot RDA spe.hel ~ env3 - scaling 1 - lc scores")

spe.sc1 <- scores(spe.rda,
		  choices = 1:2,
		  scaling = 1,
		  display = "sp")

arrows(0, 0,
       spe.sc1[, 1] * 0.92,
       spe.sc1[, 2] * 0.92,
       length = 0,
       lty = 1,
       col = "red")

#' Scaling 2
plot(spe.rda,
     display = c("sp", "lc", "cn"),
     main = "Triplot RDA spe.hel ~ env3 - scaling 2 - lc scores")

spe.sc2 <- scores(spe.rda,
		  choices = 1:2,
		  display = "sp")

arrows(0, 0,
       spe.sc2[, 1] * 0.92,
       spe.sc2[, 2] * 0.92,
       length = 0,
       lty = 1,
       col = "red")






