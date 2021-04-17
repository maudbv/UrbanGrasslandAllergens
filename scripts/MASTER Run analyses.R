# Run analyses
# on Allergenic properties of urban grasslands
# April 2021

Sys.setlocale("LC_TIME", "en_US")

# _________ Load useful R packages:  ___________________

library(data.table)  
library(sp)    
library(leaps) 
library(vegan)

library(tidyr)
library(MuMIn)
library(r2glmm)
library(performance)
library(rsq)
library(MASS)
library(lme4)
library(np)
library(DescTools)

# color and graphical packages: 
library(inlmisc)
library(ggplot2)
library(vioplot)

# Load Data:
# source('scripts/Import_all_data.R')
load("data/Urban_Grassland_Allergens_data.Rdata")
print(paste("Allergen Data from", file.info("data/Urban_Grassland_Allergens_data.Rdata")$ctime, "is loaded"))

# My own utility functions:
source('scripts/utils/add.stats.R')
source('scripts/utils/p2star.R')
source('scripts/utils/cor.print.R')

# load specific functions to fit GLMs to our data: 
source('scripts/Analyses/fit.allergen.glms.R')

# Explore distribution of species in allergen space
source("scripts/Analyses/1. ordination allergenics in allfam space.R")

# Allergenic species diversity
source('scripts/Analyses/2a. Allergenic richness vs. DIVERSITY .R')
source('scripts/Analyses/2b. Allergenic richness vs NOVELTY.R')

# PAV of communities
source('scripts/Analyses/3. Mean PAV analyses.R')

# Allergen molecule diversity
plot.all.graphs = FALSE
source('scripts/Analyses/4a. Molecule diversity vs DIVERSITY.R')
source('scripts/Analyses/4b. Molecule diversity vs. NOVELTY.R')
source('scripts/Analyses/4c. ALLFAM diversity vs. NOVELTY.R')


# Beta-dissimilarity in species and allergen composition
source('scripts/Analyses/5a. Beta dissimilarity in species.R')
source('scripts/Analyses/5b. Beta dissimilarity in molecules.R')



# Phenology of allergen production
source('scripts/Analyses/6a. Analyses temporal spectrum.R')


# Monthly averaged values based on flowering phenology:
source('scripts/Analyses/7b. Monthly SRall with novelty.R')
source('scripts/Analyses/7c. Monthly AR with novelty.R')
source('scripts/Analyses/7d. Monthly AllfamR with novelty.R')

table.monthly.models <- rbind(
  SR = glms.SRall.month$glms.table,
  AR =  glms.AR.month$glms.table,
  AllfamR =  glms.allfam.month$glms.table)
table.monthly.models$level <- as.factor(table.monthly.models$level)
table.monthly.models <- doBy::orderBy(~ group - level +type ,
                                      table.monthly.models)
write.csv(table.monthly.models,
          file = "results/table monthly models.csv")

# Check spatial autocorrelation (Moran's I)
source('scripts/Analyses/test spatial autocorrelation.R')

# Save all results
save.image(file ="results/Urban_Grassland_Allergens_results.Rdata")
