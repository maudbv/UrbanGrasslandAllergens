# Run analyses
# on Allergenic properties of urban grasslands
# April 2021

Sys.setlocale("LC_TIME", "en_US")

# _________ Load useful R packages:  ____________#### 

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

# Load my own utility functions: 
source('scripts/utils/add.stats.R')
source('scripts/utils/p2star.R')
source('scripts/utils/cor.print.R')

# load specific functions to fit GLMs to our data: 
source('scripts/Analyses/fit.allergen.glms.R')

# Load Data: ####
# source('scripts/Import_all_data.R')
load("clean data/Urban_Grassland_Allergens_data.Rdata")
print(paste("Allergen Data from", file.info("clean data/Urban_Grassland_Allergens_data.Rdata")$ctime, "is loaded"))

# 1 - Explore distribution of species in allergen space ####
source("scripts/Analyses/1. ordination allergenics in allfam space.R")

# 2 - Allergenic species diversity ####
source('scripts/Analyses/2a. Allergenic richness vs. DIVERSITY .R')
source('scripts/Analyses/2b. Allergenic richness vs NOVELTY.R')

# 3 - PAV of communities ####
source('scripts/Analyses/3. Mean PAV analyses.R')

# 4 - Allergen molecule diversity ####
plot.all.graphs = FALSE
source('scripts/Analyses/4a. Molecule diversity vs DIVERSITY.R')
source('scripts/Analyses/4b. Molecule diversity vs. NOVELTY.R')
source('scripts/Analyses/4c. ALLFAM diversity vs. NOVELTY.R')

# 5 - Beta-dissimilarity in species and allergen composition ####
source('scripts/Analyses/5a. Beta dissimilarity in species.R')
source('scripts/Analyses/5b. Beta dissimilarity in molecules.R')
source('scripts/Analyses/5c. Beta dissimilarity in Allergen families.R')

# 6 - Phenology of allergen production ####
source('scripts/Analyses/6a. Analyses flowering phenology.R')

## Monthly averaged values based on flowering phenology: 
source('scripts/Analyses/6b. Monthly SRall with novelty.R')
source('scripts/Analyses/6c. Monthly AR with novelty.R')
source('scripts/Analyses/6d. Monthly AllfamR with novelty.R')

## exported table for monthly models
table.monthly.models <- rbind(
  SR = glms.SRall.month$glms.table,
  AR =  glms.AR.month$glms.table,
  AllfamR =  glms.allfam.month$glms.table)
table.monthly.models$level <- as.factor(table.monthly.models$level)
table.monthly.models <- doBy::orderBy(~ group - level +type ,
                                      table.monthly.models)
write.csv(table.monthly.models,
          file = "results/table monthly models.csv")

# Plot figures for manuscript ####
source('scripts/Illustration/Figures manuscript.R')

# Save all results ####
save.image(file ="results/Urban_Grassland_Allergens_results.Rdata")
