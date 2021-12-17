# Master script for importing and cleaning data 
# on Allergenic properties of urban grasslands
# April 2021
# 
# ________________________ SET UP ______________________ ####

R.Version()$version.string # R version 3.6.2 (2019-12-12)
Sys.setlocale("LC_TIME", "en_US")
old.par <- par()

# _________ Load useful R packages:  ___________________####
library(data.table)  
library(vegan)
library(doBy)
library(tidyr)

#  ________________________ IMPORT DATA  _______________________####

# Import base data for cityscapelabs
source('scripts/Data import/import Cityscapelab data.R')

# Import species allergenicity data
source('scripts/Data import/import allergenicity data.R')

## Extract molecules from databases and update species allergenicity data

## Choose between definitions of allergen molecule names:
## broad == FALSE (species specific, and perfect matches with databases) 
## vs 
## broad == TRUE (genus level and additional data from the literature) 
broad <- TRUE
source('scripts/Data import/extract molecules.R')

#  ________________CALCULATE METRICS  ____________ ####

# Calculate Potential Allergenic Value (PAV)
source('scripts/Data formatting/calculate PAV.R')

# Calculate Family + Group level allergenicity data
source('scripts/Data formatting/family level comm.R')

# Molecule diversity per plot
source('scripts/Data formatting/molecules per plot.R')

# Calculate summary community metrics
source('scripts/Data formatting/community allergenicity.R')

# ADD flowering phenology community matrix 
source('scripts/Data formatting/Temporal spectrum.R')

# Clean up
rm(i, sp, spp,gen, genera, tmp, smat, pi, beg, end, x,y,z,tmp2)

## export formatted data tables:
source('scripts/Data formatting/Export formatted datatables.R')

# Remove molecule databases (too heavy)
rm(IUISallergens, AllergenOnline, SDADallergen)

## SAVE formatted data for analyses:
 save.image(file ="clean data/Urban_Grassland_Allergens_data.Rdata")
 print("Data imported and formatted in data file 'clean data/Urban_Grassland_Allergens_data.Rdata'")

