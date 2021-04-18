# Allergenic properties of Berlin grasslands

*Author*: Maud Bernard-Verdier
*Collaborators*: Birgit Seitz, Sascha Buchholz, Ingo Kowarik & Jonathan Jeschke

This repository contains the code and data to reproduce the analyses in the manuscript entitled *Allergens in grasslands increase with urbanisation and plant invasions*, currently in preparation.
This code analyses the allergenic properties of 56 plots of dry acidic grasslands in Berlin, Germany. This research work is part of the BIBS project, **Bridging in Biodiversity Science**, funded by the BMBF, Germany.

## Data
Raw data for the project are in the **data/** folder. The R script **script/Import_all_data.R** will import, clean and format all the data, and output four clean data tables in the clean data/ folder. Associated metadata are also provided as .csv files for these four tables.

**clean data/**:

    * Species traits and allergenicity (species_allergen_data.csv)
    * Allergen molecules and biochemical families (molecule_data.csv)
    * Species abundance per grassland plot (species_abundance_data.csv)
    * Environmental factors per plot (environmental_factors_data.csv)

## Analyses
The master script **script/MASTER Run analyses.R** will run all analyses sequentially to reproduce results anf figures from the article in preparation. Result tables and figures are stored in a **results/** folder.
An Rmardown file **UrbanAllergensAnalyses.Rmd** is also provided to illustrate the analyses, represent figures and a summary of results.
