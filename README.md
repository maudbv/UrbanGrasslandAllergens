# Allergenic properties of Berlin grasslands

*Author*: Maud Bernard-Verdier
*Collaborators*: Birgit Seitz, Sascha Buchholz, Ingo Kowarik & Jonathan Jeschke

This repository contains the code and data to reproduce the analyses in the manuscript entitled *Allergens in grasslands increase with urbanisation and plant invasions*, currently in preparation.
This code analyses the allergenic properties of 56 plots of dry acidic grasslands in Berlin, Germany. This research work is part of the BIBS project, **Bridging in Biodiversity Science**, funded by the BMBF, Germany.

<img src="https://user-images.githubusercontent.com/6454302/115233650-763e1800-a118-11eb-85f5-00cddf7cff4c.png" alt="Map of Berlin 56 grasslands" width="600">

## Data
Raw data for the project are in the **data/** folder. The R script **script/Import_all_data.R** will import, clean and format all the data, and output four clean data tables in the clean data/ folder. Associated metadata are also provided as .csv files for these four tables.

**clean data/**:

    * Species traits and allergenicity (species_allergen_data.csv)
    * Allergen molecules and biochemical families (molecule_data.csv)
    * Species abundance per grassland plot (species_abundance_data.csv)
    * Environmental factors per plot (environmental_factors_data.csv)

## Analyses
The master script **script/MASTER Run analyses.R** will run all analyses sequentially to reproduce results anf figures from the article in preparation. Result tables and figures are stored in a **results/** folder.
An Rmardown document **UrbanAllergensAnalyses.Rmd** is also provided to create a report summarising results with illustrations.

<img src="https://user-images.githubusercontent.com/6454302/115234901-da151080-a119-11eb-89ba-4b752414cfd6.png" alt="Figure 4" width="600">
