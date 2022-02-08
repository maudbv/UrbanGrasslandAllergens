# REmove a species to try

# species to remove:
sp_to_rm <- "Medicago_x_varia"
# sp_to_rm <- "Solidago_canadensis"


# Load data
load("~/Documents/Work/R projects/UrbanGrasslandAllergens/clean data/Urban_Grassland_Allergens_data.Rdata")

# modify the data to remove species
vegcomm <- vegcomm[ ,-which(colnames(vegcomm) == sp_to_rm)]
neophytes <- neophytes[-which(neophytes == sp_to_rm)]
allergenics <- allergenics[-which(allergenics == sp_to_rm)]
vegcomm
exotics <- exotics[-which(exotics == sp_to_rm)]

# recalculate indices
source("~/Documents/Work/R projects/UrbanGrasslandAllergens/scripts/Data formatting/community allergenicity.R", echo=TRUE)

# load packages
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


# plot new results

# Figure 1: Cover of allergenic species along novelty gradients ####
par (
  mfrow = c(3, 2),
  mar = c(2, 1, 2, 1),
  oma = c(2, 5, 1, 1),
  mgp = c(2,0.5,0),
  tcl = - 0.3,
  las = 1
)
m = 1

# graph
for (k in 1:3) {
  tmp <- allergen_summary
  y <- c("cover.all", "cover.exo.all", "all.num.exo")[k]
  mod <- c("n","n","p")[k]
  leg <- c("Allergenic\n cover ",
           "Non-native allergenic\n cover",
           "Non-native allergenic\nspecies richness")[k]
  tmp$y <- tmp[, y]
  
  for (i in 1:2) {
    xleg <- c("% Impervious surfaces",
              "Proportion of neophytes")[i]
    x <- c("Seal_500", "prop.neo")[i]
    tmp$x <- tmp[, x]
    
    plot(y ~ x,
         data = tmp,
         pch = 20, 
         col = "grey",
         ann = FALSE
    )
    
    if(mod == "n"){ 
      add.stats(f <- lm(y ~ x,
                        data = tmp), type = "lm",l.col = "black")
    }
    if(mod == "p"){ 
      add.stats(f <- glm(y ~ x,
                         data = tmp,
                         family = poisson),
                type = "glm",
                l.col = "black")
    }
    # add.stats(f <- cor.test(~ y + x,
    #                         data = tmp),
    #           adj.stats = 0, col.stats = "grey")
    
    mtext(3,text = paste(letters[m], ")",sep = ""),
          adj = 0, font = 3, cex = 0.7)
    m = m + 1
    
    # x label
    if (k ==3) {
      mtext(1, text = xleg, cex = 0.75, line = 2.2)
    }
    
    # y label
    if (i == 1) {
      mtext(
        2,
        text = leg,
        outer = FALSE,
        cex = 0.75,
        line = 2.5, las = 0
      )
    }
  }
}





