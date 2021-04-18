# Test for spatial autocorrelation in data and model residuals: 
library(sp)
library(ape)

# Calculate plot distances ####
plot_dist <- spDists(as.matrix(plot_summary[,c("Long","Lat")]), longlat = TRUE)
colnames(plot_dist) <- rownames(plot_dist) <- rownames(plot_summary)

# Calculate autocorrelation within data: ####
stopifnot(all(rownames(allergen_summary) == rownames(plot_dist)))

moranI.data <- list(
# Allergenic species richness
Allergenic_SR = Moran.I(allergen_summary$all.num,
                        plot_dist ), # NS
Allergenic_Neophyte_SR = Moran.I(allergen_summary$all.num.neo,
                                 plot_dist ), # NS

# Allergen molecule richness
Allergen_Richness = Moran.I(allergen_summary$nb.mol,
                            plot_dist ), # NS
Neophyte_AllergenRichness = Moran.I(allergen_summary$nb.mol.neo,
                                    plot_dist ), # NS
# Environmental gradients:
Proportion_Neophytes = Moran.I(allergen_summary$prop.neo, plot_dist ), # NS
Perc.Impervious =Moran.I(allergen_summary$Seal_500, plot_dist) # ***
# some spatial autocorrelation of % impervious surfaces
)

#export list in numeric dataframe:
x = purrr:::reduce(moranI.data, rbind.data.frame) 
row.names(x) = names(moranI.data)
moranI.data = x


# Calculate autocorrelation within model residuals with % sealing ####

# family richness
f <- glm(FR ~ Seal_500, allergen_summary, family = poisson)
stopifnot(all(names(f$residuals) == rownames(plot_dist)))
Moran.I(f$residuals, plot_dist ) # NS

# Allergenic cover
f <- lm(cover.all ~ Seal_500, allergen_summary)
stopifnot(all(names(f$residuals) == rownames(plot_dist)))
Moran.I(f$residuals, plot_dist ) # NS :)

# Allergenic species richness
f <-  glm(all.num ~ Seal_500, allergen_summary, family = poisson)
stopifnot(all(names(f$residuals) == rownames(plot_dist)))
Moran.I(f$residuals, plot_dist ) # NS :)

# Allergen molecule richness
f <- glm.nb(nb.mol ~ Seal_500, allergen_summary)
stopifnot(all(names(f$residuals) == rownames(plot_dist)))
Moran.I(f$residuals, plot_dist ) # NS :)

f <- glm.nb(nb.mol.neo ~ Seal_500, allergen_summary)
stopifnot(all(names(f$residuals) == rownames(plot_dist)))
Moran.I(f$residuals, plot_dist ) # NS :)

# Allergen family richness
f <- glm(nb.allfam ~ Seal_500, allergen_summary, family = poisson)
stopifnot(all(names(f$residuals) == rownames(plot_dist)))
Moran.I(f$residuals, plot_dist ) # **

# It seems like the number of allergen families is spatially autocorrelated, BUT: not the case for the subset of neophytes:

f <- glm(nb.allfam.neo ~ Seal_500, allergen_summary, family = poisson)
stopifnot(all(names(f$residuals) == rownames(plot_dist)))
Moran.I(f$residuals, plot_dist ) # NS
