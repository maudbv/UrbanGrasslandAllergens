# Test for spatial autocorrelation: 
library(sp)
library(ape)

# Calculate plot distances
plot_dist <- spDists(as.matrix(plot_summary[,c("Long","Lat")]), longlat = TRUE)
colnames(plot_dist) <- rownames(plot_dist) <- rownames(plot_summary)

# Calculate autocorrelation with  data: 
stopifnot(all(rownames(allergen_summary) == rownames(plot_dist)))
Moran.I(allergen_summary$all.num, plot_dist ) # NS
Moran.I(allergen_summary$all.num.neo, plot_dist ) # NS
Moran.I(allergen_summary$prop.neo, plot_dist ) # NS
Moran.I(allergen_summary$Seal_500, plot_dist ) # ***
# some spatial autocorrelation of % impervious surfaces

# Calculate autocorrelation with model residuals


f <- glm(FR ~ Seal_500, allergen_summary, family = poisson)
stopifnot(all(names(f$residuals) == rownames(plot_dist)))
Moran.I(f$residuals, plot_dist ) # NS

f <- lm(cover.all ~ Seal_500, allergen_summary)
stopifnot(all(names(f$residuals) == rownames(plot_dist)))
Moran.I(f$residuals, plot_dist ) # NS :)

f <-  glm(all.num ~ Seal_500, allergen_summary, family = poisson)
stopifnot(all(names(f$residuals) == rownames(plot_dist)))
Moran.I(f$residuals, plot_dist ) # NS :)

f <- glm.nb(nb.mol.neo ~ Seal_500, allergen_summary)
stopifnot(all(names(f$residuals) == rownames(plot_dist)))
Moran.I(f$residuals, plot_dist ) # NS :)


f <- glm(nb.allfam ~ Seal_500, allergen_summary, family = poisson)
stopifnot(all(names(f$residuals) == rownames(plot_dist)))
Moran.I(f$residuals, plot_dist ) # **
# It seems like the number of allergen families is spatially autocorrelated
# => only for natives and archaeophytes, not for neophytes
