# Proportion of wind pollinated species with urbanisation

temp <- allergen_summary
anemophilous <- rownames(species_allergen)[which(species_allergen$anemophilous==1)]

temp$anemo.all.rich <- sapply(1:nrow(vegcomm), FUN = function(l) {
  p <- vegcomm[l, anemophilous]
  sp <-  names(p)[p > 0]
  result <- sum(species_allergen[species_allergen$Species %in% sp , "allergenicity" ], na.rm = T)
  return(result)
})

temp$anemo.all.prop <- temp$anemo.all.rich / temp$all.num

# Trend with urbanisation: NONE
par(mar =c(1,4,2,2), mfrow = c(2,2), oma = c(4,0,0,0))
plot(anemo.all.rich ~ Seal_500, temp)
add.stats(glm(anemo.all.rich ~ Seal_500, data = temp, family = poisson))
plot(anemo.all.prop ~ Seal_500, temp)
add.stats(lm(anemo.all.prop ~ Seal_500, data = temp)) # approx, should be binomial model


# For natives
temp$anemo.all.nat.rich <- sapply(1:nrow(vegcomm), FUN = function(l) {
  p <- vegcomm[l, intersect(anemophilous,natives)]
  sp <-  names(p)[p > 0]
  result <- sum(species_allergen[species_allergen$Species %in% sp , "allergenicity" ], na.rm = T)
  return(result)
})

temp$anemo.all.nat.prop <- temp$anemo.all.nat.rich / temp$all.num.nat

# Trend with urbanisation: NONE
plot(anemo.all.nat.rich ~ Seal_500, temp)
add.stats(glm(anemo.all.nat.rich ~ Seal_500, data = temp, family = poisson))
plot(anemo.all.nat.prop ~ Seal_500, temp)
add.stats(lm(anemo.all.nat.prop ~ Seal_500, data = temp)) # approx, should be binomial model


# For exotic allergens ####
temp$anemo.all.exo.rich <- sapply(1:nrow(vegcomm), FUN = function(l) {
  p <- vegcomm[l, intersect(anemophilous,exotics)]
  sp <-  names(p)[p > 0]
  result <- sum(species_allergen[species_allergen$Species %in% sp , "allergenicity" ], na.rm = T)
  return(result)
})

temp$anemo.exo.rich <- sapply(1:nrow(vegcomm), FUN = function(l) {
  p <- vegcomm[l, intersect(anemophilous,exotics)]
  sp <-  names(p)[p > 0]
  result <- length(sp)
  return(result)
})
temp$anemo.all.exo.cover <- sapply(1:nrow(vegcomm), FUN = function(l) {
  p <- vegcomm[l, ]
  sp <-  names(p)[p > 0]
  p <- p[sp]  # select only the present species
  all <- as.numeric(sp %in% intersect(intersect(allergenics, anemophilous), exotics))
  result = sum(p*all)
  return(result)
})

temp$anemo.exo.cover <- sapply(1:nrow(vegcomm), FUN = function(l) {
  p <- vegcomm[l, ]
  sp <-  names(p)[p > 0]
  p <- p[sp]  # select only the present species
  all <- as.numeric(sp %in% intersect(anemophilous, exotics))
  result = sum(p*all)
  return(result)
})
temp$anemo.all.exo.prop <- temp$anemo.all.exo.rich / temp$all.num.exo

# Trend of allergenics with urbanisation: slight INCREASE!
par(mar =c(1,4,2,2), mfrow = c(1,2), oma = c(4,0,0,0))

plot(anemo.all.exo.rich ~ Seal_500, temp)
add.stats(glm(anemo.all.exo.rich ~ Seal_500, data = temp, family = poisson))
add.stats(cor.test( ~ Seal_500 + anemo.all.exo.rich, data = temp, method = "spearman"))
# plot(anemo.all.exo.prop ~ Seal_500, temp)
# add.stats(lm(anemo.all.exo.prop ~ Seal_500, data = temp, na.action = "na.omit"))
plot(anemo.all.exo.cover ~ Seal_500, temp, log = "y")
add.stats(lm(anemo.all.exo.cover ~ Seal_500, data = temp, na.action = "na.omit"))


# neophytes

# For neophyte allergens ####
temp$anemo.neo.rich <- sapply(1:nrow(vegcomm), FUN = function(l) {
  p <- vegcomm[l, intersect(anemophilous,neophytes)]
  sp <-  names(p)[p > 0]
  result = length(sp)
  return(result)
})

temp$anemo.all.neo.rich <- sapply(1:nrow(vegcomm), FUN = function(l) {
  p <- vegcomm[l, intersect(anemophilous,neophytes)]
  sp <-  names(p)[p > 0]
  result <- sum(species_allergen[species_allergen$Species %in% sp , "allergenicity" ], na.rm = T)
  return(result)
})



temp$anemo.all.neo.cover <- sapply(1:nrow(vegcomm), FUN = function(l) {
  p <- vegcomm[l, ]
  sp <-  names(p)[p > 0]
  p <- p[sp]  # select only the present species
  all <- as.numeric(sp %in% intersect(intersect(anemophilous, neophytes), allergenics))
  result = sum(p*all)
  return(result)
})
# Trend of allergenics with urbanisation: NS increasing trend
par(mar =c(1,4,2,2), mfrow = c(1,2), oma = c(4,0,0,0))

plot(anemo.all.neo.rich ~ Seal_500, temp)
add.stats(glm(anemo.all.neo.rich ~ Seal_500, data = temp, family = poisson),plot.abline = "always")
plot(anemo.all.neo.cover ~ Seal_500, temp, log = "y")
add.stats(lm(anemo.all.neo.cover ~ Seal_500, data = temp, na.action = "na.omit"))
