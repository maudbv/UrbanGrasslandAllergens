# Calculate community allergenicity scores and spectrums

# summary table of number of allergenic species per plot ####
# Initiate table
allergen_summary <- data.frame(plot = rownames(vegcomm))
rownames(allergen_summary) <- allergen_summary$plot

# Number of allergen species per plot 
allergen_summary$all.num <- sapply(1:nrow(vegcomm), FUN = function(l) {
  p <- vegcomm[l, ]
  sp <-  names(p)[p > 0]
  result <- sum(species_allergen[species_allergen$Species %in% sp , "allergenicity" ], na.rm = T)
  return(result)
})

# Number of allergen neophyte species per plot
allergen_summary$all.num.neo <- sapply(1:nrow(vegcomm), FUN = function(l) {
  p <- vegcomm[l, neophytes]
  sp <-  names(p)[p > 0]
  result <- sum(species_allergen[species_allergen$Species %in% sp , "allergenicity" ], na.rm = T)
  return(result)
})

# Number of allergen archaeophyte species per plot
allergen_summary$all.num.arc <- sapply(1:nrow(vegcomm), FUN = function(l) {
  p <- vegcomm[l, archaeophytes]
  sp <-  names(p)[p > 0]
  result <- sum(species_allergen[species_allergen$Species %in% sp , "allergenicity" ], na.rm = T)
  return(result)
})

# Number of allergen native species per plot
allergen_summary$all.num.nat <- sapply(1:nrow(vegcomm), FUN = function(l) {
  p <- vegcomm[l, natives]
  sp <-  names(p)[p > 0]
  result <- sum(species_allergen[species_allergen$Species %in% sp , "allergenicity" ], na.rm = T)
  return(result)
})

# Number of allergen native species per plot
allergen_summary$all.num.exo <- sapply(1:nrow(vegcomm), FUN = function(l) {
  p <- vegcomm[l, exotics]
  sp <-  names(p)[p > 0]
  result <- sum(species_allergen[species_allergen$Species %in% sp , "allergenicity" ], na.rm = T)
  return(result)
})

# Number of SEVERELY allergenic species per plot
allergen_summary$severe.all.num <- sapply(1:nrow(vegcomm), FUN = function(l) {
  p <- vegcomm[l, ]
  sp <-  names(p)[p > 0]
  result <- sum(species_allergen[species_allergen$Species %in% sp , "severe.all" ], na.rm = T)
  return(result)
})

# COVER summary per plot  #####
# Total cover of allergenic species per plot
allergen_summary$cover.all<- sapply(1:nrow(vegcomm), FUN = function(l) {
  p <- vegcomm[l, ]
  sp <-  names(p)[p > 0]
  p <- p[sp]  # select only the present species
  all <- as.numeric(sp %in% allergenics)
  result = sum(p*all)
  return(result)
})


# Relative cover of allergenic species per plot
allergen_summary$pcover.all<- sapply(1:nrow(vegcomm), FUN = function(l) {
  p <- vegcomm[l, ]
  sp <-  names(p)[p > 0]
  p <- p[sp]  # select only the present species
  all <- as.numeric(sp %in% allergenics)
  result = sum(p*all)/sum(p)
  return(result)
})

#  Total cover of SEVERELY allergenic species per plot
allergen_summary$cover.severe.all<- sapply(1:nrow(vegcomm), FUN = function(l) {
  p <- vegcomm[l, ]
  sp <-  names(p)[p > 0]
  p <- p[sp]
  all <- species_allergen[species_allergen$Species %in% sp , "severe.all" ]
  result = sum(p*all)
  return(result)
})

# Relative cover of SEVERELY allergenic species per plot
allergen_summary$pcover.severe.all<- sapply(1:nrow(vegcomm), FUN = function(l) {
  p <- vegcomm[l, ]
  sp <-  names(p)[p > 0]
  p <- p[sp]  # select only the present species
  all <- species_allergen[species_allergen$Species %in% sp , "severe.all" ]
  result = sum(p*all)/sum(p)
  return(result)
})


#  cover of ARCHAEO allergenic species per plot
allergen_summary$cover.arc.all<- sapply(1:nrow(vegcomm), FUN = function(l) {
  p <- vegcomm[l, ]
  sp <-  names(p)[p > 0]
  p <- p[sp]  # select only the present species
  all <- sp %in% allergenics * sp %in% archaeophytes
  result = sum(p*all)
  return(result)
})

#  cover of NATIVE allergenic species per plot
allergen_summary$cover.nat.all<- sapply(1:nrow(vegcomm), FUN = function(l) {
  p <- vegcomm[l, ]
  sp <-  names(p)[p > 0]
  p <- p[sp]  # select only the present species
  all <- sp %in% allergenics * sp %in% natives
  result = sum(p*all)
  return(result)
})

#  cover of NEOPHYTE allergenic species per plot
allergen_summary$cover.neo.all<- sapply(1:nrow(vegcomm), FUN = function(l) {
  p <- vegcomm[l, ]
  sp <-  names(p)[p > 0]
  p <- p[sp]  # select only the present species
  all <- sp %in% allergenics * sp %in% neophytes
  result = sum(p*all)
  return(result)
})


#  cover of EXOTIC allergenic species per plot
allergen_summary$cover.exo.all<- sapply(1:nrow(vegcomm), FUN = function(l) {
  p <- vegcomm[l, ]
  sp <-  names(p)[p > 0]
  p <- p[sp]  # select only the present species
  all <- sp %in% allergenics * sp %in% exotics
  result = sum(p*all)
  return(result)
})

# Relative cover of RESIDENTS allergenic species per plot
allergen_summary$pcover.resid.all<- sapply(1:nrow(vegcomm), FUN = function(l) {
  p <- vegcomm[l, ]
  sp <-  names(p)[p > 0]
  p <- p[sp]  # select only the present species
  all <- sp %in% allergenics * sp %in% residents
  result = sum(p*all)/sum(p)
  return(result)
})


# Relative cover of EXOTICS allergenic species per plot
allergen_summary$pcover.exo.all<- sapply(1:nrow(vegcomm), FUN = function(l) {
  p <- vegcomm[l, ]
  sp <-  names(p)[p > 0]
  p <- p[sp]  # select only the present species
  all <- sp %in% allergenics * sp %in% exotics
  result = sum(p*all)/sum(p)
  return(result)
})
# Relative cover of ARCHAEO allergenic species per plot
allergen_summary$pcover.arc.all<- sapply(1:nrow(vegcomm), FUN = function(l) {
  p <- vegcomm[l, ]
  sp <-  names(p)[p > 0]
  p <- p[sp]  # select only the present species
    all <- sp %in% allergenics * sp %in% archaeophytes
    result = sum(p*all)/sum(p)
  return(result)
})

# Relative cover of NATIVE allergenic species per plot
allergen_summary$pcover.nat.all<- sapply(1:nrow(vegcomm), FUN = function(l) {
  p <- vegcomm[l, ]
  sp <-  names(p)[p > 0]
  p <- p[sp]  # select only the present species
  all <- sp %in% allergenics * sp %in% natives
  result = sum(p*all)/sum(p)
  return(result)
})

# Relative cover of NEOPHYTE allergenic species per plot
allergen_summary$pcover.neo.all<- sapply(1:nrow(vegcomm), FUN = function(l) {
  p <- vegcomm[l, ]
  sp <-  names(p)[p > 0]
  p <- p[sp]  # select only the present species
  all <- sp %in% allergenics * sp %in% neophytes
  result = sum(p*all)/sum(p)
  return(result)
})

# MEAN ALLERGENIC SCORES per plot ####

# Overall mean allergenic score per plot
allergen_summary$mean.score <- sapply(1:nrow(vegcomm), FUN = function(l) {
  p <- vegcomm[l, ]
  sp <-  names(p)[p > 0]
  result <- mean(species_allergen[species_allergen$Species %in% sp , "allergen_score" ], na.rm = T)
  return(result)
})

# Native mean allergenic score per plot
allergen_summary$mean.nat.score <- sapply(1:nrow(vegcomm), FUN = function(l) {
  p <- vegcomm[l,natives ]
  sp <-  names(p)[p > 0]
  result <- mean(species_allergen[species_allergen$Species %in% sp , "allergen_score" ], na.rm = T)
  return(result)
})

# Archaeo mean allergenic score per plot
allergen_summary$mean.arc.score <- sapply(1:nrow(vegcomm), FUN = function(l) {
  p <- vegcomm[l,archaeophytes ]
  sp <-  names(p)[p > 0]
  result <- mean(species_allergen[species_allergen$Species %in% sp , "allergen_score" ], na.rm = T)
  return(result)
})

# Residents (nat + arc) mean allergenic score per plot
allergen_summary$mean.resid.score <- sapply(1:nrow(vegcomm), FUN = function(l) {
  p <- vegcomm[l,residents ]
  sp <-  names(p)[p > 0]
  result <- mean(species_allergen[species_allergen$Species %in% sp , "allergen_score" ], na.rm = T)
  return(result)
})

# neophyte mean allergenic score per plot
allergen_summary$mean.neo.score <- sapply(1:nrow(vegcomm), FUN = function(l) {
  p <- vegcomm[l,neophytes ]
  sp <-  names(p)[p > 0]
  result <- mean(species_allergen[species_allergen$Species %in% sp , "allergen_score" ], na.rm = T)
  return(result)
})

# exotic mean allergenic score per plot
allergen_summary$mean.exo.score <- sapply(1:nrow(vegcomm), FUN = function(l) {
  p <- vegcomm[l,exotics ]
  sp <-  names(p)[p > 0]
  result <- mean(species_allergen[species_allergen$Species %in% sp , "allergen_score" ], na.rm = T)
  return(result)
})

# maximum score per plot ####
allergen_summary$max.all.score <- sapply(1:nrow(vegcomm), FUN = function(l) {
  p <- vegcomm[l, ]
  sp <-  names(p)[p > 0]
  result <- max(species_allergen[species_allergen$Species %in% sp , "allergen_score" ], na.rm = T)
  if (is.infinite(result)) result = NA
  return(result)
})

# variance in allergenic score per plot ####
allergen_summary$variance.score <- sapply(1:nrow(vegcomm), FUN = function(l) {
  p <- vegcomm[l, ]
  sp <-  names(p)[p > 0]
  result <- var(species_allergen[species_allergen$Species %in% sp , "allergen_score" ], na.rm = T)
  return(result)
})

# CWM score per plot ####
# mean allergenic score per plot
allergen_summary$CWM.score <- sapply(1:nrow(vegcomm), FUN = function(l) {
  p <- vegcomm[l, ]
  sp <-  names(p)[p > 0]
  p <- p[sp]
  s <- species_allergen[species_allergen$Species %in% sp , "allergen_score" ]
  result <- sum(p*s, na.rm=TRUE)/sum(p)
  return(result)
})


# CWM allergenic natives score per plot
allergen_summary$CWM.nat.score <- sapply(1:nrow(vegcomm), FUN = function(l) {
  p <- vegcomm[l,natives ]
  sp <-  names(p)[p > 0]
  p <- p[sp]
  s <- species_allergen[species_allergen$Species %in% sp , "allergen_score" ]
  result <- sum(p*s, na.rm=TRUE)/sum(p)
  return(result)
})
# CWM allergenic resident score per plot
allergen_summary$CWM.resid.score <- sapply(1:nrow(vegcomm), FUN = function(l) {
  p <- vegcomm[l,residents ]
  sp <-  names(p)[p > 0]
  p <- p[sp]
  s <- species_allergen[species_allergen$Species %in% sp , "allergen_score" ]
  result <- sum(p*s, na.rm=TRUE)/sum(p)
  return(result)
})

# CWM allergenic neophyte score per plot
allergen_summary$CWM.neo.score <- sapply(1:nrow(vegcomm), FUN = function(l) {
  p <- vegcomm[l,neophytes ]
  if ( !all(p==0)) {
  sp <-  names(p)[p > 0]
  p <- p[sp]
  s <- species_allergen[species_allergen$Species %in% sp , "allergen_score" ]
 result <- sum(p*s, na.rm = T)/sum(p,na.rm = T)
 } else result = NA
  return(result)
})

# CWM allergenic archaeophytes score per plot
allergen_summary$CWM.arc.score <- sapply(1:nrow(vegcomm), FUN = function(l) {
  p <- vegcomm[l,archaeophytes]
  if ( !all(p==0)) {
    sp <-  names(p)[p > 0]
    p <- p[sp]
    s <- species_allergen[species_allergen$Species %in% sp , "allergen_score" ]
    result <- sum(p*s, na.rm = T)/sum(p,na.rm = T)
  } else result = NA
  return(result)
})


# CWM allergenic exotics score per plot
allergen_summary$CWM.exo.score <- sapply(1:nrow(vegcomm), FUN = function(l) {
  p <- vegcomm[l,exotics ]
  if ( !all(p==0)) {
    sp <-  names(p)[p > 0]
    p <- p[sp]
    s <- species_allergen[species_allergen$Species %in% sp , "allergen_score" ]
    result <- sum(p*s, na.rm = T)/sum(p,na.rm = T)
  } else result = NA
  return(result)
})




# MEAN PAV per plot ####

# Overall mean allergenic pav per plot
allergen_summary$mean.pav <- sapply(1:nrow(vegcomm), FUN = function(l) {
  p <- vegcomm[l, ]
  sp <-  names(p)[p > 0]
  result <- mean(species_allergen[species_allergen$Species %in% sp , "PAV" ], na.rm = T)
  return(result)
})

# Native mean allergenic pav per plot
allergen_summary$mean.nat.pav <- sapply(1:nrow(vegcomm), FUN = function(l) {
  p <- vegcomm[l,natives ]
  sp <-  names(p)[p > 0]
  result <- mean(species_allergen[species_allergen$Species %in% sp , "PAV" ], na.rm = T)
  return(result)
})

# Archaeo mean allergenic pav per plot
allergen_summary$mean.arc.pav <- sapply(1:nrow(vegcomm), FUN = function(l) {
  p <- vegcomm[l,archaeophytes ]
  sp <-  names(p)[p > 0]
  result <- mean(species_allergen[species_allergen$Species %in% sp , "PAV" ], na.rm = T)
  return(result)
})

# Residents (nat + arc) mean allergenic pav per plot
allergen_summary$mean.resid.pav <- sapply(1:nrow(vegcomm), FUN = function(l) {
  p <- vegcomm[l,residents ]
  sp <-  names(p)[p > 0]
  result <- mean(species_allergen[species_allergen$Species %in% sp , "PAV" ], na.rm = T)
  return(result)
})

# neophyte mean allergenic pav per plot
allergen_summary$mean.neo.pav <- sapply(1:nrow(vegcomm), FUN = function(l) {
  p <- vegcomm[l,neophytes ]
  sp <-  names(p)[p > 0]
  result <- mean(species_allergen[species_allergen$Species %in% sp , "PAV" ], na.rm = T)
  return(result)
})

# Exotic mean allergenic pav per plot
allergen_summary$mean.exo.pav <- sapply(1:nrow(vegcomm), FUN = function(l) {
  p <- vegcomm[l,exotics ]
  sp <-  names(p)[p > 0]
  result <- mean(species_allergen[species_allergen$Species %in% sp , "PAV" ], na.rm = T)
  return(result)
})
# maximum PAV per plot ####
allergen_summary$max.all.pav <- sapply(1:nrow(vegcomm), FUN = function(l) {
  p <- vegcomm[l, ]
  sp <-  names(p)[p > 0]
  result <- max(species_allergen[species_allergen$Species %in% sp , "PAV" ], na.rm = T)
  if (is.infinite(result)) result = NA
  return(result)
})

# variance in PAV per plot ####
allergen_summary$variance.pav <- sapply(1:nrow(vegcomm), FUN = function(l) {
  p <- vegcomm[l, ]
  sp <-  names(p)[p > 0]
  result <- var(species_allergen[species_allergen$Species %in% sp , "PAV" ], na.rm = T)
  return(result)
})

# CWM PAV per plot ####
# mean allergenic pav per plot
allergen_summary$CWM.pav <- sapply(1:nrow(vegcomm), FUN = function(l) {
  p <- vegcomm[l, ]
  sp <-  names(p)[p > 0]
  p <- p[sp]
  s <- species_allergen[species_allergen$Species %in% sp , "PAV" ]
  result <- sum(p*s, na.rm=TRUE)/sum(p, na.rm=TRUE)
  # Make sure we remove NA values for species with no PAV (no pollen)
  return(result)
})


# CWM allergenic natives pav per plot
allergen_summary$CWM.nat.pav <- sapply(1:nrow(vegcomm), FUN = function(l) {
  p <- vegcomm[l,natives ]
  sp <-  names(p)[p > 0]
  p <- p[sp]
  s <- species_allergen[species_allergen$Species %in% sp , "PAV" ]
  result <- sum(p*s, na.rm=TRUE)/sum(p)
  return(result)
})

# CWM allergenic resident pav per plot
allergen_summary$CWM.resid.pav <- sapply(1:nrow(vegcomm), FUN = function(l) {
  p <- vegcomm[l,residents ]
  sp <-  names(p)[p > 0]
  p <- p[sp]
  s <- species_allergen[species_allergen$Species %in% sp , "PAV" ]
  result <- sum(p*s, na.rm=TRUE)/sum(p)
  return(result)
})

# CWM allergenic neophyte pav per plot
allergen_summary$CWM.neo.pav <- sapply(1:nrow(vegcomm), FUN = function(l) {
  p <- vegcomm[l,neophytes ]
  if ( !all(p==0)) {
    sp <-  names(p)[p > 0]
    p <- p[sp]
    s <- species_allergen[species_allergen$Species %in% sp , "PAV" ]
    result <- sum(p*s, na.rm = T)/sum(p,na.rm = T)
  } else result = NA
  return(result)
})

# CWM allergenic archaeophytes pav per plot
allergen_summary$CWM.arc.pav <- sapply(1:nrow(vegcomm), FUN = function(l) {
  p <- vegcomm[l,archaeophytes]
  if ( !all(p==0)) {
    sp <-  names(p)[p > 0]
    p <- p[sp]
    s <- species_allergen[species_allergen$Species %in% sp , "PAV" ]
    result <- sum(p*s, na.rm = T)/sum(p,na.rm = T)
  } else result = NA
  return(result)
})


# CWM allergenic exotics pav per plot
allergen_summary$CWM.exo.pav <- sapply(1:nrow(vegcomm), FUN = function(l) {
  p <- vegcomm[l,exotics ]
  if ( !all(p==0)) {
    sp <-  names(p)[p > 0]
    p <- p[sp]
    s <- species_allergen[species_allergen$Species %in% sp , "PAV" ]
    result <- sum(p*s, na.rm = T)/sum(p,na.rm = T)
  } else result = NA
  return(result)
})

## Calculate diversity (number & cover) of allergenic molecules per plot ####
allergen_summary$nb.mol <- rowSums(mol_plot)[rownames(allergen_summary)]
allergen_summary$cover.mol <- rowSums(mol_plot_wtd)[rownames(allergen_summary)]

allergen_summary$nb.mol.neo <- rowSums(mol_plot_neo)[rownames(allergen_summary)]
allergen_summary$cover.mol.neo <- rowSums(mol_plot_neo_wtd)[rownames(allergen_summary)]

allergen_summary$nb.mol.nat <- rowSums(mol_plot_nat)[rownames(allergen_summary)]
allergen_summary$cover.mol.nat <- rowSums(mol_plot_nat_wtd)[rownames(allergen_summary)]

allergen_summary$nb.mol.arc <- rowSums(mol_plot_arc)[rownames(allergen_summary)]
allergen_summary$cover.mol.arc <- rowSums(mol_plot_arc_wtd)[rownames(allergen_summary)]

allergen_summary$nb.mol.exo <- rowSums(mol_plot_exo)[rownames(allergen_summary)]
allergen_summary$cover.mol.exo <- rowSums(mol_plot_exo_wtd)[rownames(allergen_summary)]

allergen_summary$nb.allfam <- rowSums(allfam_plot)[rownames(allergen_summary)]
allergen_summary$nb.allfam.wtd <- rowSums(allfam_plot_wtd)[rownames(allergen_summary)]

allergen_summary$nb.allfam.nat <- rowSums(allfam_plot_nat)[rownames(allergen_summary)]
allergen_summary$nb.allfam.neo <- rowSums(allfam_plot_neo)[rownames(allergen_summary)]
allergen_summary$nb.allfam.arc <- rowSums(allfam_plot_arc)[rownames(allergen_summary)]
allergen_summary$nb.allfam.exo <- rowSums(allfam_plot_exo)[rownames(allergen_summary)]


## Merge with plot_summary ####
# => add environmental data and general vegetation information per plot  
allergen_summary <- merge(allergen_summary,
                          plot_summary,
                          by.x = "plot", by.y = "ID_plot")
rownames(allergen_summary) <- allergen_summary$plot

# make sure the plots are in the same order: 
allergen_summary  <- allergen_summary[rownames(plot_summary),]

# Calculate the percentage vesrion of "prop.neo" for Poisson regressions:
allergen_summary$perc.neo = 100 * allergen_summary$prop.neo

# Add  number of wind pollinated species ####
allergen_summary$SR.wind <- rowSums(vegcomm[, rownames(species_allergen)[which(species_allergen$anemophilous == 1)]]>0)
allergen_summary$prop.wind <- allergen_summary$SR.wind/allergen_summary$SR


allergen_summary$SR.neo.wind <- rowSums(
  vegcomm[,
          rownames(species_allergen)[
            which(species_allergen$anemophilous == 1 &
                    species_allergen$neophyte == 1)]]>0)

allergen_summary$prop.neo.wind <- allergen_summary$SR.neo.wind/allergen_summary$SR



## Taxonomic family richness ####
allergen_summary$FR <- rowSums(vegcomm.fam>0)
allergen_summary$FR.neo <- rowSums(vegcomm.fam.neo>0)
allergen_summary$FR.nat <- rowSums(vegcomm.fam.nat>0)
allergen_summary$FR.arc <- rowSums(vegcomm.fam.arc>0)

allergen_summary$FR.all <- rowSums(vegcomm.fam.all>0)
allergen_summary$FR.all.neo <- rowSums(vegcomm.fam.all.neo>0)
allergen_summary$FR.all.nat <- rowSums(vegcomm.fam.all.nat>0)
allergen_summary$FR.all.arc <- rowSums(vegcomm.fam.all.arc>0)

# Taxonomic Family cover ####
allergen_summary$Fcover <- rowSums(vegcomm.fam)
allergen_summary$Fcover.neo <- rowSums(vegcomm.fam.neo)
allergen_summary$Fcover.nat <- rowSums(vegcomm.fam.nat)
allergen_summary$Fcover.arc <- rowSums(vegcomm.fam.arc)

allergen_summary$Fcover.all <- rowSums(vegcomm.fam.all)
allergen_summary$Fcover.all.neo <- rowSums(vegcomm.fam.all.neo)
allergen_summary$Fcover.all.nat <- rowSums(vegcomm.fam.all.nat)
allergen_summary$Fcover.all.arc <- rowSums(vegcomm.fam.all.arc)

