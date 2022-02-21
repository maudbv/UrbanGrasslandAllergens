### TEMPORAL SPECTRUM of ALLERGIES
# Take flowering phenology into account to calculate the 
# mean monthly cover and richness of allergenic species, molecules and molecule families



# Check and format the flowering phenology data:
rownames(species_allergen)[which(!apply(!is.na(cbind(species_allergen$fl.beg,species_allergen$fl.end)),1,all))]

# Missing data for one species, corrected from Bioflor database (Sept. 2020)
species_allergen["Taraxacum_sect._Ruderalia","fl.end"] <- 10
species_allergen["Taraxacum_sect._Ruderalia","fl.period"] <- 8

# Matrix of monthly flowering phenology of species  ####
months <- tolower(format(seq.Date(as.Date('2000-01-01'), by = 'month', len = 12), "%B"))


flower_month <- data.frame(matrix(0, ncol =12, nrow = nrow(species_allergen),
                                    dimnames = list(species_allergen$Species,
                                                    months)))

for (sp in rownames(species_allergen)) {
  beg = tolower(species_allergen[sp, c("Flowering_beg")])
  end = tolower(species_allergen[sp, c("Flowering_end")])
  if(  (beg!= "") | (end!= "")){
    flower_month[ sp, match(beg, months):match(end, months)] <- 1
  }
}

names(flower_month)  <-format(seq.Date(as.Date('2000-01-01'), by = 'month', len = 12), "%b")

# Create a matrix of maximum phenology for each molecule ####
# not useful at plot level as pheno is per Species 
# But givers overall pattern of seasonality
mol_month_strict <- t(sapply(colnames(mol_mat_strict), function(m) {
  sp <- rownames(mol_mat_strict)[mol_mat_strict[,m]>0]
  return(m = colSums(flower_month[sp,]))
}))

mol_month <- t(sapply(colnames(mol_mat), function(m) {
  sp <- rownames(mol_mat)[mol_mat[,m]>0]
  return(m = colSums(flower_month[sp,]))
}))

allfam_month <- t(sapply(colnames(allfam_mat), function(m) {
  sp <- rownames(allfam_mat)[allfam_mat[,m]>0]
  return(m = colSums(flower_month[sp,]))
}))

allfam_month_nat <- t(sapply(colnames(allfam_mat), function(m) {
  sp <- rownames(allfam_mat)[allfam_mat[,m]>0]
  sp <- intersect(sp,natives)
  return(m = colSums(flower_month[sp,]))
}))
allfam_month_exo <- t(sapply(colnames(allfam_mat), function(m) {
  sp <- rownames(allfam_mat)[allfam_mat[,m]>0]
  sp <- intersect(sp, exotics)
  return(m = colSums(flower_month[sp,]))
}))

allfam_month_neo <- t(sapply(colnames(allfam_mat), function(m) {
  sp <- rownames(allfam_mat)[allfam_mat[,m]>0]
  sp <- intersect(sp,neophytes)
  return(m = colSums(flower_month[sp,]))
}))


allfam_month_arc <- t(sapply(colnames(allfam_mat), function(m) {
  sp <- rownames(allfam_mat)[allfam_mat[,m]>0]
  sp <- intersect(sp, archaeophytes)
  return(m = colSums(flower_month[sp,]))
}))



# Matrices of species flowering per month in each plot ####
fl_plot <- as.matrix((vegcomm)>0) %*% as.matrix(flower_month)
fl_neo_plot <- as.matrix((vegcomm[, (neophytes)])>0) %*%
  as.matrix(flower_month[(neophytes),])
fl_nat_plot <- as.matrix((vegcomm[, (natives)])>0) %*%
  as.matrix(flower_month[(natives),])
fl_arc_plot <- as.matrix((vegcomm[, (archaeophytes)])>0) %*%
  as.matrix(flower_month[(archaeophytes),])
fl_exo_plot <- as.matrix((vegcomm[, (exotics)])>0) %*%
  as.matrix(flower_month[(exotics),])


fl_plot_wtd <- as.matrix(vegcomm) %*% as.matrix(flower_month)
fl_neo_plot_wtd <- as.matrix(vegcomm[, (neophytes)]) %*%
  as.matrix(flower_month[(neophytes),])
fl_nat_plot_wtd <- as.matrix(vegcomm[, (natives)]) %*%
  as.matrix(flower_month[(natives),])
fl_arc_plot_wtd <- as.matrix(vegcomm[, (archaeophytes)]) %*%
  as.matrix(flower_month[(archaeophytes),])
fl_exo_plot_wtd <- as.matrix(vegcomm[, (exotics)]) %*%
  as.matrix(flower_month[(exotics),])

# Matrices of allergenic species flowering per month in each plot ####
fl_all_plot <- as.matrix((vegcomm[, allergenics])>0) %*% as.matrix(flower_month[allergenics,])
fl_all_neo_plot <- as.matrix((vegcomm[, intersect(allergenics,neophytes)])>0) %*%
  as.matrix(flower_month[intersect(allergenics,neophytes),])
fl_all_nat_plot <- as.matrix((vegcomm[, intersect(allergenics,natives)])>0) %*%
  as.matrix(flower_month[intersect(allergenics,natives),])
fl_all_arc_plot <- as.matrix((vegcomm[, intersect(allergenics,archaeophytes)])>0) %*%
  as.matrix(flower_month[intersect(allergenics,archaeophytes),])
fl_all_exo_plot <- as.matrix((vegcomm[, intersect(allergenics,exotics)])>0) %*%
  as.matrix(flower_month[intersect(allergenics,exotics),])


fl_all_plot_wtd <- as.matrix(vegcomm[, allergenics]) %*% as.matrix(flower_month[allergenics,])
fl_all_neo_plot_wtd <- as.matrix(vegcomm[, intersect(allergenics,neophytes)]) %*%
  as.matrix(flower_month[intersect(allergenics,neophytes),])
fl_all_nat_plot_wtd <- as.matrix(vegcomm[, intersect(allergenics,natives)]) %*%
  as.matrix(flower_month[intersect(allergenics,natives),])
fl_all_arc_plot_wtd <- as.matrix(vegcomm[, intersect(allergenics,archaeophytes)]) %*%
  as.matrix(flower_month[intersect(allergenics,archaeophytes),])

fl_all_exo_plot_wtd <- as.matrix(vegcomm[, intersect(allergenics,exotics)]) %*%
  as.matrix(flower_month[intersect(allergenics,exotics),])

# Matrices of number of molecules produced per month in each plot ####
fl_mol_plot <- t(sapply(rownames(vegcomm), function(m) {
  sp <- colnames(vegcomm)[vegcomm[m,]>0]
  mols <- as.matrix(t(mol_mat[sp,]))
  fls <- as.matrix(flower_month[sp,])
  out <- colSums(mols %*%fls)
  return(out)
  }))
fl_mol_neo_plot <- t(sapply(rownames(vegcomm), function(m) {
  sp <- colnames(vegcomm)[vegcomm[m,]>0]
  sp <- intersect(sp, neophytes) # select only neophytes
  mols <- as.matrix(t(mol_mat[sp,]))
  fls <- as.matrix(flower_month[sp,])
  out <- colSums(mols %*%fls)
  return(out)
}))
fl_mol_nat_plot <- t(sapply(rownames(vegcomm), function(m) {
  sp <- colnames(vegcomm)[vegcomm[m,]>0]
  sp <- intersect(sp, natives) 
  mols <- as.matrix(t(mol_mat[sp,]))
  fls <- as.matrix(flower_month[sp,])
  out <- colSums(mols %*%fls)
  return(out)
}))
fl_mol_arc_plot <- t(sapply(rownames(vegcomm), function(m) {
  sp <- colnames(vegcomm)[vegcomm[m,]>0]
  sp <- intersect(sp, archaeophytes) 
  mols <- as.matrix(t(mol_mat[sp,]))
  fls <- as.matrix(flower_month[sp,])
  out <- colSums(mols %*%fls)
  return(out)
}))
fl_mol_exo_plot <- t(sapply(rownames(vegcomm), function(m) {
  sp <- colnames(vegcomm)[vegcomm[m,]>0]
  sp <- intersect(sp, exotics) 
  mols <- as.matrix(t(mol_mat[sp,]))
  fls <- as.matrix(flower_month[sp,])
  out <- colSums(mols %*%fls)
  return(out)
}))

# Weighted by percent cover: 
fl_mol_plot_wtd <- t(sapply(rownames(vegcomm), function(m) {
  sp <- colnames(vegcomm)[vegcomm[m,]>0]
  cov <- as.numeric(vegcomm[m,sp])
  mols <- as.matrix(t(mol_mat[sp,]))
  fls <- as.matrix(flower_month[sp,])
  out <- colSums(mols %*% (as.numeric(cov)*fls))
  return(out)
}))
fl_mol_neo_plot_wtd <- t(sapply(rownames(vegcomm), function(m) {
  sp <- colnames(vegcomm)[vegcomm[m,]>0]
  sp <- intersect(sp, neophytes) # select only neophytes
  cov <- as.numeric(vegcomm[m,sp])
  mols <- as.matrix(t(mol_mat[sp,]))
  fls <- as.matrix(flower_month[sp,])
  out <- colSums(mols %*% (as.numeric(cov)*fls))
  return(out)
}))
fl_mol_arc_plot_wtd <- t(sapply(rownames(vegcomm), function(m) {
  sp <- colnames(vegcomm)[vegcomm[m,]>0]
  sp <- intersect(sp, archaeophytes) # select only archaeophytes
  cov <- as.numeric(vegcomm[m,sp])
  mols <- as.matrix(t(mol_mat[sp,]))
  fls <- as.matrix(flower_month[sp,])
  out <- colSums(mols %*% (as.numeric(cov)*fls))
  return(out)
}))
fl_mol_nat_plot_wtd <- t(sapply(rownames(vegcomm), function(m) {
  sp <- colnames(vegcomm)[vegcomm[m,]>0]
  sp <- intersect(sp, natives) # select only natives
  cov <- as.numeric(vegcomm[m,sp])
  mols <- as.matrix(t(mol_mat[sp,]))
  fls <- as.matrix(flower_month[sp,])
  out <- colSums(mols %*% (as.numeric(cov)*fls))
  return(out)
}))
fl_mol_exo_plot_wtd <- t(sapply(rownames(vegcomm), function(m) {
  sp <- colnames(vegcomm)[vegcomm[m,]>0]
  sp <- intersect(sp, exotics) # select only exotics
  cov <- as.numeric(vegcomm[m,sp])
  mols <- as.matrix(t(mol_mat[sp,]))
  fls <- as.matrix(flower_month[sp,])
  out <- colSums(mols %*% (as.numeric(cov)*fls))
  return(out)
}))


# Number of UNIQUE allfam per month per plot ####
## Subsets along gradient: 
z = plot_summary$Seal_levels
table(z)


fl_allfam_plot  <- t(sapply(rownames(vegcomm), function(m) {
  presence = colSums(vegcomm[m,])>0
  sp <- names(presence)[presence]
  fls <- as.matrix(flower_month[sp,])
  allfams <- allfam_mat[sp,]
  fam.month <- apply( fls, 2, function(x) colSums(allfams[ x == 1,])>0)
  out <- colSums(fam.month)
  return(out)
}))


fl_allfam_neo_plot<- t(sapply(rownames(vegcomm), function(m) {
  presence = colSums(vegcomm[m,])>0
  sp <- names(presence )[presence]
  sp <- intersect(sp, neophytes) # select only neophytes
  
  fls <- as.matrix(flower_month[sp,])
  allfams <- allfam_mat[sp,]
  fam.month <- apply( fls, 2, function(x) colSums(allfams[ x == 1,])>0)
  out <- colSums(fam.month)
  return(out)
}))

fl_allfam_nat_plot <- t(sapply(rownames(vegcomm), function(m) {
  presence = colSums(vegcomm[m,])>0
  sp <- names(presence )[presence]
  sp <- intersect(sp, natives) # select only natives
  
  fls <- as.matrix(flower_month[sp,])
  allfams <- allfam_mat[sp,]
  fam.month <- apply( fls, 2, function(x) colSums(allfams[ x == 1,])>0)
  out <- colSums(fam.month)
  return(out)
}))


fl_allfam_arc_plot <- t(sapply(rownames(vegcomm), function(m) {
  presence = colSums(vegcomm[m,])>0
  sp <- names(presence )[presence]
  sp <- intersect(sp, archaeophytes) # select only archaeo
  
  fls <- as.matrix(flower_month[sp,])
  allfams <- allfam_mat[sp,]
  fam.month <- apply( fls, 2, function(x) colSums(allfams[ x == 1,])>0)
  out <- colSums(fam.month)
  return(out)
}))


fl_allfam_exo_plot <- t(sapply(rownames(vegcomm), function(m) {
  presence = colSums(vegcomm[m,])>0
  sp <- names(presence )[presence]
  sp <- intersect(sp, exotics) # select only exotics
  
  fls <- as.matrix(flower_month[sp,])
  allfams <- allfam_mat[sp,]
  fam.month <- apply( fls, 2, function(x) colSums(allfams[ x == 1,])>0)
  out <- colSums(fam.month)
  return(out)
}))


# Weighted by percent cover: 
fl_allfam_plot_wtd <- t(sapply(rownames(vegcomm), function(m) {
  sp <- colnames(vegcomm)[vegcomm[m,]>0]
  cov <- as.numeric(vegcomm[m,sp])
  allfams <- as.matrix(t(allfam_mat[sp,]))
  fls <- as.matrix(flower_month[sp,])
  out <- colSums(allfams %*% (as.numeric(cov)*fls))
  return(out)
}))


fl_allfam_neo_plot_wtd <- t(sapply(rownames(vegcomm), function(m) {
  sp <- colnames(vegcomm)[vegcomm[m,]>0]
  sp <- intersect(sp, neophytes) # select only neophytes
  cov <- as.numeric(vegcomm[m,sp])
  allfams <- as.matrix(t(allfam_mat[sp,]))
  fls <- as.matrix(flower_month[sp,])
  out <- colSums(allfams %*% (as.numeric(cov)*fls))
  return(out)
}))



fl_allfam_arc_plot_wtd <- t(sapply(rownames(vegcomm), function(m) {
  sp <- colnames(vegcomm)[vegcomm[m,]>0]
  sp <- intersect(sp, archaeophytes) # select only archaeophytes
  cov <- as.numeric(vegcomm[m,sp])
  allfams <- as.matrix(t(allfam_mat[sp,]))
  fls <- as.matrix(flower_month[sp,])
  out <- colSums(allfams %*% (as.numeric(cov)*fls))
  return(out)
}))


fl_allfam_nat_plot_wtd <- t(sapply(rownames(vegcomm), function(m) {
  sp <- colnames(vegcomm)[vegcomm[m,]>0]
  sp <- intersect(sp, natives) # select only natives
  cov <- as.numeric(vegcomm[m,sp])
  allfams <- as.matrix(t(allfam_mat[sp,]))
  fls <- as.matrix(flower_month[sp,])
  out <- colSums(allfams %*% (as.numeric(cov)*fls))
  return(out)
}))

fl_allfam_exo_plot_wtd <- t(sapply(rownames(vegcomm), function(m) {
  sp <- colnames(vegcomm)[vegcomm[m,]>0]
  sp <- intersect(sp, exotics) # select only exotics
  cov <- as.numeric(vegcomm[m,sp])
  allfams <- as.matrix(t(allfam_mat[sp,]))
  fls <- as.matrix(flower_month[sp,])
  out <- colSums(allfams %*% (as.numeric(cov)*fls))
  return(out)
}))



# Number of UNIQUE allfam per month by subset of plots ####
## Subsets along gradient: 
z = plot_summary$Seal_levels
table(z)


SRallfam_subset_month <- t(sapply(levels(z), function(m) {
  presence = colSums(vegcomm[plot_summary$Seal_levels == m,])>0
  sp <- names(presence )[presence]
  fls <- as.matrix(flower_month[sp,])
  allfams <- allfam_mat[sp,]
  fam.month <- apply( fls, 2, function(x) colSums(allfams[ x == 1,])>0)
  out <- colSums(fam.month)
  return(out)
}))


SRallfam_subset_month_neo <- t(sapply(levels(z), function(m) {
  presence = colSums(vegcomm[plot_summary$Seal_levels == m,])>0
  sp <- names(presence )[presence]
  sp <- intersect(sp, neophytes) # select only neophytes
  
  fls <- as.matrix(flower_month[sp,])
  allfams <- allfam_mat[sp,]
  fam.month <- apply( fls, 2, function(x) colSums(allfams[ x == 1,])>0)
  out <- colSums(fam.month)
  return(out)
}))

SRallfam_subset_month_nat <- t(sapply(levels(z), function(m) {
  presence = colSums(vegcomm[plot_summary$Seal_levels == m,])>0
  sp <- names(presence )[presence]
  sp <- intersect(sp, natives) # select only natives
  
  fls <- as.matrix(flower_month[sp,])
  allfams <- allfam_mat[sp,]
  fam.month <- apply( fls, 2, function(x) colSums(allfams[ x == 1,])>0)
  out <- colSums(fam.month)
  return(out)
}))


SRallfam_subset_month_arc <- t(sapply(levels(z), function(m) {
  presence = colSums(vegcomm[plot_summary$Seal_levels == m,])>0
  sp <- names(presence )[presence]
  sp <- intersect(sp, archaeophytes) # select only archaeo
  
  fls <- as.matrix(flower_month[sp,])
  allfams <- allfam_mat[sp,]
  fam.month <- apply( fls, 2, function(x) colSums(allfams[ x == 1,])>0)
  out <- colSums(fam.month)
  return(out)
}))




# MOLPHENO : Molecule pheno metrics for all plots: ####
molpheno <-  data.frame(
    beg.fl =  apply(fl_mol_plot  , 1, function(x) min(which(x>0))),
    end.fl =   apply(fl_mol_plot , 1, function(x) max(which(x>0))),
    period.fl =  apply(fl_mol_plot , 1,
                       function(x) max(which(x>0)) - min(which(x>0)) +1),
    peak.fl =apply(fl_mol_plot  , 1, function(x) if (sum(x)!=0) min(which(x == max(x))) else NA),
    cum.fl = apply(fl_mol_plot  , 1, function(x) sum(x)/12), # mean monthly value
    
    beg.fl.5pc =  apply(fl_mol_plot_wtd  , 1, function(x) min(which(x>=5))),
    end.fl.5pc =   apply(fl_mol_plot_wtd , 1, function(x) max(which(x>=5))),
    period.fl.5pc =  apply(fl_mol_plot_wtd , 1,
                           function(x) max(which(x>=5)) - min(which(x>0)) +1),
    peak.fl.wtd =apply(fl_mol_plot_wtd  , 1,
                       function(x) if (sum(x)!=0) min(which(x == max(x))) else NA),
    cum.fl.wtd = apply(fl_mol_plot_wtd  , 1, function(x) sum(x)/12), # mean monthly value
    
    # Neophytes
    beg.fl.neo =  apply(fl_mol_neo_plot, 1, function(x) min(which(x>0))),
    end.fl.neo =   apply(fl_mol_neo_plot , 1, function(x) max(which(x>0))),
    period.fl.neo =  apply(fl_mol_neo_plot , 1,
                       function(x) max(which(x>0)) - min(which(x>0)) +1),
    peak.fl.neo =apply(fl_mol_neo_plot  , 1,
                       function(x) if (sum(x)!=0) min(which(x == max(x))) else NA),
    cum.fl.neo = apply(fl_mol_neo_plot  , 1, function(x) sum(x)/12), # mean monthly value
    beg.fl.neo.5pc =  apply(fl_mol_neo_plot_wtd  , 1,
                            function(x) min(which(x>=5))),
    end.fl.neo.5pc =   apply(fl_mol_neo_plot_wtd , 1,
                             function(x) max(which(x>=5))),
    period.fl.neo.5pc =  apply(fl_mol_neo_plot_wtd , 1,
                           function(x) max(which(x>=5)) - min(which(x>0)) +1),
    peak.fl.neo.wtd =apply(fl_mol_neo_plot_wtd  , 1,
                       function(x) if (sum(x)!=0) min(which(x == max(x))) else NA),
    cum.fl.neo.wtd = apply(fl_mol_neo_plot_wtd  , 1, function(x) sum(x)/12), # mean monthly value
    
    
    # archaeophytes
    beg.fl.arc =  apply(fl_mol_arc_plot, 1, function(x) min(which(x>0))),
    end.fl.arc =   apply(fl_mol_arc_plot , 1, function(x) max(which(x>0))),
    period.fl.arc =  apply(fl_mol_arc_plot , 1,
                           function(x) max(which(x>0)) - min(which(x>0)) +1),
    peak.fl.arc =apply(fl_mol_arc_plot  , 1,
                       function(x) if (sum(x)!=0) min(which(x == max(x))) else NA),
    cum.fl.arc = apply(fl_mol_arc_plot  , 1, function(x) sum(x)/12), # mean monthly value
    beg.fl.arc.5pc =  apply(fl_mol_arc_plot_wtd  , 1,
                            function(x) min(which(x>=5))),
    end.fl.arc.5pc =   apply(fl_mol_arc_plot_wtd , 1,
                             function(x) max(which(x>=5))),
    period.fl.arc.5pc =  apply(fl_mol_arc_plot_wtd , 1,
                               function(x) max(which(x>=5)) - min(which(x>0)) +1),
    peak.fl.arc.wtd =apply(fl_mol_arc_plot_wtd  , 1,
                           function(x) if (sum(x)!=0) min(which(x == max(x))) else NA),
    cum.fl.arc.wtd = apply(fl_mol_arc_plot_wtd  , 1, function(x) sum(x)/12), # mean monthly value
    
    
    # natives
    beg.fl.nat =  apply(fl_mol_nat_plot, 1, function(x) min(which(x>0))),
    end.fl.nat =   apply(fl_mol_nat_plot , 1, function(x) max(which(x>0))),
    period.fl.nat =  apply(fl_mol_nat_plot , 1,
                           function(x) max(which(x>0)) - min(which(x>0)) +1),
    peak.fl.nat =apply(fl_mol_nat_plot  , 1,
                       function(x) if (sum(x)!=0) min(which(x == max(x))) else NA),
    cum.fl.nat = apply(fl_mol_nat_plot  , 1, function(x) sum(x)/12), # mean monthly value
    beg.fl.nat.5pc =  apply(fl_mol_nat_plot_wtd  , 1,
                            function(x) min(which(x>=5))),
    end.fl.nat.5pc =   apply(fl_mol_nat_plot_wtd , 1,
                             function(x) max(which(x>=5))),
    period.fl.nat.5pc =  apply(fl_mol_nat_plot_wtd , 1,
                               function(x) max(which(x>=5)) - min(which(x>0)) +1),
    peak.fl.nat.wtd =apply(fl_mol_nat_plot_wtd  , 1,
                           function(x) if (sum(x)!=0) min(which(x == max(x))) else NA),
    cum.fl.nat.wtd = apply(fl_mol_nat_plot_wtd  , 1, function(x) sum(x)/12), # mean monthly value
    
    # exotics
    beg.fl.exo =  apply(fl_mol_exo_plot, 1, function(x) min(which(x>0))),
    end.fl.exo =   apply(fl_mol_exo_plot , 1, function(x) max(which(x>0))),
    period.fl.exo =  apply(fl_mol_exo_plot , 1,
                           function(x) max(which(x>0)) - min(which(x>0)) +1),
    peak.fl.exo =apply(fl_mol_exo_plot  , 1,
                       function(x) if (sum(x)!=0) min(which(x == max(x))) else NA),
    cum.fl.exo = apply(fl_mol_exo_plot  , 1, function(x) sum(x)/12), # mean monthly value
    beg.fl.exo.5pc =  apply(fl_mol_exo_plot_wtd  , 1,
                            function(x) min(which(x>=5))),
    end.fl.exo.5pc =   apply(fl_mol_exo_plot_wtd , 1,
                             function(x) max(which(x>=5))),
    period.fl.exo.5pc =  apply(fl_mol_exo_plot_wtd , 1,
                               function(x) max(which(x>=5)) - min(which(x>0)) +1),
    peak.fl.exo.wtd =apply(fl_mol_exo_plot_wtd  , 1,
                           function(x) if (sum(x)!=0) min(which(x == max(x))) else NA),
    cum.fl.exo.wtd = apply(fl_mol_exo_plot_wtd  , 1, function(x) sum(x)/12) # mean monthly value
  )   


rownames(molpheno) <- allergen_summary$plot
molpheno[is.infinite(as.matrix(molpheno))] <- NA
molpheno <- as.data.frame(molpheno, as.is = T)


# ALLFAMPHENO : ALLFAM pheno metrics for all plots: ####
allfampheno <-  data.frame(
  beg.fl =  apply(fl_allfam_plot  , 1, function(x) min(which(x>0))),
  end.fl =   apply(fl_allfam_plot , 1, function(x) max(which(x>0))),
  period.fl =  apply(fl_allfam_plot , 1,
                     function(x) max(which(x>0)) - min(which(x>0)) +1),
  peak.fl =apply(fl_allfam_plot  , 1, function(x) if (sum(x)!=0) min(which(x == max(x))) else NA),
  cum.fl = apply(fl_allfam_plot  , 1, function(x) sum(x)/12/12),
  
  beg.fl.5pc =  apply(fl_allfam_plot_wtd  , 1, function(x) min(which(x>=5))),
  end.fl.5pc =   apply(fl_allfam_plot_wtd , 1, function(x) max(which(x>=5))),
  period.fl.5pc =  apply(fl_allfam_plot_wtd , 1,
                         function(x) max(which(x>=5)) - min(which(x>0)) +1),
  peak.fl.wtd =apply(fl_allfam_plot_wtd  , 1,
                     function(x) if (sum(x)!=0) min(which(x == max(x))) else NA),
  cum.fl.wtd = apply(fl_allfam_plot_wtd  , 1, function(x) sum(x)/12),
  
  # Neophytes
  beg.fl.neo =  apply(fl_allfam_neo_plot, 1, function(x) min(which(x>0))),
  end.fl.neo =   apply(fl_allfam_neo_plot , 1, function(x) max(which(x>0))),
  period.fl.neo =  apply(fl_allfam_neo_plot , 1,
                         function(x) max(which(x>0)) - min(which(x>0)) +1),
  peak.fl.neo =apply(fl_allfam_neo_plot  , 1,
                     function(x) if (sum(x)!=0) min(which(x == max(x))) else NA),
  cum.fl.neo = apply(fl_allfam_neo_plot  , 1, function(x) sum(x)/12),
  beg.fl.neo.5pc =  apply(fl_allfam_neo_plot_wtd  , 1,
                          function(x) min(which(x>=5))),
  end.fl.neo.5pc =   apply(fl_allfam_neo_plot_wtd , 1,
                           function(x) max(which(x>=5))),
  period.fl.neo.5pc =  apply(fl_allfam_neo_plot_wtd , 1,
                             function(x) max(which(x>=5)) - min(which(x>0)) +1),
  peak.fl.neo.wtd =apply(fl_allfam_neo_plot_wtd  , 1,
                         function(x) if (sum(x)!=0) min(which(x == max(x))) else NA),
  cum.fl.neo.wtd = apply(fl_allfam_neo_plot_wtd  , 1, function(x) sum(x)/12),
  
  
  # archaeophytes
  beg.fl.arc =  apply(fl_allfam_arc_plot, 1, function(x) min(which(x>0))),
  end.fl.arc =   apply(fl_allfam_arc_plot , 1, function(x) max(which(x>0))),
  period.fl.arc =  apply(fl_allfam_arc_plot , 1,
                         function(x) max(which(x>0)) - min(which(x>0)) +1),
  peak.fl.arc =apply(fl_allfam_arc_plot  , 1,
                     function(x) if (sum(x)!=0) min(which(x == max(x))) else NA),
  cum.fl.arc = apply(fl_allfam_arc_plot  , 1, function(x) sum(x)/12),
  beg.fl.arc.5pc =  apply(fl_allfam_arc_plot_wtd  , 1,
                          function(x) min(which(x>=5))),
  end.fl.arc.5pc =   apply(fl_allfam_arc_plot_wtd , 1,
                           function(x) max(which(x>=5))),
  period.fl.arc.5pc =  apply(fl_allfam_arc_plot_wtd , 1,
                             function(x) max(which(x>=5)) - min(which(x>0)) +1),
  peak.fl.arc.wtd =apply(fl_allfam_arc_plot_wtd  , 1,
                         function(x) if (sum(x)!=0) min(which(x == max(x))) else NA),
  cum.fl.arc.wtd = apply(fl_allfam_arc_plot_wtd  , 1, function(x) sum(x)/12),
  
  
  # natives
  beg.fl.nat =  apply(fl_allfam_nat_plot, 1, function(x) min(which(x>0))),
  end.fl.nat =   apply(fl_allfam_nat_plot , 1, function(x) max(which(x>0))),
  period.fl.nat =  apply(fl_allfam_nat_plot , 1,
                         function(x) max(which(x>0)) - min(which(x>0)) +1),
  peak.fl.nat =apply(fl_allfam_nat_plot  , 1,
                     function(x) if (sum(x)!=0) min(which(x == max(x))) else NA),
  cum.fl.nat = apply(fl_allfam_nat_plot  , 1, function(x) sum(x)/12),
  beg.fl.nat.5pc =  apply(fl_allfam_nat_plot_wtd  , 1,
                          function(x) min(which(x>=5))),
  end.fl.nat.5pc =   apply(fl_allfam_nat_plot_wtd , 1,
                           function(x) max(which(x>=5))),
  period.fl.nat.5pc =  apply(fl_allfam_nat_plot_wtd , 1,
                             function(x) max(which(x>=5)) - min(which(x>0)) +1),
  peak.fl.nat.wtd =apply(fl_allfam_nat_plot_wtd  , 1,
                         function(x) if (sum(x)!=0) min(which(x == max(x))) else NA),
  cum.fl.nat.wtd = apply(fl_allfam_nat_plot_wtd  , 1, function(x) sum(x)/12),
  
  # exotics
  beg.fl.exo =  apply(fl_allfam_exo_plot, 1, function(x) min(which(x>0))),
  end.fl.exo =   apply(fl_allfam_exo_plot , 1, function(x) max(which(x>0))),
  period.fl.exo =  apply(fl_allfam_exo_plot , 1,
                         function(x) max(which(x>0)) - min(which(x>0)) +1),
  peak.fl.exo =apply(fl_allfam_exo_plot  , 1,
                     function(x) if (sum(x)!=0) min(which(x == max(x))) else NA),
  cum.fl.exo = apply(fl_allfam_exo_plot  , 1, function(x) sum(x)/12),
  beg.fl.exo.5pc =  apply(fl_allfam_exo_plot_wtd  , 1,
                          function(x) min(which(x>=5))),
  end.fl.exo.5pc =   apply(fl_allfam_exo_plot_wtd , 1,
                           function(x) max(which(x>=5))),
  period.fl.exo.5pc =  apply(fl_allfam_exo_plot_wtd , 1,
                             function(x) max(which(x>=5)) - min(which(x>0)) +1),
  peak.fl.exo.wtd =apply(fl_allfam_exo_plot_wtd  , 1,
                         function(x) if (sum(x)!=0) min(which(x == max(x))) else NA),
  cum.fl.exo.wtd = apply(fl_allfam_exo_plot_wtd  , 1, function(x) sum(x)/12)
)   


rownames(allfampheno) <- allergen_summary$plot
allfampheno[is.infinite(as.matrix(allfampheno))] <- NA
allfampheno <- as.data.frame(allfampheno, as.is = T)




# COMPHENO : Allergenic species pheno metrics for all plots  ####
compheno <-  data.frame(
  beg.fl =  apply(fl_all_plot  , 1, function(x) min(which(x>0))),
  end.fl =   apply(fl_all_plot , 1, function(x) max(which(x>0))),
  period.fl =  apply(fl_all_plot , 1,
                     function(x) max(which(x>0)) - min(which(x>0)) +1),
  peak.fl =apply(fl_all_plot  , 1, function(x) min(which(x == max(x)))) ,
  cum.fl = apply(fl_all_plot  , 1, function(x) sum(x)/12),
  
  beg.fl.5pc =  apply(fl_all_plot_wtd  , 1, function(x) min(which(x>=5))),
  end.fl.5pc =   apply(fl_all_plot_wtd , 1, function(x) max(which(x>=5))),
  period.fl.5pc =  apply(fl_all_plot_wtd , 1,
                         function(x) max(which(x>=5)) - min(which(x>0)) +1),
  peak.fl.wtd =apply(fl_all_plot_wtd  , 1,
                     function(x) if (sum(x)!=0) min(which(x == max(x))) else NA),
  cum.fl.wtd = apply(fl_all_plot_wtd  , 1, function(x) sum(x)/12),
  
  # Neophytes
  beg.fl.neo =  apply(fl_all_neo_plot, 1, function(x) min(which(x>0))),
  end.fl.neo =   apply(fl_all_neo_plot , 1, function(x) max(which(x>0))),
  period.fl.neo =  apply(fl_all_neo_plot , 1,
                         function(x) max(which(x>0)) - min(which(x>0)) +1),
  peak.fl.neo =apply(
    fl_all_neo_plot,1,
    function(x) if (sum(x)!=0) min(which(x == max(x))) else NA),
  
  cum.fl.neo = apply(fl_all_neo_plot  , 1, function(x) sum(x)/12),
  beg.fl.neo.5pc =  apply(fl_all_neo_plot_wtd  , 1,
                          function(x) min(which(x>=5))),
  end.fl.neo.5pc =   apply(fl_all_neo_plot_wtd , 1,
                           function(x) max(which(x>=5))),
  period.fl.neo.5pc =  apply(fl_all_neo_plot_wtd , 1,
                             function(x) max(which(x>=5)) - min(which(x>0)) +1),
  peak.fl.neo.wtd =apply(fl_all_neo_plot_wtd  , 1,
                         function(x) if (sum(x)!=0) min(which(x == max(x))) else NA),
  cum.fl.neo.wtd = apply(fl_all_neo_plot_wtd  , 1, function(x) sum(x)/12),
  
  
  # archaeophytes
  beg.fl.arc =  apply(fl_all_arc_plot, 1, function(x) min(which(x>0))),
  end.fl.arc =   apply(fl_all_arc_plot , 1, function(x) max(which(x>0))),
  period.fl.arc =  apply(fl_all_arc_plot , 1,
                         function(x) max(which(x>0)) - min(which(x>0)) +1),
  peak.fl.arc =apply(fl_all_arc_plot  , 1,
                     function(x) if (sum(x)!=0) min(which(x == max(x))) else NA) ,
  cum.fl.arc = apply(fl_all_arc_plot  , 1, function(x) sum(x)/12),
  beg.fl.arc.5pc =  apply(fl_all_arc_plot_wtd  , 1,
                          function(x) min(which(x>=5))),
  end.fl.arc.5pc =   apply(fl_all_arc_plot_wtd , 1,
                           function(x) max(which(x>=5))),
  period.fl.arc.5pc =  apply(fl_all_arc_plot_wtd , 1,
                             function(x) max(which(x>=5)) - min(which(x>0)) +1),
  peak.fl.arc.wtd =apply(fl_all_arc_plot_wtd,1,
                         function(x) if (sum(x)!=0) min(which(x == max(x))) else NA), 
  cum.fl.arc.wtd = apply(fl_all_arc_plot_wtd  , 1, function(x) sum(x)/12),
  
  
  # natives
  beg.fl.nat =  apply(fl_all_nat_plot, 1, function(x) min(which(x>0))),
  end.fl.nat =   apply(fl_all_nat_plot , 1, function(x) max(which(x>0))),
  period.fl.nat =  apply(fl_all_nat_plot , 1,
                         function(x) max(which(x>0)) - min(which(x>0)) +1),
  peak.fl.nat =apply(fl_all_nat_plot  , 1,
                     function(x) if (sum(x)!=0) min(which(x == max(x))) else NA),
  cum.fl.nat = apply(fl_all_nat_plot, 1, function(x) sum(x)/12),
  beg.fl.nat.5pc =  apply(fl_all_nat_plot_wtd  , 1,
                          function(x) min(which(x>=5))),
  end.fl.nat.5pc =   apply(fl_all_nat_plot_wtd , 1,
                           function(x) max(which(x>=5))),
  period.fl.nat.5pc =  apply(fl_all_nat_plot_wtd , 1,
                             function(x) max(which(x>=5)) - min(which(x>0)) +1),
  peak.fl.nat.wtd =apply(fl_all_nat_plot_wtd  , 1,
                         function(x) if (sum(x)!=0) min(which(x == max(x))) else NA),
  cum.fl.nat.wtd = apply(fl_all_nat_plot_wtd  , 1, function(x) sum(x)/12),
  
  
  # exotics
  beg.fl.exo =  apply(fl_all_exo_plot, 1, function(x) min(which(x>0))),
  end.fl.exo =   apply(fl_all_exo_plot , 1, function(x) max(which(x>0))),
  period.fl.exo =  apply(fl_all_exo_plot , 1,
                         function(x) max(which(x>0)) - min(which(x>0)) +1),
  peak.fl.exo =apply(fl_all_exo_plot  , 1,
                     function(x) if (sum(x)!=0) min(which(x == max(x))) else NA),
  cum.fl.exo = apply(fl_all_exo_plot  , 1, function(x) sum(x)/12),
  beg.fl.exo.5pc =  apply(fl_all_exo_plot_wtd  , 1,
                          function(x) min(which(x>=5))),
  end.fl.exo.5pc =   apply(fl_all_exo_plot_wtd , 1,
                           function(x) max(which(x>=5))),
  period.fl.exo.5pc =  apply(fl_all_exo_plot_wtd , 1,
                             function(x) max(which(x>=5)) - min(which(x>0)) +1),
  peak.fl.exo.wtd =apply(fl_all_exo_plot_wtd  , 1,
                         function(x) if (sum(x)!=0) min(which(x == max(x))) else NA),
  cum.fl.exo.wtd = apply(fl_all_exo_plot_wtd  , 1, function(x) sum(x)/12)
)   


rownames(compheno) <- allergen_summary$plot
compheno[is.infinite(as.matrix(compheno))] <- NA
compheno <- as.data.frame(compheno, as.is = T)


## Calculate community mean phenological traits - allergenic SPECIES
tmp <- ceiling(vegcomm[,allergenics]/100)
tmp[tmp==0]<- NA
tmp2 <- species_allergen[allergenics,]

compheno$mean.fl.beg <- rowSums(tmp* tmp2[, "fl.beg"],na.rm = T)/
  rowSums(tmp,na.rm = T)
compheno$sd.fl.beg <- apply(tmp* tmp2[, "fl.beg"],1, sd, na.rm = T)

compheno$mean.fl.end <- rowSums(tmp* tmp2[, "fl.end"], na.rm = T)/
  rowSums(tmp,na.rm = T)
compheno$sd.fl.end <- apply(tmp* tmp2[, "fl.end"],1, sd, na.rm = T)

compheno$mean.fl.period <- rowSums(tmp* tmp2[, "fl.period"], na.rm = T)/
  rowSums(tmp,na.rm = T)
compheno$sd.fl.period <- apply(tmp* tmp2[, "fl.period"],1, sd, na.rm = T)

# intersect(allergenics, neophytes) only
compheno$mean.neo.fl.beg <- apply((tmp* tmp2[, "fl.beg"])[,intersect(allergenics, neophytes)],1, mean, na.rm = T)
compheno$sd.neo.fl.beg <- apply((tmp* tmp2[, "fl.beg"])[ , intersect(allergenics, neophytes)],1, sd, na.rm = T)

compheno$mean.neo.fl.end <- apply((tmp* tmp2[, "fl.end"])[ , intersect(allergenics, neophytes)],1, mean, na.rm = T)
compheno$sd.neo.fl.end <- apply((tmp* tmp2[, "fl.end"])[ , intersect(allergenics, neophytes)],1, sd, na.rm = T)

compheno$mean.neo.fl.period <- apply((tmp* tmp2[, "fl.period"])[ , intersect(allergenics, neophytes)],1, mean, na.rm = T)
compheno$sd.neo.fl.period <- apply((tmp* tmp2[, "fl.period"])[ , intersect(allergenics, neophytes)],1, sd, na.rm = T)

# intersect(allergenics, archaeophytes) only
compheno$mean.arc.fl.beg <- apply((tmp* tmp2[, "fl.beg"])[ , intersect(allergenics, archaeophytes)],1, mean, na.rm = T)
compheno$sd.arc.fl.beg <- apply((tmp* tmp2[, "fl.beg"])[ , intersect(allergenics, archaeophytes)],1, sd, na.rm = T)

compheno$mean.arc.fl.end <- apply((tmp* tmp2[, "fl.end"])[ , intersect(allergenics, archaeophytes)],1, mean, na.rm = T)
compheno$sd.arc.fl.end <- apply((tmp* tmp2[, "fl.end"])[ , intersect(allergenics, archaeophytes)],1, sd, na.rm = T)

compheno$mean.arc.fl.period <- apply((tmp* tmp2[, "fl.period"])[ , intersect(allergenics, archaeophytes)],1, mean, na.rm = T)
compheno$sd.arc.fl.period <- apply((tmp* tmp2[, "fl.period"])[ , intersect(allergenics, archaeophytes)],1, sd, na.rm = T)

# intersect(allergenics, natives) only
compheno$mean.nat.fl.beg <- apply((tmp* tmp2[, "fl.beg"])[ , intersect(allergenics, natives)],1, mean, na.rm = T)
compheno$sd.nat.fl.beg <- apply((tmp* tmp2[, "fl.beg"])[ , intersect(allergenics, natives)],1, sd, na.rm = T)

compheno$mean.nat.fl.end <- apply((tmp* tmp2[, "fl.end"])[ , intersect(allergenics, natives)],1, mean, na.rm = T)
compheno$sd.nat.fl.end <- apply((tmp* tmp2[, "fl.end"])[ , intersect(allergenics, natives)],1, sd, na.rm = T)

compheno$mean.nat.fl.period <- apply((tmp* tmp2[, "fl.period"])[ , intersect(allergenics, natives)],1, mean, na.rm = T)
compheno$sd.nat.fl.period <- apply((tmp* tmp2[, "fl.period"])[ , intersect(allergenics, natives)],1, sd, na.rm = T)


## CALCULATE AN ABUNDANCE WEIGHTED VERSION  - ALL SPECIES
tmp <- vegcomm[,allergenics]
tmp2 <- tmp2[allergenics,]

compheno$CWM.fl.beg <- rowSums((tmp/rowSums(tmp)) * tmp2[, "fl.beg"], na.rm = T)
compheno$CWM.fl.end <- rowSums((tmp/rowSums(tmp)) * tmp2[, "fl.end"], na.rm = T)
compheno$CWM.fl.period <- rowSums((tmp/rowSums(tmp)) * tmp2[, "fl.period"], na.rm = T)

# intersect(allergenics, natives) only
compheno$CWM.nat.fl.beg <- rowSums((tmp* tmp2[, "fl.beg"])[,intersect(allergenics, natives)], na.rm = T)/ rowSums(tmp[,intersect(allergenics, natives)])
compheno$CWM.nat.fl.end <-  rowSums((tmp* tmp2[, "fl.end"])[,intersect(allergenics, natives)], na.rm = T)/ rowSums(tmp[,intersect(allergenics, natives)])
compheno$CWM.nat.fl.period <-  rowSums((tmp* tmp2[, "fl.period"])[,intersect(allergenics, natives)], na.rm = T)/ rowSums(tmp[,intersect(allergenics, natives)])

# intersect(allergenics, archaeophytes) only
compheno$CWM.arc.fl.beg <- rowSums((tmp* tmp2[, "fl.beg"])[,intersect(allergenics, archaeophytes)], na.rm = T)/ rowSums(tmp[,intersect(allergenics, archaeophytes)])
compheno$CWM.arc.fl.end <-  rowSums((tmp* tmp2[, "fl.end"])[,intersect(allergenics, archaeophytes)], na.rm = T)/ rowSums(tmp[,intersect(allergenics, archaeophytes)])
compheno$CWM.arc.fl.period <-  rowSums((tmp* tmp2[, "fl.period"])[,intersect(allergenics, archaeophytes)], na.rm = T)/ rowSums(tmp[,intersect(allergenics, archaeophytes)])

# intersect(allergenics, neophytes) only
compheno$CWM.neo.fl.beg <- rowSums((tmp* tmp2[, "fl.beg"])[,intersect(allergenics, neophytes)], na.rm = T)/ rowSums(tmp[,intersect(allergenics, neophytes)])
compheno$CWM.neo.fl.end <-  rowSums((tmp* tmp2[, "fl.end"])[,intersect(allergenics, neophytes)], na.rm = T)/ rowSums(tmp[,intersect(allergenics, neophytes)])
compheno$CWM.neo.fl.period <-  rowSums((tmp* tmp2[, "fl.period"])[,intersect(allergenics, neophytes)], na.rm = T)/ rowSums(tmp[,intersect(allergenics, neophytes)])




# COMPHENO_total : Species pheno metrics for all plots - including non allergenics ####
compheno_total<-  data.frame(
  beg.fl =  apply(fl_plot  , 1, function(x) min(which(x>0))),
  end.fl =   apply(fl_plot , 1, function(x) max(which(x>0))),
  period.fl =  apply(fl_plot , 1,
                     function(x) max(which(x>0)) - min(which(x>0)) +1),
  peak.fl =apply(fl_plot  , 1, 
                 function(x) if (sum(x)!=0) min(which(x == max(x))) else NA),
  cum.fl = apply(fl_plot  , 1, function(x) sum(x)/12),
  
  beg.fl.5pc =  apply(fl_plot_wtd  , 1, function(x) min(which(x>=5))),
  end.fl.5pc =   apply(fl_plot_wtd , 1, function(x) max(which(x>=5))),
  period.fl.5pc =  apply(fl_plot_wtd , 1,
                         function(x) max(which(x>=5)) - min(which(x>0)) +1),
  peak.fl.wtd =apply(fl_plot_wtd  , 1,
                     function(x) if (sum(x)!=0) min(which(x == max(x))) else NA),
  cum.fl.wtd = apply(fl_plot_wtd  , 1, function(x) sum(x)/12),
  
  # Neophytes
  beg.fl.neo =  apply(fl_neo_plot, 1, function(x) min(which(x>0))),
  end.fl.neo =   apply(fl_neo_plot , 1, function(x) max(which(x>0))),
  period.fl.neo =  apply(fl_neo_plot , 1,
                         function(x) max(which(x>0)) - min(which(x>0)) +1),
  peak.fl.neo =apply(fl_neo_plot  , 1,
                     function(x) if (sum(x)!=0) min(which(x == max(x))) else NA),
  cum.fl.neo = apply(fl_neo_plot  , 1, function(x) sum(x)/12),
  beg.fl.neo.5pc =  apply(fl_neo_plot_wtd  , 1,
                          function(x) min(which(x>=5))),
  end.fl.neo.5pc =   apply(fl_neo_plot_wtd , 1,
                           function(x) max(which(x>=5))),
  period.fl.neo.5pc =  apply(fl_neo_plot_wtd , 1,
                             function(x) max(which(x>=5)) - min(which(x>0)) +1),
  peak.fl.neo.wtd =apply(fl_neo_plot_wtd  , 1,
                         function(x) if (sum(x)!=0) min(which(x == max(x))) else NA),
  cum.fl.neo.wtd = apply(fl_neo_plot_wtd  , 1, function(x) sum(x)/12),
  
  
  # archaeophytes
  beg.fl.arc =  apply(fl_arc_plot, 1, function(x) min(which(x>0))),
  end.fl.arc =   apply(fl_arc_plot , 1, function(x) max(which(x>0))),
  period.fl.arc =  apply(fl_arc_plot , 1,
                         function(x) max(which(x>0)) - min(which(x>0)) +1),
  peak.fl.arc =apply(fl_arc_plot  , 1,
                     function(x) if (sum(x)!=0) min(which(x == max(x))) else NA),
  cum.fl.arc = apply(fl_arc_plot  , 1, function(x) sum(x)/12),
  beg.fl.arc.5pc =  apply(fl_arc_plot_wtd  , 1,
                          function(x) min(which(x>=5))),
  end.fl.arc.5pc =   apply(fl_arc_plot_wtd , 1,
                           function(x) max(which(x>=5))),
  period.fl.arc.5pc =  apply(fl_arc_plot_wtd , 1,
                             function(x) max(which(x>=5)) - min(which(x>0)) +1),
  peak.fl.arc.wtd =apply(fl_arc_plot_wtd  , 1,
                         function(x) if (sum(x)!=0) min(which(x == max(x))) else NA),
  cum.fl.arc.wtd = apply(fl_arc_plot_wtd  , 1, function(x) sum(x)/12),
  
  
  # natives
  beg.fl.nat =  apply(fl_nat_plot, 1, function(x) min(which(x>0))),
  end.fl.nat =   apply(fl_nat_plot , 1, function(x) max(which(x>0))),
  period.fl.nat =  apply(fl_nat_plot , 1,
                         function(x) max(which(x>0)) - min(which(x>0)) +1),
  peak.fl.nat =apply(fl_nat_plot  , 1,
                     function(x) if (sum(x)!=0) min(which(x == max(x))) else NA),
  cum.fl.nat = apply(fl_nat_plot  , 1, function(x) sum(x)/12),
  beg.fl.nat.5pc =  apply(fl_nat_plot_wtd  , 1,
                          function(x) min(which(x>=5))),
  end.fl.nat.5pc =   apply(fl_nat_plot_wtd , 1,
                           function(x) max(which(x>=5))),
  period.fl.nat.5pc =  apply(fl_nat_plot_wtd , 1,
                             function(x) max(which(x>=5)) - min(which(x>0)) +1),
  peak.fl.nat.wtd =apply(fl_nat_plot_wtd  , 1,
                         function(x) if (sum(x)!=0) min(which(x == max(x))) else NA),
  cum.fl.nat.wtd = apply(fl_nat_plot_wtd  , 1, function(x) sum(x)/12)
)   


rownames(compheno_total) <- allergen_summary$plot
compheno_total[is.infinite(as.matrix(compheno_total))] <- NA
compheno_total <- as.data.frame(compheno_total, as.is = T)

## Calculate community mean phenological traits - ALL SPECIES
compheno_total$mean.fl.beg <- rowSums((ceiling(vegcomm/100)) * species_allergen[, "fl.beg"],
                                na.rm = T)/rowSums(ceiling(vegcomm/100))
compheno_total$sd.fl.beg <- apply((ceiling(vegcomm/100)) * species_allergen[, "fl.beg"],1, sd, na.rm = T)

compheno_total$mean.fl.end <- rowSums((ceiling(vegcomm/100)) * species_allergen[, "fl.end"], na.rm = T)/rowSums(ceiling(vegcomm/100))
compheno_total$sd.fl.end <- apply((ceiling(vegcomm/100)) * species_allergen[, "fl.end"],1, sd, na.rm = T)

compheno_total$mean.fl.period <- rowSums((ceiling(vegcomm/100)) * species_allergen[, "fl.period"], na.rm = T)/rowSums(ceiling(vegcomm/100))
compheno_total$sd.fl.period <- apply((ceiling(vegcomm/100)) * species_allergen[, "fl.period"],1, sd, na.rm = T)

# Neophytes only
compheno_total$mean.neo.fl.beg <- apply(((ceiling(vegcomm/100)) * species_allergen[, "fl.beg"])[,neophytes],1, mean, na.rm = T)
compheno_total$sd.neo.fl.beg <- apply(((ceiling(vegcomm/100)) * species_allergen[, "fl.beg"])[ , neophytes],1, sd, na.rm = T)

compheno_total$mean.neo.fl.end <- apply(((ceiling(vegcomm/100)) * species_allergen[, "fl.end"])[ , neophytes],1, mean, na.rm = T)
compheno_total$sd.neo.fl.end <- apply(((ceiling(vegcomm/100)) * species_allergen[, "fl.end"])[ , neophytes],1, sd, na.rm = T)

compheno_total$mean.neo.fl.period <- apply(((ceiling(vegcomm/100)) * species_allergen[, "fl.period"])[ , neophytes],1, mean, na.rm = T)
compheno_total$sd.neo.fl.period <- apply(((ceiling(vegcomm/100)) * species_allergen[, "fl.period"])[ , neophytes],1, sd, na.rm = T)

# archaeophytes only
compheno_total$mean.arc.fl.beg <- apply(((ceiling(vegcomm/100)) * species_allergen[, "fl.beg"])[ , archaeophytes],1, mean, na.rm = T)
compheno_total$sd.arc.fl.beg <- apply(((ceiling(vegcomm/100)) * species_allergen[, "fl.beg"])[ , archaeophytes],1, sd, na.rm = T)

compheno_total$mean.arc.fl.end <- apply(((ceiling(vegcomm/100)) * species_allergen[, "fl.end"])[ , archaeophytes],1, mean, na.rm = T)
compheno_total$sd.arc.fl.end <- apply(((ceiling(vegcomm/100)) * species_allergen[, "fl.end"])[ , archaeophytes],1, sd, na.rm = T)

compheno_total$mean.arc.fl.period <- apply(((ceiling(vegcomm/100)) * species_allergen[, "fl.period"])[ , archaeophytes],1, mean, na.rm = T)
compheno_total$sd.arc.fl.period <- apply(((ceiling(vegcomm/100)) * species_allergen[, "fl.period"])[ , archaeophytes],1, sd, na.rm = T)

# natives only
compheno_total$mean.nat.fl.beg <- apply(((ceiling(vegcomm/100)) * species_allergen[, "fl.beg"])[ , natives],1, mean, na.rm = T)
compheno_total$sd.nat.fl.beg <- apply(((ceiling(vegcomm/100)) * species_allergen[, "fl.beg"])[ , natives],1, sd, na.rm = T)

compheno_total$mean.nat.fl.end <- apply(((ceiling(vegcomm/100)) * species_allergen[, "fl.end"])[ , natives],1, mean, na.rm = T)
compheno_total$sd.nat.fl.end <- apply(((ceiling(vegcomm/100)) * species_allergen[, "fl.end"])[ , natives],1, sd, na.rm = T)

compheno_total$mean.nat.fl.period <- apply(((ceiling(vegcomm/100)) * species_allergen[, "fl.period"])[ , natives],1, mean, na.rm = T)
compheno_total$sd.nat.fl.period <- apply(((ceiling(vegcomm/100)) * species_allergen[, "fl.period"])[ , natives],1, sd, na.rm = T)


## CALCULATE AN ABUNDANCE WEIGHTED VERSION  - ALL SPECIES

compheno_total$CWM.fl.beg <- apply((vegcomm/rowSums(vegcomm)) * species_allergen[, "fl.beg"],1, mean, na.rm = T)
compheno_total$CWM.fl.end <- rowSums((vegcomm/rowSums(vegcomm)) * species_allergen[, "fl.end"], na.rm = T)
compheno_total$CWM.fl.period <- rowSums((vegcomm/rowSums(vegcomm)) * species_allergen[, "fl.period"], na.rm = T)

# Natives only
compheno_total$CWM.nat.fl.beg <- rowSums((vegcomm* species_allergen[, "fl.beg"])[,natives], na.rm = T)/ rowSums(vegcomm[,natives])
compheno_total$CWM.nat.fl.end <-  rowSums((vegcomm* species_allergen[, "fl.end"])[,natives], na.rm = T)/ rowSums(vegcomm[,natives])
compheno_total$CWM.nat.fl.period <-  rowSums((vegcomm* species_allergen[, "fl.period"])[,natives], na.rm = T)/ rowSums(vegcomm[,natives])

# archaeophytes only
compheno_total$CWM.arc.fl.beg <- rowSums((vegcomm* species_allergen[, "fl.beg"])[,archaeophytes], na.rm = T)/ rowSums(vegcomm[,archaeophytes])
compheno_total$CWM.arc.fl.end <-  rowSums((vegcomm* species_allergen[, "fl.end"])[,archaeophytes], na.rm = T)/ rowSums(vegcomm[,archaeophytes])
compheno_total$CWM.arc.fl.period <-  rowSums((vegcomm* species_allergen[, "fl.period"])[,archaeophytes], na.rm = T)/ rowSums(vegcomm[,archaeophytes])

# neophytes only
compheno_total$CWM.neo.fl.beg <- rowSums((vegcomm* species_allergen[, "fl.beg"])[,neophytes], na.rm = T)/ rowSums(vegcomm[,neophytes])
compheno_total$CWM.neo.fl.end <-  rowSums((vegcomm* species_allergen[, "fl.end"])[,neophytes], na.rm = T)/ rowSums(vegcomm[,neophytes])
compheno_total$CWM.neo.fl.period <-  rowSums((vegcomm* species_allergen[, "fl.period"])[,neophytes], na.rm = T)/ rowSums(vegcomm[,neophytes])

