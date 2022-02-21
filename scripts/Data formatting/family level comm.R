## Taxonomic family level community analyses

## Groups of species by status:
all.nat<- intersect(natives, allergenics)
all.neo <-intersect(neophytes, allergenics)
all.arc <- intersect(archaeophytes, allergenics)
all.exo <- intersect(exotics, allergenics)

vegcomm <- vegcomm[rownames(plot_summary),]

## Modified community matrices
vegcomm.fam <- t(apply(vegcomm[,species_allergen$Species], 1, FUN = function(x) tapply(x, species_allergen$family, FUN = sum)))

vegcomm.fam.neo <- t(apply(vegcomm[, neophytes], 1, FUN = function(x) tapply(x, species_allergen[neophytes,]$family, FUN = sum)))

vegcomm.fam.nat <- t(apply(vegcomm[, natives], 1, FUN = function(x) tapply(x, species_allergen[natives,]$family, FUN = sum)))

vegcomm.fam.arc <- t(apply(vegcomm[, archaeophytes], 1, FUN = function(x) tapply(x, species_allergen[archaeophytes,]$family, FUN = sum)))



## allergenic neophytes:
vegcomm.fam.all <- t(apply(vegcomm[, allergenics], 1, FUN = function(x) tapply(x, species_allergen[allergenics,]$family, FUN = sum)))
vegcomm.fam.all.neo <- t(apply(vegcomm[, all.neo], 1, FUN = function(x) tapply(x, species_allergen[ all.neo,]$family, FUN = sum)))
vegcomm.fam.all.nat <- t(apply(vegcomm[,  all.nat], 1, FUN = function(x) tapply(x, species_allergen[ all.nat,]$family, FUN = sum)))
vegcomm.fam.all.arc <- t(apply(vegcomm[,  all.arc], 1, FUN = function(x) tapply(x, species_allergen[ all.arc,]$family, FUN = sum)))


