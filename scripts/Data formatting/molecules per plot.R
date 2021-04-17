# Matrix of allergenic molecules per plot

# Presence of allergenic molecules per plot ###
mol_plot  <- data.frame(matrix(NA,
                              nrow= nrow(vegcomm),
                              ncol = ncol(mol_mat),
                              dimnames = list(rownames(vegcomm),
                                              names(mol_mat)
                                              )
                              )
                       )

for (i in rownames(mol_plot)) {
  mol_plot[i,] <- colSums(mol_mat[colnames(vegcomm)[which(vegcomm[i,]>0)],])
}


# STRICT : Presence of allergenic molecules per plot ###
mol_plot_strict  <- data.frame(matrix(NA,
                               nrow= nrow(vegcomm),
                               ncol = ncol(mol_mat_strict),
                               dimnames = list(rownames(vegcomm),
                                               names(mol_mat_strict)
                               )
)
)

for (i in rownames(mol_plot)) {
  mol_plot_strict[i,] <- colSums(mol_mat_strict[colnames(vegcomm)[which(vegcomm[i,]>0)],])
}

# Presence of allergenic molecules at genus level ###
mol_plot_genus<- data.frame(matrix(NA,
                                       nrow= nrow(vegcomm),
                                       ncol = ncol(mol_mat_gen),
                                       dimnames = list(rownames(vegcomm),
                                                       colnames(mol_mat_gen)
                                       )
)
)

for (i in rownames(mol_plot_genus)) {
  spp <- colnames(vegcomm)[which(vegcomm[i,]>0)]
  genera <- species_allergen$Genus[species_allergen$Species %in% spp]
  mol_plot_genus[i,] <- colSums(mol_mat_gen[genera,])
}


# Cover of different allergenic molecules per plot ####
mol_plot_wtd<- data.frame(matrix(NA, nrow= nrow(vegcomm), ncol = ncol(mol_mat),
                                 dimnames = list(rownames(vegcomm), names(mol_mat))))

for (i in rownames(mol_plot)) {
  sp = colnames(vegcomm)[which(vegcomm[i,]>0)]
  pi = vegcomm[i,sp]
  mol_plot_wtd[i,] <- colSums(mol_mat[sp,] * t(pi))
}

# Neophyte molecules only ####
mol_plot_neo <-data.frame(matrix(NA, nrow= nrow(vegcomm), ncol = ncol(mol_mat),
                                 dimnames = list(rownames(vegcomm), names(mol_mat))))

for (i in rownames(mol_plot_neo)) {
  mol_plot_neo[i,] <- colSums(mol_mat[neophytes[which(vegcomm[i,neophytes]>0)],])
}


mol_plot_neo_wtd <-data.frame(matrix(NA, nrow= nrow(vegcomm), ncol = ncol(mol_mat),
                                 dimnames = list(rownames(vegcomm), names(mol_mat))))

for (i in rownames(mol_plot_neo_wtd)) {
  sp = neophytes[which(vegcomm[i,neophytes]>0)]
  pi = vegcomm[i,sp]
  mol_plot_neo_wtd[i,] <- colSums(mol_mat[sp,] * t(pi))
}

# Native molecules only ###
mol_plot_nat <-data.frame(matrix(NA, nrow= nrow(vegcomm), ncol = ncol(mol_mat),
                                 dimnames = list(rownames(vegcomm), names(mol_mat))))

for (i in rownames(mol_plot_nat)) {
  mol_plot_nat[i,] <- colSums(mol_mat[natives[which(vegcomm[i,natives]>0)],])
}


mol_plot_nat_wtd <-data.frame(matrix(NA, nrow= nrow(vegcomm), ncol = ncol(mol_mat),
                                     dimnames = list(rownames(vegcomm), names(mol_mat))))

for (i in rownames(mol_plot_nat_wtd)) {
  sp = natives[which(vegcomm[i,natives]>0)]
  pi = vegcomm[i,sp]
  mol_plot_nat_wtd[i,] <- colSums(mol_mat[sp,] * t(pi))
}



# Exotic molecules only ####
mol_plot_exo <-data.frame(matrix(NA, nrow= nrow(vegcomm), ncol = ncol(mol_mat),
                                 dimnames = list(rownames(vegcomm), names(mol_mat))))

for (i in rownames(mol_plot_neo)) {
  mol_plot_exo[i,] <- colSums(mol_mat[exotics[which(vegcomm[i,exotics]>0)],])
}


mol_plot_exo_wtd <-data.frame(matrix(NA, nrow= nrow(vegcomm), ncol = ncol(mol_mat),
                                     dimnames = list(rownames(vegcomm), names(mol_mat))))

for (i in rownames(mol_plot_neo_wtd)) {
  sp = exotics[which(vegcomm[i,exotics]>0)]
  pi = vegcomm[i,sp]
  mol_plot_exo_wtd[i,] <- colSums(mol_mat[sp,] * t(pi))
}

## archaeophytes only ####
mol_plot_arc <-data.frame(matrix(NA, nrow= nrow(vegcomm), ncol = ncol(mol_mat),
                                 dimnames = list(rownames(vegcomm), names(mol_mat))))

for (i in rownames(mol_plot_arc)) {
  mol_plot_arc[i,] <- colSums(mol_mat[archaeophytes[which(vegcomm[i,archaeophytes]>0)],])
}


mol_plot_arc_wtd <-data.frame(matrix(NA, nrow= nrow(vegcomm), ncol = ncol(mol_mat),
                                     dimnames = list(rownames(vegcomm), names(mol_mat))))

for (i in rownames(mol_plot_arc_wtd)) {
  sp = archaeophytes[which(vegcomm[i,archaeophytes]>0)]
  pi = vegcomm[i,sp]
  mol_plot_arc_wtd[i,] <- colSums(mol_mat[sp,] * t(pi))
}



### ALLERGEN FAMILIES ######

# Presence of allergenic molecules families per plot ####
allfam_plot <- matrix(NA, nrow= nrow(vegcomm), ncol = ncol(allfam_mat),
                              dimnames = list(rownames(vegcomm), names(allfam_mat)))

for (i in rownames(mol_plot)) {
  allfam_plot[i,] <- as.numeric(
    colSums(
      allfam_mat[colnames(vegcomm)[which(vegcomm[i,]>0)],])>0)
}

# Cover of allergenic molecules families per plot ####
allfam_plot_wtd <- matrix(NA, nrow= nrow(vegcomm), ncol = ncol(allfam_mat),
                      dimnames = list(rownames(vegcomm), names(allfam_mat)))

for (i in rownames(mol_plot)) {
    sp = colnames(vegcomm)[which(vegcomm[i,]>0)]
    pi = vegcomm[i,sp]
    allfam_plot_wtd[i,] <- colSums(allfam_mat[sp,] * t(pi))
  }



# Presence of allergenic molecules families per plot represented by residents
allfam_plot_nat <- matrix(NA, nrow= nrow(vegcomm), ncol = ncol(allfam_mat),
                                 dimnames = list(rownames(vegcomm), names(allfam_mat)))

for (i in rownames(mol_plot)) {
  allfam_plot_nat[i,] <- as.numeric(
    colSums(
      allfam_mat[natives[which(vegcomm[i,natives]>0)],
        ]
    )>0)
}

# Presence of allergenic molecules families per plot represented by neophytes
allfam_plot_neo <- matrix(NA, nrow= nrow(vegcomm), ncol = ncol(allfam_mat),
                                       dimnames = list(rownames(vegcomm), names(allfam_mat)))

for (i in rownames(mol_plot)) {
  allfam_plot_neo[i,] <- as.numeric(
    colSums(
      allfam_mat[neophytes[which(vegcomm[i,neophytes]>0)],
        ]
    )>0)
}

# Presence of allergenic molecules families per plot represented by Archaeophytes
allfam_plot_arc <- matrix(NA, nrow= nrow(vegcomm), ncol = ncol(allfam_mat),
                          dimnames = list(rownames(vegcomm), names(allfam_mat)))

for (i in rownames(mol_plot)) {
  allfam_plot_arc[i,] <- as.numeric(
    colSums(
      allfam_mat[archaeophytes[which(vegcomm[i,archaeophytes]>0)],
        ]
    )>0)
}

# Presence of allergenic molecules families per plot represented by EXOTICS
allfam_plot_exo <- matrix(NA,
                          nrow= nrow(vegcomm),
                          ncol = ncol(allfam_mat),
                          dimnames = list(rownames(vegcomm), names(allfam_mat)))

for (i in rownames(mol_plot)) {
  allfam_plot_exo[i,] <- as.numeric(
    colSums(
      allfam_mat[exotics[which(vegcomm[i,exotics]>0)],
        ]
    )>0)
}
