# Import data from the Cityscapelab Berlin

# Import vegetation and environmental data from Cityscapelabs ####

## Abiotic plot data
plot_summary <- fread(file = "data/environmental_factors_data.csv",
                      data.table = FALSE)
row.names(plot_summary) <- plot_summary$ID_plot

## vegetation survey data (community matrix)
vegcomm <- fread(file = "data/species_abundance_data.csv",
                 data.table = FALSE, )
row.names(vegcomm) <- vegcomm$V1  
vegcomm <- vegcomm[,-1]


## Species taxonomic and trait data
species_data <- fread(file = "data/species_data.csv",
                      data.table = FALSE )
row.names(species_data) <- species_data$Species_col

# Remove tree/large shrub species from vegetation dataset (if any): ####

## list of species:
tree.species <- c('Acer', 'Acer campestre', 'Acer negundo', 'Acer platanoides', 'Acer pseudoplatanus','Crataegus monogyna' ,'Prunus', 'Prunus serotina','Prunus spinosa', 'Prunus domestica', 'Pinus sylvestris', 'Quercus', 'Tilia', 'Populus', 'Populus tremula', 'Quercus robur', 'Robinia pseudoacacia',  "Pyrus_communis_agg.")

## Check out how abundant the tree seedlings were:
tree_comm <- vegcomm[,which(colnames(vegcomm) %in%
                 sub(" " , "_", tree.species))]

sum(rowSums(tree_comm>0)>0) # 32/56 plots with at least 1 tree
summary(rowSums(tree_comm>0)) 
sd(rowSums(tree_comm>0))
# max = 5 tree species, median = 1, mean  = 0.94, sd = 1.09
range(rowSums(tree_comm)) # max = 2.4 % cumulated cover for trees
median(rowSums(tree_comm)) # median = 0.1 % cover
mean(rowSums(tree_comm)) # mean = 0.26 % cover
sd(rowSums(tree_comm)) # sd = 0.47 % cover

# Remove trees from vegetation matrix:
if (sum(colnames(vegcomm)%in% sub(" " , "_", tree.species)) >0) {
  vegcomm <- vegcomm[,-which(colnames(vegcomm) %in%
                               sub(" " , "_", tree.species))]
}

# Reorder dataframes: ####
## Order plots like in the community matrix:
plot_summary <- plot_summary[row.names(vegcomm),]

## Order species like in vegcomm:
#species_data <- species_data[colnames(vegcomm),]


# Extract species introduction status ( I = Indigenous, A = Archaeophyte, N = Neophyte) ####
natives = row.names(species_data)[which(species_data$Introduction_status_Seitz2012 %in% c("I"))]
neophytes = row.names(species_data)[which(species_data$Introduction_status_Seitz2012 %in% c("N"))]
archaeophytes = row.names(species_data)[which(species_data$Introduction_status_Seitz2012 %in% c("A"))]
exotics = row.names(species_data)[which(species_data$Introduction_status_Seitz2012 %in% c("N", "A"))]
residents = row.names(species_data)[which(species_data$Introduction_status_Seitz2012 %in% c("I", "A"))]


# Calculate community composition, diversity and cover indices ####

## Calculate species richness:
plot_summary$SR <- rowSums(vegcomm>0)
plot_summary$SR.nat <- rowSums(vegcomm[, !colnames(vegcomm) %in% exotics]>0)
plot_summary$SR.neo <- rowSums(vegcomm[, colnames(vegcomm) %in% neophytes]>0)
plot_summary$SR.arch <- rowSums(vegcomm[, colnames(vegcomm) %in% archaeophytes]>0)
plot_summary$SR.exo <- rowSums(vegcomm[, colnames(vegcomm) %in% exotics]>0)
plot_summary$SR.resid <- rowSums(vegcomm[, colnames(vegcomm) %in% residents]>0)

## Calculate cumulated percent cover:
plot_summary$cover.nat <- rowSums(vegcomm[, colnames(vegcomm) %in% natives])
plot_summary$cover.neo <- rowSums(vegcomm[, colnames(vegcomm) %in% neophytes])
plot_summary$cover.arch <- rowSums(vegcomm[, colnames(vegcomm) %in% archaeophytes])
plot_summary$cover.exo <- rowSums(vegcomm[, colnames(vegcomm) %in% exotics])
plot_summary$cover.resid <- rowSums(vegcomm[, colnames(vegcomm) %in% residents])

## Calculate cumulated relative percent cover:
plot_summary$relcover.nat <- rowSums(vegcomm[, !colnames(vegcomm) %in% exotics])/rowSums(vegcomm)
plot_summary$relcover.neo <- rowSums(vegcomm[, colnames(vegcomm) %in% neophytes])/rowSums(vegcomm)
plot_summary$relcover.arch <- rowSums(vegcomm[, colnames(vegcomm) %in% archaeophytes])/rowSums(vegcomm)
plot_summary$relcover.exo <- rowSums(vegcomm[, colnames(vegcomm) %in% exotics])/rowSums(vegcomm)
plot_summary$relcover.resid <- rowSums(vegcomm[, colnames(vegcomm) %in% residents])/rowSums(vegcomm)

# Calculate proportion of species:
plot_summary$prop.nat <- as.numeric(plot_summary$SR.nat)/as.numeric(plot_summary$SR)
plot_summary$prop.neo <- as.numeric(plot_summary$SR.neo)/as.numeric(plot_summary$SR)
plot_summary$prop.arch <- as.numeric(plot_summary$SR.arch)/ as.numeric(plot_summary$SR)
plot_summary$prop.exo <- as.numeric(plot_summary$SR.exo)/ as.numeric(plot_summary$SR)
plot_summary$prop.resid <- as.numeric(plot_summary$SR.resid)/ as.numeric(plot_summary$SR)

