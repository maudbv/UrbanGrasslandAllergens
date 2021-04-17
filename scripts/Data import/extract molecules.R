## Extract information about allergenic molecules per species

# FUNCTION Assemble all the molecular data collected from databases ####
# We try matching species names and species Genus
# Matches for species name are match = 1
# Matches at the genus level are match = 3
# fuzzy matching is match = 2 (deprecated)

complete_molecule_data <- (function () {
  
# Create template to store database matches  ####
molecule_template <- data.frame(Species = NA, Matched.name = NA,
                                Allergen.mol=NA, Type=NA, Group=NA,
                                row.index = NA, Family=NA, 
                                Allergen.name = NA,
                                Ref=NA, match = NA)

#  Allergen Online (also includes all IUIS) 
AllergenOnline_data <- molecule_template
## Test each species in a loop:
for (sp in species_allergen$Species_acceptedname) {
  #extract the genus
gen <- species_allergen[species_allergen$Species_acceptedname == sp, "Genus"]

# match species and genus in the database:
sp.e <-  grep(sp, x =  AllergenOnline$Species)
gen.e <-  grep(gen, x =  AllergenOnline$Species)

# collect all unique matched names:
ind <- unique(c(sp.e, gen.e))

# Assign match level (1 or 3) to each name: 
match = rep(3, length(ind))
match[ind %in% sp.e] <- 1

# if at least one match: 
if (length(match)>0) {
  # extract the information in the order of the template: 
     tmp <-  data.frame(
       sp, 
       AllergenOnline[ind,c("Species", "IUIS.Allergen", "Type", "Group")],
       ind,
       Family=NA, 
       Allergen.name = NA,
       Ref="AllergenOnline",
       match)
    # correct col names to match template
     names(tmp) <- names(molecule_template)
     
    # Add new data as additional rows for the output:  
     AllergenOnline_data = rbind(AllergenOnline_data, tmp)
    }
}

# IUIS allergens
IUISallergens_data <- molecule_template
for (sp in species_allergen$Species_acceptedname) {
  gen <- species_allergen[species_allergen$Species_acceptedname == sp, "Genus"]
  
  sp.e <-  grep(sp, x =   IUISallergens$Species)
  gen.e <-  grep(gen, x =   IUISallergens$Species)
  
  ind <- unique(c(sp.e, gen.e))
  match = rep(3, length(ind))
  match[ ind %in% sp.e ] <- 1
  
  if (length(match)>0) {
    tmp <-  data.frame(
      sp,
      IUISallergens[ind,c("Species", "Allergen", "Route.of.Allergen.Exposure")],
      Group = NA,
      ind,
      Family=NA,
      Allergen.name = IUISallergens[ind,"Biochemical.name"],
      Ref="IUISallergens",
      match)
    
    names(tmp) <- names(molecule_template)
    IUISallergens_data = rbind(IUISallergens_data, tmp)
  }
}

# SDAPallergen
SDAPallergen_data <- molecule_template
for (sp in species_allergen$Species_acceptedname) {
  gen <- species_allergen[species_allergen$Species_acceptedname == sp, "Genus"]
  
  sp.e <-  grep(sp, x =  SDAPallergen$Species)
  gen.e <-  grep(gen, x =  SDAPallergen$Species)
  
  ind <- unique(c(sp.e, gen.e))
  match = rep(3, length(ind))
  match[ind %in% sp.e] <- 1
  
  if (length(match)>0) {
    tmp <-  data.frame(
      sp, 
      SDAPallergen[ind,c("Species", "Allergen", "Allergen.Type")],
      Group = NA,
      ind,
      Family=NA, 
      Allergen.name = SDAPallergen[ind,"Allergen.Description"],
      Ref="SDAPallergen",
      match)
    names(tmp) <- names(molecule_template)
    SDAPallergen_data = rbind( SDAPallergen_data, tmp)
  }
}

# Allergome
AllergomeData_data <- molecule_template
AllergomeData$Source_mod <- gsub(x = AllergomeData$Source, pattern = " " , replacement= "_")

for (sp in species_allergen$Species_acceptedname) {
  gen <- species_allergen[species_allergen$Species_acceptedname == sp, "Genus"]
  gen <- paste(gen, "_", sep = "")
  
  sp.e <-  grep(sp, x =  AllergomeData$Source_mod)
  gen.e <-  grep(gen, x =  AllergomeData$Source_mod, fixed = TRUE)
  
  ind <- unique(c(sp.e, gen.e))
  match = rep(3, length(ind))
  match[ind %in% sp.e] <- 1
  
  if (length(match)>0) {
    tmp <-  data.frame(sp, AllergomeData[ind,c("Source_mod", "Name")],
                       Type = NA, Group = NA, ind , Family=NA, 
                       Allergen.name = NA, Ref="AllergomeData",
                       match)
    
    names(tmp) <- names(molecule_template)
    AllergomeData_data = rbind(AllergomeData_data, tmp)
  }
}

# Allfam (NO LONGER MATCHED - cf below)
AllFamData_data <- molecule_template

# for (sp in species_allergen$Species_acceptedname) {
#   gen <- species_allergen[species_allergen$Species_acceptedname == sp, "Genus"]
#   
#   sp.e <-  grep(sp, x =  AllFamData$Species)
#   gen.e <-  grep(gen, x =  AllFamData$Species)
#   
#   ind <- unique(c(sp.e, gen.e))
#   match = rep(3, length(ind))
#   match[ind %in% sp.e] <- 1
#   
#   if (length(match)>0) {
#     tmp <-  data.frame(sp, AllFamData[ind,c("Species", "Name", "Routes.of.exposure")],Group = NA,
#                        ind, Family=  AllFamData[ind,c("Family")],
#                        Allergen.name = NA, Ref="AllFam", match)
#     names(tmp) <- names(molecule_template)
#     AllFamData_data = rbind( AllFamData_data, tmp)
#   }
# }


## Assemble all data: ####

complete_molecule_data <- rbind( IUISallergens_data[-1,], AllergenOnline_data[-1,],SDAPallergen_data[-1,], AllergomeData_data[-1,])

return(complete_molecule_data)
})()


# FUNCTION Clean up data and add additional data manually: ####
molecule_data_raw <- (function () {
# Create table to clean up and store original complete data
molecule_data <- complete_molecule_data

# Replace species names by the names in the vegetation survey 
# Matching between datatables done via "accepted names":
molecule_data$Species_acceptedname <- molecule_data$Species
molecule_data$Species <- species_allergen$Species[
  match(molecule_data$Species_acceptedname, 
        species_allergen$Species_acceptedname)
  ]


## correct manually the Unassigned allergen names :
# write.csv( molecule_data[which(molecule_data$Allergen.mol %in% c("Unassigned", "")),], file = "Unassigned molecules new.csv")
unassigned_molecules <- fread(
  file = "data/manual updates/Unassigned molecules.csv",
  data.table = FALSE)

# Extract manually corrected allergen molecules in a new column:
# Match by group = original long version of the molecule name, 
# which I corrected manually in some cases
molecule_data$Allergen.mol.manual <- unassigned_molecules$Allergen.mol.manual[
  match(molecule_data$Group, unassigned_molecules$Group)]

# Replace allergen names by the manually corrected names
molecule_data$Allergen.mol[!is.na(molecule_data$Allergen.mol.manual)]<- molecule_data$Allergen.mol.manual[!is.na(molecule_data$Allergen.mol.manual)]

# visually check that corrected names match the original "group" column:
# View(molecule_data[molecule_data$Group %in% unassigned_molecules$Group,c("Allergen.mol","Group")])

 
## Manually remove erroneous allergen molecules :
allergens_to_be_removed <- c(
  "Jun c 1",   # erroneously matches Pear tree to Juniperus communis
  "Tar o RAP", # it is a root allergen
  "Art v 60kD", # it is the same as Art v1
  "Art v 47kD", # it is the same as Art v1
  "cyclophilin", "Dau c 1", "Dau c 2", # Remove Daucus carota
  "Unassigned", # useless unsassigned entry for genus Ambrosia 
  "Poa p b"  # Only mentioned once in SDAP database, no other support found anywhere
)
rows_to_remove <- which(molecule_data$Allergen.mol %in% allergens_to_be_removed )
if (length(rows_to_remove) >0) molecule_data <- molecule_data[-rows_to_remove,]


# Simplify molecule names: 
molecule_data$Allergen.simple <- make.names(
  sapply(molecule_data$Allergen.mol,
         FUN = function(s) {
           paste(strsplit(strsplit(x = s , split = "\\.")[[1]][1],
                          split = " ")[[1]], collapse = "_")
           }
         )
  )


# remove food allergens (= not pollen)
rows_food_allergens <- grep("food", molecule_data$Type, ignore.case = TRUE)
if (length(rows_food_allergens)>0){
  molecule_data <- molecule_data[-rows_food_allergens,]
}

## Assign IUIS or non-IUIS criteria
molecule_data$IUIS <- molecule_data$Allergen.simple %in% molecule_data$Allergen.simple[molecule_data$Ref == "IUISallergens"]

# remove most obvious repeated entries: 
molecule_data <- molecule_data[- which(duplicated(molecule_data[,c("Species", "Matched.name","Allergen.simple","Group", "Ref")])),]
print(nrow(molecule_data))
#755 rows left

## Extract allergen family from the ALLFAM database ##
# AT FIRST: Simple match with allergen names to the AllFam database
molecule_data$AllFam <- AllFamData[
  match(molecule_data$Allergen.simple, 
        gsub(AllFamData$Name, pattern = " ", replacement = "_")),
  "Family"]
molecule_data$AllFam[is.na(molecule_data$AllFam)] <- "Unknown"

# Unknown Allergen Families will be updated progressively in a new column:
molecule_data$Allergen.Family <- molecule_data$AllFam

# Add genus level allergens:  #####
# add genus column
molecule_data$genus <- species_data[
  match(molecule_data$Species,species_data$Species_col),
  "genus"]


## Homogenize molecule names within genus  for case when match == 3 :
# Only 27 unique genera in total:

molecule_data$Allergen.genus <- NA

for (g in unique(molecule_data$genus)) {
  mat <-  molecule_data[molecule_data$genus == g,]
  mols <- unique(mat$Allergen.simple)
  
  # Create simplified allergen names, 
  # assuming numbers within genera are equivalent molecules
  mat$Allergen.genus <- paste(mat$genus, 
                              sapply(mat$Allergen.simple,
                                     function(x) {
                                       y = strsplit(x, split= "_")[[1]]
                                       y<- y[length(y)]
                                       return(y)
                                     }),
                              sep = "_")
  
  # Re-assign the allergen family: 
  for (i in unique(mat$Allergen.genus)) {
    allf <- mat$Allergen.Family[mat$Allergen.genus == i]
    if ( "Unknown" %in% allf)  allf <- allf[-which(allf == "Unknown")]
    if (length(unique(allf)>1)) allf <- names(which(table(allf) == max(table(allf))))
    if (length(unique(allf))==0) allf = "Unknown"
    mat$Allergen.Family[mat$Allergen.genus == i] <- unique(allf)[1]
  }
  molecule_data[molecule_data$genus == g,] <-  mat
}



## Add simplified GRASS allergens PER TRIBE  ######
# add subfamilymy and tribe column
molecule_data$subfamily <- species_data[
  match(molecule_data$Species,species_data$Species_col),
  "subfamily"]

molecule_data$tribe <- species_data[
  match(molecule_data$Species,species_data$Species_col),
  "tribe"]

grass_allergens <- fread('data/grass allergens.csv', 
                         data.table = FALSE,
                         check.names = TRUE)
grasses<- species_data[which(species_data$family == "Poaceae"),"Species_col"]
pooideae<- species_data[which(species_data$subfamily == "Pooideae"),"Species_col"]
poas<- species_data[which(species_data$genus == "Poa"),"Species_col"]


tmp <- molecule_data[molecule_data$Species %in% grasses,]
# 15 grasses in th emolecule data, out of 36 in the vegetatino survey
# We are going to add new species and neww molecules to the molecule dataframe:
# 31 of which are pooideae...
# Molecules assigned per genus or tribe
dim(tmp)

grass.matrix <- sapply(grasses, function(sp) {
  alls = c("1","4","7","12")  # Allergen groups shared by most grasses 
  if (sp %in% pooideae) alls = c(alls,"5") #(5 only for pooideae)
  if (sp %in% poas) alls = c(alls,"6") # 6 (5 only for genus Poa)
  
  tribe <- species_data[which(species_data$Species_col == sp), "tribe"]
  Allergen.created <- paste(substr(tribe, 1,3), alls, sep = "_")
  
  o.mat = NULL
  mat = NULL
  
  # Case when there are already some known allergens for the species
  if (sp %in% tmp$Species) {
    o.mat <- tmp[tmp$Species == sp,]
    o.mat$grass.group <- sapply(o.mat$Allergen.simple,function(x) strsplit(x, split ="_")[[1]][3])
    o.mat$Allergen.created [o.mat$grass.group %in% alls] <-  paste(substr(tribe, 1,3),
                                                                   o.mat$grass.group[o.mat$grass.group %in% alls],
                                                                   sep = "_")
    
    # remove the grass allergens which are already matched
    alls <- alls [- which(alls %in% o.mat$grass.group)]
  }
  
  if (length(alls) >0) {
    # For the molecules which did not exist yet: 
    mat = as.data.frame(matrix(NA, nrow= length(alls), ncol = (ncol( tmp))))
    names(mat) <- names(tmp)
    mat$Species <- sp
    mat$grass.group <- alls
    mat$Allergen.created <- paste(substr(tribe, 1,3), alls, sep = "_")
  }
  
  # Bind the original and the new data:
  if (!is.null(o.mat)) mat <- rbind(o.mat, mat)
  
  # return the matrix:
  return(mat)
},simplify =FALSE)
grass.matrix<- do.call(rbind.data.frame, grass.matrix)

#correct allergen family
grass.matrix$Allergen.Family.grass <- grass_allergens$Protein.family[
  match(grass.matrix$grass.group, grass_allergens$Allergen.group)]

# Merge columns with grass allergens to molecule_data
tmp <- merge(molecule_data, grass.matrix, all = TRUE)

# Replace the name of grass allergens by the tribe level name :
tmp$Allergen.simple.old <- tmp$Allergen.simple
tmp$Allergen.simple[!is.na(tmp$Allergen.created)] <- tmp$Allergen.created[!is.na(tmp$Allergen.created)] 

# These grasses will not be counted as "genus level" matches (ie match = 3), but as almost species level matches  (ie match = 1)
tmp$init.match <- tmp$match # Store records of initial matching
tmp$match[!is.na(tmp$Allergen.created)] <- "1b"

# Correct allergen family for the grasses: 
tmp$Allergen.Family[!is.na(tmp$Allergen.Family.grass)] <- 
  tmp$Allergen.Family.grass[!is.na(tmp$Allergen.Family.grass)]

# update molecule_data:
molecule_data <- tmp

## Remove extra row of Fes p 4w :
molecule_data <- molecule_data[-which(molecule_data$Allergen.simple %in% c("Fes_p_4w","Fes_p_4")),] 

# correct genus allergens to include tribe level molecules:
molecule_data$Allergen.genus[!is.na(molecule_data$Allergen.created)] <-
  molecule_data$Allergen.created[!is.na(molecule_data$Allergen.created)]

### Manually correct the rest of the missing Allergen Families ####
## Identify missing families:  98 molecules out of 804
missing.families <- molecule_data[  molecule_data$Allergen.Family %in% c("Unknown","") |
                                    is.na(molecule_data$Allergen.Family) ,]
dim(missing.families)
# remove duplicated allergen molecules: 
missing.families <- missing.families[-which(duplicated(missing.families$Allergen.mol)),]
dim(missing.families)
# Only 29 unique allergen molecules are concerned

# Import manually corrected Allergen Families: 
molecule_update <- fread(
  "data/manual updates/manual_updates_allergen_families_Oct2020.csv",
  data.table = FALSE, header = TRUE)

# Test if the updates cover all the missing data: OK
which(!missing.families$Allergen.simple %in% molecule_data$Allergen.simple)


# One allergen is not matched (because removed)
molecule_update$Allergen.simple[which(!molecule_update$Allergen.simple %in% molecule_data$Allergen.simple)]

# Merge molecule data with the manually updated entries:
molecule_data_tmp <-cbind(
  molecule_data, 
  molecule_update[match(molecule_data$Allergen.simple, molecule_update$Allergen.simple),
                  c("AllFam.manual","manual.source")])

# molecule_data_tmp$Allergen.Family <- molecule_data_tmp$AllFam
molecule_data_tmp$Allergen.Family[which(molecule_data_tmp$Allergen.Family == "Unknown")] <-
  molecule_data_tmp$AllFam.manual[which(molecule_data_tmp$Allergen.Family == "Unknown")]
      
# replace by corrected version:
molecule_data <- molecule_data_tmp



# One allergen is not matched (because removed)
molecule_update$Allergen.simple[which(!molecule_update$Allergen.simple %in% molecule_data$Allergen.simple)]

# Merge molecule data with the manually updated entries:
molecule_data_tmp <-cbind(
  molecule_data, 
  molecule_update[match(molecule_data$Allergen.simple, molecule_update$Allergen.simple), c("AllFam.manual","manual.source")])

# molecule_data_tmp$Allergen.Family <- molecule_data_tmp$AllFam
molecule_data_tmp$Allergen.Family[which(molecule_data_tmp$AllFamily == "Unknown")] <-
  molecule_data_tmp$AllFam.manual[which(molecule_data_tmp$AllFamily == "Unknown")]

# Check for missing families : 
which(molecule_data_tmp$AllFamily == "Unknown")

# replace by corrected version:
molecule_data <- molecule_data_tmp


# Add manually some additional molecule data  ####
# = data only suggested in the literature, no serious evidence yet

# Taraxacum sp. in Yong et al. 2007: Show that T. is cross reactive with Artemisia.
# "Chrysanthemum and dandelion showed several protein fractions similar in molecular weight to Art v 1, Art v 4, or Art v 6. "``
sp <- species_allergen$Species[grep("Taraxacum", species_allergen$Species)]

tmp <- unique(molecule_data[
  which(molecule_data$Species == "Artemisia_vulgaris" &
          molecule_data$Allergen.simple %in% c("Art_v_1","Art_v_4","Art_v_6")&
          molecule_data$Ref ==  "IUISallergens"),])
tmp<- rbind(tmp, tmp) # There are 2 taraxacum species to update
tmp$Ref <- "cross reactivity with Artemisia vulgaris in Yong et al. 2007"
tmp$Species <- c(rep(sp[1],3), rep(sp[2],3))
tmp$genus <- "Taraxacum" 
tmp$match <- "4"  # Even less reliable matching
molecule_data <- rbind(molecule_data,tmp)


# FOR RUMEX: perhaps there is data in: 
# (NOT DONE) https://www.sciencedirect.com/science/article/abs/pii/S0091674907013930?via%3Dihub
# (DONE) https://www.uniprot.org/uniprot/Q9SWD5
# At least one Beta-expansin allergen for Rumex:
sp <- species_allergen$Species[grep("Rumex", species_allergen$Species)]
tmp <- unique(molecule_data[
  grep("expansin",molecule_data$Allergen.Family),])[1,]
tmp[, which(colnames(tmp) != "Allergen.Family")] <- NA
tmp<- rbind(tmp, tmp) # There are 2 taraxacum species to update
tmp$Species <- sp
tmp$Allergen.simple <- "Rum_expansin"
tmp$match <- 4
tmp$Ref <- "Uniprot Q9SWD5"
tmp$Allergen.genus <- "Rum_expansin"
molecule_data <- rbind(molecule_data,tmp)


# update the genus column for the new rows of data:
molecule_data$genus <- species_data[
  match(molecule_data$Species,species_data$Species_col),
  "genus"]


# Remove duplicates ####
molecule_data <- molecule_data[- which(duplicated(molecule_data[,c("Species", "Matched.name","Allergen.simple", "Allergen.Family", "Ref")])),]
dim(molecule_data) # 751


return(molecule_data)
})()


molecule_data <- molecule_data_raw

# Create a simplified table per species x allergen.simple:  ####
library(plyr)
molecule_data <- ddply(.data = molecule_data_raw,
                      .variables = .(Species, genus, Allergen.simple, Allergen.Family),
                      .fun = summarise,
                      Databases = paste(unique(Ref), collapse = ";"),
                      Ref.manual = paste(unique(manual.source), collapse = ";"),
                      Mol.alt = paste(unique(Allergen.mol), collapse = ";"),
                      Allergen.genus = paste(unique(Allergen.genus), collapse = ";"),
                      Types = paste(unique(Type), collapse = ";"),
                      Matched.names = paste(unique(Matched.name), collapse = ";"),
                      Species.in.data = paste(unique(Species), collapse = ";"),
                      Accepted.names.in.data = paste(unique(Species), collapse = ";"),
                      
                      Matches.with.database =  paste(unique(match), collapse = ";")

)
dim(molecule_data) ## 416 unique entries per species
unique(molecule_data$Allergen.simple) # 122 unique molecules

# Create a simplified table per simple allergen name:  ####
molecule_sum <- ddply(.data = molecule_data_raw, 
                      .variables = .(Allergen.simple, Allergen.Family),
                      .fun = summarise,
                      Databases = paste(unique(Ref), collapse = ";"),
                      Ref.manual = paste(unique(manual.source), collapse = ";"),
                      Mol.alt = paste(unique(Allergen.mol), collapse = ";"),
                      Types = paste(unique(Type), collapse = ";"),
                      Matched.names = paste(unique(Matched.name), collapse = ";"),
                      Species.in.data = paste(unique(Species), collapse = ";"),
                      Accepted.names.in.data = paste(unique(Species), collapse = ";"),
                      Matches.with.database =  paste(unique(match), collapse = ";")
                      
)
dim(molecule_sum) ## 122 unique entries 
rownames(molecule_sum) <- molecule_sum$Allergen.simple



# MATRICES of allergen molecule composition ####

## Species x  Molecule matrix 
mol_mat_strict <- data.frame(
  matrix(0,
         nrow = length(species_allergen$Species),
         ncol = length(unique(molecule_data$Allergen.simple)),
         dimnames = list(species_allergen$Species,
                         sort(unique(molecule_data$Allergen.simple)))))

dim(mol_mat_strict)                            
for (sp in species_allergen$Species) {
  smat <- molecule_data[intersect(which(molecule_data$Species == sp),
                                # Only perfect match:
                                grep("1",molecule_data$Matches.with.database)
                                ), ] 
  
  if (nrow(smat) >0) {
    mol_mat_strict[sp, unique(smat$Allergen.simple)] <- 1
  }
}



## Add allergen as "Unknowned" if the species is allergenic but does not show up in molecule database
mol_mat_strict$Unknown <- 0
mol_mat_strict[names(which(rowSums(mol_mat_strict[allergenics,])==0)),
        "Unknown"] <- 1

# Check that the number of allergenic species has not changed:
if (length(allergenics) != sum(species_allergen$allergenicity)) {
e <- simpleError(message = "Inconsistent number of allergenic species between species data and strict molecule data" )
stop(e)
} 
          
## Add info concerning missing molecule info
 species_allergen$molecule.unknown <- 0
 species_allergen$molecule.unknown[
   (rownames(species_allergen) %in% allergenics) &
     (rowSums(mol_mat_strict) == mol_mat_strict$Unknown)
   ] <- 1
 sum(species_allergen$molecule.unknown)  
 # 31 missing = 41 % of the allergenic species
  

## Matrix including genus matches = "Broad" :
mol_mat_broad <- data.frame(
  matrix(0,
         nrow = length(species_allergen$Species),
         ncol = length(unique(molecule_data$Allergen.genus)),
         dimnames = list(species_allergen$Species,
                         sort(unique(molecule_data$Allergen.genus)))))
 
for (sp in species_allergen$Species) {
  gen <- species_allergen$Genus[species_allergen$Species ==sp]
  smat <- molecule_data[which(molecule_data$genus == gen), ]
  
  if (nrow(smat) >0) {
    mol_mat_broad[sp, unique(smat$Allergen.genus)] <- 1
  }
}
 
 rowSums(mol_mat_broad)

## Add allergen as "Unassigned" if the species is allergenic but does not show up in molecule database
 mol_mat_broad$Unknown <- 0
 mol_mat_broad[names(which(rowSums(mol_mat_broad[allergenics,]) ==0)),
         "Unknown"] <- 1
 # 21 unknown molecules even at genus level =>10 were added compared to strict
 if (broad) print(paste("OPTION = BROAD added molecule info for",
       sum(mol_mat_strict$Unknown) - sum(mol_mat_broad$Unknown),
       "additional species, bringing the number of species with unknown molecules to",
       sum(mol_mat_broad$Unknown),
       sep = " "))

###  Add allergenicity in the BROAD sense to the species_allergen table: 
species_allergen$allergenicity.broad <- 0
species_allergen$allergenicity.broad[which(rowSums(mol_mat_broad)>0)] <- 1

# Check that this did not add new allergenic species, 
# only new molecule information:
  if (length(allergenics) != sum(species_allergen$allergenicity.broad)) {
   e <- simpleError(message = "Inconsistent number of allergenic species between species data and broad molecule data" )
   stop(e)
  } 

## Add info concerning missing molecule info
species_allergen$molecule.unknown <- 0
species_allergen$molecule.unknown[
  (rownames(species_allergen) %in% allergenics) &
    (rowSums(mol_mat_broad) == mol_mat_broad$Unknown)
  ] <- 1

missing.mols <- species_allergen[species_allergen$molecule.unknown ==1, ]

# OPTION : Choose to use only the broad version?  ####
if (broad == TRUE) {
  mol_mat <- mol_mat_broad
} else {
  mol_mat <- mol_mat_strict
}

## Matrix at Genus level : ####

mol_mat_gen <- t(sapply( na.omit(unique(species_allergen$Genus)), function(x) {
  mat <- species_allergen[species_allergen$Genus == x,]
  return( as.numeric(colSums(mol_mat_broad[unique(mat$Species),], na.rm = TRUE) >0))
}))
colnames(mol_mat_gen) <- colnames(mol_mat_broad)
rownames(mol_mat_gen) <- na.omit(unique(species_allergen$Genus))

## Taxonomic Family X Molecule matrix ####
mol_mat_fam <- t(sapply( na.omit(unique(species_allergen$family)), function(x) {
  mat <- species_allergen[species_allergen$family == x,]
  return( as.numeric(colSums(mol_mat_broad[unique(mat$Species),], na.rm = TRUE) >0))
}))
colnames(mol_mat_fam) = colnames(mol_mat_broad)



## Matrix by Allergen Family per species matrix ######

allfam_mat <-
  data.frame(matrix(
    0,
    nrow = length(species_allergen$Species),
    ncol = length(unique(molecule_data$Allergen.Family)),
    dimnames = list(species_allergen$Species,
                    sort(unique(
                      molecule_data$Allergen.Family
                    )))
  ),
  check.names = FALSE)

dim(allfam_mat)

for (sp in species_allergen$Species) {
  smat <- molecule_data[molecule_data$Species == sp ,]
  smat <- smat[ which(smat$Matches.with.database %in% c("1", "1b", "3","4" )),]
  if (nrow(smat) > 0) {
    allfam_mat[sp, unique(smat$Allergen.Family)] <- 1
  }
}
colSums(allfam_mat) # No unknown families

## Add allergen family as "Unknown" if the species is allergenic but does not show up
allfam_mat$Unknown<- 0
allfam_mat[names(which(! rowSums(allfam_mat[allergenics, ])>0)), "Unknown"] <- 1

## Taxonomic Family X Allergen family matrix
allfam_mat_fam <- as.data.frame(t(sapply( unique(species_allergen$family), function(x) {
  mat <- species_allergen[species_allergen$family == x,]
  y = 0
  if (nrow(mat) >0) {
    y = ceiling(colSums(allfam_mat[rownames(mat),], na.rm = T)/100)
  }
  return(y)
})))

colSums(allfam_mat_fam)
rowSums(allfam_mat_fam)



## Add molecule counts per species in species_allergen ####
species_allergen$known.mols <- rowSums(mol_mat_broad[
  rownames(species_allergen), -which(names(mol_mat_broad) == "Unknown")])
                                                     
species_allergen$known.allfam <- rowSums(allfam_mat[
  rownames(species_allergen), -which(names(allfam_mat) == "Unknown")])


# clean up
#rm(AllergenOnline_data, IUISallergens_data, SDAPallergen_data, AllergomeData_data, molecule_template)

