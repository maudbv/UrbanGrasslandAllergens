# Import data on species allergenic properties and molecules

# _______  Import allergenic properties of species ___________####

## Import raw flowering phenology data:
species_flowering <-fread("data/species_flowering_months_April 2020.csv",
                          check.names = TRUE, data.table = FALSE)

# transform Species_names in a simplified version:
species_flowering$Species <- sapply(as.character(species_flowering$Species_name), FUN = function(x) paste(strsplit(x, split = " ")[[1]], collapse = "_"))
species_flowering <- species_flowering[,-which(colnames(species_flowering)=="Species_name")]

## Import allergen severity datasheet 
# (only for the herbaceaous vegetation = 216 species) 
species_allergen_severity <- fread("data/Species allergenicity scores.csv", data.table = FALSE, check.names = TRUE)
rownames(species_allergen_severity) <- species_allergen_severity$Species

# make sure we chose a numeric allergenic score
species_allergen_severity$allergen_score <- 
  as.numeric(species_allergen_severity$allergen_score_corrected) 

# Merge datasets 
species_allergen <- merge(species_allergen_severity,species_flowering,
                          by = "Species", all.x = TRUE)


# Correct Taxonomy
species_allergen$Species.match <- species_data[
  match(sub("-", "\\.", species_allergen$Species_name),
        species_data$Species_name),
  "Species"]

species_allergen$family <- 
  species_data[species_allergen$Species.match,]$family

species_allergen$subfamily <- 
  species_data[species_allergen$Species.match,]$subfamily

species_allergen <- separate(species_allergen,
                             col =  "Species.match",
                             into = "Genus", 
                             sep = "_", remove = FALSE, extra = "drop")


# Add missing family for Carex praecox
species_allergen$family[
  species_allergen$Species=="Carex_praecox"] <- "Cyperaceae"

# add accepted name columns to species allergen to match with other databases
species_allergen$Species_acceptedname <-  gsub(
  pattern=" ", replacement= "_",
  species_data$acceptedname[match(species_allergen$Species.match,
                                  species_data$Species)])



# _____ Import Online databases _________________________ ####

# Allergen Online :
AllergenOnline<-fread(file = "data/AllergenOnline_v19- September 2019.csv", check.names = TRUE, data.table = FALSE ) 
AllergenOnline$Species <- gsub(AllergenOnline[,1], pattern = " ", replacement = "_")
names(AllergenOnline) <- sub("\\..", "\\.", names(AllergenOnline))

unique(match(AllergenOnline$Species, species_allergen$Species_acceptedname))
unique(match(AllergenOnline$Species, species_allergen$Species))

# SDAP allergen database: 
SDAPallergen<-fread(file = "data/SDAP allergen database.csv", check.names = T,  data.table = FALSE ) 
SDAPallergen$Species <- sub(SDAPallergen$Species...Scientific.Name, pattern = " ", replacement = "_")
SDAPallergen$Species <- gsub(SDAPallergen$Species, pattern = " ", replacement = "")
names(SDAPallergen) <- sub("\\.\\.\\.", "\\.",names(SDAPallergen))

unique(match (SDAPallergen$Species, species_allergen$Species_acceptedname))
unique(match (SDAPallergen$Species, species_allergen$Species))


# WHO_IUIS allergen database
IUISallergens <- (function(){
  
  IUISallergens <- fread(file = "data/WHO_IUIS allergen database 2019.csv",  check.names = T,data.table = FALSE) 
  
  # Reorganize the IUIS data and complete species names column:
  tmp <- IUISallergens[,1]
  i = 1
  while ( i < length(tmp) ) {
    test = FALSE
    while (test == FALSE &  (i < length(tmp)) ) {
      if (tmp[i+1] == "" ) {
        tmp[i+1] = tmp[i]
      } else { test <- TRUE }
      i = i + 1
      print (tmp[i])
    }
  }
  IUISallergens$SpeciesName <- tmp
  
  IUISallergens$Species <- sapply(tmp, FUN = function(x) {
    paste(strsplit(x, " ")[[1]][1:2], collapse = "_")
  })
  
  IUISallergens <- IUISallergens[- which(IUISallergens$Allergen == ""),]
  
  unique(match (IUISallergens$Species, species_allergen$Species_acceptedname))
  unique(match (IUISallergens$Species, species_allergen$Species))
  
  return(IUISallergens)
})()

#  Allergome : protein-centered dataset 
AllergomeData <- fread(file = "data/Allergome molecules_2019.csv",  check.names = T,data.table = FALSE) 


# ALLFAMM family of molecules database
AllFamData <- fread(file = "data/Allfam molecule database.csv", data.table = FALSE,  check.names = T,header = TRUE) 
AllFamData  <- separate(AllFamData, "Source", into = c("Species","common_name"), sep = "_\\(", extra = "merge")
AllFamData <- AllFamData[ grep("Inhalation", AllFamData$Routes.of.exposure ),]
AllFamData$common_name <- gsub(pattern = "\\)", replacement = "", AllFamData$common_name)
AllFamData$common_name <- gsub(pattern = "\\(", replacement = "", AllFamData$common_name)

colnames(AllFamData)[c(2,3)] <- c("IUIS","AllergenOnline")
AllFamData[,"IUIS"] <- as.numeric(AllFamData[,"IUIS"] == "IUIS")
AllFamData[,"AllergenOnline"] <- as.numeric(AllFamData[,"AllergenOnline"] == "AllergenOnline")

unique(match (AllFamData$Species, species_allergen$Species_acceptedname))
unique(match (AllFamData$Species, species_allergen$Species))


## __________Merge allergen data in one table____________________ ####

# add online records as allergenicity binary criterion
species_allergen$IUIS.record <- 
  species_allergen$SDAP.record <-
  species_allergen$Allergenonline.record <-
  species_allergen$AllFam.record <-  0

# If species appears in the database, mark as "1"
species_allergen$IUIS.record [
  species_allergen$Species_acceptedname %in% IUISallergens$Species] <- 1
species_allergen$SDAP.record [
  species_allergen$Species_acceptedname  %in% SDAPallergen$Species] <- 1
species_allergen$Allergenonline.record [
  species_allergen$Species_acceptedname  %in% AllergenOnline$Species] <- 1
species_allergen$AllFam.record [
  species_allergen$Species_acceptedname  %in% AllFamData$Species] <- 1

# Different method needed for Allergome: 
# extract names of species from Allergome database
species_allergen$Allergome.record <-  sapply(
  1:length(species_allergen$Species_acceptedname),
  FUN = function(k) {
    nam <-
      sub(
        species_allergen$Species_acceptedname[k],
        pattern = "_",
        replacement = " "
      )
    rec <- grep(AllergomeData$Source, pattern = nam)
    if (sum(rec > 0) > 0) {
      x <- 1
      print(paste("We were looking for", nam))
      print(paste("We found", unlist(AllergomeData[rec, "Source"])))
    } else {
      x <- 0
    }
    return(x)
  }
)


#_____________________ Update allerginicity classification _________ ####
species_allergen$allergenicity <- as.numeric(!(rowSums(
  data.frame(
    species_allergen$allergen_score,
    species_allergen$IUIS.record,
    species_allergen$SDAP.record,
    species_allergen$Allergenonline.record,
    species_allergen$Allergome.record
  ),
  na.rm = T
) == 0))

# Summarise whether species appears in one of the online databases: 
species_allergen$database.record <- rowSums(
  species_allergen[, c("IUIS.record","Allergenonline.record",
                       "Allergome.record","AllFam.record")])>0
# total: 22 species

# Remove allergenicity entry for "Daucus carota",
# which only creates known food allergies (root), not pollen: 
which((species_allergen$allergen_score==0) & 
        species_allergen$database.record ==1)

species_allergen[species_allergen$Species.match == "Daucus_carota",
                 "allergenicity"] <- 0

# Identify the severely allergenics (score >3)
species_allergen$severe.all <- 0
species_allergen$severe.all[which(species_allergen$allergen_score > 2)] <- 1

# Add row names
row.names(species_allergen) <- species_allergen$Species.match


# add frequency in vegetation surveys from 2017
species_allergen <- species_allergen[match(colnames(vegcomm),
                                           rownames(species_allergen)),]
species_allergen$Frequency <- colSums(vegcomm>0)

# Check that all species do occur in the plot: 
stopifnot(all(!which(species_allergen$Frequency == 0)))

# Add frequency along urbanisation gradient: 
# three levels of urbanisation: 
plot_summary$Seal_levels <- cut(
    plot_summary$Seal_500,
    breaks = c(0,quantile(plot_summary$Seal_500,c(1/3,2/3)),100))
table(plot_summary$Seal_levels)
levels(plot_summary$Seal_levels) <- c(
  "near-rural","low-urban", "high-urban")

vegcomm <- vegcomm[rownames(plot_summary),]
species_allergen$Freq.rural <- colSums(
vegcomm[which(plot_summary$Seal_levels == "near-rural"), ]>0)
species_allergen$Freq.suburb <- colSums(
  vegcomm[which(plot_summary$Seal_levels == "low-urban"), ]>0)
species_allergen$Freq.urban <- colSums( 
  vegcomm[which(plot_summary$Seal_levels == "high-urban"), ]>0)

## ADD species trait data and florsitic status #### 
species_allergen <- merge(
  species_allergen,
  species_data[,c("Species","Introduction_status_Seitz2012",
                  "Origin","status_num","first_evidence",
                  "SLA_mean", "PlantHeight_mean",
                  "SeedDryMass_mean", 
                  "poll_vect_B", "flw_muell",      
                  "FunGroup","LifeForm"                       
  )],
  by.x = "Species", by.y ="Species", all.x = TRUE)

species_allergen <- as.data.frame(species_allergen)

# Format species names into "column names" to match with vegcomm 
# replaces "-" with "." :
species_allergen$Species <- make.names(species_allergen$Species)
rownames(species_allergen) <- species_allergen$Species

## Correct Flowergin Phenology ####
# Add numerical versions of flowering months
species_allergen$fl.beg <- match(species_allergen$Flowering_beg,
                                 month.name) 
species_allergen$fl.end <- match(species_allergen$Flowering_end,
                                 month.name) 
species_allergen$fl.period <- 
  species_allergen$fl.end - species_allergen$fl.beg

# Correct missing pheno data from Bioflor (Sept. 2020) :
# Taraxacum_sect._Ruderalia flowers from months 3 to 10 (= 8 months)
species_allergen$fl.end [which(species_allergen$Species =="Taraxacum_sect._Ruderalia")] <- 10
species_allergen$fl.period [which(species_allergen$Species =="Taraxacum_sect._Ruderalia")] <- 8

# Equisetaceae: no pollen produced
species_allergen$fl.beg[grep("Equisetum",species_allergen$Species)] <- 0
species_allergen$fl.end[grep("Equisetum",species_allergen$Species)] <- 0
species_allergen$fl.period[grep("Equisetum",species_allergen$Species)] <- 0


## Correct Introduction Status ####
## Manually correct floristic status from BioFlor (Sept. 2020): 
# "Koeleria_pyramidata","Medicago_falcata","Papaver","Tragopogon"     
species_allergen[is.na(species_allergen$Introduction_status_Seitz2012),]

species_allergen["Koeleria_pyramidata","Introduction_status_Seitz2012" ] <- "I"
species_allergen["Medicago_falcata","Introduction_status_Seitz2012" ] <- "I"

# Most papaver in Bioflor are archaeophytes except Papaver alpinum (I) and Papaver somniferum (N)
species_allergen["Papaver","Introduction_status_Seitz2012" ] <- "A" 

# All Tragopogon species are indigenous in Bioflor except Tragopogon porrifolius (N)
species_allergen["Tragopogon","Introduction_status_Seitz2012" ] <- "I"  

# Correct numerical status: 
species_allergen$status_num <- NA
species_allergen$status_num[which(species_allergen$Introduction_status_Seitz2012 == "I")] <- 1
species_allergen$status_num[which(species_allergen$Introduction_status_Seitz2012 == "A")] <- 2
species_allergen$status_num[which(species_allergen$Introduction_status_Seitz2012 == "N")] <- 3


# Create groups of species by name
natives <- rownames(species_allergen)[which(species_allergen$Introduction_status_Seitz2012 == "I")]
archaeophytes <- rownames(species_allergen)[which(species_allergen$Introduction_status_Seitz2012  == "A")]
neophytes <- rownames(species_allergen)[which(species_allergen$Introduction_status_Seitz2012 == "N")]

exotics <-  union(neophytes, archaeophytes)
residents  <-  union(natives, archaeophytes)

# Vector or allergenic species (n = 74)
allergenics <- species_allergen[species_allergen$allergenicity == 1 , "Species"]

# merge natives with archaeophytes as "residents"
species_allergen$resident <- 0
species_allergen$resident[which(species_allergen$Introduction_status_Seitz2012 %in% c("I","A"))] <- 1

species_allergen$native <- 0
species_allergen$native[which(species_allergen$Introduction_status_Seitz2012 %in% c("I"))] <- 1

species_allergen$exotic <- 0
species_allergen$exotic[which(species_allergen$Introduction_status_Seitz2012 %in% c("A", "N"))] <- 1

species_allergen$neophyte <- 0
species_allergen$neophyte[which(species_allergen$Introduction_status_Seitz2012 %in% c("N"))] <- 1

# clean up objects
rm(species_flowering, species_allergen_severity)