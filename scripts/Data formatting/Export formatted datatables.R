# Export clean data tables

## Clean allergen data tables for re-analyses ####
write.csv( species_allergen ,file = "clean data/species_allergenicity.csv")
write.csv( molecule_data,file = "clean data/allergen_molecule_data.csv")
write.csv( molecule_sum,file = "clean data/allergen_molecule_summary.csv")
write.csv( mol_mat,file = "clean data/allergen_species_by_molecule_matrix.csv")

# Simplified data tables in the "results" folder: ####
# Species allergen and trait data 

# clean up the columns
tmp <- species_allergen
tmp$tribe <- species_data[match(species_allergen$Species.match, species_data$Species.simple),
             "tribe"]
tmp <- tmp[, c(
  "Species_name","Species_acceptedname","family","tribe",
  
  # Frequency in Berlin grasslands:
  "Frequency","Freq.rural","Freq.suburb", "Freq.urban",
  
  # Floristic status: 
  "Introduction_status_Seitz2012","first_evidence" ,
  
  # Phenology
  "Flowering_beg" , "Flowering_end" ,"fl.period", "phenology_score" ,
  
  #Pollination
  "main.pollination","flw_muell", "pollen_score",  "Source",

  # Allergenicity:
  "allergenicity_strict","allergenicity","allergen_score",
  
  # PAV
  "PAV",
  
  # allergenicity references
  "IUIS.record","SDAP.record","Allergenonline.record", 
  "Allergome.record",
  "AllFam.record" , "additional_ref","react.europe.Allergome", "FDA",
  "ref_severity", "Comment_severity",
  
  # Molecules: 
  "molecule.unknown"
  )
  ]

colnames(tmp) <- c("Species_name", "Accepted_name","Family", "Tribe",
 "Frequency","Freq_rural","Freq_suburb","Freq_urban",
 "Status_Berlin","first_evidence",
 "Flowering_beg", "Flowering_end","fl_period", "phenology_score",
 "main_pollination", "Mueller_flower_type", "pollen_score","trait_source",
 "allergenicity_strict","allergenicity","allergen_score",
 "PAV",
 "IUIS","SDAP","Allergenonline","Allergome",
 "AllFam", "additional_ref",
 "react_europe_Allergome", "FDA_allergen_code",
 "ref_severity", "Comment_severity",
 "molecule_unknown"
)

# Export csv: 
write.csv( tmp,  na = "" ,row.names = FALSE,
           file = "results/species_allergenicity_data.csv")

# allergen molecule data:
tmp  <- molecule_data[, c( "Species","genus",
  "Allergen.simple","Allergen.genus","Allergen.Family",
  "Databases","Ref.manual","Mol.alt","Matched.names", "Matches.with.database"
  )]

colnames(tmp) <- c("Species","Genus",
  "Allergen.name","Allergen.genus.name","Allergen.family",
  "Databases","Manual.reference",
  "Alternative.molecule.names","Matched.species.names", "Match.types"
) 

tmp[tmp == "NA"] <- ""

# Export csv: 
write.csv( tmp, na = "" ,row.names = FALSE,
           file = "results/molecule_data.csv")
rm(tmp)

# Table of subset of the most allergenic species

# clean up the columns
tmp<- species_allergen

# Collate status data: 
tmp$Status <- as.character(tmp$Introduction_status_Seitz2012)
tmp$Status[tmp$Status == "N"] <- paste(tmp$Status[tmp$Status == "N"],
                                       " (",
                                       tmp$first_evidence[tmp$Status == "N"],
                                       ")",
                                       sep = ""
)

# Pretty "flowering" column
tmp$Flowering <- paste(substr(tmp$Flowering_beg, 1,3),
                       substr(tmp$Flowering_end,1,3),
                       sep = "-")

# Pretty "frequencies" column
tmp$Frequencies <- paste(tmp$Freq.rural,tmp$Freq.suburb,tmp$Freq.urban,
                         sep = "_")


tmp <- do.call(rbind,
               sapply( c("N","A","I"), 
                       function(i) {
                         x <-  tmp[tmp$Introduction_status_Seitz2012 == i,]
                         x <- x[order(x$PAV, decreasing = TRUE),][1:10,]
                         x <- x[!x$PAV == 0,]
                         return(x)
                       },
                       simplify = FALSE))

tmp <- tmp[, c(
  "Species_name","family",
  # Floristic status: 
  "Status" ,
  # Frequency in Berlin grasslands:
  "Frequencies",
  # Phenology
  "Flowering",
  #Pollination
  "main.pollination", 
  # Allergenicity:
  "allergen_score",
  # PAV
  "PAV"
)]


# Export csv: 
write.csv( tmp,  na = "" ,row.names = FALSE,
           file = "results/Table_5_select_species_allergen.csv")
rm(tmp)



