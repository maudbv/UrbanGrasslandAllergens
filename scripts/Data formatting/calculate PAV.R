# Calculate PAV = Potential Allergenic Value
# For each species 

# Flowering period: phenology score  ####
species_allergen$phenology_score <- (function(){
  
tmp <- species_allergen

# 1 month, 2 months, over 2 months
tmp$phenology_score <- NA
tmp$phenology_score [which(species_allergen$fl.period ==1)] <- 1
tmp$phenology_score [which(species_allergen$fl.period ==2)] <- 2
tmp$phenology_score [which(species_allergen$fl.period >2)] <- 3

# Bioflor: Taraxacum_sect._Ruderalia flowers from months 3 to 10 (= 8 months)
tmp$phenology_score [which(species_allergen$Species =="Taraxacum_sect._Ruderalia")] <- 3

# Equisetaceae: no pollen produced
tmp$phenology_score[grep("Equisetum",species_allergen$Species)] <- 0

return(tmp$phenology_score)
})()

# Exctract main pollination syndrome:  ####
species_allergen$main.pollination<- (function(){
  
tmp <- species_allergen

tmp$main.pollination <- sapply(
  tmp$poll_vect_B, 
  FUN = function(x) {
    l <- data.frame(
      vector = as.character(strsplit(x, split = " - ")[[1]]),
                  stringsAsFactors = FALSE)
   
  # Create the "type" score corresponding to syndrome frequency: 
  # if (length(l) <= 1) {  (# What is this condition for ???)
  l$type <- 0 
  l$type[grep("always", l$vector)] <- 1
  l$type[grep("rule", l$vector)] <- 2
  l$type[grep("often", l$vector)] <- 3
  l$type[grep("possible", l$vector)] <- 4
  l$type[grep("rare", l$vector)] <- 5
  l$type[grep("at failure of outcrossing", l$vector)] <- 6
  l$type[grep("unknown", l$vector)] <- 3  # Average ranking for 'unknown'
  # }
  
  #  extract the main syndrome, which is the minimum "type score"
  main = l$vector[l$type == min(l$type)]
  
  # extract the second main syndrome
  second = l$vector[l$type %in% c((min(l$type)  +1):5)]
  
  # Prepare the output: 
  out = NA
  
  
  # If selfing is the main syndrome: 
  if (length(grep("selfing",main)) == 1){
    # if wind pollination is one of the two main syndromes
    if (length(grep("insect", second)) == 1 | length(grep("insect", main)) == 1 ) {
      out = "insects-selfing"
    } else {
    out = "selfing"
    }
  }
  
  # If cleistogamy is the main syndrome:
  if (length(grep("ogamy",main)) == 1){
    if (length(grep("insect", second)) == 1 | length(grep("insect", main)) == 1 ) {
      out = "insects-selfing"
    } else {
    out = "cleistogamy"
    }
  }
  
  # Insect pollination is the main syndrome: 
  if (length(grep("insect",main)) == 1) {
    # if wind pollination is one of the two main syndromes
    if (length(grep("wind", second)) == 1 | length(grep("wind", main)) == 1 ) {
      out = "wind-insects"
    } else {
      out = "insects"
    }
  }
  
  # Wind pollination is the main syndrome:
  if (length(grep("wind",main)) == 1) {
    # If there is also insect pollination (ie not a "pure" wind pollinated"): 
    if (length(grep("insect",main)) == 1 | length(grep("insect",second)) == 1 ) {
      out = "wind-insects"
    } else {
      out = "wind"
    }
  }
return(out) })

# Add info from column "pollen dispersal" for the missing data: 
missing_poll <- which(is.na(tmp$main.pollination))

tmp$main.pollination[missing_poll] <- sapply(
  tmp$Pollen_dispersal[missing_poll], 
  FUN = function(x) {
    l <- as.character(strsplit(x, split = ",")[[1]])
    
    # Prepare the output: 
    main = l[1]
    out = main
    
    # Wind pollination is the main syndrome but also insect pollination 
    # (ie not a "pure" wind pollinated"): 
      if (length(grep("wind",main)) == 1 & length(grep("insect",l)) == 1 ) {
        out = "wind-insects"
      }
    
    # If cleistogamy or geitonogamy or pseudocleistogamy is the main syndrome:
    if (length(grep("ogamy",main)) == 1){
      out = "cleistogamy"
    }
    return(out) })

tmp$main.pollination <- tolower(tmp$main.pollination)

# Correct manually one missing value of pollination syndrome: 
tmp["Rubus_sect._Rubus","main.pollination"] <- "insects"

return(tmp$main.pollination)
})()


# Calculate pollen score based on main pollination syndrome:  ####
species_allergen$pollen_score <- NA
species_allergen$pollen_score[species_allergen$main.pollination =="wind"] <- 3
species_allergen$pollen_score[species_allergen$main.pollination =="wind-insects"] <- 2
species_allergen$pollen_score[species_allergen$main.pollination %in% c("insects","insects-selfing" )] <- 1
species_allergen$pollen_score[species_allergen$main.pollination %in% c( "cleistogamy","selfing")] <- 0


# Correct some of the "selfing" species as mild pollen producers
# According to Mueller classification in Biolflor 2018:
species_allergen$pollen_score[
  intersect(which(species_allergen$pollen_score %in% c(0,1)),
            grep("W ",species_allergen$flw_muell))] <- 1

# Add a non-zero pollen score to an allergenic species which is mostly selfing, but does produce pollens for insects (hymenopters) (
species_allergen$pollen_score[grep("H", species_allergen$flw_muell)] <- 1

# Calculate PAV  ####
# PAV = pheno_score * pollination_score * allergenicity score
species_allergen$PAV <-  
  species_allergen$allergen_score *
  species_allergen$pollen_score *
  species_allergen$phenology_score

# PAV = 0 for the two Equisetum species
species_allergen$PAV[is.na(species_allergen$PAV) ] <- 0


# Correct anemophilous species: 
species_allergen$anemophilous <- 0
species_allergen$anemophilous[grep("wind",
                                   species_allergen$main.pollination)] <- 1
