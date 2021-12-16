# focus on severe allergens

# list "severe" ( non-native allergen
species_allergen[ species_allergen$severe.all ==1,]
hist(species_allergen[ species_allergen$PAV>0, "PAV"])

hist(species_allergen[ species_allergen$severe.all ==1 & species_allergen$PAV>0, "PAV"], 
     breaks  = 20)

# Create a threshold for PAV ###
thresh = 10
species_allergen[  species_allergen$PAV>thresh,]


# recalculate severity based on PAV
# Number of SEVERELY allergenic species per plot
allergen_summary$thresh.all.num <- sapply(1:nrow(vegcomm), FUN = function(l) {
  p <- vegcomm[l, ]
  sp <-  names(p)[p > 0]
  result <- sum(species_allergen[species_allergen$Species %in% sp , "PAV" ]>thresh, na.rm = T)
  return(result)
})


#  Total cover of SEVERELY allergenic species per plot
allergen_summary$cover.thresh.all<- sapply(1:nrow(vegcomm), FUN = function(l) {
  p <- vegcomm[l, ]
  sp <-  names(p)[p > 0]
  p <- p[sp]
  all <- as.numeric(species_allergen[species_allergen$Species %in% sp , "PAV" ]> thresh)
  result = sum(p*all)
  return(result)
})



# relationship with total allergenics

plot(jitter(thresh.all.num) ~ jitter(all.num), allergen_summary)
add.stats(glm(thresh.all.num ~ all.num, allergen_summary, family = poisson))

plot(jitter(severe.all.num) ~ jitter(all.num), allergen_summary)
add.stats(glm(severe.all.num ~ all.num, allergen_summary, family = poisson))

plot(jitter(cover.thresh.all) ~ jitter(cover.all), allergen_summary)
add.stats(lm(cover.thresh.all ~ cover.all, allergen_summary))

# trends in cover and richness
par (
  mfrow = c(4, 2),
  mar = c(2, 1, 2, 1),
  oma = c(2, 5, 1, 1),
  mgp = c(2,0.5,0),
  tcl = - 0.3,
  las = 1
)
m = 1

# graph
for (k in 1:4) {
  tmp <- allergen_summary
  y <- c("cover.severe.all","cover.thresh.all", "severe.all.num","thresh.all.num")[k]
  mod <- c("n","n", "p", "p")[k]
  leg <- c("Severe\nAllergenic cover ",
           "Threshold\nAllergenic cover ",
           "Severe\nAllergenic richness",
           "Threshold\nAllergenic richness")[k]
  tmp$y <- tmp[, y]
  
  for (i in 1:2) {
    xleg <- c("% Impervious surfaces",
              "Proportion of neophytes")[i]
    x <- c("Seal_500", "prop.neo")[i]
    tmp$x <- tmp[, x]
    
    plot(y ~ x,
         data = tmp,
         pch = 20, 
         col = "grey",
         ann = FALSE
    )
    
    if(mod == "n"){ 
      add.stats(f <- lm(y ~ x,
                        data = tmp), type = "lm",l.col = "black")
    }
    if(mod == "p"){ 
      add.stats(f <- glm(y ~ x,
                         data = tmp,
                         family = poisson),
                type = "glm",
                l.col = "black")
    }
    # add.stats(f <- cor.test(~ y + x,
    #                         data = tmp),
    #           adj.stats = 0, col.stats = "grey")
    
    mtext(3,text = paste(letters[m], ")",sep = ""),
          adj = 0, font = 3, cex = 0.7)
    m = m + 1
    
    # x label
    if (k ==3) {
      mtext(1, text = xleg, cex = 0.75, line = 2.2)
    }
    
    # y label
    if (i == 1) {
      mtext(
        2,
        text = leg,
        outer = FALSE,
        cex = 0.75,
        line = 2.5, las = 0
      )
    }
  }
}
