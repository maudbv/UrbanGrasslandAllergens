# identify high non native allergenic cover ####
hi_exo_all <- rownames(allergen_summary)[which(allergen_summary$cover.exo.all>20)]
plot_summary[hi_exo_all,]
allergen_summary[hi_exo_all,]

# what neophytes are leading this?
vegcomm[hi_exo_all,intersect(allergenics,exotics)]
sort(colSums(vegcomm[hi_exo_all,intersect(allergenics,exotics)]))

## Main neophyte species with high cover in these 4 plots:
# Medicago x varia,
# Arrhenatherum elatius,  => neophyte to Berlin, but often considered indigenous to Europe/germany
# Solidago_canadensis
# Plantago_arenaria

## Archaeotphytes:
# Digitaria_ischaemum
# Setaria_viridis
# Plantago_lanceolata





# effect of ignoring plots on fig 1: we loose correlation of cover with urbanisation, but increase in richness of non native allergenics remains marginally significant. Relationships with proportion of neophyte remain unchanged.

par (
  mfrow = c(3, 2),
  mar = c(2, 1, 2, 1),
  oma = c(2, 5, 1, 1),
  mgp = c(2,0.5,0),
  tcl = - 0.3,
  las = 1
)
m = 1

# graph
for (k in 1:3) {
  tmp <- allergen_summary[hi_exo_all,]
  y <- c("cover.all", "cover.exo.all", "all.num.exo")[k]
  mod <- c("n","n","p")[k]
  leg <- c("Allergenic\ncumulative cover ",
           "Non-native allergenic\ncumulative cover",
           "Non-native allergenic\nspecies richness")[k]
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

# Identifying high PAV plots and neophytes
allergen_summary[hi_exo_all,"CWM.neo.pav"] # => not very high except for one plot with A. elatius
allergen_summary[hi_exo_all,"mean.neo.pav"] 
