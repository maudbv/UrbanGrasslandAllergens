# PERMANOVA on strictly defined molecular data

# allergen biochemical families per species, according to strict species matching, not genus- or family-level matching.
tmp <- strict_allfam_mat 

# remove column of "Unknown" allergen family
# to avoid meaningless grouping of species with missing data
tmp <- tmp[,- which(colnames(tmp) == "Unknown") ] 

# Remove species which do not have any known allergen families
# These are both non-allergenics and allergenics without known molecule
tmp <- tmp[-which(rowSums(tmp) == 0), ]

# Create dataframe with species variables:
expl.data= species_allergen[rownames(tmp),
                            c("allergenicity", 
                              "allergen_score" ,
                              "status_num",
                              "family")]

# Make sure introduction status is a numeric:
expl.data$status_num <- as.numeric(expl.data$status_num)

# Remove potentially missing data: 
expl.data <- na.omit(expl.data)

# Make sure the order of rownames is the same: 
tmp <- tmp[rownames(expl.data),]
stopifnot(nrow(tmp) == nrow(expl.data))

# Calculate jaccard dissimilarities between species : 
jacc.mol.fam <- vegdist(tmp, method = "jaccard")

### PERMANOVA: function adonis from package vegan ######

## Significant differences between families : R = 0.56 , P < 0.0001
(adonis.mol.fam <-  adonis2(
  jacc.mol.fam ~  as.factor(family),
  data = expl.data,
  sqrt.dist = TRUE,
  permutations = 10000,
  by = "margin"
))
# Very Signif: Family R2 = 0.56 , P < 0.0001


## TEST EFFECT OF ALLERGEN SCORE + status num :
# Status only marginally signif once family taken into account
(adonis.mol.2 <-  adonis2(
  jacc.mol.fam ~  as.factor(family) + status_num ,
  data = expl.data,
  sqrt.dist = TRUE,
  by = "margin",
  permutations = 10000
)) 


# return data
write.csv( rbind(
                 as.matrix(adonis.mol.fam), "",
                 as.matrix(adonis.mol.2)),
           "results/PERMANOVA STRICT results in allergen family space.csv")

# Figure S2A: allergen families by floristic status ####
png(width = 14, height = 15, unit ="cm", res = 600,
    file = "results/figure S2A_STRICT.png")

par (
  mfrow = c(1, 1),
  bty = "l",
  mar = c(4, 16, 2, 2),
  oma = c(0, 0, 0, 0),
  mgp = c(2,0.4,0),
  tcl = - 0.2,
  las = 1,
  cex.axis = 0.9,
  xpd = TRUE
)
col.line <- c(GetColors(1, alpha = 1,stop = c(0.4,0.41),
                        scheme = "iridescent"),
              GetColors(1, alpha = 1,stop = c(0.65,0.66),
                        scheme = "iridescent"),
              GetColors(1, alpha = 1,stop = c(0.8,0.81),
                        scheme = "iridescent"))

# shared allergen families

tmp<- strict_allfam_mat[,c(29:26, 1, 25:2 )]

barplot(rbind(colSums(tmp[natives,]),
              colSums(tmp[archaeophytes,]),
              colSums(tmp[neophytes,])),
        beside = T, las = 2, horiz = T, cex.names = 0.8,
        col = col.line, border = col.line)

mtext(1, text = "Number of species", line = 1.5, cex = 1)

legend("bottomright", 
       bty = "n", inset = c(-0.05,0.08),
       legend = c( "neophytes", "archaeophytes","natives"),
       fill = col.line[3:1], border = col.line[3:1],
       cex = 0.6)

## AF118: Group 5 ragweed allergen is only one unique to exotics (for now)
# many unknown family groups make it difficult to judge.

dev.off()
# Figure S2B: Allergen family accumulation curve #####
png(width = 11, height = 12, unit ="cm", res = 600,
    file = "results/figure S2B_STRICT.png")

par (
  mfrow = c(1, 1),
  bty = "l",
  mar = c(4, 3, 2, 1),
  oma = c(0, 0, 0, 0),
  mgp = c(2,0.4,0),
  tcl = - 0.2,
  las = 1,
  cex.axis = 0.8,
  xpd = TRUE
)

col.area<- c(GetColors(1, alpha = 0.6,stop = c(0.4,0.41),
                       scheme = "iridescent"),
             GetColors(1, alpha = 0.6,stop = c(0.65,0.66),
                       scheme = "iridescent"),
             GetColors(1, alpha = 0.6,stop = c(0.8,0.81),
                       scheme = "iridescent"))
col.line <- c(GetColors(1, alpha = 1,stop = c(0.4,0.41),
                        scheme = "iridescent"),
              GetColors(1, alpha = 1,stop = c(0.65,0.66),
                        scheme = "iridescent"),
              GetColors(1, alpha = 1,stop = c(0.8,0.81),
                        scheme = "iridescent"))

rare.est <- data.frame(matrix(NA,nrow = 3, ncol = 4,
                              dimnames = list(1:3,
                                              c("type", "SR", "mean","sd"))))

for (i in 1:3) {
  pool <- list(natives, archaeophytes, neophytes)
  type <- c("native", "archaeophyte", "neophyte")[i]
  leg <- c("Natives","Archaeophytes","Neophytes")[i]
  sac <- specaccum(strict_allfam_mat[intersect(pool[[i]],allergenics),], 
                   method = "random", permutations = 10000)
  
  print(paste("For 10", type ,"species:",
              round(sac$richness[10],2),"(sd =",
              round(sac$sd[10],2),")", "families"))
  
  rare.est[i,] <- c(type,10,round(sac$richness[10],2),round(sac$sd[10],2))
  
  if( i == 1) {
    plot(sac$sites,sac$richness, type = "l", lwd = 2,
         ylim = c(0,45),col = col.line[i],
         xlab = "Species", ylab = "Allergen Families")
  } else {
    lines(sac$sites,sac$richness, type = "l", lwd = 2, 
          col = col.line[i])
  }
  
  lo <- smooth.spline(sac$sites, apply(sac$perm, 1, quantile, 0.025),
                      spar = 0.7)
  hi <-smooth.spline(sac$sites, apply(sac$perm, 1, quantile, 0.975),
                     spar = 0.7)
  polygon(rbind(data.frame(predict(lo)),
                data.frame(predict(hi))[length(hi$x):1,]),
          col = col.area[i], lty = 0)
  
  # # add prediction: no longer enough data points for that
  # sac <- specaccum(strict_allfam_mat[intersect(pool[[i]],allergenics),], 
  #                  method = "exact")
  # mod1 <- fitspecaccum(sac, model = "lomolino", method = "random")
  # pred.model <- predict(mod1, newdata = 1:50)
  # lines(pred.model,
  #       col=col.line[i], lty=2)
  # if (i ==1) y.pos<- max(pred.model) -3.5 else y.pos<- max(pred.model)+1
  # text(x = 52, y =y.pos, label = leg, 
  #      cex = 0.7, col=col.line[i],adj = 1)
}

abline(v = 10, lty = 3, xpd = FALSE)
print(rare.est)
# 
# legend("bottomright", legend = c("Natives","Archaeophytes","Neophytes"),
#        fill= col.line, bty ="n",
#        cex = 0.7)


dev.off()

# we loose half of the neophyte allergen molecules in the strict version. The difference with broad version concerns the following neophytes:
# Plantago_arenaria                           
# Senecio_inaequidens                        
# Senecio_vernalis
# extrapolating molecules at genus level for Senecio, Solidago and Plantago
# + Ambrosia coronopifolia also gets all the ambrosia allergens.

# Figure S2C: NMDS of species in allergen family space ####
png(width = 17, height = 13, unit ="cm", res = 600,
    file = "results/figure S2C_STRICT.png")

# allergen biochemical families per species:
tmp <- strict_allfam_mat 
# remove column of "Unknown" allergen family
# to avoid meaningless grouping of species with missing data
tmp <- tmp[,- which(colnames(tmp) == "Unknown") ] 
# Remove species which do not have any known allergen families
# These are both non-allergenics and allergenics without known molecule
tmp <- tmp[-which(rowSums(tmp) == 0), ]

# Create dataframe with species variables:
expl.data= species_allergen[rownames(tmp),
                            c("status_num",
                              "Introduction_status_Seitz2012",
                              "family")]

# Make sure status is numeric:
expl.data$status_num <- as.numeric(expl.data$status_num)

# Remove potentially missing data: 
expl.data <- na.omit(expl.data)

# Make sure the order of rownames is the same: 
tmp <- tmp[rownames(expl.data),]
stopifnot(nrow(tmp) == nrow(expl.data))

# Run NMDS
(nmds.molfam.species<- metaMDS(tmp ,
                               k = 2,
                               try = 200,
                               distance = "jaccard",
                               autotransform = TRUE, 
                               trymax = 40,
                               trace = FALSE
))

mol.scores <- as.data.frame(nmds.molfam.species$species)
mol.scores$names <- sapply(rownames(mol.scores),
                           FUN = function(x){
                             strsplit(x ,":")[[1]][1]
                           })
site.scores <- as.data.frame(nmds.molfam.species$points)

## GRAPH NMDS with STATUS + FAMILY:

par (
  mfrow = c(1, 1),
  bty = "l",
  mar = c(3, 3, 2, 1),
  oma = c(0, 0, 0, 0),
  mgp = c(2,0.4,0),
  tcl = - 0.2,
  las = 1,
  cex.axis = 0.8,
  xpd = TRUE
)
plot(site.scores[,c(1,2)],
     ylim = c(-3.5,3),
     xlim =  c(-3,4),
     type = "n")
arrows(0,0,
       mol.scores$MDS1*1.3,mol.scores$MDS2*1.3,
       col = "lightgrey",length = 0.1 )


# Group labels for duplicate molecule scores:
mol.scores$label <- NA
for (i in 1:nrow(mol.scores)) {
  x <- mol.scores[i,]
  n = which(round(mol.scores[,1],2) == round(x[1,1],2) &
              round(mol.scores[,2],2) == round(x[1,2],2))
  mol.scores$label[i] <- paste(mol.scores$names[n], collapse  = "\n")
}
redux.scores <- mol.scores[- which(duplicated(mol.scores$label)),]
redux.scores$lab.x <- redux.scores$MDS1 + 0.05*sign(redux.scores$MDS1)

# identify the position best for each label on the plot:
# id<- identify(mol.scores$MDS1*1.3,mol.scores$MDS2*1.3,pos = TRUE,
#          label = mol.scores$names, cex = 0.5,col = "grey")

text(redux.scores[,1:2]*1.3,
     label = redux.scores$label,
     pos = c(4,4,1,1,4,4,1,4,2,3,3,2,4,4,4,2,2,2,3),
     cex = 0.5,
     col = "grey")

# Choose colors for florsitic status: 
# GetColors(n = 3,alpha = 0.1,stops = c(0.4,0.8))
soft.col.status <- c(GetColors(1, alpha = 0.6,stop = c(0.4,0.41),
                               scheme = "iridescent"),
                     GetColors(1, alpha = 0.6,stop = c(0.65,0.66),
                               scheme = "iridescent"),
                     GetColors(1, alpha = 0.6,stop = c(0.8,0.81),
                               scheme = "iridescent"))
col.status <- c(GetColors(1, alpha = 1,stop = c(0.4,0.41),
                          scheme = "iridescent"),
                GetColors(1, alpha = 1,stop = c(0.65,0.66),
                          scheme = "iridescent"),
                GetColors(1, alpha = 1,stop = c(0.8,0.81),
                          scheme = "iridescent"))

status.ellipse <- with(expl.data, 
                       ordiellipse(nmds.molfam.species, status_num,
                                   kind = "ehull",
                                   label = FALSE,
                                   draw= "polygon",
                                   col =col.status ,alpha = 0.1,
                                   border =soft.col.status ))

# points(jitter(site.scores[,1],factor = 20),
#        jitter(site.scores[,2],factor = 20),
#        cex = 0.7,
#        pch = 20,
#        col = col.status[expl.data$status_num])

# represent families: 
col.fam <-  GetColors(n = 7, scheme = "smooth rainbow",
                      alpha = NULL, 
                      start = 0.55, end = 1, bias = 1,
                      reverse = FALSE, 
                      gray = FALSE)
names(col.fam) <- unique(expl.data$family)

points(jitter(site.scores[,1],factor = 20),
       jitter(site.scores[,2],factor = 20),
       cex = 0.7,
       pch = 20,
       col = col.fam[expl.data$family])

# family centroids: 
fam.spid <- as.data.frame(t(summary(
  with(expl.data, ordihull(nmds.molfam.species,
                           family,
                           label =FALSE,
                           lwd =0,))
)))

text(x = fam.spid$NMDS1 +0.1 ,
     y = jitter(fam.spid$NMDS2,factor = 1),
     adj = 0,
     labels  = rownames(fam.spid),
     # col = col.fam[rownames(fam.spid)],
     cex = 0.6)

# Legend
legend(x = -3, y =3.2, xjust = 0, yjust = 1,
       legend= c("Native","Archaeophyte","Neophyte"),
       fill = col.status, border= col.status,  bty = "n",
       cex = 0.6)

# Add stress data
mtext(1, text = paste("stress =", round( nmds.molfam.species$stress,3)),
      font = 3, cex  = 0.8, adj = 0.99, line = -1,col = "darkgrey")

dev.off()



