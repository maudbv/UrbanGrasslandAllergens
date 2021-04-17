## ALLERGENIC SPECIES in Molecule FAMILY space ordination

library(vegan)
library(ggplot2)
library(inlmisc)

# # Heatmap of species in ALLFAM: ####
# y <- allfam_mat[allergenics,]
# y <- y[,-which(colnames(y) == "Unknown")]
# y <- y[rowSums(y)>0,]
# heatmap(x = as.matrix(y),Colv = NA,
#         margins =c(11,5), cexRow = 0.5, cexCol = 0.6 )

# Prepare oridnation of allergenic species in molecule family space  ####

# allergen biochemical families per species:
tmp <- allfam_mat 

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
                              "Introduction_status_Seitz2012",
                              "native",
                              "anemophilous",
                              "family")]

# Make sure status is numeric:
expl.data$status_num <- as.numeric(expl.data$status_num)

# Remove potentially missing data: 
# 
expl.data <- na.omit(expl.data)

# Make sure the order of rownames is the same: 
tmp <- tmp[rownames(expl.data),]
stopifnot(nrow(tmp) == nrow(expl.data))

# Calculate jaccard dissimilarities between species : 
jacc.mol.fam <- vegdist(tmp, method = "jaccard")

### PERMANOVA: adonis ######

## TEST EFFECT OF ALLERGEN SCORE : SIGNIF R2 = 0.083
# (adonis.mol <-  adonis2(
#   jacc.mol.fam ~  allergen_score,
#   data = expl.data,
#   sqrt.dist = TRUE,
#   permutations = 1000
# ))

## # EFFECT OF wind pollination # Signif R2 = 0.074
# (adonis.mol <-  adonis2(
#   jacc.mol.fam ~  anemophilous,
#   data = expl.data,
#   sqrt.dist = TRUE,
#   permutations = 1000
# ))

## Significant differences by status: signif marginal (P= 0.0340, R2 = 0.04)
(adonis.mol.status <- adonis2(
  jacc.mol.fam ~  as.numeric(status_num),
  data = expl.data,
  sqrt.dist = TRUE,
  permutations = 1000
))

## SIgnificant differences between families : 
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


# NMDS: representation ########

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

col.status <- GetColors(n = 3,alpha = 0.8,stops = c(0.4,0.8))


### GRAPH NMDS with STATUS + FAMILY:

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
 
soft.col.status <- GetColors(n = 3,alpha = 0.1,stops = c(0.4,0.8))
with(expl.data, ordiellipse(nmds.molfam.species, status_num,
                            kind = "ehull", label = FALSE,
                            draw= "polygon",
                            col =col.status ,alpha = 0.1,
                            border =soft.col.status ))

points(jitter(site.scores[,1],factor = 20),
       jitter(site.scores[,2],factor = 20),
       cex = 0.7,
       pch = 20,
       col = col.status[expl.data$status_num])

fam.spid <- with(expl.data, ordihull(nmds.molfam.species, family,
                          label = TRUE,
                          lwd =0, cex = 0.7))

legend("topleft",
       legend= c("Native","Archaeophyte","Neophyte"),
       fill = col.status, border= col.status,  bty = "n",
       cex = 0.6)
       
# represent families: 
col.fam <-  GetColors(n = 7, scheme = "smooth rainbow", alpha = NULL, 
 start = 0.1, end = 1, bias = 1, reverse = FALSE, blind = NULL, gray = FALSE)
names(col.fam) <- unique(expl.data$family)

# Add stress data
mtext(1, text = paste("stress =", round( nmds.molfam.species$stress,3)),
      font = 3, cex  = 0.6, adj = 0.99, line = -1,col = "darkgrey")


# ### dbRDA:  Distance-based redundancy analysis = SAME AS adonis #####
# anova(dbrda.mol <- dbrda(tmp ~  as.numeric(status_num),
#                          distance = "jaccard",
#                          sqrt.dist = TRUE,
#                          data = expl.data))
# 
# 
# anova(dbrda.mol.2 <- dbrda(tmp ~   as.numeric(status_num) + as.factor(family),
#                            distance = "jaccard",
#                            sqrt.dist = TRUE,
#                            data = expl.data), by = "margin")
# 
# anova(dbrda.mol, dbrda.mol.2,  permutations = 1000) 
# anova(dbrda.mol, permutations = 1000)
# anova(dbrda.mol.2, permutations = 1000)
# anova(dbrda.mol, by="axis", permutations = 1000)
# anova(dbrda.mol.2, by="axis", permutations = 1000)
# anova(dbrda.mol.2, by="margin", permutations = 1000)
# 
# dbrda.terms <- anova(dbrda.mol.2, by="margin", permutations = 1000)
# dbrda.terms$R2 <- dbrda.terms$SumOfSqs/sum(dbrda.terms$SumOfSqs)
# 
