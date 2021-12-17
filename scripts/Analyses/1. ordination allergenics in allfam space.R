## PERMUTATIONAL MANOVA of ALLERGENIC SPECIES in ALLERGEN MOLECULE FAMILY space
library(vegan)
library(ggplot2)
library(inlmisc)

(function() {
# Prepare ordination of allergenic species in molecule family space  ####

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
expl.data <- na.omit(expl.data)

# Make sure the order of rownames is the same: 
tmp <- tmp[rownames(expl.data),]
stopifnot(nrow(tmp) == nrow(expl.data))

# Calculate jaccard dissimilarities between species : 
jacc.mol.fam <- vegdist(tmp, method = "jaccard")

### PERMANOVA: function adonis from package vegan ######

## Significant differences by introduction status: signif (P= 0.03, R2 = 0.04)
(adonis.mol.status <- adonis2(
  jacc.mol.fam ~  as.numeric(status_num),
  data = expl.data,
  sqrt.dist = TRUE,
  permutations = 1000
))

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
write.csv( rbind(as.matrix(adonis.mol.status), "", 
                 as.matrix(adonis.mol.fam), "",
                 as.matrix(adonis.mol.2)),
           "results/PERMANOVA results in allergen family space.csv")



# ## NMDS: representation of the permanova effect with families and status ########
# ## Final figure for paper is redrawn in the "scripts/illustration/figures manuscript.R" script
# 
# (nmds.molfam.species<- metaMDS(tmp ,
#                                k = 2,
#                                try = 200,
#                                distance = "jaccard",
#                                autotransform = TRUE,
#                                trymax = 40,
#                                trace = FALSE
#                                 ))
# 
# mol.scores <- as.data.frame(nmds.molfam.species$species)
# mol.scores$names <- sapply(rownames(mol.scores),
#                            FUN = function(x){
#                              strsplit(x ,":")[[1]][1]
#                            })
# site.scores <- as.data.frame(nmds.molfam.species$points)
# 
# col.status <- GetColors(n = 3,alpha = 0.8,stops = c(0.4,0.8))
# 
# 
# ### GRAPH NMDS with STATUS + FAMILY:
# 
# plot(site.scores[,c(1,2)],
#        ylim = c(-3.5,3),
#        xlim =  c(-3,4),
#         type = "n")
# arrows(0,0,
#        mol.scores$MDS1*1.3,mol.scores$MDS2*1.3,
#        col = "lightgrey",length = 0.1 )
# 
# 
# # Group labels for duplicate molecule scores:
# mol.scores$label <- NA
# for (i in 1:nrow(mol.scores)) {
#   x <- mol.scores[i,]
#   n = which(round(mol.scores[,1],2) == round(x[1,1],2) &
#           round(mol.scores[,2],2) == round(x[1,2],2))
#   mol.scores$label[i] <- paste(mol.scores$names[n], collapse  = "\n")
# }
# redux.scores <- mol.scores[- which(duplicated(mol.scores$label)),]
# redux.scores$lab.x <- redux.scores$MDS1 + 0.05*sign(redux.scores$MDS1)
# 
# # identify the position best for each label on the plot:
# # id<- identify(mol.scores$MDS1*1.3,mol.scores$MDS2*1.3,pos = TRUE,
# #          label = mol.scores$names, cex = 0.5,col = "grey")
# 
# text(redux.scores[,1:2]*1.3,
#      label = redux.scores$label,
#      pos = c(4,4,1,1,4,4,1,4,2,3,3,2,4,4,4,2,2,2,3),
#      cex = 0.5,
#      col = "grey")
# 
# soft.col.status <- GetColors(n = 3,alpha = 0.1,stops = c(0.4,0.8))
# with(expl.data, ordiellipse(nmds.molfam.species, status_num,
#                             kind = "ehull", label = FALSE,
#                             draw= "polygon",
#                             col =col.status ,alpha = 0.1,
#                             border =soft.col.status ))
# 
# points(jitter(site.scores[,1],factor = 20),
#        jitter(site.scores[,2],factor = 20),
#        cex = 0.7,
#        pch = 20,
#        col = col.status[expl.data$status_num])
# 
# fam.spid <- with(expl.data, ordihull(nmds.molfam.species, family,
#                           label = TRUE,
#                           lwd =0, cex = 0.7))
# 
# legend("topleft",
#        legend= c("Native","Archaeophyte","Neophyte"),
#        fill = col.status, border= col.status,  bty = "n",
#        cex = 0.6)
# 
# # represent families:
# col.fam <-  GetColors(n = 7, scheme = "smooth rainbow", alpha = NULL,
#  start = 0.1, end = 1, bias = 1, reverse = FALSE, blind = NULL, gray = FALSE)
# names(col.fam) <- unique(expl.data$family)
# 
# # Add stress data
# mtext(1, text = paste("stress =", round( nmds.molfam.species$stress,3)),
#       font = 3, cex  = 0.6, adj = 0.99, line = -1,col = "darkgrey")
# 
})()