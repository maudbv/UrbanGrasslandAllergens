## Ordination of PLOTS in Molecule space ordination

library(vegan)
library(ggplot2)

### DATA for ordination of plots in individual molecule space  ####

tmp <- mol_plot # allergen present per plot
tmp <- tmp[, -which(names(tmp) == "Unknown")]

expl.data= cbind(plot_summary[rownames(tmp),
                              c("Seal_500" , "prop.neo", "BNIs",
                                "SR","Rao")],
                 allergen_summary[rownames(tmp),
                                  c("all.num" , "all.num.neo",
                                    "all.num.arc", "all.num.nat")]
)

expl.data$Seal_cut <- cut_interval(expl.data$Seal_500,3 )
expl.data <- na.omit(expl.data)

tmp <- tmp[rownames(expl.data),]

# Check we have the exact same plot name rows
stopifnot(all(rownames(tmp) == rownames(expl.data)))


# Jaccard Beta-dissimilarities: 
jacc.mol.dist <- vegdist(tmp, method = "jaccard")


# # NMDS ####
# nmds.mol.plot <- metaMDS(tmp,distance = "jacc",
#                          k = 3,
#                          try = 200,
#                          trymax = 40)
# 
# 
# plot(nmds.mol.plot,
#      ylim = c(-1,1),
#      xlim = c(-2,2),
#      display="sites",
#      type = "t")
# text(nmds.mol.plot,
#      display="species",col = "red", cex = 0.5)

# # PERMANOVA in individual molecule space ####
# (adonis.mol.plot <- adonis2(
#   jacc.mol.dist ~  Seal_500 + prop.neo,
#   data = expl.data, sqrt.dist = TRUE))
# 
# betadiver(tmp)
# (betadisper(jacc.mol.dist,group = expl.data$Seal_cut))
# 
# 
# # With all factors combined: 
# (adonis.mol.plot <- adonis2(
#   jacc.mol.dist ~  SR + BNIs + Rao + Seal_500 + prop.neo,
#   data = expl.data, sqrt.dist = TRUE))
# 


# DBRDA = plots in molecule space - SEALING + PROP.NEO ####

dbrda.mol <-  dbrda(tmp ~  Seal_500 + prop.neo ,
                    distance = "jaccard",
                    sqrt.dist = TRUE,
                    data = expl.data)
dbrda.mol.two.jacc  <- dbrda.mol

# Extract stats
dbrda.model <- anova(dbrda.mol, by=NULL, perm.max=999)
dbrda.model$R2 <-dbrda.model$SumOfSqs/
    sum(dbrda.model$SumOfSqs)
dbrda.terms <- anova(dbrda.mol, by="margin", perm.max=500)
dbrda.terms$R2 <- dbrda.terms$SumOfSqs/sum(dbrda.terms$SumOfSqs)
dbrda.terms$names <- c("% Impervious","Prop.Neophytes", "resid")

# Graph:
quartz(width = 5, height = 4.5,
       file = "results/DBRDA plot in MOLECULE_only seal+PropNeo.pdf",
       type = "pdf")

par(mar = c(3,3,2,1))
ordi <- plot(dbrda.mol, type = "n",
             mgp = c(1.6,0.4,0), las = 1, tck = -0.01, 
             cex.axis = 0.8, )
points(dbrda.mol, type = "p" , pch = 20,  col = "#48876280")
#dbrda.mol, display = "cn", cex = 0.7,adj = 0)
labels = c("% Impervious", "Prop.Neo")

arrows(0,0, dbrda.mol$CCA$biplot[,1]*2, dbrda.mol$CCA$biplot[,2]*2,
       length = 0.1, angle =30, col = "darkgrey")

text(dbrda.mol$CCA$biplot[,1]*2.2,
     dbrda.mol$CCA$biplot[,2]*2.2,
     labels = c("% Impervious", "Prop.Neo"), cex = 0.8)

mtext(3, text= "b) Allergen molecule composition",
      cex = 1, font = 2, line =1)
mtext(3, adj = 0, line = 0.1, cex= 0.7,
      text = "Jaccard dissimilarities")

mtext(3, adj = 1, line = 0.1, cex= 0.7,
      text = as.expression(
        substitute("R"^2*" = "*a*b,
                   list(  a = round(dbrda.model[1,"R2"],3),
                          b = p2star(dbrda.model[1,"Pr(>F)"], marginal = TRUE)))))
#Legend 
list.expr = sapply(1:(nrow(dbrda.terms)-1), function(i){
  x <-  dbrda.terms$names[i]
  a <- round(dbrda.terms[i,"R2"],3)
  b <- p2star(dbrda.terms[i,"Pr(>F)"], marginal = TRUE)
  return(as.expression(substitute(x*": R"^2*" = "*a*b)))
})


legend("bottomright", cex = 0.6,adj = 0,bty = "n",
       legend = list.expr )
dev.off()





# DBRDA = plots in molecule space - all predictors ####
# WIth all predictors: 
dbrda.mol <-  dbrda(tmp ~  Seal_500 + prop.neo + SR + BNIs + Rao,
                    distance = "jaccard",
                    sqrt.dist = TRUE,add = TRUE,
                    data = expl.data)

dbrda.mol.all.jacc <- dbrda.mol

# Extract stats: 
dbrda.model <- anova(dbrda.mol, by=NULL, perm.max=999)
dbrda.model$R2 <-dbrda.model$SumOfSqs/
    sum(dbrda.model$SumOfSqs)

dbrda.terms <- anova(dbrda.mol, by="margin", perm.max=9999)
dbrda.terms$R2 <- dbrda.terms$SumOfSqs/sum(dbrda.terms$SumOfSqs)
dbrda.terms$names <- c("% Impervious","Prop.Neophytes","SR", "BNIs","Rao's Q", "resid")

# GRaph
quartz(width = 5, height = 4.5,
       file = "results/DBRDA plot in MOLECULE.pdf",
       type = "pdf")
par(mar = c(3,3,2,1))
ordi <- plot(dbrda.mol, type = "n",
             mgp = c(1.6,0.4,0), las = 1, tck = -0.01, 
             cex.axis = 0.8, )
points(dbrda.mol, type = "p" , pch = 20,  col = "#48876280")
#dbrda.mol, display = "cn", cex = 0.7,adj = 0)
labels = c("% Impervious", "Prop.Neo","SR","BNIs","Rao")

arrows(0,0, dbrda.mol$CCA$biplot[,1]*3, dbrda.mol$CCA$biplot[,2]*3,
       length = 0.1, angle =30, col = "darkgrey")

text(dbrda.mol$CCA$biplot[,1]*3.2,
     dbrda.mol$CCA$biplot[,2]*3.2,
     labels = c("% Impervious", "Prop.Neo","SR"
                ,"BNIs","Rao"), cex = 0.8)

mtext(3, text= "b) Allergen molecule composition",
      cex = 1, font = 2, line =1)
mtext(3, adj = 0, line = 0.1, cex= 0.7,
      text = "Jaccard dissimilarities")
mtext(3, adj = 1, line = 0.1, cex= 0.7,
      text = as.expression(
        substitute("R"^2*" = "*a*b,
                   list(  a = round(dbrda.model[1,"R2"],3),
                          b = p2star(dbrda.model[1,"Pr(>F)"], marginal = TRUE)))))
#Legend 
list.expr = sapply(1:(nrow(dbrda.terms)-1), function(i){
  x <-  dbrda.terms$names[i]
  a <- round(dbrda.terms[i,"R2"],3)
  b <- p2star(dbrda.terms[i,"Pr(>F)"], marginal = TRUE)
  return(as.expression(substitute(x*": R"^2*" = "*a*b)))
})


legend("bottomright", cex = 0.6,adj = 0,bty = "n",
       legend = list.expr )
dev.off()


# DBRDA = plots in molecule space - SEALING + PROP.NEO - BRAY ####

dbrda.mol <-  dbrda(tmp ~  Seal_500 + prop.neo,
                    distance = "bray",
                    sqrt.dist = TRUE,
                    data = expl.data)
dbrda.mol.two.bray <- dbrda.mol 

# Extract stats
dbrda.model <- anova(dbrda.mol, by=NULL, perm.max=999)
dbrda.model$R2 <-dbrda.model$SumOfSqs/
    sum(dbrda.model$SumOfSqs)
dbrda.terms <- anova(dbrda.mol, by="margin", perm.max=500)
dbrda.terms$R2 <- dbrda.terms$SumOfSqs/sum(dbrda.terms$SumOfSqs)
dbrda.terms$names <- c("% Impervious","Prop.Neophytes","resid")

# Graph:
quartz(width = 5, height = 4.5,
       file = "results/DBRDA plot in MOLECULE_only seal+PropNeo_BRAY.pdf",
       type = "pdf")

par(mar = c(3,3,2,1))
ordi <- plot(dbrda.mol, type = "n",
             mgp = c(1.6,0.4,0), las = 1, tck = -0.01, 
             cex.axis = 0.8, )
points(dbrda.mol, type = "p" , pch = 20,  col = "#48876280")
#dbrda.mol, display = "cn", cex = 0.7,adj = 0)
labels = c("% Impervious", "Prop.Neo")

arrows(0,0, dbrda.mol$CCA$biplot[,1]*2, dbrda.mol$CCA$biplot[,2]*2,
       length = 0.1, angle =30, col = "darkgrey")

text(dbrda.mol$CCA$biplot[,1]*2.2,
     dbrda.mol$CCA$biplot[,2]*2.2,
     labels = c("% Impervious", "Prop.Neo"), cex = 0.8)

mtext(3, text= "b) Allergen molecule composition",
      cex = 1, font = 2, line =1)
mtext(3, adj = 0, line = 0.1, cex= 0.7,
      text = "Bray-Curtis dissimilarities")

mtext(3, adj = 1, line = 0.1, cex= 0.7,
      text = as.expression(
        substitute("R"^2*" = "*a*b,
                   list(  a = round(dbrda.model[1,"R2"],3),
                          b = p2star(dbrda.model[1,"Pr(>F)"], marginal = TRUE)))))
#Legend 
list.expr = sapply(1:(nrow(dbrda.terms)-1), function(i){
  x <-  dbrda.terms$names[i]
  a <- round(dbrda.terms[i,"R2"],3)
  b <- p2star(dbrda.terms[i,"Pr(>F)"], marginal = TRUE)
  return(as.expression(substitute(x*": R"^2*" = "*a*b)))
})


legend("bottomright", cex = 0.6,adj = 0,bty = "n",
       legend = list.expr )
dev.off()





# DBRDA = plots in molecule space - all predictors - BRAY ####
# WIth all predictors: 
dbrda.mol <-  dbrda(tmp ~  Seal_500 + prop.neo + SR + BNIs + Rao,
                    distance = "bray",
                    sqrt.dist = TRUE,add = TRUE,
                    data = expl.data)
dbrda.mol.all.bray <- dbrda.mol

#extract stats: 
dbrda.model <- anova(dbrda.mol, by=NULL, perm.max=999)
dbrda.model$R2 <-dbrda.model$SumOfSqs/
    sum(dbrda.model$SumOfSqs)

dbrda.terms <- anova(dbrda.mol, by="margin", perm.max=9999)
dbrda.terms$R2 <- dbrda.terms$SumOfSqs/sum(dbrda.terms$SumOfSqs)
dbrda.terms$names <- c("% Impervious","Prop.Neophytes","SR", "BNIs","Rao's Q", "resid")

## both still signif

quartz(width = 5, height = 4.5,
       file = "results/DBRDA plot in MOLECULE_BRAY.pdf",
       type = "pdf")
par(mar = c(3,3,2,1))
ordi <- plot(dbrda.mol, type = "n",
             mgp = c(1.6,0.4,0), las = 1, tck = -0.01, 
             cex.axis = 0.8, )
points(dbrda.mol, type = "p" , pch = 20,  col = "#48876280")
#dbrda.mol, display = "cn", cex = 0.7,adj = 0)
labels = c("% Impervious", "Prop.Neo","SR","BNIs","Rao")

arrows(0,0, dbrda.mol$CCA$biplot[,1]*3, dbrda.mol$CCA$biplot[,2]*3,
       length = 0.1, angle =30, col = "darkgrey")

text(dbrda.mol$CCA$biplot[,1]*3.2,
     dbrda.mol$CCA$biplot[,2]*3.2,
     labels = c("% Impervious", "Prop.Neo","SR"
                ,"BNIs","Rao"), cex = 0.8)

mtext(3, text= "b) Allergen molecule composition",
      cex = 1, font = 2, line =1)
mtext(3, adj = 0, line = 0.1, cex= 0.7,
      text = "Bray-curtis dissimilarities")

mtext(3, adj = 1, line = 0.1, cex= 0.7,
      text = as.expression(
        substitute("R"^2*" = "*a*b,
                   list(  a = round(dbrda.model[1,"R2"],3),
                          b = p2star(dbrda.model[1,"Pr(>F)"], marginal = TRUE)))))
#Legend 
list.expr = sapply(1:(nrow(dbrda.terms)-1), function(i){
  x <-  dbrda.terms$names[i]
  a <- round(dbrda.terms[i,"R2"],3)
  b <- p2star(dbrda.terms[i,"Pr(>F)"], marginal = TRUE)
  return(as.expression(substitute(x*": R"^2*" = "*a*b)))
})


legend("bottomright", cex = 0.6,adj = 0,bty = "n",
       legend = list.expr )
dev.off()





### DATA forOrdination of plots in molecule family space  ####

tmp <- allfam_plot # allergen families present per plot
tmp <- tmp[,- which(colnames(tmp) %in% c("Unknown"))]
#tmp <- tmp[-which(rowSums(tmp) == 0), ]

expl.data= cbind(plot_summary[rownames(tmp),
                              c("Seal_500" , "prop.neo", "BNIs",
                                "SR","Rao")],
                 allergen_summary[rownames(tmp),
                                  c("all.num" , "all.num.neo",
                                    "all.num.arc", "all.num.nat")]
)


expl.data$Seal_cut <- cut_interval(expl.data$Seal_500, 6 )
expl.data <- na.omit(expl.data)

tmp <- tmp[rownames(expl.data),]

## Jaccard distance
jacc.mol.fam <- vegdist(tmp, method = "jaccard")

## Bray-Curtis distance ## SAME AS FOR NOW NO ABUNDANCE
bray.mol.fam <- vegdist(tmp, method = "bray")


# # NMDS PLOTS in ALLERGEN FAMILY SPACE ####
# (nmds.mol.fam<- metaMDS(tmp ,
#                         k = 2,
#                         try = 200,
#                         distance = "jaccard",
#                         autotransform = TRUE, 
#                         trymax = 40))
# 
# mol.scores <- as.data.frame(nmds.mol.fam$species)
# mol.scores$names <- sapply(rownames(mol.scores),
#                            FUN = function(x){
#                              strsplit(x ," ")[[1]][1]
#                            })
# site.scores <- as.data.frame(nmds.mol.fam$points)
# 
# plot(site.scores[,1:2],
#      ylim = c(-2,2),
#      xlim = c(-2,2),
#      pch = 20, col = heat.colors(6)[expl.data$Seal_cut])
# text(mol.scores[,1:2]*1.5,
#      label = mol.scores$names, cex = 0.6,
#      col = "black")
# arrows(0,0,
#        mol.scores$MDS1*1.5,mol.scores$MDS2*1.5,
#        col = "black",length = 0.1 )
# 
# # ### Ellipsoid hulls show treatment
# with(expl.data, ordiellipse(nmds.mol.fam, Seal_cut,
#                             kind = "ehull", label = FALSE,
#                             col=heat.colors(6)))
# with(expl.data, ordispider(nmds.mol.fam, Seal_cut,
#                            lty=3,
#                            col=heat.colors(6)))
# legend("topright",
#        legend= c(levels(expl.data$Seal_cut)),
#        fill = heat.colors(6),
#        cex = 0.5
# )    


# # TEST PERMANOVA ADONIS - ALLERGEN FAMILY SPACE ####
# (adonis.mol <-  adonis2(
#   jacc.mol.fam  ~   Seal_500  + prop.neo ,
#   data = expl.data,
#   sqrt.dist = TRUE,
#   permutations = 1000,
#   by = "margin"
# ))
# 
# (adonis.mol <-  adonis2(
#   jacc.mol.fam ~ Seal_500 +  prop.neo  + SR + all.num +  BNIs + Rao,
#   data = expl.data,
#   sqrt.dist = TRUE,
#   permutations = 999,
#   by = "margin"
# ))

#  DBRDA = plots in ALLFAM - Sealing + Prop.neo ####
(dbrda.mol <- dbrda(tmp ~   Seal_500 + prop.neo  ,
                    distance = "jaccard",
                    sqrt.dist = TRUE,
                    data = expl.data))
dbrda.allfam.two.jacc  <- dbrda.mol

(dbrda.model <- anova(dbrda.mol, by=NULL, perm.max=999))
(dbrda.model$R2 <-dbrda.model$SumOfSqs/
    sum(dbrda.model$SumOfSqs))
anova(dbrda.mol, by="axis", perm.max=500)
(dbrda.terms <- anova(dbrda.mol, by="margin", perm.max=9999))
dbrda.terms$R2 <- dbrda.terms$SumOfSqs/sum(dbrda.terms$SumOfSqs)
(dbrda.terms)
dbrda.terms$names <- c("% Impervious","Prop.Neophytes", "resid")

# GRaph
quartz(width = 5, height = 4.5,
       file = "results/DBRDA plot in ALLFAM.pdf",
       type = "pdf")
par(mar = c(3,3,2,1))
ordi <- plot(dbrda.mol, type = "n",
             mgp = c(1.6,0.4,0), las = 1, tck = -0.01, 
             cex.axis = 0.8, )
points(dbrda.mol, type = "p" , pch = 20,  col = "#48876280")
#dbrda.mol, display = "cn", cex = 0.7,adj = 0)

arrows(0,0, dbrda.mol$CCA$biplot[,1]*2, dbrda.mol$CCA$biplot[,2]*2,
       length = 0.1, angle =30, col = "darkgrey")

text(dbrda.mol$CCA$biplot[,1]*2.2,
     dbrda.mol$CCA$biplot[,2]*2.2,
     labels = c("% Impervious", "Prop.Neo"), cex = 0.8)

mtext(3, text= "c) Allergen Biochemical Family composition",
      cex = 1, font = 2, line =1)
mtext(3, adj = 0, line = 0.1, cex= 0.7,
      text = "Jaccard dissimilarities")
mtext(3, adj = 1, line = 0.1, cex= 0.7,
      text = as.expression(
        substitute("R"^2*" = "*a*b,
                   list(  a = round(dbrda.model[1,"R2"],3),
                          b = p2star(dbrda.model[1,"Pr(>F)"], marginal = TRUE)))))

#Legend 
list.expr = sapply(1:(nrow(dbrda.terms)-1), function(i){
  x <-  dbrda.terms$names[i]
  a <- round(dbrda.terms[i,"R2"],3)
  b <- p2star(dbrda.terms[i,"Pr(>F)"], marginal = TRUE)
  return(as.expression(substitute(x*": R"^2*" = "*a*b)))
})

legend("bottomleft", cex = 0.6, adj = 0,bty = "n",
       legend = list.expr )
dev.off()


#  DBRDA = plots in ALLFAM - ALL PREDICTORS ####
(dbrda.mol <- dbrda(tmp ~  Seal_500 + prop.neo + SR +  BNIs + Rao,
                    distance = "jaccard",
                    sqrt.dist = TRUE,
                    data = expl.data))
dbrda.allfam.all.jacc  <- dbrda.mol


(dbrda.model <- anova(dbrda.mol, by=NULL, perm.max=999))
(dbrda.model$R2 <-dbrda.model$SumOfSqs/
    sum(dbrda.model$SumOfSqs))
anova(dbrda.mol, by="axis", perm.max=500)
(dbrda.terms <- anova(dbrda.mol, by="margin", perm.max=9999))
dbrda.terms$R2 <- dbrda.terms$SumOfSqs/sum(dbrda.terms$SumOfSqs)
(dbrda.terms)
dbrda.terms$names <- c("% Impervious","Prop.Neophytes","SR", "BNIs","Rao's Q", "resid")

#GRaph

quartz(width = 5, height = 4.5,
       file = "results/DBRDA plot in ALLFAM_all preds.pdf",
       type = "pdf")
par(mar = c(3,3,2,1))
ordi <- plot(dbrda.mol, type = "n",
             mgp = c(1.6,0.4,0), las = 1, tck = -0.01, 
             cex.axis = 0.8, )
points(dbrda.mol, type = "p" , pch = 20,  col = "#48876280")
#dbrda.mol, display = "cn", cex = 0.7,adj = 0)

arrows(0,0, dbrda.mol$CCA$biplot[,1]*2, dbrda.mol$CCA$biplot[,2]*2,
       length = 0.1, angle =30, col = "darkgrey")

text(dbrda.mol$CCA$biplot[,1]*2.2,
     dbrda.mol$CCA$biplot[,2]*2.2,
     labels = c("% Impervious", "Prop.Neo","SR",
                "BNIs","Rao"), cex = 0.8)

mtext(3, text= "c) Allergen Biochemical Family composition",
      cex = 1, font = 2, line =1)
mtext(3, adj = 0, line = 0.1, cex= 0.7,
      text = "Jaccard dissimilarities")
mtext(3, adj = 1, line = 0.1, cex= 0.7,
      text = as.expression(
        substitute("R"^2*" = "*a*b,
                   list(  a = round(dbrda.model[1,"R2"],3),
                          b = p2star(dbrda.model[1,"Pr(>F)"], marginal = TRUE)))))
#Legend 
list.expr = sapply(1:(nrow(dbrda.terms)-1), function(i){
  x <-  dbrda.terms$names[i]
  a <- round(dbrda.terms[i,"R2"],3)
  b <- p2star(dbrda.terms[i,"Pr(>F)"], marginal = TRUE)
  return(as.expression(substitute(x*": R"^2*" = "*a*b)))
})

legend("bottomright", cex = 0.6,adj = 0,bty = "n",
       legend = list.expr )
dev.off()

