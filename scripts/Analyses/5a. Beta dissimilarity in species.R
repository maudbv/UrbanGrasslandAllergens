## Ordination of PLOTS in Molecule space ordination

library(vegan)
library(ggplot2)

# DATA for Ordination of plots in ALLERGENIC SPECIES space  ####

# plant community matrix of allergenic species:
tmp <- vegcomm[,allergenics] 

# Explanatory variables: 
expl.data= cbind(plot_summary[rownames(tmp),
                        c("Seal_500" , "prop.neo", "BNIs",
                        "SR","Rao")],
                 allergen_summary[rownames(tmp),
                          c("all.num" , "all.num.neo",
                            "all.num.arc", "all.num.nat")]
                 )

# Create 5 categories ot Sealing (quartiles) :
expl.data$Seal_cut <- cut(
  expl.data$Seal_500, 
  breaks = c(0,
             quantile(expl.data$Seal_500,probs = c(0.25,0.50,0.75)),
             100))

# Create 2 categories ot Sealing (quartiles) :
expl.data$Seal_cut <- cut(
  expl.data$Seal_500, 
  breaks = c(0,
             quantile(expl.data$Seal_500,probs = c(1/3, 2/3)),
             100))
expl.data <- na.omit(expl.data)

tmp <- tmp[rownames(expl.data),]

# Check we have the exact same plot name rows
stopifnot(all(rownames(tmp) == rownames(expl.data)))

# vector of colors for the levels of Sealing: 
col.vec <- colorRampPalette(
  c("forestgreen","yellow","firebrick"))(
    length(levels(expl.data$Seal_cut)))

# Jaccard dissimilarities
jacc.veg.dist <- vegdist(tmp, method = "jaccard")

# ## NMDS  for Allergenic species space ####
# nmds.veg.plot <- metaMDS(tmp, distance = "jacc",
#                          k = 3,
#                          try = 200,
#                          trymax = 40)
# 
# 
# plot(nmds.veg.plot$points[, c(1,2)],
#      ylim = c(-1,1),
#      xlim = c(-2,2),
#      display="sites",
#      type = "p", pch = 20,
#      col = col.vec[as.numeric(expl.data$Seal_cut)])
# text(nmds.veg.plot,
#      display="species",col = "lightgrey", cex = 0.5)
# 
# 
# # Test of effect of gradients on molecule space:
# (adonis.veg.plot <- adonis2(
#   jacc.veg.dist ~  Seal_500 + prop.neo,
#   data = expl.data, sqrt.dist = TRUE,
#   by = "margin"))
# 
# 
# redundancy analysis Seal + Prop.neo #####

dbrda.sp <-  dbrda(tmp ~   Seal_500 + prop.neo  ,
                    distance = "jaccard",
                    sqrt.dist = TRUE,
                    data = expl.data)

dbrda.sp.two.jacc <- dbrda.sp

# extract stats: 
dbrda.model <- anova(dbrda.sp, by=NULL, perm.max=999)
dbrda.model$R2 <-dbrda.model$SumOfSqs/
  sum(dbrda.model$SumOfSqs)

dbrda.terms <- anova(dbrda.sp, by="margin", perm.max=999)
dbrda.terms$R2 <- dbrda.terms$SumOfSqs/sum(dbrda.terms$SumOfSqs)
dbrda.terms$names <- c("% Impervious","Prop.Neophytes","Resid")

# GRAPH
quartz(width = 5, height = 4.5,
       file = "results/DBRDA plot in Allergen SR_seal+prop.neo.pdf",
       type = "pdf")
par(mar = c(3,3,2,1))
ordi <- plot(dbrda.sp, type = "n",
             mgp = c(1.6,0.4,0), las = 1, tck = -0.01, 
             cex.axis = 0.8, )
points(dbrda.sp, type = "p" , pch = 20,  col = "#B2266280")
#dbrda.sp, display = "cn", cex = 0.7,adj = 0)

arrows(0,0, dbrda.sp$CCA$biplot[,1]*2, dbrda.sp$CCA$biplot[,2]*2,
       length = 0.1, angle =30, col = "darkgrey")

text(dbrda.sp$CCA$biplot[,1]*2.2,
     dbrda.sp$CCA$biplot[,2]*2.2,
     labels = c("% Impervious", "Prop.Neo"), cex = 0.8)

mtext(3, text= "a) Allergenic species composition",
      cex = 1, font = 2, line =1)
mtext(3, adj = 0, line = 0.1, cex= 0.7,
      text = "Jaccard dissimilarities")

mtext(3, adj = 1, line = 0.1, cex= 0.7,
      text = as.expression(substitute("R"^2*" = "*a*b,
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



## db redundancy analysis ALL predictors ####
dbrda.sp <- dbrda(tmp ~  Seal_500  + prop.neo + SR + BNIs + Rao ,
                    distance = "jaccard",
                    sqrt.dist = TRUE,
                    data = expl.data)

dbrda.sp.all.jacc <- dbrda.sp

# extract stats: 
dbrda.model <- anova(dbrda.sp, by=NULL, perm.max=999)
dbrda.model$R2 <-dbrda.model$SumOfSqs/
    sum(dbrda.model$SumOfSqs)

dbrda.terms <- anova(dbrda.sp, by="margin", perm.max=999)
dbrda.terms$R2 <- dbrda.terms$SumOfSqs/sum(dbrda.terms$SumOfSqs)
dbrda.terms$names <- c("% Impervious","Prop.Neophytes","SR", "BNIs", "Rao","resid")

# GRAPH
quartz(width = 5, height = 4.5,
       file = "results/DBRDA plot in Allergen SR_all preds.pdf",
       type = "pdf")
par(mar = c(3,3,2,1))
ordi <- plot(dbrda.sp, type = "n",
             mgp = c(1.6,0.4,0), las = 1, tck = -0.01, 
             cex.axis = 0.8, )
points(dbrda.sp, type = "p" , pch = 20,  col = "#B2266280")
#dbrda.sp, display = "cn", cex = 0.7,adj = 0)

arrows(0,0, dbrda.sp$CCA$biplot[,1]*1.8, dbrda.sp$CCA$biplot[,2]*1.8,
       length = 0.1, angle =30, col = "darkgrey")

text(dbrda.sp$CCA$biplot[,1]*1.9,
     dbrda.sp$CCA$biplot[,2]*1.9,
     labels = c("% Impervious", "Prop.Neo","SR",
                "BNIs","Rao"), cex = 0.8, adj = 0)

mtext(3, text= "a) Allergenic species composition",
      cex = 1, font = 2, line =1)
mtext(3, adj = 0, line = 0.1, cex= 0.7,
      text = "Jaccard dissimilarities")

mtext(3, adj = 1, line = 0.1, cex= 0.7,
      text = as.expression(substitute("R"^2*" = "*a*b,
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


# redundancy analysis Seal + Prop.neo - BRAY #####

dbrda.sp <-  dbrda(tmp ~   Seal_500 + prop.neo  ,
                    distance = "bray",
                    sqrt.dist = TRUE,
                    data = expl.data)
dbrda.sp.two.bray <-dbrda.sp

# extract stats: 
dbrda.model <- anova(dbrda.sp, by=NULL, perm.max=999)
dbrda.model$R2 <-dbrda.model$SumOfSqs/
    sum(dbrda.model$SumOfSqs)

dbrda.terms <- anova(dbrda.sp, by="margin", perm.max=999)
dbrda.terms$R2 <- dbrda.terms$SumOfSqs/sum(dbrda.terms$SumOfSqs)
dbrda.terms$names <- c("% Impervious","Prop.Neophytes","Resid")


# GRAPH
quartz(width = 5, height = 4.5,
       file = "results/DBRDA plot in Allergen SR_seal+prop.neo_BRAY.pdf",
       type = "pdf")
par(mar = c(3,3,2,1))
ordi <- plot(dbrda.sp, type = "n",
             mgp = c(1.6,0.4,0), las = 1, tck = -0.01, 
             cex.axis = 0.8, )
points(dbrda.sp, type = "p" , pch = 20,  col = "#B2266280")
#dbrda.sp, display = "cn", cex = 0.7,adj = 0)

arrows(0,0, dbrda.sp$CCA$biplot[,1]*2, dbrda.sp$CCA$biplot[,2]*2,
       length = 0.1, angle =30, col = "darkgrey")

text(dbrda.sp$CCA$biplot[,1]*2.2,
     dbrda.sp$CCA$biplot[,2]*2.2,
     labels = c("% Impervious", "Prop.Neo"), cex = 0.8)

mtext(3, text= "a) Allergenic species composition",
      cex = 1, font = 2, line =1)
mtext(3, adj = 0, line = 0.1, cex= 0.7,
      text = "Bray-Curtis dissimilarities")

mtext(3, adj = 1, line = 0.1, cex= 0.7,
      text = as.expression(substitute("R"^2*" = "*a*b,
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



## db redundancy analysis ALL predictors - BRAY ####
dbrda.sp <- dbrda(tmp ~  Seal_500 + prop.neo + SR + BNIs + Rao,
                    distance = "bray",
                    sqrt.dist = TRUE,
                    data = expl.data)

dbrda.sp.all.bray <-dbrda.sp

# Extract stats:
dbrda.model <- anova(dbrda.sp, by=NULL, perm.max=999)
dbrda.model$R2 <-dbrda.model$SumOfSqs/
    sum(dbrda.model$SumOfSqs)

dbrda.terms <- anova(dbrda.sp, by="margin", perm.max=999)
dbrda.terms$R2 <- dbrda.terms$SumOfSqs/sum(dbrda.terms$SumOfSqs)
dbrda.terms$names <- c("% Impervious","Prop.Neophytes","SR", "BNIs", "Rao", "resid")

# GRAPH
quartz(width = 5, height = 4.5,
       file = "results/DBRDA plot in Allergen SR_all preds_BRAY.pdf",
       type = "pdf")
par(mar = c(3,3,2,1))
ordi <- plot(dbrda.sp, type = "n",
             mgp = c(1.6,0.4,0), las = 1, tck = -0.01, 
             cex.axis = 0.8, )
points(dbrda.sp, type = "p" , pch = 20,  col = "#B2266280")
#dbrda.sp, display = "cn", cex = 0.7,adj = 0)

arrows(0,0, dbrda.sp$CCA$biplot[,1]*1.8, dbrda.sp$CCA$biplot[,2]*1.8,
       length = 0.1, angle =30, col = "darkgrey")

text(dbrda.sp$CCA$biplot[,1]*1.9,
     dbrda.sp$CCA$biplot[,2]*1.9,
     labels = c("% Impervious", "Prop.Neo","SR",
                "BNIs","Rao"), cex = 0.8, adj = 0)

mtext(3, text= "a) Allergenic species composition",
      cex = 1, font = 2, line =1)
mtext(3, adj = 0, line = 0.1, cex= 0.7,
      text = "Bray-Curtis dissimilarities")

mtext(3, adj = 1, line = 0.1, cex= 0.7,
      text = as.expression(substitute("R"^2*" = "*a*b,
                                      list(  a = round(dbrda.model[1,"R2"],3),
                                             b = p2star(dbrda.model[1,"Pr(>F)"], marginal = TRUE)))))
#Legend 
list.expr = sapply(1:(nrow(dbrda.terms)-1), function(i){
  x <-  dbrda.terms$names[i]
  a <- round(dbrda.terms[i,"R2"],3)
  b <- p2star(dbrda.terms[i,"Pr(>F)"], marginal = TRUE)
  return(as.expression(substitute(x*": R"^2*" = "*a*b)))
})

legend("topright", cex = 0.6,adj = 0,bty = "n",
       legend = list.expr )
dev.off()
