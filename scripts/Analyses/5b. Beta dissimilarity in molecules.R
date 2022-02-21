## # Analyse turnover in allergen molecules along the two novelty gradients


library(vegan)
library(ggplot2)

### DATA for ordination of plots in allergen molecule space  ####

tmp <- mol_plot # allergen molecules expected per plot
# abundance values correspond to the number of different species contributing this allergen, not species cover

tmp <- tmp[, -which(names(tmp) == "Unknown")]

expl.data= cbind(plot_summary[rownames(tmp),
                              c("Seal_500" , "prop.neo")],
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
       file = "results/DBRDA plot in MOLECULE.pdf",
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


legend("bottomleft", cex = 0.6,adj = 0,bty = "n",
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
       file = "results/DBRDA plot in MOLECULE_BRAY.pdf",
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


legend("bottomleft", cex = 0.6,adj = 0,bty = "n",
       legend = list.expr )
dev.off()




