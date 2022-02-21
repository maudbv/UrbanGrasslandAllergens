# Analyse turnover in allergen families along the two novelty gradients

### DATA for Ordination of plots in molecule family space  ####
tmp <- allfam_plot # allergen families present per plot/ no abundance calculated
tmp <- tmp[,- which(colnames(tmp) %in% c("Unknown"))]
#tmp <- tmp[-which(rowSums(tmp) == 0), ]

expl.data= cbind(plot_summary[rownames(tmp),
                              c("Seal_500" , "prop.neo")],
                 allergen_summary[rownames(tmp),
                                  c("all.num" , "all.num.neo",
                                    "all.num.arc", "all.num.nat")])


expl.data$Seal_cut <- cut_interval(expl.data$Seal_500, 6 )
expl.data <- na.omit(expl.data)

tmp <- tmp[rownames(expl.data),]

## Jaccard distance
jacc.mol.fam <- vegdist(tmp, method = "jaccard")

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


