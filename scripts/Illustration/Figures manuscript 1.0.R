## Figures for manuscript

(function() {
# Figure 1: Cover of allergenic species along novelty gradients ####
png(file = "results/figure 1.png",width = 13, height = 18,
    unit = "cm", res = 600)
# pdf(file = "results/figure 1.pdf",width = 7, height = 6)
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
   tmp <- allergen_summary
   y <- c("cover.all", "cover.exo.all", "all.num.exo")[k]
   mod <- c("n","n","p")[k]
   leg <- c("Allergenic\n cover ",
            "Non-native allergenic\n cover",
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
dev.off()




# Figure 2: Mean PAV along gradients ####
# pdf(file = "results/figure 2.pdf",
#     width = 7, height = 4)
 png(file = "results/figure 2.png",
     width = 18, height = 10,unit ="cm", res = 600)
par (
   mfrow = c(1, 2),
   bty = "l",
   mar = c(4, 3, 2, 1),
   oma = c(0, 0, 1, 0),
   mgp = c(2,0.4,0),
   tcl = - 0.2,
   las = 1,
   cex.axis = 0.8,
   xpd = TRUE
)


# Panel 1: interaction plot for mean.PAV
tmp <- allergen_summary
f <- lm(mean.pav ~ Seal_500 * prop.neo, tmp)
library(interactions)
p <- interactions::interact_plot(f, pred = Seal_500,  modx = prop.neo,
                                 legend.main = "Prop. Neophytes",
                                 plot.points = TRUE , interval = FALSE, 
                                 x.label = "% Impervious surfaces",
                                 y.label = "mean.PAV",
                                 colors = "seagreen",
                                 point.size = 1,line.thickness = 0.8,
                                 partial.residuals = TRUE,
                                 facet.modx = FALSE
                                 )
# Extract representative cut points for prop.neo :
(unique(p$data$prop.neo)[1:2] + unique(p$data$prop.neo)[2:3])/2
tmp$cut_prop.neo <- cut(allergen_summary$prop.neo, 
                        breaks = c(0, 0.0405,0.0886,1))
cols <- colorRampPalette(c("peachpuff3","firebrick"))(3)

plot(mean.pav ~ Seal_500 ,data = tmp,
     pch = 20,
     col = cols[tmp$cut_prop.neo],
     ann = FALSE)

add.stats(f, plot.abline = "no")
#mtext(3, text = "mean.PAV ~ %Imp.Surf * Prop.Neo",
#      adj = 0, cex = 0.7)

mtext( 1,text = "% Impervious surfaces", outer = FALSE,
       cex = 0.9,line = 2, las = 0)
mtext(2,text = substitute("Mean" [a],list( a="PAV")),outer = FALSE,
      cex = 0.9, line = 1.5, las = 0)

mtext(3,text = "(a)", adj = -0.2)

# Add illustrative regression lines: 
ltys <- c("dotted","dashed","solid")
lines(mean.pav ~ Seal_500 ,
      data = p$data[p$data$modx_group == "+ 1 SD",],
      col = cols[3], lty = ltys[3], lwd = 1.2)
lines(mean.pav ~ Seal_500 ,
      data = p$data[p$data$modx_group == "Mean",],
      col = cols[2], lty = ltys[2], lwd = 1.2)
lines(mean.pav ~ Seal_500 ,
      data = p$data[p$data$modx_group == "- 1 SD",],
      col = cols[1] , lty = ltys[1], lwd = 1.2)

legend(2.5,8.5,
       legend = c("-SD","Mean","+SD"),
       title = "Proportion of Neophytes",
       lty = ltys, col = cols, lwd = 1.2,
       text.col = "grey10",
       bty =  "n",
       cex = 0.6,
       inset = -0.07
)

# Second panel : mean archaeophyte PAVs
plot(mean.arc.pav ~ Seal_500, data = allergen_summary,
     pch = 20,col = "grey",ann = FALSE)

add.stats(f <- lm(mean.arc.pav ~ Seal_500, data = allergen_summary, na.action = "na.omit"),
          type = "lm",l.col = "black")
mtext( 1,text = "% Impervious surfaces", outer = FALSE,
       cex = 0.9,line = 2, las = 0)
mtext(2,text = substitute("Archaeophyte Mean" [a],list( a="PAV")),
      outer = FALSE,
      cex = 0.9, line = 1.5, las = 0)
mtext(3,text = "(b)", adj = -0.2)

dev.off()



# Figure 3: Violin plots of pheno traits ####

png(file = "results/figure 3.png", 
    width = 24, height = 11, unit = "cm", res = 600)

species_allergen$Introduction_status_Seitz2012<-factor(species_allergen$Introduction_status_Seitz2012, levels = c("I","A","N"))

# GRAPH:
par(mfrow = c(1,3), oma = c(0,3,0,0), mar = c(3,2,4,1),
    mgp = c(2,0.7,0),tcl = - 0.3,
    cex = 0.9)
months.labels <- months
substr(months.labels, 1, 1) <- toupper(substr(months.labels, 1, 1))

# Onset of Flowering
plot(as.numeric(fl.beg) ~ jitter(as.numeric(status_num),factor = 1),
     data =species_allergen, axes =FALSE,ann = FALSE,
     ylim = c(1,12),xlim = c(0.5, 3.5))
axis(1, at = c(1,2,3), labels= c("Natives","Archaeo.","Neophytes"),
     cex.axis = 0.9)
axis(2, at = 1:12, labels= months.labels, las = 1, cex.axis=0.9)
abline(h = 1:12, col = "grey90", lty = "dotted")
title(main ="a) Onset of Flowering")
vioplot(as.numeric(fl.beg) ~ Introduction_status_Seitz2012, 
        data =species_allergen,
        outlier = FALSE,
        col ="#ff7cc250" , border = "#ff7cc2",
        ann = FALSE,axes = FALSE, add = TRUE)
mtext( 3, adj = 1, cex = 0.7,
   text =substitute("R"^2*"="*a*b,
                 list( a = round(pheno.test$fl.beg$r2,2),
                       b =p2star(pheno.test$fl.beg$P, marginal = T)))
)
# Neophytes start flowering later
text(x = 1:3, y = 11.5, labels = c("a", "a", "b"))

# End of flowering
par(mar = c(3,1,4,2))
plot(as.numeric(fl.end) ~ jitter(as.numeric(status_num),factor = 1),
     data =species_allergen, axes =FALSE,ann = FALSE,
     ylim = c(1,12),xlim = c(0.5, 3.5))
axis(1, at = c(1,2,3), labels= c("Natives","Archaeo.","Neophytes"),
     cex.axis = 0.9)
title(main ="b) End of Flowering")
abline(h = 1:12, col = "grey90", lty = "dotted")
vioplot(as.numeric(fl.end) ~ Introduction_status_Seitz2012, 
        data =species_allergen,
        outlier = FALSE, col ="#ff7cc250" , border = "#ff7cc2",
        ann = FALSE,axes = FALSE, add = TRUE)
axis(2, at = 1:12, labels = NA, las = 1, cex.axis=0.8)
mtext( 3, adj = 1, cex = 0.7,
       text =substitute("R"^2*"="*a*b,
                        list( a = round(pheno.test$fl.end$r2,2),
                              b =p2star(pheno.test$fl.end$P, marginal = T)))
)
text(x = 1:3, y = 3, labels = c("a", "ab", "b"))

# Length of flowering
par(mar = c(3,2,4,1))
plot(as.numeric(fl.period) ~ jitter(as.numeric(status_num),factor = 1),
     data =species_allergen, axes =FALSE,ann = FALSE,
     ylim = c(1,12),xlim = c(0.5, 3.5))
axis(1, at = c(1,2,3), labels= c("Natives","Archaeo.","Neophytes"),
     cex.axis = 0.9)
title(main ="c) Length of Flowering")
#abline(h = 1:12, col = "grey90", lty = "dotted")
vioplot(as.numeric(fl.period) ~ Introduction_status_Seitz2012, 
        data =species_allergen,
        outlier = FALSE, col ="#fc5a5a50" , border = "#fc5a5a",
        ann = FALSE,axes = FALSE, add = TRUE)
axis(2, at = 1:12, labels=1:12, las = 1, cex.axis=0.9)
mtext(2, text = "Number of months", cex = 0.8, line = 1.8)
mtext( 3, adj = 1, cex = 0.7,
       text =substitute(
          "R"^2*"="*a*b,
          list( a = round(pheno.test$fl.period$r2,2),
                b =p2star(pheno.test$fl.period$P, marginal = T)))
)
text(x = 1:3, y = 12, labels = c("a", "b", "a"))

dev.off()

# Figure 4: Flowering species per months in 3 clusters ####

# col.status <- c(GetColors(1, alpha = 0.6,stop = c(0.5,0.51)),
#                 GetColors(1, alpha = 0.6,stop = c(0.8,0.81)),
#                 GetColors(1, alpha = 0.6,stop = c(0.9,0.91)))

col.status <- c(GetColors(1, alpha = 0.6,stop = c(0.3,0.31),
                          scheme = "iridescent"),
                GetColors(1, alpha = 0.6,stop = c(0.65,0.66),
                          scheme = "iridescent"),
                GetColors(1, alpha = 0.6,stop = c(0.8,0.81),
                          scheme = "iridescent"))

# GRAPH POLYGONS - NUMBER OF ALLERGEN SPECIES & ALLFAM PER MONTH

png(file = "results/figure 4.png", 
     width = 24, height = 18, unit = "cm", res = 600)

## Subsets along gradient: 
z = plot_summary$Seal_levels
table(z)
par(mfrow = c(2,3), mar = c(3,2,1,0), oma = c(0,2,2,1),
    mgp = c(2,0.7,0),tcl = - 0.3,
    cex = 0.8)

# Number of allergenic species
for (p in c(1:3)) {
   
   # Select allergenic species in the plot
   mat <- vegcomm[,colSums(vegcomm[as.numeric(z) == p,])>0]
   sp <- intersect(names(mat), allergenics)
   
   # flowering species richness per month
   
   #indigenous species
   yi <- colSums(flower_month[intersect(sp, natives),]>0)
   
   # Neophytes
   if (length(intersect(sp, neophytes))>0) {
      yn <- colSums(flower_month[intersect(sp, neophytes),]>0)
   } else { yn <- NA}
   
   # Archaeophytes
   if (length(intersect(sp, archaeophytes))>0) {
      ya <- colSums(flower_month[intersect(sp, archaeophytes),]>0)
   } else { ya <- NA}
   
   # Plot
   plot(yi, type = "n",
        ann = F, axes = F, bty = "l", ylim = c(0, 35))
   abline(v = 1:12, col = "grey90", lty = "dotted")
   
   lines(1:12, yi, col = col.status[1],type = "l")
   polygon( x= c(1:12,12:1),
            y =c(yi, rep(0,12)),
            border=col.status[1], 
            col = col.status[1])

      polygon( x= c(1:12,12:1),
               y =  c(ya, rep(0,12)),
               border = col.status[2],
               col = col.status[2])

      polygon( x= c(1:12,12:1),
               y =  c(yn, rep(0,12)),
               border= col.status[3],
               col = col.status[3])
   axis(1, at = 1:12,
        labels = substr(months.labels,1,3),
        las = 1,
        cex.axis = 1,
        las = 2
        )
   axis(2, las = 1)
   # Letter
   mtext(3, 
         text = paste("(", letters[p], ")",
                      sep = "" ),
         line= 0 , adj =0,cex = 0.8, font = 3)
   
   # Level
   mtext(3, 
         text = levels(z)[p],
         line= 0.8, cex = 1, font = 2)
   mtext(3, 
         text = paste("(n = ",table(z)[p],")", sep = ""),
         line= -0.2, cex = 0.7, font = 0)                  

   if (p ==1) {
      mtext(2, text = "Flowering allergenic species",
            outer = FALSE, cex = 1, line= 2.3)     
   }
}

# Allergen families: 

for (p in c(1:3)) {
   yt <- SRallfam_subset_month[p,]
   yn <- SRallfam_subset_month_neo[p,]
   ya <- SRallfam_subset_month_arc[p,]
   yi <- SRallfam_subset_month_nat[p,]
   
   
   # Plot
   plot(yt, type = "n",
        ann = F, axes = F, bty = "l", ylim = c(0, 35))
   # vertical grid lines
   abline(v = 1:12, col = "grey90", lty = "dotted")
   
   # total allergen family richness
   lines( x= 1:12, yt, col= "grey")
   
   # For each floristic group:
   for (f in 1:3) {
   y <- list(yi, ya,yn)[[f]]
    polygon( x= c(1:12,12:1),
             y =c(y, rep(0,12)),
             border=col.status[f],
             col = col.status[f])
   }   
   
   # y-axis label
   if (p == 1) {
      mtext(2, text = "Allergen families",
            outer = FALSE, cex = 1,line= 2.3)
   }
   
   #x-axis
   axis(1, at = 1:12,
        labels = substr(months.labels,1,3),
        cex.axis = 1,las = 2)
   # y-axis
   axis(2, las = 1)
   
   # letters:
   mtext(3, 
         text = paste("(", letters[3+p], ")",
                      sep = "" ),
         line= 0 , adj = 0,cex = 0.8, font = 3)
}

# Legend
legend("topleft",
       legend = c("All", "Natives","Archaeophytes","Neophytes"),
       fill = c(NA,col.status), border = c("grey",NA,NA,NA),
       bty ="n",  cex = 0.9)



dev.off()




# Figure S1: Broad vs strict definitions of allergenicity ####

png(width = 14, height = 15, unit ="cm", res = 600,
    file = "results/figure S1.png")
# Represent change between strict and broad methods :
# x <- rowSums(mol_plot_strict[, -ncol(mol_plot_strict)]) # removing "unknown"
# y <- rowSums(mol_plot_genus[, -ncol(mol_plot_genus)])# removing "unknown"
x <- rowSums(mol_plot_strict) # with "unknown"
y <- rowSums(mol_plot_genus)# with "unknown"

plot(x, y,   
     pch = 20,
     xlim = c(10,110),ylim = c(10,110),
     xlab = "Strict Molecule Richness",
     ylab = "Broad Molecule Richness")
abline(0,1)
add.stats(cor.test(rowSums(mol_plot_strict),
                   rowSums(mol_plot_genus), 
                   method = "spearman"))
dev.off()


# Figure S2A: allergen families by floristic status ####
png(width = 14, height = 15, unit ="cm", res = 600,
    file = "results/figure S2A.png")

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

tmp<- allfam_mat[,c(29:26, 1, 25:2 )]

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
# Figure S2B: Allergen family accumulatino curve #####
png(width = 11, height = 12, unit ="cm", res = 600,
    file = "results/figure S2B.png")

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
   sac <- specaccum(allfam_mat[intersect(pool[[i]],allergenics),], 
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
   
   # add prediction
   sac <- specaccum(allfam_mat[intersect(pool[[i]],allergenics),], 
                    method = "exact")
   mod1 <- fitspecaccum(sac, model = "lomolino", method = "random")
   pred.model <- predict(mod1, newdata = 1:50)
   lines(pred.model,
         col=col.line[i], lty=2)
   if (i ==1) y.pos<- max(pred.model) -3.5 else y.pos<- max(pred.model)+1
   text(x = 52, y =y.pos, label = leg, 
        cex = 0.7, col=col.line[i],adj = 1)
}

abline(v = 10, lty = 3, xpd = FALSE)
print(rare.est)
# 
# legend("bottomright", legend = c("Natives","Archaeophytes","Neophytes"),
#        fill= col.line, bty ="n",
#        cex = 0.7)


dev.off()



# Figure S2C: NMDS of species in allergen family space ####
png(width = 17, height = 13, unit ="cm", res = 600,
    file = "results/figure S2C.png")

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




# Figure S3: temporal distribution of allergenics by family ####
#  FAMILY allergenic flowering
png(width = 26, height = 11, unit ="cm", res = 600,
    file = "results/figure S3.png")

par(mfrow = c(1,2), mar = c(2,4,2,1),
    mgp = c(2,0.5,0),
    tcl = - 0.3,
    las = 1,
    cex.axis = 0.8)

# All allergenics
y <- flower_month
y$family <- species_allergen[rownames(species_allergen),"family"]
y <- y[allergenics,]
y <- doBy::summaryBy(.~ family,data = y, FUN = sum)
rownames(y) <- y$family
y<- y[-1]

# Vector of colors per family:
col.vec <- GetColors(nrow(y), stops = c(0.1,1))
names(col.vec) <- rownames(y)

plot(1:12, seq(0, max(y)+2, length.out = 12),  type = "n",
     axes = FALSE, xlab = "", ylab = "Nb. Flowering allergenic species")
axis(1, at = 1:12, labels = substr(months, 1,3))
axis(2, las = 1)
abline(v = 1:12, col = "grey90", lty = "dotted")

for (i in 1:nrow(y)){
   rf <- range(which(y[i,]>0))+c(-1,+1)
   f<- smooth.spline(x = 1:12, y = as.numeric(y[i,]), spar = 0.2)
   pf <- predict(f, x = seq(rf[1], rf[2], length.out = 100))
   lines(pf,col = col.vec[rownames(y)[i]])
}
legend( "topleft",legend = rownames(y),
        lty = "solid", lwd = 1.4,
        col = col.vec, cex = 0.5, bty = "n")


# For neophytes:  FAMILY allergenic flowering
y <- flower_month
y$family <- species_allergen[rownames(species_allergen),"family"]
y <- y[intersect(allergenics,exotics),]
y <- doBy::summaryBy(.~ family,data = y, FUN = sum)
rownames(y) <- y$family
y<- y[-1]

plot(1:12, seq(0, max(y)+2, length.out = 12),  type = "n",
     axes = FALSE, xlab = "", ylab = "Nb. Flowering non-native allergenic species")
axis(1, at = 1:12, labels = substr(months, 1,3))
axis(2, las = 1)
abline(v = 1:12, col = "grey90", lty = "dotted")

for (i in 1:nrow(y)){
   rf <- range(which(y[i,]>0))+c(-1,+1)
   f<- smooth.spline(x = 1:12, y = as.numeric(y[i,]), spar = 0.2)
   pf <- predict(f, x = seq(rf[1], rf[2], length.out = 100))
   lines(pf,col = col.vec[rownames(y)[i]])
}
dev.off()

# Figure SB.1: Beta-dissimilarity and dbRDA ####
# GRAPH
png(width = 16, height = 22,unit ="cm", res = 600,
    file = "results/figure SB.1.png")


par(mfcol = c(3,2), mar = c(3,3,2,1),
    mgp = c(1.6,0.4,0), las = 1, tck = -0.01)

for (i in 1:3) {
   
   #Select dbrda:
   dbrda.sp <- list(dbrda.sp.all.jacc,
                    dbrda.mol.all.jacc, 
                    dbrda.allfam.all.jacc)[[i]]
   
   tit <- c("a)",
            "c)",
            "e)")[i]
   col.code <- c("#B2266270", "#48876270", "#3681af70")[i]
   
   # Calculate overall statistics:
   (dbrda.model <- anova(dbrda.sp, by=NULL, perm.max=999))
   (dbrda.model$R2 <-dbrda.model$SumOfSqs/
         sum(dbrda.model$SumOfSqs))
   
   # Calculate partial R2:
   (dbrda.terms <- anova(dbrda.sp, by="margin", perm.max=999))
   dbrda.terms$R2 <- dbrda.terms$SumOfSqs/sum(dbrda.terms$SumOfSqs)
   (dbrda.terms)
   dbrda.terms$names <- c("% Impervious","Prop.Neophytes","SR", "BNIs","Rao's Q", "resid")
   
   # Plot :
   ordi <- plot(dbrda.sp, type = "n",
                cex.axis = 0.8, 
                xlim = c(-3,3), ylim = c(-3,3))
   points(dbrda.sp, type = "p" , pch = 20,  col = col.code)
   
   arrows(0,0, dbrda.sp$CCA$biplot[,1]*2.5, dbrda.sp$CCA$biplot[,2]*2.5,
          length = 0.08,lwd = 1.5, angle = 35, col = "grey20")
   
   text(dbrda.sp$CCA$biplot[,1]*3.1,
        dbrda.sp$CCA$biplot[,2]*3,
        labels = c("% Impervious", "Prop.Neo","SR",
                   "BNIs","Rao"), cex = 0.8, adj = 0)
   
   mtext(3, text= tit,adj = 0,
         cex = 0.8, font = 2, line =1)
   
   # mtext(3, adj = 0, line = 0.1, cex= 0.7,
   #       text = "Jaccard dissimilarities")
   
   mtext(3, adj = 1, line = 0.1, cex= 0.65,
         text = as.expression(
            substitute("R"^2*" = "*a*b,
                       list(  a = round(dbrda.model[1,"R2"],2),
                              b = p2star(dbrda.model[1,"Pr(>F)"],
                                         marginal = TRUE)))))
   #Legend 
   list.expr = sapply(1:(nrow(dbrda.terms)-1), function(i){
      x <-  dbrda.terms$names[i]
      a <- round(dbrda.terms[i,"R2"],3)
      b <- p2star(dbrda.terms[i,"Pr(>F)"], marginal = TRUE)
      return(as.expression(substitute(x*": R"^2*" = "*a*b)))
   })
   
   legend("topleft", cex = 0.65,adj = 0,bty = "n",
          legend = list.expr )
}

# rows <- c("Allergenic\nspecies",
#           "Allergen\nmolecules",
#           "Allergen\nfamilies")
# mtext(2, at = c(0.25,0.5,0.75), text= rows, 
#       las = 1, adj = 1,
#       outer = TRUE, cex = 0.8, font = 2, line =1)
# 
# 

for (i in 1:3) {
   
   #Select dbrda:
   dbrda.sp <- list(dbrda.sp.two.jacc,
                    dbrda.mol.two.jacc, 
                    dbrda.allfam.two.jacc)[[i]]
   tit <- c("b)","d)","f)")[i]
   
   col.code <- c("#B2266270", "#48876270", "#3681af70")[i]
   
   # Calculate overall statistics:
   (dbrda.model <- anova(dbrda.sp, by=NULL, perm.max=999))
   (dbrda.model$R2 <-dbrda.model$SumOfSqs/
         sum(dbrda.model$SumOfSqs))
   
   # Calculate partial R2:
   (dbrda.terms <- anova(dbrda.sp, by="margin", perm.max=999))
   dbrda.terms$R2 <- dbrda.terms$SumOfSqs/sum(dbrda.terms$SumOfSqs)
   (dbrda.terms)
   dbrda.terms$names <- c("% Impervious","Prop.Neophytes","resid")
   
   # Plot :
   ordi <- plot(dbrda.sp, type = "n",
                cex.axis = 0.8, 
                xlim = c(-3,3), ylim = c(-3,3))
   points(dbrda.sp, type = "p" , pch = 20,  col = col.code)
   
   arrows(0,0, dbrda.sp$CCA$biplot[,1]*2.5, dbrda.sp$CCA$biplot[,2]*2.5,
          length = 0.08,lwd = 1.5, angle = 35, col = "grey20")
   
   text(dbrda.sp$CCA$biplot[,1]*3.1,
        dbrda.sp$CCA$biplot[,2]*3,
        labels = c("% Impervious", "Prop.Neo"), cex = 0.8, adj = 0)
   
   mtext(3, text= tit,adj = 0,
         cex = 0.8, font = 2, line =1)
   # mtext(3, adj = 0, line = 0.1, cex= 0.7,
   #       text = "Jaccard dissimilarities")
   
   mtext(3, adj = 1, line = 0.1, cex= 0.65,
         text = as.expression(
            substitute("R"^2*" = "*a*b,
                       list(  a = round(dbrda.model[1,"R2"],2),
                              b = p2star(dbrda.model[1,"Pr(>F)"],
                                         marginal = TRUE)))))
   #Legend 
   list.expr = sapply(1:(nrow(dbrda.terms)-1), function(i){
      x <-  dbrda.terms$names[i]
      a <- round(dbrda.terms[i,"R2"],3)
      b <- p2star(dbrda.terms[i,"Pr(>F)"], marginal = TRUE)
      return(as.expression(substitute(x*": R"^2*" = "*a*b)))
   })
   
   legend("bottomleft", cex = 0.65,adj = 0,bty = "n",
          legend = list.expr )
}

dev.off()
# Figure 2 alt:

# alt
png(file = "results/figure 2 alt.png",
    width = 18, height = 10,unit ="cm", res = 600)
par (
   mfrow = c(1, 2),
   bty = "l",
   mar = c(4, 3, 2, 1),
   oma = c(0, 0, 1, 0),
   mgp = c(2,0.4,0),
   tcl = - 0.2,
   las = 1,
   cex.axis = 0.8,
   xpd = TRUE
)

# Panel 1: interaction plot for CWM.PAV
tmp <- allergen_summary
f <- lm(CWM.pav ~ Seal_500 * prop.neo, tmp)
library(interactions)
p <- interactions::interact_plot(f, pred = Seal_500,  modx = prop.neo,
                                 legend.main = "Prop. Neophytes",
                                 plot.points = TRUE , interval = FALSE, 
                                 x.label = "% Impervious surfaces",
                                 y.label = "CWM.PAV",
                                 colors = "seagreen",
                                 point.size = 1,line.thickness = 0.8,
                                 partial.residuals = TRUE,
                                 facet.modx = FALSE
)
# Extract representative cut points for prop.neo :
(unique(p$data$prop.neo)[1:2] + unique(p$data$prop.neo)[2:3])/2
tmp$cut_prop.neo <- cut(allergen_summary$prop.neo, 
                        breaks = c(0, 0.0405,0.0886,1),include.lowest = TRUE)
cols <- colorRampPalette(c("peachpuff3","firebrick"))(3)

plot(CWM.pav ~ Seal_500 ,data = tmp,
     pch = 20,
     col = cols[tmp$cut_prop.neo],
     ann = FALSE)

add.stats(f, plot.abline = "no")
#mtext(3, text = "mean.PAV ~ %Imp.Surf * Prop.Neo",
#      adj = 0, cex = 0.7)

mtext( 1,text = "% Impervious surfaces", outer = FALSE,
       cex = 0.9,line = 2, las = 0)
mtext(2,text = substitute("CWM" [a],list( a="PAV")),outer = FALSE,
      cex = 0.9, line = 1.5, las = 0)

mtext(3,text = "(a)", adj = -0.2)

# Add illustrative regression lines: 
ltys <- c("dotted","dashed","solid")
lines(CWM.pav ~ Seal_500 ,
      data = p$data[p$data$modx_group == "+ 1 SD",],
      col = cols[3], lty = ltys[3], lwd = 1.2)
lines(CWM.pav ~ Seal_500 ,
      data = p$data[p$data$modx_group == "Mean",],
      col = cols[2], lty = ltys[2], lwd = 1.2)
lines(CWM.pav ~ Seal_500 ,
      data = p$data[p$data$modx_group == "- 1 SD",],
      col = cols[1] , lty = ltys[1], lwd = 1.2)

legend(2.5,20,
       legend = c("-SD","Mean","+SD"),
       title = "Proportion of Neophytes",
       lty = ltys, col = cols, lwd = 1.2,
       text.col = "grey10",
       bty =  "n",
       cex = 0.6,
       inset = -0.07
)

# Second panel : CWM archaeophyte PAVs
plot(CWM.arc.pav ~ Seal_500, data = allergen_summary,
     pch = 20,col = "grey",ann = FALSE)

add.stats(f <- lm(CWM.arc.pav ~ Seal_500, data = allergen_summary, na.action = "na.omit"),
          type = "lm",l.col = "black")
mtext( 1,text = "% Impervious surfaces", outer = FALSE,
       cex = 0.9,line = 2, las = 0)
mtext(2,text = substitute("Archaeophyte CWM" [a],list( a="PAV")),
      outer = FALSE,
      cex = 0.9, line = 1.5, las = 0)
mtext(3,text = "(b)", adj = -0.2)

dev.off()


})()