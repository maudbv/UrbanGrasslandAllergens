### Analyses of temporal spectrum of allergen production
library(vioplot)
library(np)

species_allergen$Introduction_status_Seitz2012<-factor(species_allergen$Introduction_status_Seitz2012, levels = c("I","A","N"))

## PHENO TRAIT Differences with status categories ####
## RUN ONCE with BOOTSTRAP (takes a few minutes)
### Statistical tests of difference between status
## NON PARAMETRIC TESTS: similar R2 found as for the linear model ANOVA
## But not normal distributions, so use non-parametric:
# 
# pheno.test <- NULL
# set.seed(50)
# 
# # Non-parametric regression:
# bw <-  np::npregbw(fl.beg ~ Introduction_status_Seitz2012,
#                data = species_allergen)
# summary(f.np<-npreg(bw))
# sig <- npsigtest(f.np,  boot.num = 999) #kernel
# # Get bootstrap confidence intervals:
# b <- plot(bw,
#           plot.errors.method="bootstrap",
#           plot.behavior = "plot-data",
#           plot.errors.type = "quantiles",
#           plot.errors.quantiles =c(0.05,0.95),
#           plot.errors.boot.num = 999,
#           plot.bxp = TRUE,
#           ylim = c(4,6.5)
#           )
# # Sometimes no longer working: reload all data with only 'np' package
# 
# pheno.test[["fl.beg"]] <- list(
#         stats = data.frame(mean = b$r1$mean,
#                            ci.low = b$r1$mean +b$r1$merr [,1],
#                            ci.high = b$r1$mean +b$r1$merr [,2]),
#         bws = f.np$bws,
#         r2 = f.np$R2,
#         P.boot = sig$P,
#         n.boot = sig$boot.num
# )
# 
# ### End of flowering:
# 
# # Non-parametric regression:
# bw <-  np::npregbw(fl.end ~ Introduction_status_Seitz2012,
#                    data = species_allergen)
# summary(f.np<-npreg(bw))
# sig <- npsigtest(f.np,  boot.num = 999, random.seed = 50 ) #kernel
# # Get bootstrap confidence intervals:
# b <- plot(bw,
#           plot.errors.method="bootstrap",
#           plot.behavior = "plot-data",
#           plot.errors.type = "quantiles",
#           plot.errors.quantiles =c(0.05,0.95),
#           plot.errors.boot.num = 999,
#           plot.bxp = TRUE,
#           ylim = c(7,10)
# )
# # Sometimes no longer working: reload all data with only 'np' package
# 
# pheno.test[["fl.end"]] <- list(
#         stats = data.frame(mean = b$r1$mean,
#                            ci.low = b$r1$mean +b$r1$merr [,1],
#                            ci.high = b$r1$mean +b$r1$merr [,2]),
#         bws = f.np$bws,
#         r2 = f.np$R2,
#         P.boot = sig$P,
#         n.boot = sig$boot.num
# )
# 
# ### Length of flowering:
# # Non-parametric regression:
# bw <-  np::npregbw(fl.period ~ Introduction_status_Seitz2012,
#                    data = species_allergen)
# summary(f.np<-npreg(bw))
# sig <- npsigtest(f.np,  boot.num = 999, random.seed = 50 ) #kernel
# b <- plot(bw,
#           plot.errors.method="bootstrap",
#           plot.behavior = "plot-data",
#           plot.errors.type = "quantiles",
#           plot.errors.quantiles =c(0.05,0.95),
#           plot.errors.boot.num = 999,
#           plot.bxp = TRUE,
#           ylim = c(2,5)
# )
# # Sometimes no longer working: reload all data with only 'np' package
# 
# pheno.test[["fl.period"]] <- list(
#         stats = data.frame(mean = b$r1$mean,
#                            ci.low = b$r1$mean +b$r1$merr [,1],
#                            ci.high = b$r1$mean +b$r1$merr [,2]),
#         bws = f.np$bws,
#         r2 = f.np$R2,
#         P.boot = sig$P,
#         n.boot = sig$boot.num
# )
# 
# 
# save(pheno.test, file = "data/phenology.bootstrap.test.Rdata")

load(file = "data/phenology.bootstrap.test.Rdata")

## NOT USED: SPECIES PHENO TRAIT Differences with allergenicity ####
# # Non-parametric regression: 
# summary(f.np<-npreg(fl.beg ~ as.factor(allergenicity), 
#                data =species_allergen))
# 
# npsigtest(f.np) #kernel regression significance test (bootstrap)
# # signif:  P < 2.22e-16 ***
# summary(f.np)  #R2 = 0.035 = NOT A BETTER FIT
# # plot(f.np, plot.errors.method="bootstrap")
# 
# # SPECIES Differences with allergenicity * status ####
# 
# tmp<-factor(species_allergen$exotic, levels = c(0,1), labels = c("Native", "Exotic"))
# tmp2<-factor(species_allergen$allergenicity, levels = c(0,1), labels = c("non-Allergenic", "Allergenic"))
# 
# par(mfrow = c(1,3), las = 2, mar = c(10,3,4,1))
# summary(lm(species_allergen$fl.beg ~ tmp2 * tmp2)) 
# TukeyHSD(aov(species_allergen$fl.beg ~ tmp * tmp2))
# boxplot(species_allergen$fl.beg ~ tmp + tmp2, range = 0, varwidth = T, main = "Onset of flowering (r2 = 0.04 *)")
# text(x = 1:4, y = 6.5,pos = 4, labels = c("a", "b", "b", "b"))
# 
# 
# 
# summary(lm(species_allergen$fl.end ~ tmp * tmp2)) # 
# TukeyHSD(aov(species_allergen$fl.end ~ tmp * tmp2))
# boxplot(species_allergen$fl.end ~ tmp * tmp2, range = 0, varwidth = T, main = "End of flowering (r2 = 0.04 *)")
# text(x = 1:4, y = 10.5,pos = 4, labels = c("a", "b", "a", "b"))
# 
# 
# summary(lm(species_allergen$fl.period ~ tmp : tmp2)) 
# TukeyHSD(aov(species_allergen$fl.period ~ tmp : tmp2))
# boxplot(species_allergen$fl.period ~ tmp * tmp2, range = 0, varwidth = T, main = "Period of flowering (r2 = 0.05 *)")
# text(x = 1:4, y = 6.5,pos = 4, labels = c("a", "b", "ac", "bc"))
# 


# ### NOT USED: Correlation matrices for community phenology :  ####
# 
# # Add environmental variables to compheno
# tmp <- cbind(compheno[rownames(allergen_summary),], allergen_summary)
# 
# # CWM of PHENO traits  # NOTHING HAPPENING WITH GRADIENT ####
# 
# cor.print(data = tmp[,c("CWM.fl.beg","CWM.fl.end", "CWM.fl.period",
#                         "CWM.nat.fl.beg","CWM.nat.fl.end","CWM.nat.fl.period",
#                         "CWM.neo.fl.beg","CWM.neo.fl.end","CWM.neo.fl.period",
#                         "Seal_500","prop.neo","CWM.pav","SR")],
#           return.print = TRUE,plotting = TRUE,method = "spearman")
#               
# 
# #  MEAN Flowering dates # NOTHING HAPPENING WITH GRADIENT ####
# cor.print(data = tmp[,c("mean.fl.beg","mean.fl.end", "mean.fl.period",
#                         "mean.nat.fl.beg","mean.nat.fl.end","mean.nat.fl.period",
#                         "mean.neo.fl.beg","mean.neo.fl.end","mean.neo.fl.period",
#                         "Seal_500","prop.neo","mean.pav", "SR")],
#           return.print = TRUE,plotting = TRUE,method = "spearman")
# 
# 
# # Overall flowering dates : OVERALL end of flowering later with prop.neo  ####
# cor.print(data = tmp[,c("beg.fl","end.fl","beg.fl.5pc","end.fl.5pc",
#                         "period.fl","period.fl.5pc","peak.fl",
#                         "cum.fl","cum.fl.wtd",
#                         "Seal_500","prop.neo","BNIs","Rao","mean.pav", "SR")],
#           return.print = TRUE,plotting = TRUE,method = "spearman")
# 
# 
# # Neophyte flowering dates: later peak of flowering with prop.neo ####
# cor.print(data = tmp[,c("beg.fl.neo","end.fl.neo",
#                         "period.fl.neo","peak.fl.neo","cum.fl.neo",
#                         "cum.fl.neo.wtd",
#                         "Seal_500","prop.neo","BNIs","Rao","mean.pav", "SR")],
#           return.print = TRUE,plotting = TRUE,method = "spearman")
# 
# 
# # natives flowering dates:  ####
# cor.print(data = tmp[,c("beg.fl.nat","end.fl.nat",
#                         "period.fl.nat","peak.fl.nat","cum.fl.nat",
#                         "cum.fl.nat.wtd",
#                         "Seal_500","prop.neo","BNIs","Rao","mean.pav", "SR")],
#           return.print = TRUE,plotting = TRUE,method = "spearman")
# 
# # archaeo flowering dates:  ####
# cor.print(data = tmp[,c("beg.fl.arc","end.fl.arc",
#                         "period.fl.arc","peak.fl.arc","cum.fl.arc",
#                         "cum.fl.arc.wtd",
#                         "Seal_500","prop.neo","BNIs","Rao","mean.pav", "SR")],
#           return.print = TRUE,plotting = TRUE,method = "spearman")
# 
# 
# # exotic flowering dates: peak is later in urban + prop.neo  ####
# cor.print(data = tmp[,c("beg.fl.exo","end.fl.exo",
#                         "period.fl.exo","peak.fl.exo","cum.fl.exo",
#                         "cum.fl.exo.wtd",
#                         "Seal_500","prop.neo","BNIs","Rao","mean.pav", "SR")],
#           return.print = TRUE, plotting = TRUE,method = "spearman")

## Illustrate correlations Main results:  ####
tmp <- cbind(compheno[rownames(allergen_summary),], allergen_summary)

# Remove entries for plots with no allergenic species from a group:
tmp[tmp$cum.fl.neo == 0, 
    c("beg.fl.neo","end.fl.neo","period.fl.neo","peak.fl.neo")] <- NA
tmp[tmp$cum.fl.arc == 0, 
    c("beg.fl.arc","end.fl.arc","period.fl.arc","peak.fl.arc")] <- NA
tmp[tmp$cum.fl.exo == 0, 
    c("beg.fl.exo","end.fl.exo","period.fl.exo","peak.fl.exo")] <- NA
tmp[tmp$cum.fl.nat == 0, 
    c("beg.fl.nat","end.fl.nat","period.fl.nat","peak.fl.nat")] <- NA

tmp <- tmp [, c("Seal_500","prop.neo",
                "beg.fl","end.fl","period.fl","peak.fl", 
                "beg.fl.neo","end.fl.neo","period.fl.neo","peak.fl.neo",
                "beg.fl.arc","end.fl.arc","period.fl.arc","peak.fl.arc",
                "beg.fl.exo","end.fl.exo","period.fl.exo","peak.fl.exo",
                "beg.fl.nat","end.fl.nat","period.fl.nat","peak.fl.nat")]
cor.pheno <- cor.print(tmp, return.print = FALSE,
                       method = "spearman", plotting = TRUE)

# cor.pheno <- data.frame(matrix(NA, nrow = 6, ncol = 6))
# names(cor.pheno) <- c("trait","var","nobs", "rho","P","df")
# cor.pheno$trait <- c(rep("end.fl",2),
#                      rep("peak.fl.neo",2),
#                      rep("peak.fl.exo",2))
# cor.pheno$var <- rep(c("Seal_500","prop.neo"),3)
# 
# par(mfrow = c(3,2), mar = c(4,4,1,1))
# 
# for (i in 1:6) {
# tmp2 <- tmp               
# 
# # Make sure the phenological dates make sense, 
# # ie there are neophyte or exotic allergenics are present
# if (i %in% 3:4) tmp2 <- tmp[which(tmp$cum.fl.neo >0),]
# if (i %in% 5:6) tmp2 <- tmp[which(tmp$cum.fl.exo >0),]
# 
# plot(as.formula(paste("jitter(", cor.pheno$trait[i],")",
#                       "~ ", cor.pheno$var[i])),
#      tmp2,
#      pch = 20,col = "darkgrey", cex = 1)
# f <- cor.test(tmp2[,cor.pheno$trait[i]], tmp2[,cor.pheno$var[i]],
#               method = "spearman")
# cor.pheno[i,3:6] <- c(nrow(tmp2),unlist(add.stats(f)))
# }
