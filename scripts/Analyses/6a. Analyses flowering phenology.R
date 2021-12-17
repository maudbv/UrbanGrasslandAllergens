### Analyses of temporal spectrum of allergen production
library(vioplot)
library(np)

# make sure the order of factor levels is Indigenous<Archaeophyte<Neophytes
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


## Illustrate correlations main trends:  ####
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
                       method = "spearman", plotting = FALSE)


