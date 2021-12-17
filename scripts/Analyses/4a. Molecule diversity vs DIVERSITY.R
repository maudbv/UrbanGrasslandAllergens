## Trends in Allergen richness with DIVERSITY 
library(Hmisc)
library(corrplot)


## CORRELATION MATRIX for Allergen Molecule Richness ######
tmp <- allergen_summary[,c("SR","SR.nat", "SR.arch", "SR.neo","all.num",
                           "FR","FR.nat", "FR.arc", "FR.neo",
                           "nb.mol","nb.mol.neo","nb.mol.arc","nb.mol.nat", "nb.mol.exo",
                           "nb.allfam","nb.allfam.neo","nb.allfam.arc","nb.allfam.nat", "nb.allfam.exo")]

cor.AR.print <- cor.print (data = tmp, method = "spearman",
                           plotting = FALSE)
cor.AR <- cor.print(tmp,return.print = FALSE, plotting = FALSE)

write.csv(cor.AR.print,
          "results/Spearman Cor table diversity vs. AR.csv")

# Summary statistics for Allergen Richness :

# Mean statistics to report: ####

library(doBy)
AR.stats <- data.frame(
  
  prop <- data.frame(
    prop.nb.mol = c(median = NA,
                    mean = NA,
                     sd = NA,
                     min = NA,
                     max = NA
    ),
    prop = sapply(c("nb.mol.nat","nb.mol.arc", "nb.mol.neo","nb.mol.exo"),
                  function(x) {
                    out = data.frame(
                      median = median(allergen_summary[,x]/allergen_summary$nb.mol,na.rm = TRUE),
                      mean = mean(allergen_summary[,x]/allergen_summary$nb.mol,
                                  na.rm = TRUE),
                      sd = sd(allergen_summary[,x]/allergen_summary$nb.mol,
                              na.rm = TRUE),
                      min = min(allergen_summary[,x]/allergen_summary$nb.mol,
                                na.rm = TRUE),
                      max = max(allergen_summary[,x]/allergen_summary$nb.mol,
                                na.rm = TRUE)
                    )
                    return(out)
                  })),
  
  SR <- sapply(c("nb.mol","nb.mol.nat","nb.mol.arc", "nb.mol.neo","nb.mol.exo"),
               function(x) {
                 out = data.frame(
                   median = median(allergen_summary[,x],na.rm = TRUE),
                   mean = mean(allergen_summary[,x],na.rm = TRUE),
                   sd = sd(allergen_summary[,x],
                           na.rm = TRUE),
                   min = min(allergen_summary[,x],
                             na.rm = TRUE),
                   max = max(allergen_summary[,x],
                             na.rm = TRUE)
                 )
                 return(out)
               }),
  
  # Number of Allergen families 
  #(!counting all the zeros = underestimate Neophytes)
  SR <- sapply(c("nb.allfam","nb.allfam.nat","nb.allfam.arc", "nb.allfam.neo","nb.allfam.exo"),
               function(x) {
                 out = data.frame(
                   median = median(allergen_summary[,x],na.rm = TRUE),
                   mean = mean(allergen_summary[,x],na.rm = TRUE),
                   sd = sd(allergen_summary[,x],
                           na.rm = TRUE),
                   min = min(allergen_summary[,x],
                             na.rm = TRUE),
                   max = max(allergen_summary[,x],
                             na.rm = TRUE)
                 )
                 return(out)
               })
  
  # Cover <- sapply(c("cover.mol","cover.mol.nat",
  #                   "cover.mol.arc", "cover.mol.neo","cover.mol.exo"),
  #                 function(x) {
  #                   out = data.frame(
  #                     mean = mean(allergen_summary[,x],na.rm = TRUE),
  #                     sd = sd(allergen_summary[,x],
  #                             na.rm = TRUE),
  #                     min = min(allergen_summary[,x],
  #                               na.rm = TRUE),
  #                     max = max(allergen_summary[,x],
  #                               na.rm = TRUE)
  #                   )
  #                   return(out)
  #                 }),
  # 
  # Cover.rural <- sapply(c("cover.mol","cover.mol.nat",
  #                         "cover.mol.arc", "cover.mol.neo","cover.mol.exo"),
  #                       function(x) {
  #                         tmp <- allergen_summary[
  #                           (allergen_summary$Seal_500)< 2,]
  #                         out = data.frame(
  #                           mean = mean( tmp[,x],na.rm = TRUE),
  #                           sd = sd( tmp[,x],na.rm = TRUE),
  #                           min = min( tmp[,x],na.rm = TRUE),
  #                           max = max( tmp[,x],na.rm = TRUE)
  #                         )
  #                         return(out)
  #                       }),
  # 
  # Cover.urban <- sapply(c("cover.mol","cover.mol.nat",
  #                         "cover.mol.arc", "cover.mol.neo","cover.mol.exo"),
  #                       function(x) {
  #                         tmp <- allergen_summary[
  #                           (allergen_summary$Seal_500)> 50,]
  #                         out = data.frame(
  #                           mean = mean( tmp[,x],na.rm = TRUE),
  #                           sd = sd( tmp[,x],na.rm = TRUE),
  #                           min = min( tmp[,x],na.rm = TRUE),
  #                           max = max( tmp[,x],na.rm = TRUE)
  #                         )
  #                         return(out)
  #                       })
)
