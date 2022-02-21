# Allergenic Diversity vs. Community diversity
# Creates matrix of spearman correlation for Allergen richness
# with species richness and Rao.
library(Hmisc) 
library(corrplot)
library(MuMIn)
library(r2glmm)

## CORRELATION MATRICES with allergenic SPECIES RICHNESS #####
tmp <- allergen_summary[,c("SR","SR.nat", "SR.arch", "SR.neo",
                           "all.num","severe.all.num",
                           "all.num.neo","all.num.arc","all.num.nat",
                           "all.num.exo",
                           "cover.all","cover.nat.all",
                           "cover.arc.all","cover.neo.all","cover.exo.all")]

cor.SRall.print <- cor.print(tmp, plotting = FALSE)
cor.SRall <- cor.print(tmp,return.print = FALSE, plotting = FALSE)

write.csv(cor.SRall.print,
          "results/Spearman Cor table diversity vs. allergenic SR.csv")



# Mean statistics to report: ####

library(doBy)
SRall.stats <- data.frame(
  
prop <- data.frame(
  prop.all.num = c(mean = mean(allergen_summary$all.num/allergen_summary$SR,
                   na.rm = TRUE),
              sd = sd(allergen_summary$all.num/allergen_summary$SR,
               na.rm = TRUE),
              min = min(allergen_summary$all.num/allergen_summary$SR,
                 na.rm = TRUE),
              max = max(allergen_summary$all.num/allergen_summary$SR,
                 na.rm = TRUE)
  ),
  prop = sapply(c("all.num.nat","all.num.arc", "all.num.neo","all.num.exo"),
       function(x) {
         out = data.frame(
          mean = mean(allergen_summary[,x]/allergen_summary$all.num,
              na.rm = TRUE),
          sd = sd(allergen_summary[,x]/allergen_summary$all.num,
                    na.rm = TRUE),
          min = min(allergen_summary[,x]/allergen_summary$all.num,
                  na.rm = TRUE),
          max = max(allergen_summary[,x]/allergen_summary$all.num,
                    na.rm = TRUE)
         )
         return(out)
       })),

SR <- sapply(c("all.num","all.num.nat","all.num.arc", "all.num.neo","all.num.exo"),
       function(x) {
         out = data.frame(
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


Cover <- sapply(c("cover.all","cover.nat.all",
               "cover.arc.all", "cover.neo.all","cover.exo.all"),
             function(x) {
               out = data.frame(
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

Cover.rural <- sapply(c("cover.all","cover.nat.all",
                  "cover.arc.all", "cover.neo.all","cover.exo.all"),
                function(x) {
                  tmp <- allergen_summary[
                    (allergen_summary$Seal_500)< 2,]
                  out = data.frame(
                    mean = mean( tmp[,x],na.rm = TRUE),
                    sd = sd( tmp[,x],na.rm = TRUE),
                    min = min( tmp[,x],na.rm = TRUE),
                    max = max( tmp[,x],na.rm = TRUE)
                  )
                  return(out)
                }),

Cover.urban <- sapply(c("cover.all","cover.nat.all",
                        "cover.arc.all", "cover.neo.all","cover.exo.all"),
                      function(x) {
                        tmp <- allergen_summary[
                          (allergen_summary$Seal_500)> 50,]
                        out = data.frame(
                          mean = mean( tmp[,x],na.rm = TRUE),
                          sd = sd( tmp[,x],na.rm = TRUE),
                          min = min( tmp[,x],na.rm = TRUE),
                          max = max( tmp[,x],na.rm = TRUE)
                        )
                        return(out)
                      })
)
