## Analysing trends in Allergenic molecule diversity along gradients

library(vegan)
library(jtools)
library(interactions)
library(MuMIn)
library(MASS) # for negative binomials
library(performance)
library(rsq)
# Community Allergenicity along the gradient
# Fitting generalized linear models to allergenic species richness, proportion and cover
# Extract both the "best model" with interaction (prop.neo*seal_500)
# and three independant models for each of the three predictors: Seal_500, prop.neo, BNIs

glms.AR <- (function(exclude.absences = TRUE,
                     show.plots = FALSE){
  
  require(Hmisc) 
  require(corrplot)
  require(MuMIn)
  require(r2glmm)
  require(performance)
  options(na.action = "na.fail") 
  stopifnot(!all(c(is.null(fit.binom.glms),
                   is.null(fit.negbin.glms))))
  
  # Create a table of results:   ####
  glms.table <- data.frame(matrix(NA, nrow = 0,ncol =24) )
  colnames(glms.table) <- c("group", "type","n.obs","class","var",
                            "Best.model", "df.resid","P.lrt","R2",
                            "PropNeo.coef", "PropNeo.se",
                            "PropNeo.df","PropNeo.P","PropNeo.R2",
                            "Seal.coef", "Seal.se","Seal.df","Seal.P","Seal.R2",
                            "BNIs.coef", "BNIs.se","BNIs.df","BNIs.P","BNIs.R2")
  
  # Create a list to store the models: 
  all.models <- list()
  
  # Create table for best model predictors: 
  glms.pred<- data.frame(matrix(NA, nrow = 0,ncol =12 ) )
  colnames(glms.pred) <- c("group", "type","n.obs",
                           "class", "var","pred","est","se","z","P","R2beta","Pchi")
  
  # Trends in ALLERGEN RICHNESS = negative binomial ####
  for (i in 1:5) {
    tmp <- allergen_summary
    y <- c("nb.mol", "nb.mol.nat","nb.mol.arc","nb.mol.neo","nb.mol.exo")[i]
    g <- c("all", "nat","arc","neo","exo")[i]
    type <- "richness"
    
    if (exclude.absences) {
      z <- c("SR","SR.nat","SR.arch","SR.neo","SR.exo")[i]
      tmp <- tmp[which(tmp[,z]>0 & !is.na(tmp[,y])),]
    }
    n = length(tmp[,y])
    fit <- fit.negbin.glms(dataset = tmp, var = y,
                           BNI.include = F,plot.graphs = show.plots)
    
    glms.table[i,] <- c(g,type, n, fit$result.table)
    rownames(glms.table)[i] <- y
    
    glms.pred <- rbind(glms.pred, 
                       cbind(group = g, type = type, n = n, 
                             fit$predictor.table))
    
    all.models [[i]] <- fit$models
    names(all.models )[i]<- y
    rm(tmp)
  }
  
  # trends in allergenic species proportions  ####
  for (i in 1:4) {
    tmp <- allergen_summary
    x <- c("nb.mol.nat","nb.mol.arc","nb.mol.neo","nb.mol.exo")[i]
    y <-  "nb.mol"
    nam <- paste("prop",x, sep="_")
    g <- c("nat","arc","neo","exo")[i]
    type <- "proportion"
    
    
    if (exclude.absences) {
      z <- c("SR.nat","SR.arch","SR.neo","SR.exo")[i]
      tmp <- tmp[which(tmp[,z]>0 & !is.na(tmp[,x])),]
    }
    
    n = length(tmp[,y])
    
    # Fit Binomials for proportions
    fit <- fit.binom.glms(tmp, x, y, BNI.include = F,plot.graphs = show.plots)
    glms.table[i + 5,] <- c(g,type, n, fit$result.table)
    rownames(glms.table)[i + 5] <- nam
    
    
    glms.pred <- rbind(glms.pred, 
                       cbind(group = g, type = type, n = n, 
                             fit$predictor.table))
    
    all.models [[i + 5]] <- fit$models
    names(all.models )[i + 5]<- nam
    rm(tmp)
  }
  
  
  # Format main result table ####
  glms.table$group <- factor(glms.table$group,
                             levels = c("all","nat","arc","neo","exo"))
  glms.table$type <- factor(glms.table$type,
                            levels = c("richness","proportion","cover"))
  glms.table <- glms.table[order(glms.table$group, glms.table$type),]
  
  # # Output predictors
  
  # REturn: 
  return(list(glms.table = glms.table,
              glms.best.pred = glms.pred
  ))
  
})()

# Export result table: ####
write.csv(glms.AR$glms.table,
          "results/allergen molecule richness models table.csv")
write.csv(glms.AR$glms.best.pred,
          "results/allergen molecule richness best predictor table.csv")

