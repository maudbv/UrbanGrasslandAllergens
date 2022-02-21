## Analyzing trends in Allergenic FAMILY molecule diversity along NOVELTY gradients

glms.allfam <- (function(exclude.absences = TRUE,
                        show.plots = FALSE){
  
  require(MuMIn)
  require(r2glmm)
  require(MASS)
  require(performance)
  
  options(na.action = "na.fail") 
  stopifnot(!all(c(is.null(fit.poisson.glms),
                   is.null(fit.binom.glms),
                   is.null(fit.lms))))
  
  # Create a table of results:   ####
  glms.table <- data.frame(matrix(NA, nrow = 0,ncol =19) )
  colnames(glms.table) <- c("group", "type","n.obs","class","var",
                            "Best.model", "df.resid","P.lrt","R2",
                            "PropNeo.coef", "PropNeo.se",
                            "PropNeo.df","PropNeo.P","PropNeo.R2",
                            "Seal.coef", "Seal.se","Seal.df","Seal.P","Seal.R2")
  
  # Create a list to store the models: 
  all.models <- list()
  
  # Create table for best models predictors: 
  glms.pred<- data.frame(matrix(NA, nrow = 0,ncol =12 ) )
  colnames(glms.pred) <- c("group", "type","n.obs",
                           "class", "var","pred","est","se","z","P","R2beta","Pchi")
  
  # Trends in allfam RICHNESS - Poisson GLMs ####
  for (i in 1:5) {
    tmp <- allergen_summary
    y <- c("nb.allfam", "nb.allfam.nat","nb.allfam.arc","nb.allfam.neo","nb.allfam.exo")[i]
    g <- c("all", "nat","arc","neo","exo")[i]
    type <- "richness"
    
    if (exclude.absences) {
      z <- c("SR","SR.nat","SR.arch","SR.neo","SR.exo")[i]
      tmp <- tmp[which(tmp[,z]>0 & !is.na(tmp[,y])),]
    }
    n = length(tmp[,y])
    fit <- fit.poisson.glms(tmp, y,plot.graphs = show.plots)
    
    glms.table[i,] <- c(g,type, n, fit$result.table)
    rownames(glms.table)[i] <- y
    
    glms.pred <- rbind(glms.pred, 
                       cbind(group = g, type = type, n = n, 
                             fit$predictor.table))
    
    all.models [[i]] <- fit$models
    names(all.models )[i]<- y
  }
  
   # trends in allfam proportions - Binomial GLM  ####
  for (i in 1:4) {
    tmp <- allergen_summary
    x <- c( "nb.allfam.nat","nb.allfam.arc","nb.allfam.neo","nb.allfam.exo")[i]
    y <-  "nb.allfam"
    nam <- paste("prop",x, sep="_")
    g <- c("nat","arc","neo","exo")[i]
    type <- "proportion"
    
    
    if (exclude.absences) {
      z <- c("SR.nat","SR.arch","SR.neo","SR.exo")[i]
      tmp <- tmp[which(tmp[,z]>0 & !is.na(tmp[,x])),]
    }
    
    n = length(tmp[,y])
    
    # Fit Binomials for proportions
    fit <- fit.binom.glms(tmp, x, y, plot.graphs = show.plots)
    glms.table[i + 10,] <- c(g,type, n, fit$result.table)
    rownames(glms.table)[i + 10] <- nam
    
    
    glms.pred <- rbind(glms.pred, 
                       cbind(group = g, type = type, n = n, 
                             fit$predictor.table))
    
    all.models [[i + 10]] <- fit$models
    names(all.models )[i + 10]<- nam
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
write.csv(glms.allfam $glms.table,
          "results/allergenic families models table.csv")
write.csv(glms.allfam $glms.best.pred,
          "results/allergenic families best predictor table.csv")


