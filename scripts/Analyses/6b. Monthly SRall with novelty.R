# Analyse trends in mean monthly flowering allergenic species 
### unweighted cumulated molecules: only Neophyte significant

glms.SRall.month <- (function(exclude.absences = TRUE,
                        show.plots = TRUE){
  
  require(Hmisc) 
  require(corrplot)
  require(MuMIn)
  require(r2glmm)
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
  
  # Create main data set:
  dat <- cbind(compheno[rownames(plot_summary),], plot_summary)
  
  
  # Trends in mean monthly SRall - LM ####
  for (i in 1:5) {
    tmp <- dat
    y <- c("cum.fl", "cum.fl.nat","cum.fl.arc","cum.fl.neo","cum.fl.exo")[i]
    g <- c("all", "nat","arc","neo","exo")[i]
    type <- "richness"
    
    if (exclude.absences) {
      z <- c("SR","SR.nat","SR.arch","SR.neo","SR.exo")[i]
      tmp <- tmp[which(tmp[,z]>0 & !is.na(tmp[,y])),]
    }
    n = length(tmp[,y])
    fit <- fit.lms(dataset = tmp, var = y,
                            plot.graphs = show.plots)
    
    glms.table[i,] <- c(g,type, n, fit$result.table)
    rownames(glms.table)[i] <- y
    
    glms.pred <- rbind(glms.pred, 
                       cbind(group = g, type = type, n = n, 
                             fit$predictor.table))
    
    all.models [[i]] <- fit$models
    names(all.models )[i]<- y
    rm(tmp)
  }
  
  # trends in wtd.mean monthly SRall - LM ####
  for (i in 1:5) {
    tmp <- dat
    y <- c("cum.fl.wtd", "cum.fl.nat.wtd","cum.fl.arc.wtd","cum.fl.neo.wtd","cum.fl.exo.wtd")[i]
    g <- c("all", "nat","arc","neo","exo")[i]
    type <- "cover"
  
    
    # transform the variable:
    name <- paste("trans.",y, sep = "")
    tmp$trans.y <- sqrt(tmp[,y])
    names(tmp)[names(tmp) == "trans.y"] <- name
    y <- name

    if (exclude.absences) {
      z <- c("SR","SR.nat","SR.arch","SR.neo","SR.exo")[i]
      tmp <- tmp[which(tmp[,z]>0 & !is.na(tmp[,y])),]
    }
    n = length(tmp[,y])
    
    fit <- fit.lms(dataset = tmp, var = y,
                   plot.graphs = show.plots)
    
    glms.table[i + 5,] <- c(g,type, n, fit$result.table)
    rownames(glms.table)[i + 5] <- y
    
    glms.pred <- rbind(glms.pred, 
                       cbind(group = g, type = type, n = n, 
                             fit$predictor.table))
    
    all.models [[i + 5]] <- fit$models
    names(all.models )[i + 5]<- y
    rm(tmp)
  }
  
  # Format main result table ####
  glms.table$group <- factor(glms.table$group,
                             levels = c("all","nat","arc","neo","exo"))
  glms.table$type <- factor(glms.table$type,
                            levels = c("richness","proportion","cover"))
  glms.table$level <- "Species"
  glms.table <- glms.table[order(glms.table$group, glms.table$type),]
  
  # # Output predictors
  
  # REturn: 
  return(list(glms.table = glms.table,
              glms.best.pred = glms.pred
  ))
  
})()

# Export result table: ####
write.csv(glms.SRall.month$glms.table,
          "results/Monthly SRall models table.csv")
write.csv(glms.SRall.month$glms.best.pred,
          "results/Monthly SRall  best predictor table.csv")


