# Functions to fit LM and GLMs on allergen diversity data for Berlin grasslands
# with three possible predictors: 
# % Impervious surfaces (Seal_500)
# Proportion of neophytes (prop.neo)


require(Hmisc) 
require(corrplot)
require(MuMIn)
require(r2glmm)
require(performance)
options(na.action = "na.fail") 

# Utility FUNCTION extract predictors from a glm
extract.preds <- function(mod = f.best) {
  
  require (r2glmm)
  
  coef <- summary(mod)$coef
  anv <- anova(mod,test = "LRT")
  anv <- anv[rownames(coef),]
  
  # If model is the null model : stats = NA
  if (nrow(coef) == 1) {
    output <-  cbind(pred = NA,
                     t(rep(NA,4)),
                     r2b = NA,
                     P.chi = NA
    )
    rownames(output) <- "null"
  } else { 
  # If it is a negative binomial model: 
  if ("negbin" %in% class(mod)) {
    # refit as a glm with the same theta, to be able to apply r2beta
    f2 <-  glm(as.formula(formula(mod)),
               data = mod$model,
               family = negative.binomial(theta = mod$theta,
                                          link = "log"))
    source('/Users/maud/Dropbox/Work/doc boulot/postdoc Berlin/R projects/utility functions/r2beta.glm.modif.R')
    r2b <-   as.data.frame(r2beta.glm.local(f2))
  } else {
    # If it is a non-null glm or lm models: 
    f2 <- mod
    r2b <-   as.data.frame(r2beta(f2))
  }
    
    # Reorder R2 beta table: 
    r2b <-  r2b[-1,]
    r2b <-  r2b[match(rownames(coef)[-1],as.character(r2b$Effect)),]
    
    # Calculate statistics if 2 or more predictors (#coef = predictors + 1)
  if (nrow(coef)>2) {
    output <- cbind(pred = as.character(r2b$Effect),
                    coef[2:nrow(coef),],
                    r2b = r2b$Rsq,
                    P.chi = anv[2:nrow(coef),5]
    )
  }
    # If only one predictor: 
  if (nrow(coef) == 2) {
    output <-  cbind(pred = as.character(r2b$Effect),
                     t(coef[2:nrow(coef),]),
                     r2b = r2b$Rsq,
                     P.chi = anv[2:nrow(coef),5]
    )
    rownames(output) <- rownames(anv)[2:nrow(coef)]
  }
  } 
  # store info abou the model class in the output:
  output <- cbind(class(mod)[1], output)
  
  # rename columns (z: for glms, t: for lms)
  colnames(output) <- c("class", "pred","Estimate","Std. Error", "z_or_t",
                        "Pr(z_or_t)","r2b","P.chi")
  return(output)
}


# Poisson GLMs
fit.poisson.glms <- function(dataset, var,
                             plot.graphs = FALSE) {
  require(MuMIn)
  require(performance)
  
  dataset$y <- dataset[, var]
  
  # Poisson GLM : 
  mods.list <- list(
    f0 = glm(y ~ 1, data = dataset, family = poisson),
    f1a = glm(y ~ prop.neo, data = dataset, family = poisson),
    f1b = glm(y ~ Seal_500 , data = dataset, family = poisson),
    f2a = glm(y ~ Seal_500  + prop.neo  , data = dataset, family = poisson),
    f4a = glm(y ~ Seal_500 * prop.neo , data = dataset, family = poisson)
  )
  
  # Extract best models: 
 
  sel = model.sel(mods.list)
  mods <- rownames(subset(sel, delta <2))
  
  f.best <- mods.list[[mods[1]]]

  if ("f0" %in% mods) {
    if (mods[1] == "f0") {
      f.best <-  mods.list$f0
    } else {
      # check that best model is significantly better than null:
      p.fbest <- anova(mods.list$f0,f.best, test = "LRT")$P[2] 
      if (p.fbest> 0.05) {
        f.best <-  mods.list$f0
      }
    }
  }
  
  # AVONA tables with Likelihood ratio test: 
  aov.mods <- sapply(names(mods.list), function(x) {
    anova(mods.list$f0, mods.list[[x]])[2,]
  }, simplify = FALSE)
  
  # Summary of coefficients:
  coef.mods <- sapply(names(mods.list), function(x) {
    summary(mods.list[[x]])$coefficient
  }, simplify = FALSE)
  
  # Add f.best to the list: 
  mods.list$f.best <- f.best
  
  # AVONA tables with Likelihood ratio test: 
  aov.mods <- sapply(names(mods.list), function(x) {
    anova(mods.list$f0, mods.list[[x]],
          test = "LRT")[2,]
  }, simplify = FALSE)
  
  # Summary of coefficients:
  coef.mods <- sapply(names(mods.list), function(x) {
    summary(mods.list[[x]])$coefficient
  }, simplify = FALSE)
  
  # Store statistics
  out.indiv <- round(c(r2_nagelkerke(f.best)[[1]],
                       coef.mods$f1a[2,1:2],
                       aov.mods$f1a$`Resid. Df`, aov.mods$f1a$`Pr(>Chi)`,
                       r2_nagelkerke(mods.list$f1a)[[1]],
                       coef.mods$f1b[2,1:2],
                       aov.mods$f1b$`Resid. Df`, aov.mods$f1b$`Pr(>Chi)`,
                       r2_nagelkerke(mods.list$f1b)[[1]]),
                     digits = 4)
  
  out.indiv <-  c(class(f.best)[1], var,
                                as.character(formula(f.best))[3],
                                aov.mods$f.best$`Resid. Df`,
                                aov.mods$f.best$`Pr(>Chi)`,
                                out.indiv)
  
  # other vector for predictors in best model
  out.best.pred <- cbind(var, extract.preds(f.best))
  
  # Plot graphs:
  if (plot.graphs) {
    par (mfrow = c(1,2),
         mar = c(4,3,2,2),
         oma =c(1,2,0,0))
    
    plot(y ~ Seal_500,
         dataset,
         pch = 20,
         ylab = "",
         xlab = "% Impervious surfaces")
    add.stats(mods.list$f1b, type = "glm")
    
    plot( y ~ prop.neo,
          data= dataset,
          pch = 20,
          ylab = "",
          xlab = "Proportion of Neophytes")
    add.stats(mods.list$f1a, type = "glm")
    
    mtext(2, text = var,
          outer = TRUE, cex = 0.7)
    }
  return(list(result.table = out.indiv,
               predictor.table = out.best.pred,
               models = mods.list))
}

# Binomial GLMs: used for proportion data
fit.binom.glms <- function(dataset, var, tot,
                           plot.graphs = FALSE) {
  
  dataset$x <- dataset[, var]
  dataset$y <- dataset[, tot]
  
  # GLM : 
  quasi.mods.list <- list(
    f0 = glm(cbind(x,y-x) ~ 1, data = dataset, family = quasibinomial),
    f1a = glm(cbind(x,y-x) ~ prop.neo, data = dataset, family = quasibinomial),
    f1b = glm(cbind(x,y-x) ~ Seal_500 , data = dataset, family = quasibinomial),
    f2a = glm(cbind(x,y-x) ~ Seal_500  + prop.neo  , data = dataset, family = quasibinomial),
    f4a = glm(cbind(x,y-x) ~ Seal_500 * prop.neo , data = dataset, family = quasibinomial)
   )
  
  disp.param <- lapply( quasi.mods.list, function(x) summary(x)$dispersion)
  
  mods.list <- list(
    f0 = glm(cbind(x,y-x) ~ 1, data = dataset, family = binomial),
    f1a = glm(cbind(x,y-x) ~ prop.neo, data = dataset, family = binomial),
    f1b = glm(cbind(x,y-x) ~ Seal_500 , data = dataset, family = binomial),
    
    f2a = glm(cbind(x,y-x) ~ Seal_500  + prop.neo  , data = dataset, family = binomial),
    f4a = glm(cbind(x,y-x) ~ Seal_500 * prop.neo , data = dataset, family = binomial)
   )
  disp.param <- lapply( mods.list, function(x) summary(x)$dispersion)
  
  # Extract best models: 
  sel = model.sel(mods.list)
  
  mods <- rownames(subset(sel, delta <2))

  f.best <- mods.list[[mods[1]]]

  if ("f0" %in% mods) {
    if (mods[1] == "f0") {
      f.best <-  mods.list$f0
    } else {
      # if the best model is still significant
      if (anv$P [2] > 0.05) {
        f.best <-  mods.list$f0
      }
    }
  }
  
 
  # Add f.best to the list: 
  mods.list$f.best <- f.best
  
  # AVONA tables with Likelihood ratio test: 
  aov.mods <- sapply(names(mods.list), function(x) {
    anova(mods.list$f0, mods.list[[x]] ,
          dispersion = disp.param[[x]],
          test = "LRT")[2,]
  }, simplify = FALSE)
  
  # Summary of coefficients:
  coef.mods <- sapply(names(mods.list), function(x) {
    summary(mods.list[[x]] ,
            dispersion = disp.param[[x]])$coefficient
  }, simplify = FALSE)
  
  
  # R2 values
  # the Nagelkerke tends to give really high values
  r2.best <- r2_tjur(mods.list$f.best)[[1]]
  r2.f1a <-  r2_tjur(mods.list$f1a)[[1]]
  r2.f1b <- r2_tjur(mods.list$f1b)[[1]]
 

  # alternative based on r2beta
  # r2.best <- r2glmm::r2beta(f.best)[1,"Rsq"]
  # r2.f1a <- r2beta(mods.list$f1a)[1,"Rsq"]
  # r2.f1b <- r2beta(mods.list$f1b)[1,"Rsq"]

  
  # Store statistics
  out.indiv <- round(c(r2.best,
                       coef.mods$f1a[2,1:2],
                       aov.mods$f1a$`Resid. Df`, aov.mods$f1a$`Pr(>Chi)`,
                       r2.f1a,
                       coef.mods$f1b[2,1:2],
                       aov.mods$f1b$`Resid. Df`, aov.mods$f1b$`Pr(>Chi)`,
                       r2.f1b),
                     digits = 4)
  
  out.indiv <-  c(class(f.best)[1], var,
                  as.character(formula(f.best))[3],
                  aov.mods$f.best$`Resid. Df`,
                  aov.mods$f.best$`Pr(>Chi)`,
                  out.indiv)
  
  
  # other vector for predictors in best model
  out.best.pred <- cbind(var, extract.preds(f.best))
  
  # Plot graphs: (# wrong dispersion parameter is being used = overly optimistic)
  if (plot.graphs) {
    par (mfrow = c(1,2),
         mar = c(4,3,2,2),
         oma =c(1,2,0,0))
    
    plot( (x/y) ~ Seal_500,
         dataset,
         pch = 20,
         ylab = "",
         xlab = "% Impervious surfaces")
    add.stats(mods.list$f1b, type = "glm")
    
    plot( (x/y) ~ prop.neo,
          data= dataset,
          pch = 20,
          ylab = "",
          xlab = "Proportion of Neophytes")
    
    add.stats(mods.list$f1a, type = "glm")

    mtext(2, text = var,
          outer = TRUE, cex = 0.7)
  }
  return(list( result.table = out.indiv,
               predictor.table = out.best.pred,
               models = mods.list))
}

# Negative binomial GLMs
fit.negbin.glms <- function(dataset, var,
                             plot.graphs = FALSE) {
  
  dataset$y <- dataset[, var]
  require(MASS)
  # Negative binomial : 
  mods.list <- list(
    f0 = glm.nb(y ~ 1, data = dataset),
    f1a = glm.nb(y ~ prop.neo, data = dataset),
    f1b = glm.nb(y ~ Seal_500 , data = dataset),
    f2a = glm.nb(y ~ Seal_500  + prop.neo  , data = dataset),
    f4a = glm.nb(y ~ Seal_500 * prop.neo , data = dataset)
  )
  
  # Extract best models: 
  sel = model.sel(mods.list)
  mods <- rownames(subset(sel, delta <2))
  
  f.best <- mods.list[[mods[1]]]
  if ("f0" %in% mods) {
    if (mods[1] == "f0") {
      f.best <-  mods.list$f0
    } else {
      # if the best model is still significant
      if (anova(mods.list$f0,f.best)[2,6] > 0.05) {
        f.best <-  mods.list$f0
      }
    }
  }
  
  
  # Add f.best to the list: 
  mods.list$f.best <- f.best
  
  # AVONA tables with Likelihood ratio test: 
  aov.mods <- sapply(names(mods.list), function(x) {
    anova(mods.list$f0, mods.list[[x]])[2,]
  }, simplify = FALSE)
  
  # Summary of coefficients:
  coef.mods <- sapply(names(mods.list), function(x) {
    summary(mods.list[[x]])$coefficient
  }, simplify = FALSE)
  
  # Store statistics
  out.indiv <- round(c(r2_nagelkerke(f.best)[[1]],
                       coef.mods$f1a[2,1:2],
                       aov.mods$f1a$`Resid. df`, aov.mods$f1a$P,
                       r2_nagelkerke(mods.list$f1a)[[1]],
                       coef.mods$f1b[2,1:2],
                       aov.mods$f1b$`Resid. df`, aov.mods$f1b$P,
                       r2_nagelkerke(mods.list$f1b)[[1]]),
                     digits = 4)
  
  out.indiv <-  c(class(f.best)[1], var,
                  as.character(formula(f.best))[3],
                  aov.mods$f.best$`Resid. df`,
                  aov.mods$f.best$P,
                  out.indiv)
  
  
  # other vector for predictors in best model
  out.best.pred <- cbind(var, extract.preds(f.best))
  
  # Plot graphs:
  if (plot.graphs) {
    par (mfrow = c(1,2),
         mar = c(4,3,2,2),
         oma =c(1,2,0,0))
    
    plot(y ~ Seal_500,
         dataset,
         pch = 20,
         ylab = "",
         xlab = "% Impervious surfaces")
    add.stats(mods.list$f1b, type = "negbin")
    
    plot( y ~ prop.neo,
          data= dataset,
          pch = 20,
          ylab = "",
          xlab = "Proportion of Neophytes")
    
    add.stats(mods.list$f1a, type = "negbin")
   
    mtext(2, text = var,
          outer = TRUE, cex = 0.7)
  }
  return(list( result.table = out.indiv,
               predictor.table = out.best.pred,
               models = mods.list))
}


# Normal linear models
fit.lms <- function(dataset, var ,
                    plot.graphs = FALSE) {
  
  dataset$y <- dataset[, var]
  
  # LM : 
  mods.list <- list(
    f0 = lm(y ~ 1, data = dataset),
    f1a = lm(y ~ prop.neo, data =  dataset),
    f1b = lm(y ~ Seal_500 , data =  dataset),
    f2a = lm(y ~  Seal_500 + prop.neo , data =  dataset),
    f4a = lm(y ~  prop.neo * Seal_500  , data =  dataset)
  )
  
  # Extract Best model
  sel = model.sel(mods.list)
  mods <- rownames(subset(sel, delta <2))
  
  f.best <- mods.list[[mods[1]]]
  if ("f0" %in% mods) {
    if (mods[1] == "f0") {
      f.best <-  mods.list$f0
    } else {
      # if the best model is still significant
      if (anova(mods.list$f0,f.best)[2,6] > 0.05) {
        f.best <-  mods.list$f0
      }
    }
  }
  
  
  # Add f.best to the list: 
  mods.list$f.best <- f.best
  
  # AVONA tables with Likelihood ratio test: 
  aov.mods <- sapply(names(mods.list), function(x) {
    anova(mods.list$f0, mods.list[[x]])
  }, simplify = FALSE)
  
  # Summary of coefficients:
  coef.mods <- sapply(names(mods.list), function(x) {
    summary(mods.list[[x]])$coefficient
  }, simplify = FALSE) 
  # Store statistics
  out.indiv <- round(c(r2(f.best)[[1]],
                       coef.mods$f1a[2,1:2],
                       aov.mods$f1a$Res.Df[2], aov.mods$f1a$P[2],
                       r2(mods.list$f1a)[[1]],
                       coef.mods$f1b[2,1:2],
                       aov.mods$f1b$Res.Df[2], aov.mods$f1b$P[2],
                       r2(mods.list$f1b)[[1]]),
                     digits = 4)
  
  out.indiv <-  c(class(f.best)[1], var,
                  as.character(formula(f.best))[3],
                  aov.mods$f.best$Res.Df[2], aov.mods$f.best$P[2],
                  out.indiv)

  # other vector for predictors in best model
  out.best.pred <- cbind(var, extract.preds(f.best))
  
  # Plots graphs 
  if (plot.graphs) {
    par (mfrow = c(1,2),
         mar = c(4,3,2,2),
         oma =c(1,2,0,0))
    
    plot(y ~ Seal_500,
         dataset,
         pch = 20,
         ylab = "",
         xlab = "% Impervious surfaces")
    add.stats(mods.list$f1b, type = "lm")
    
    plot( y ~ prop.neo,
          data= dataset,
          pch = 20,
          ylab = "",
          xlab = "Proportion of Neophytes")
    
    add.stats(mods.list$f1a, type = "lm")
    
    mtext(2, text = var,
          outer = TRUE, cex = 0.7)
  }
  
  return(list( result.table = out.indiv,
               predictor.table = out.best.pred,
               models = mods.list))
}
