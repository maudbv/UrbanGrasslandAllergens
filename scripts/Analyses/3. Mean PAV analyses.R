# Community mean allergergenic severity pav along the gradient

pav.analyses <- (function(plot.results = FALSE) {
## Plotting function ####
plot.trends.pav <-function(df = tmp,
                             type.model = "lm",
                             y.title = "CWM pav",
                             plot.title = "All species",
                           interaction = TRUE,
                           plotting = plot.results){
   
   df <- na.omit(df[, c("y", "Seal_500","prop.neo")])
   
   # function to fit null model
   f.null <- function(dat = df) {
     
      # Model fitting = spearman correlation
      if (type.model== "spearman.cor") {
         f <- 0
      }
      
      # Model fitting = spearman correlation
      if (type.model== "lm") {
         f <- lm( y ~ 1, dat)
      }
      
      if (type.model== "glm") {
         f <- glm( y ~ 1, dat,
                   family = gaussian(log),
                   start = c(1,0))
      }
      return(f)
   }
    f0 <- f.null()
   
   # Function to fit one model:
   fn <- function(dat = df,
                  x.name ="Seal_500" , 
                  x.label= "% impervious surfaces",
                  plot = plotting) {
      dat$x <- dat[,x.name]
      # Model fitting = spearman correlation
      if (type.model== "spearman.cor") {
         f <- cor.test(dat$y, dat$x, method = "spearman")
      }
      
      # Model fitting = spearman correlation
      if (type.model== "lm") {
         f <- lm( y ~ x, dat)
      }
      
      if (type.model== "glm") {
         f <- glm( y ~ x, dat,
                   family = gaussian(log),
                   start = c(1,0))
      }
      
   if (plot) {
      plot( y ~ x, dat,
         pch = 20, ylab = y.title,
         xlab = x.label)
      add.stats(f, type = type.model)
   }
      
   return(f)
   }
   
   # Format plotting device:
   if (plotting) par (mfrow = c(1,3),mar = c(4,4, 2, 1), oma = c(0,0, 2,0))
   
   # Apply to each variable:
   f1 <- fn(df, "Seal_500","% impervious surfaces")
   f2 <-  fn(df, "prop.neo", "Prop. Neophytes")
   if (plotting) mtext(3, text = plot.title, outer = TRUE)
   
   # Prepare results for testing interaction model
   f.int.a <- NULL

   # If testing interaction model
   if (interaction) {
      # temporary plotting function: 
   f.int <- function(dat = df,
                     x1.name ="Seal_500", 
                     x2.name = "prop.neo") {
         dat$x1 <- dat[,x1.name]
         dat$x2 <- dat[,x2.name]
         # Model fitting = spearman correlation
         if (type.model== "spearman.cor") {
            f <- NA
         }
         
         # Model fitting = lm
         if (type.model== "lm") {
            f <- lm( y ~ x1 * x2, dat)
         }
         
         if (type.model== "glm") {
            f <- glm( y ~ x1 * x2, dat,
                      family = gaussian(log),
                      start = c(1,0,0,0))
         }
         return(f)
   }
   
   f.int.a <- f.int(df, "Seal_500","prop.neo")
 
   }
   
   # Formatted output:
   mean.pav <- c(f1$df.residual,
                 summary(f1)$coefficients[2, 1:4],
                 r2(f1)[[1]][1],
                 summary(f2)$coefficients[2, 1:4],
                 r2(f2)[[1]][1])
   
   mean.pav.int.a <- c(f.int.a$df.residual,
                                  r2(f.int.a)[[1]][1],
                                  anova(f0, f.int.a)[2,6],
                                  summary(f.int.a)$coefficients[2, 1:4],
                                  summary(f.int.a)$coefficients[3, 1:4],
                                  summary(f.int.a)$coefficients[4, 1:4]
   )
   
 
   

   return(output <- list(
      models = list(f1 = f1,
                    f2 = f2,
                    f0 = f0,
                    f.int.a = f.int.a),
      results = list(mean.pav = mean.pav,
                     mean.pav.int.a = mean.pav.int.a)
   )
   )
}



#create result table: ####
pav.models <- data.frame(matrix(NA, nrow = 0,ncol =11))
colnames(pav.models) <- c("df.resid",
                          "Seal.coef","Seal.sd", "Seal.t","Seal.P","Seal.r2",
                          "p.neo.coef","p.neo.sd", "p.neo.t","p.neo.P","p.neo.r2")


pav.models.inter.pneo <- data.frame(matrix(NA, nrow = 0,ncol = 15))
colnames(pav.models.inter.pneo ) <- c("df.resid", "r2","P",
                          "Seal.coef","Seal.sd", "Seal.t","Seal.P",
                          "p.neo.coef","p.neo.sd", "p.neo.t","p.neo.P",
                          "p.int.coef","p.int.sd", "p.int.t","p.int.P"
                          )



# Run the analyses for each types of mean PAV: ####
for (i in 1:10) {
var <- c("mean.pav","mean.nat.pav","mean.arc.pav",
         "mean.neo.pav","mean.exo.pav",
         "CWM.pav","CWM.nat.pav","CWM.arc.pav",
         "CWM.neo.pav","CWM.exo.pav")[i]
         
tmp <- allergen_summary
tmp$y <- tmp[,var]
tmp <- tmp[!is.na(tmp$y),]
 m <- plot.trends.pav(df = tmp,
                  type.model = "lm",
                  y.title = "Mean pav",
                  plot.title = "All species")
 
 pav.models[var,] <- m$results$mean.pav
 pav.models.inter.pneo[var,] <- m$results$mean.pav.int.a
}

# Add mean and SD values ####
 tmp <- cbind(mean = round(apply(allergen_summary[, c("mean.pav","CWM.pav",
                                              "mean.nat.pav", "CWM.nat.pav",
                                              "mean.arc.pav","CWM.arc.pav",
                                              "mean.neo.pav", "CWM.neo.pav",
                                              "mean.exo.pav", "CWM.exo.pav")],
                         2, mean, na.rm = TRUE), 2),
      
      sd = round(apply(allergen_summary[, c("mean.pav","CWM.pav",
                                            "mean.nat.pav", "CWM.nat.pav",
                                            "mean.arc.pav","CWM.arc.pav",
                                            "mean.neo.pav", "CWM.neo.pav",
                                            "mean.exo.pav", "CWM.exo.pav")],
                       2, sd, na.rm = TRUE), 2)
)


# Reformat and add means and sd: 
pav.models <- as.data.frame(pav.models)
pav.models <- cbind(tmp, pav.models[rownames(tmp),])
pav.models.inter.pneo <- as.data.frame(pav.models.inter.pneo)[rownames(tmp),]

# return results ####
return(pav.analyses = list(pav.models = pav.models,
                            pav.models.inter.pneo = pav.models.inter.pneo))
})()

# Write table : ####
write.csv(pav.analyses$pav.models , "results/PAV.models.csv")
write.csv(pav.analyses$pav.models.inter.pneo , "results/PAV.models.inter.pneo.csv")



