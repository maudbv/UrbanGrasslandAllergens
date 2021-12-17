# Community Allergernicity along the gradient
# Fitting generalized linear models to allergenic species richness, proportion and cover
# Extract both the "best model" with interaction (prop.neo*seal_500)
# and three independant models for each of the two predictors: Seal_500 & prop.neo

glms.SRall <- (function(exclude.absences = TRUE,
                        show.plots = FALSE){

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
                          "Seal.coef", "Seal.se","Seal.df","Seal.P","Seal.R2"
                          )

# Create a list to store the models: 
all.models <- list()

# Create table for best models predictors: 
glms.pred<- data.frame(matrix(NA, nrow = 0,ncol =12 ) )
colnames(glms.pred) <- c("group", "type","n.obs",
                         "class", "var","pred","est","se","z","P","R2beta","Pchi")

# Trends in RICHNESS - Poisson GLMs ####
for (i in 1:5) {
tmp <- allergen_summary
y <- c("all.num", "all.num.nat","all.num.arc","all.num.neo","all.num.exo")[i]
g <- c("all", "nat","arc","neo","exo")[i]
type <- "richness"

if (exclude.absences) {
   z <- c("SR","SR.nat","SR.arch","SR.neo","SR.exo")[i]
   tmp <- tmp[which(tmp[,z]>0 & !is.na(tmp[,y])),]
}
n = length(tmp[,y])
fit <- fit.poisson.glms(dataset = tmp, var = y,
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

# trends in allergenic species COVER - LM ####
for (i in 1:5) {
   tmp <- allergen_summary
   y <- c("cover.all", "cover.nat.all","cover.arc.all",
          "cover.neo.all","cover.exo.all")[i]
   g <- c("all", "nat","arc","neo","exo")[i]
   type <- "cover"
   
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

# trends in allergenic species proportions - Binomial GLM  ####
for (i in 1:4) {
   tmp <- allergen_summary
   x <- c("all.num.nat","all.num.arc","all.num.neo","all.num.exo")[i]
   y <-  "all.num"
   nam <- paste("prop",x, sep="_")
   g <- c("nat","arc","neo","exo")[i]
   type <- "proportion"
   

   if (exclude.absences) {
      z <- c("SR.nat","SR.arch","SR.neo","SR.exo")[i]
      tmp <- tmp[which(tmp[,z]>0 & !is.na(tmp[,x])),]
   }
   
   n = length(tmp[,y])
   
   # Fit Binomials for proportions
   fit <- fit.binom.glms(dataset = tmp, var = x,tot = y,
                         plot.graphs = show.plots)
   glms.table[i + 10,] <- c(g,type, n, fit$result.table)
   rownames(glms.table)[i + 10] <- nam
   
   
   glms.pred <- rbind(glms.pred, 
                      cbind(group = g, type = type, n = n, 
                            fit$predictor.table))
   
   all.models [[i + 10]] <- fit$models
   names(all.models )[i + 10]<- nam
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
write.csv(glms.SRall$glms.table,
          "results/allergenic diversity models table.csv")
write.csv(glms.SRall$glms.best.pred,
          "results/allergenic diversity best predictor table.csv")

