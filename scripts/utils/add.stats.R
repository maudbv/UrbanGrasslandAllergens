# function to quickly add stats to a plot : a simple linear regression

add.stats <- function(f=NULL, formula = NULL, data= NULL, 
                      type = "", 
                      col.stats = "black",
                      adj.stats = 0.90,
                      line.stats = 0,
                      l.col= "darkgrey",
                      show.marginal = TRUE,
                      plot.abline = "yes",
                      se.show = TRUE,
                      add.text2plot = TRUE,
                      test.disp = FALSE,
                      output = TRUE) {
  
  require(performance)
  
  # lm ####
  if ((type == "lm") | (class(f)[1] == "lm")) {
    
    if (is.null(f)) {
      f <- lm(as.formula(formula), data)
    }
    
    
    est <-  summary(f)$r.squared
    fstat <- summary(f)$fstatistic
    df <- f$df.residual
    
    p <- pf(q = fstat[1],df1 = fstat[2],df2 =fstat[3],lower.tail = FALSE)
    
    stts <- substitute("R"^2*"="*a*b,
                       list(
                         a = round(est,2),
                         b =p2star(p, marginal = show.marginal)
                       )
    )
    
    
    # plot line
    if (plot.abline != "no") {
      l.typ = NULL
      if ((show.marginal == TRUE & p < 0.1)| plot.abline == "always") {
        l.typ = "dashed"
      }
      if (p < 0.05)  l.typ = "solid"
      
      if (! is.null(l.typ) ) {
        x <- as.character(formula(f))[3]
        nwdat <-  f$model
        nwdat <-  nwdat[order(nwdat[,x]),]
        pr <- predict.lm(f,
                         newdata = nwdat,
                         se.fit = TRUE)
        nwdat <- cbind(nwdat, pr)
        
        lines(nwdat[,x], nwdat$fit,
              col = l.col, lty= l.typ)
        if (se.show == TRUE) {
          lines(nwdat[,x], nwdat$fit + 1.96 * nwdat$se.fit,
                col = l.col, lty= "dotted")
          lines(nwdat[,x], nwdat$fit - 1.96 * nwdat$se.fit,
                col = l.col,lty= "dotted")
        }
      }
    }
  }
  
  # lme
  if (type == "lme" | class(f)[1]  == "lme") {
    est = MuMIn::r.squaredGLMM(f)[1]
    p = anova(f)[2,4]
    df = NA
    # write stats
    stts <- substitute("R"^2*".marginal ="*a*b,
                       list(
                         a = round(est,2),
                         b =p2star(p, marginal = show.marginal)
                       )
    )
    # plot line
    if (plot.abline != "no" ) {
      l.typ = NULL
      if (p < 0.05) {
        l.typ = "solid"
      }
      if (show.marginal == TRUE & p < 0.1 | plot.abline == "always"){
        l.typ = "dashed"
      }
      abline(lm(as.formula(as.character(f$call)[2]), data = f$model),
             col = col.stats, lty = l.typ)
    }
  }
  
  # lmer ####
  if (type == "lmer" | class(f)[1]  == "lmerModLmerTest") {
    require(lmerTest)
    require(performance)
    
    f <- lmerTest::as_lmerModLmerTest(f)
    
    est = MuMIn::r.squaredGLMM(f)[1]
    p = anova(f)[1,6]
    df = summary(f)$coef[2,"df"]
    
    # write stats
    stts <- substitute("R"^2*".marginal ="*a*b,
                       list(
                         a = round(est,2),
                         b =p2star(p, marginal = show.marginal)
                       )
    )
    
    # plot line
    if (plot.abline  != "no") {
      
      #extract formula for fixed effects
      form <- as.character(f@call$formula)
      form.pred <- strsplit(form[3],split = "\\+")[[1]]
      form.pred <- form.pred[-grep("\\|", form.pred)]
      form <- paste(form[2], form[1], paste(form.pred , collapse = "+"))
      
      if (p < 0.05) abline(lm(as.formula(form), data = f@frame),
                           col = col.stats)
      if ((show.marginal == TRUE & p < 0.1) | plot.abline == "always") {
        abline(lm(as.formula(form), data = f@frame),
               col = col.stats,
               lty = "dashed"
        )
      }
    }
  }
  
  # negbin ####
  if (type == "negbin" | class(f)[1]  == "negbin") {
    require(performance)
    fit <-summary(f)
    p <- fit$coefficients[2, 4]
    est <- as.numeric(r2_nagelkerke(f)[1])
    df = f$df.residual
    
    # write stats
    stts <-  substitute("R" ^2 * "=" * rtwo * ps,
                        list( rtwo=round(est, 2),
                              ps = p2star(p,marginal = show.marginal)))
    
    # plot line
    if (plot.abline != "no" ) {
      
      l.typ = NULL
      if ((show.marginal == TRUE & p < 0.1) | plot.abline == "always") l.typ = "dashed"
      if (p < 0.05)  l.typ = "solid"
      
      if ( ! is.null(l.typ) ) {
        
        x <- as.character(formula(f))[3]
        nwdat <-  f$model
        nwdat <-  nwdat[order(nwdat[,x]),]
        pr <- predict.glm(f,
                          newdata = nwdat,
                          se.fit = TRUE)
        nwdat <- cbind(nwdat, pr)
        
        lines(nwdat[,x], exp(nwdat$fit),
              col = l.col, lty= l.typ)
        if (se.show == TRUE) {
          lines(nwdat[,x], exp(nwdat$fit + 1.96 * nwdat$se.fit),
                col = l.col, lty= "dotted")
          lines(nwdat[,x], exp(nwdat$fit - 1.96 * nwdat$se.fit),
                col = l.col,lty= "dotted")
        }
      }
    }
  }
  
  # glm ####
  if (type == "glm" | class(f)[1]  == "glm") {
    
    # fit <-summary(f)
    # p <- fit$coefficients[2, 4]
    # dev <- 1 - fit$deviance/fit$null.deviance
    #   
    # 
    # # write stats
    #  stts <- paste("%dev: ", round(dev,2), p2star(p, marginal = show.marginal), sep = "")
    # mtext(3, text = stts, cex = 0.7,col = col.stats, adj = adj.stats)
    # 
    # # plot line
    # if (plot.abline ) {
    # if (p < 0.05) abline(lm(f$formula, data = f$data), col = col.stats)
    # if (show.marginal == TRUE & p < 0.1) abline(lm(f$formula, data = f$data), col = col.stats, lty = "dashed")
    # }
    
    require(performance)
    fit <-summary(f)
    # p <- fit$coefficients[2, 4]
    p <- anova(f, test = "LRT")[2, 5]
    est <- as.numeric(r2_nagelkerke(f)[1])
    df <- f$df.residual
    
    # write stats
    stts <-  substitute("R" ^2 * "=" * rtwo * ps,
                        list( rtwo=round(est, 2),
                              ps = p2star(p,marginal = show.marginal)))
    
    # plot line
    if (plot.abline  != "no") {
      
      l.typ = NULL
      if ((show.marginal == TRUE & p < 0.1) | plot.abline == "always") l.typ = "dashed"
      if (p < 0.05)  l.typ = "solid"
      
      if ( ! is.null(l.typ) ) {
        
        x <- as.character(f$formula)[3]
        nwdat <-  f$data
        nwdat <-  nwdat[order(nwdat[,x]),]
        pr <- predict.glm(f,
                          newdata = nwdat,
                          se.fit = TRUE,
                          type = "response")
        nwdat <- cbind(nwdat, pr)
        lines(nwdat[,x], nwdat$fit)
        
        if (se.show == TRUE) {
          lines(nwdat[,x], nwdat$fit + 1.96 *  nwdat$se.fit,
                col = l.col, lty= "dotted")
          lines(nwdat[,x], nwdat$fit - 1.96 * nwdat$se.fit,
                col = l.col,lty= "dotted")
        }
      }
      
      
      
      # Disperson for poisson glm
      if (f$family$family == "poisson" & test.disp) {
        require (AER)
        disp.test <- dispersiontest(f)
        disp <- dispersiontest(f)$estimate
        p.disp <- dispersiontest(f)$p.value
        
        dsp.stts <- paste("Disp = ",round(disp,1),
                          p2star(p.disp, marginal = show.marginal), sep = "")
        if (add.text2plot){
          mtext(4, text = dsp.stts,
                cex = 0.7, col = col.stats, 
                adj = 0, line = line.stats)
        }
      }
      
      
    }
  }
  
  # spearman correlation ####
  if (type == "cor" | class(f)[1]  == "htest") {
    
    est = f$estimate
    p = f$p.value
    df = NA
    # write stats
    stts <- substitute("rho"*"="*a*b,
                       list(
                         a = round(est,2),
                         b =p2star(p, marginal = show.marginal)
                       )
    )
    
  }
  
  # Add stats to an existing plot ####
  if (add.text2plot) {
    mtext(3, text = stts, cex = 0.7,
          col = col.stats, adj = adj.stats,
          line = line.stats)
  }
  
  # return the formatted stats: ####
  
  if (output) return( data.frame(R2 = est, P = p, df.resid = df))
  
}


