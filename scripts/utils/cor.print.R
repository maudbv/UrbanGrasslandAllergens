## Printing and saving formatted correlation tables

cor.print <- function(data = NULL, cor.mat = NULL,
                      method = "spearman",
                      return.print = TRUE,
                      lower.triangle = TRUE,
                      plotting = FALSE) {

  require(Hmisc)
  
  stopifnot (is.null(data) | is.null(cor.mat)) 
  
if (!is.null(data)) {
  stopifnot (all(apply(data,2,is.numeric))) 
  cor.mat <-rcorr(as.matrix(data), type = method)
}
  
if ( (!is.null(cor.mat)) & (class(cor.mat) != "rcorr") ) {
  print("Error: cor.mat must be of class 'rcorr'")
  stop()
}

  if (plotting) {
    par(mfrow= c(1,1))
    corrplot(corr = cor.mat$r, method = "ellipse", 
             p.mat = cor.mat$P, insig = "label_sig",
             sig.level = c(0.001, 0.01,0.05 ),
             pch.cex = 0.7, tl.pos = "lt",type = "lower")
    corrplot(corr = cor.mat$r, method =  "number", 
             cex = 0.5, type = "upper",diag = FALSE, add=TRUE, tl.pos = "n" , cl.pos = "n")
  }
  
  # If we want a nice printable table with stars for significance levels
  if (return.print) {
    
  # function to transform pvalues into stars:
  p2star <- function(x, marginal=F, quantiles=NULL) {
    # x : a p value 
    s<-NA
    x <- as.numeric(x)
    if(!is.na(x)){
      if (!is.null(quantiles)) q<-quantiles
      else q<-c(0.1,0.05,0.01,0.001)
      if (x>q[2]) s<-'ns'
      if(marginal) if (x<q[1]) s <-'a'
      if (x<=q[2]) s<-'*'
      if (x<=q[3]) s<-'**'
      if (x<=q[4]) s<-'***'
    }
    return(s)
  } 
  
  # Create a pretty csv table for printing:
  cor.print <- sapply(colnames(cor.mat$r), function(i) {
    y <- paste(round(cor.mat$r[,i],2), 
               sapply(cor.mat$P[,i],p2star)
    )
  })
  
  rownames(cor.print) <- colnames(cor.print)
  diag(cor.print) <- ""
  if (lower.triangle) cor.print[ upper.tri(cor.print,diag = TRUE) ]<- ""
  
} else {
  cor.print <- cor.mat
}
  return(cor.print)
}
