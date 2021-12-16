# Identifying high PAV plots

# Extract high PAV plots
hi_pav_plots <- unique(c(
  rownames(allergen_summary)[allergen_summary$CWM.pav >10], # weighted
  rownames(allergen_summary)[allergen_summary$mean.pav >6])) # unweighted

hi_pav_plots <-hi_pav_plots[order(allergen_summary[hi_pav_plots,"CWM.pav"])]


# Extract communities
hipav_comm <- vegcomm[hi_pav_plots,]
hipav_comm <- hipav_comm[,which(colnames(hipav_comm) %in% allergenics &
                                      colSums(hipav_comm>1)>0)
]
hipav_comm <- hipav_comm[order(species_allergen[colnames(hipav_comm),"PAV"])]


par(mar = c(2,8,1,1))
image( 1:nrow(hipav_comm),1:ncol(hipav_comm), 
       z = as.matrix(hipav_comm), 
      xaxt = "n", yaxt = "n", ann = FALSE)
axis(2, at = 1:ncol(hipav_comm), labels = colnames(hipav_comm), 
     cex.axis = 0.5, las = 1)
axis(1, at = 1:nrow(hipav_comm), labels = rownames(hipav_comm), 
     cex.axis = 0.6, las = 1)


# Extract high allergen cover plots
hi_cover_plots <-  rownames(allergen_summary)[allergen_summary$cover.all >0]
hi_cover_plots <-hi_cover_plots[order(allergen_summary[hi_cover_plots,"cover.all"])]

# Extract communities
hicover_comm <- vegcomm[hi_cover_plots,]
hicover_comm <- hicover_comm[,which(colnames(hicover_comm) %in% allergenics &
                                      colSums(hicover_comm>1)>0)
                             ]
hicover_comm <- hicover_comm[order(species_allergen[colnames(hicover_comm),"PAV"])]


par(mar = c(2,8,1,1))
image( 1:nrow(hicover_comm),1:ncol(hicover_comm), 
       z = as.matrix(hicover_comm), 
       xaxt = "n", yaxt = "n", ann = FALSE)
axis(2, at = 1:ncol(hicover_comm), 
     labels = colnames(hicover_comm), 
     cex.axis = 0.5, las = 1)
axis(1, at = 1:nrow(hicover_comm), labels = rownames(hicover_comm), 
     cex.axis = 0.6, las = 1)


# Extract high allergen SRplots
hi_SR_plots <-  rownames(allergen_summary)[allergen_summary$all.num > 0]
hi_SR_plots <-hi_SR_plots[order(allergen_summary[hi_SR_plots,"all.num"])]

# Extract communities
hiSR_comm <- vegcomm[hi_SR_plots,]
hiSR_comm <- hiSR_comm[,which(colnames(hiSR_comm) %in% allergenics &
                                      colSums(hiSR_comm>0)>0)]
hiSR_comm <- hiSR_comm[order(species_allergen[colnames(hiSR_comm),"PAV"])]


par(mar = c(2,8,1,1))
image( 1:nrow(hiSR_comm),1:ncol(hiSR_comm), 
       z = log(as.matrix(hiSR_comm)), 
       xaxt = "n", yaxt = "n", ann = FALSE)
axis(2, at = 1:ncol(hiSR_comm), 
     labels = colnames(hiSR_comm), 
     cex.axis = 0.5, las = 1)
axis(1, at = 1:nrow(hiSR_comm), labels = rownames(hiSR_comm), 
     cex.axis = 0.6, las = 1)

