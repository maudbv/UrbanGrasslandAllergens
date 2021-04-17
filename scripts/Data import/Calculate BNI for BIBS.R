## Calculate BNI on Berlin BIBS data
## only for herbaceous species

BNI.berlin <- (function(){

# Load function: 
source('data/Schittko et al. 2020/BNI_function.R')

# Load data ####

### load abundance table; sites as rows, species as cloumns
abund<-read.csv("data/Schittko et al. 2020/abundance.csv", sep=";", row.names=1, check.names = FALSE)
abund<-as.matrix(abund)

# Remove trees
tree.species <- sub(" " , "_", 
                    c('Acer', 'Acer campestre', 'Acer negundo', 'Acer platanoides', 'Acer pseudoplatanus','Crataegus monogyna' ,'Prunus', 'Prunus serotina','Prunus spinosa', 'Prunus domestica', 'Pinus sylvestris', 'Quercus', 'Tilia', 'Populus', 'Populus tremula', 'Quercus robur', 'Robinia pseudoacacia',  "Pyrus_communis_agg."))

abund <- abund[,which(! colnames(abund) %in% tree.species)]

### calculate relative abundance
abund <- abund / rowSums(abund)
head(abund)

### load status table; species as rows, status as column
status.data<-read.csv("data/Schittko et al. 2020/factor.csv", sep=";", row.names=1)
# select species:
status.data <- status.data[which(rownames(status.data)%in%colnames(abund) ),]
status.data <- status.data[match(colnames(abund), rownames(status.data)),]

####load trait table; species as rows, traits as columns
trait.mat<-read.csv("data/Schittko et al. 2020/trait4.csv", sep=";", row.names=1)
#select species:
trait.mat <- trait.mat[which(rownames(trait.mat)%in%colnames(abund)),]   
trait.mat <-trait.mat[match(colnames(abund), rownames(trait.mat)),]

# Make sure the order of species is the same in all tables:
stopifnot(all(rownames(trait.mat) == rownames(status.data)))
stopifnot(all(rownames(trait.mat) == colnames(abund)))


## Calculate BNI on the BIBS data  ####
BNI <- BNI.calc(com = abund,
                       trait = trait.mat,
                       YSI = status.data$years_since_introduction,
                       dist.method = "gower")

neos <- rownames(status.data)[status.data$status_category == "N"]
BNI$index$pneo <- rowSums(abund[, colnames(abund) %in% neos]>0)/rowSums(abund>0)

return(BNI$index)
})()

