#############################################################################
#Read in csv file with shape coordinates
tmp <- read.csv("caran_coords.csv",header=T, row.names =1, stringsAsFactors = FALSE) 
# here we are reading in a tab-delimited text file from MorphoJ, but this can be an issue with data from any outside program. The stringsAsFactors = FALSE is VERY important here

shape <- as.matrix(sapply(tmp[,-(1)], as.numeric)) 
# here we say, use all columns except the first three.

is.numeric(shape)
[1] TRUE # now it's numeric. Ready to go!
#Gonna Need names and Classifiers
names <- tmp[,1]
coords <- arrayspecs(shape[,1:ncol(shape)], 136,3)
#arrayspecs- substitute the column number where the coordinates begin for the 2, p=number of landmarks, k=number of dimensions
dimnames(coords)[[3]] <- names
################################################
###remove landmarks here if neccessary#####
omit <- c(29:30)
skull <- coords[-omit,,]
###########################################################
###Read in Sliding semi-landmark file

slide<-read.csv("flat_slide.csv",header=T)
######################################################
##Perfrom Procrustes Superimposition and incorporate sliding sem-landmarks
Y.gpa <- gpagen(skull,curves=slide,ProcD=FALSE)   

## Run a principal components analysis
# Use a PCA to 
car.skull.pca <- gm.prcomp(Y.gpa$coords)
########################################
##Seperate PC extremes
PC1.min <- car.skull.pca$shapes$shapes.comp1$min
PC1.max <- car.skull.pca$shapes$shapes.comp1$max
PC2.min <- car.skull.pca$shapes$shapes.comp2$min
PC2.max <- car.skull.pca$shapes$shapes.comp2$max
###################################################
##Warp PC extremes to view axes of shape variation
# PC 1 extremes
plotRefToTarget(PC1.min,PC1.max,method="vector")
plotRefToTarget(PC2.min,PC2.max,method="vector")
##########################################
##Warp the mean shape to PC extremes
M <- mshape(Y.gpa$coords)
plotRefToTarget(M,PC1.max,method="vector")
plotRefToTarget(M,PC1.min,method="vector")