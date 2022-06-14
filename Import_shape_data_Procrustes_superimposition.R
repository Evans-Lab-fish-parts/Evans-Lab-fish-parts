require(geomorph)


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
summary(Y.gpa)
plot(Y.gpa)