#######################################
###   Example of mvMORPH analysis   ###
###       by Olivier Larouche       ###
#######################################

#This is an example of an analysis used for the Systematic Biology Submission. It performs model-fitting and extracts model parameters for two very simple models of trait evolution: a single rate BM model and a single peak OU model.

#Before performing these analyses, you should consult the vignettes which explain what the various parameters do in mvMORPH. These vignettes can be accessed via the command below:
browseVignettes("mvMORPH")

library(mvMORPH)
library(geomorph)

#Loading the data. In my case, this is a file that contains geomorph dataframes including the superimposed coordinates, the phylogeny pruned to the morphological data, and categorical traits (e.g., diet, taxonomy)
load("final_all3LMdatasetsandtree_procD.Rdata")#minimize chord distance 


#Extracting the coordinates for each of the partitions of shape
premax.lm<-skull_GPA_coords[1:27,,]
ang.lm<-skull_GPA_coords[28:39,,]
dent.lm<-skull_GPA_coords[40:47,,]
nas.lm<-skull_GPA_coords[48:57,,]
neuro.lm<-skull_GPA_coords[58:109,,]
pjaws.lm<-skull_GPA_coords[110:125,,]
max.lm<-skull_GPA_coords[126:132,,]
hyo.lm<-skull_GPA_coords[133:149,,]
uro.lm<-skull_GPA_coords[150:163,,]
cer.lm<-skull_GPA_coords[164:182,,]
oper.lm<-skull_GPA_coords[183:187,,]
palat.lm<-skull_GPA_coords[188:191,,]
uppjaws.lm<-skull_GPA_coords[192:197,,]


#Running mvMORPH on individual shape partitions. For each partition, first we need to perform PCA analysis and extract the PC scores (numbers in brackets) explaining up to a target amount of variation (in this case I chose 70% of variance explained). The number of PCs required will vary according to the shape complexity of the various structures. The model-fitting in mvmorph is computationally extensive and will take several days to complete. This is why we use a few PC's instead of the raw landmark data: computational time increases exponentially with the number of traits. 
PCA_GPA<-gm.prcomp(skull_GPA_coords)
model_fit.GPA.BM1<- mvBM(pruned_tree, PCA_GPA$x[,1:8], model="BM1")
model_fit.GPA.OU1<- mvOU(pruned_tree, PCA_GPA$x[,1:8], model="OU1")

PCA_premax<-gm.prcomp(premax.lm)
model_fit.GPA.BM1.premax<- mvBM(pruned_tree, PCA_premax$x[,1:4], model="BM1")
model_fit.GPA.OU1.premax<- mvOU(pruned_tree, PCA_premax$x[,1:4], model="OU1")

PCA_ang<-gm.prcomp(ang.lm)
model_fit.GPA.BM1.ang<- mvBM(pruned_tree, PCA_ang$x[,1:3], model="BM1")
model_fit.GPA.OU1.ang<- mvOU(pruned_tree, PCA_ang$x[,1:3], model="OU1")

PCA_dent<-gm.prcomp(dent.lm)
model_fit.GPA.BM1.dent<- mvBM(pruned_tree, PCA_dent$x[,1:3], model="BM1")
model_fit.GPA.OU1.dent<- mvOU(pruned_tree, PCA_dent$x[,1:3], model="OU1")

PCA_nas<-gm.prcomp(nas.lm)
model_fit.GPA.BM1.nas<- mvBM(pruned_tree, PCA_nas$x[,1:3], model="BM1")
model_fit.GPA.OU1.nas<- mvOU(pruned_tree, PCA_nas$x[,1:3], model="OU1")

PCA_neuro<-gm.prcomp(neuro.lm)
model_fit.GPA.BM1.neuro<- mvBM(pruned_tree, PCA_neuro$x[,1:7], model="BM1")
model_fit.GPA.OU1.neuro<- mvOU(pruned_tree, PCA_neuro$x[,1:7], model="OU1")

PCA_pjaws<-gm.prcomp(pjaws.lm)
model_fit.GPA.BM1.pjaws<- mvBM(pruned_tree, PCA_pjaws$x[,1:4], model="BM1")
model_fit.GPA.OU1.pjaws<- mvOU(pruned_tree, PCA_pjaws$x[,1:4], model="OU1")

PCA_max<-gm.prcomp(max.lm)
model_fit.GPA.BM1.max<- mvBM(pruned_tree, PCA_max$x[,1:3], model="BM1")
model_fit.GPA.OU1.max<- mvOU(pruned_tree, PCA_max$x[,1:3], model="OU1")

PCA_hyo<-gm.prcomp(hyo.lm)
model_fit.GPA.BM1.hyo<- mvBM(pruned_tree, PCA_hyo$x[,1:4], model="BM1")
model_fit.GPA.OU1.hyo<- mvOU(pruned_tree, PCA_hyo$x[,1:4], model="OU1")

PCA_uro<-gm.prcomp(uro.lm)
model_fit.GPA.BM1.uro<- mvBM(pruned_tree, PCA_uro$x[,1:3], model="BM1")
model_fit.GPA.OU1.uro<- mvOU(pruned_tree, PCA_uro$x[,1:3], model="OU1")

PCA_cer<-gm.prcomp(cer.lm)
model_fit.GPA.BM1.cer<- mvBM(pruned_tree, PCA_cer$x[,1:2], model="BM1")
model_fit.GPA.OU1.cer<- mvOU(pruned_tree, PCA_cer$x[,1:2], model="OU1")

PCA_oper<-gm.prcomp(oper.lm)
model_fit.GPA.BM1.oper<- mvBM(pruned_tree, PCA_oper$x[,1:3], model="BM1")
model_fit.GPA.OU1.oper<- mvOU(pruned_tree, PCA_oper$x[,1:3], model="OU1")

PCA_palat<-gm.prcomp(palat.lm)
model_fit.GPA.BM1.palat<- mvBM(pruned_tree, PCA_neuro$x[,1:2], model="BM1")
model_fit.GPA.OU1.palat<- mvOU(pruned_tree, PCA_neuro$x[,1:2], model="OU1")

PCA_uppjaws<-gm.prcomp(uppjaws.lm)
model_fit.GPA.BM1.uppjaws<- mvBM(pruned_tree, PCA_uppjaws$x[,1:2], model="BM1")
model_fit.GPA.OU1.uppjaws<- mvOU(pruned_tree, PCA_uppjaws$x[,1:2], model="OU1")

#Exporting the results for post-processing
save(model_fit.GPA.BM1, model_fit.GPA.OU1, model_fit.GPA.BM1.premax, model_fit.GPA.OU1.premax, model_fit.GPA.BM1.ang, model_fit.GPA.OU1.ang, model_fit.GPA.BM1.dent, model_fit.GPA.OU1.dent, model_fit.GPA.BM1.nas, model_fit.GPA.OU1.nas, model_fit.GPA.BM1.neuro, model_fit.GPA.OU1.neuro, model_fit.GPA.BM1.pjaws,  model_fit.GPA.OU1.pjaws, model_fit.GPA.BM1.max, model_fit.GPA.OU1.max, model_fit.GPA.BM1.hyo,  model_fit.GPA.OU1.hyo,  model_fit.GPA.BM1.uro, model_fit.GPA.OU1.uro, model_fit.GPA.BM1.cer,  model_fit.GPA.OU1.cer, model_fit.GPA.BM1.oper, model_fit.GPA.OU1.oper, model_fit.GPA.BM1.palat, model_fit.GPA.OU1.palat,  model_fit.GPA.BM1.uppjaws, model_fit.GPA.OU1.uppjaws,  file="final_mvMORPH_GPA.Rdata")


#This is the same analysis as above, except that the variance-covariance matrix is reduced to just its diagonal elements. This means that it only incorporates the variance of individual traits, and ignores their covariances in model-fitting. For the downstream analyses, you need to compare the AICc between the models including and excluding the covariances, and proceed with the models that yield lower AICc values.

PCA_GPA<-gm.prcomp(skull_GPA_coords)
model_fit.GPA.BM1.diag<- mvBM(pruned_tree, PCA_GPA$x[,1:8], model="BM1",param=list(decomp="diagonal"))
model_fit.GPA.OU1.diag<- mvOU(pruned_tree, PCA_GPA$x[,1:8], model="OU1",param=list(decomp="diagonal"))

PCA_premax<-gm.prcomp(premax.lm)
model_fit.GPA.BM1.premax.diag<- mvBM(pruned_tree, PCA_premax$x[,1:4], model="BM1",param=list(decomp="diagonal"))
model_fit.GPA.OU1.premax.diag<- mvOU(pruned_tree, PCA_premax$x[,1:4], model="OU1",param=list(decomp="diagonal"))

PCA_ang<-gm.prcomp(ang.lm)
model_fit.GPA.BM1.ang.diag<- mvBM(pruned_tree, PCA_ang$x[,1:3], model="BM1",param=list(decomp="diagonal"))
model_fit.GPA.OU1.ang.diag<- mvOU(pruned_tree, PCA_ang$x[,1:3], model="OU1",param=list(decomp="diagonal"))

PCA_dent<-gm.prcomp(dent.lm)
model_fit.GPA.BM1.dent.diag<- mvBM(pruned_tree, PCA_dent$x[,1:3], model="BM1",param=list(decomp="diagonal"))
model_fit.GPA.OU1.dent.diag<- mvOU(pruned_tree, PCA_dent$x[,1:3], model="OU1",param=list(decomp="diagonal"))

PCA_nas<-gm.prcomp(nas.lm)
model_fit.GPA.BM1.nas.diag<- mvBM(pruned_tree, PCA_nas$x[,1:3], model="BM1",param=list(decomp="diagonal"))
model_fit.GPA.OU1.nas.diag<- mvOU(pruned_tree, PCA_nas$x[,1:3], model="OU1",param=list(decomp="diagonal"))

PCA_neuro<-gm.prcomp(neuro.lm)
model_fit.GPA.BM1.neuro.diag<- mvBM(pruned_tree, PCA_neuro$x[,1:7], model="BM1",param=list(decomp="diagonal"))
model_fit.GPA.OU1.neuro.diag<- mvOU(pruned_tree, PCA_neuro$x[,1:7], model="OU1",param=list(decomp="diagonal"))

PCA_pjaws<-gm.prcomp(pjaws.lm)
model_fit.GPA.BM1.pjaws.diag<- mvBM(pruned_tree, PCA_pjaws$x[,1:4], model="BM1",param=list(decomp="diagonal"))
model_fit.GPA.OU1.pjaws.diag<- mvOU(pruned_tree, PCA_pjaws$x[,1:4], model="OU1",param=list(decomp="diagonal"))

PCA_max<-gm.prcomp(max.lm)
model_fit.GPA.BM1.max.diag<- mvBM(pruned_tree, PCA_max$x[,1:3], model="BM1",param=list(decomp="diagonal"))
model_fit.GPA.OU1.max.diag<- mvOU(pruned_tree, PCA_max$x[,1:3], model="OU1",param=list(decomp="diagonal"))

PCA_hyo<-gm.prcomp(hyo.lm)
model_fit.GPA.BM1.hyo.diag<- mvBM(pruned_tree, PCA_hyo$x[,1:4], model="BM1",param=list(decomp="diagonal"))
model_fit.GPA.OU1.hyo.diag<- mvOU(pruned_tree, PCA_hyo$x[,1:4], model="OU1",param=list(decomp="diagonal"))

PCA_uro<-gm.prcomp(uro.lm)
model_fit.GPA.BM1.uro.diag<- mvBM(pruned_tree, PCA_uro$x[,1:3], model="BM1",param=list(decomp="diagonal"))
model_fit.GPA.OU1.uro.diag<- mvOU(pruned_tree, PCA_uro$x[,1:3], model="OU1",param=list(decomp="diagonal"))

PCA_cer<-gm.prcomp(cer.lm)
model_fit.GPA.BM1.cer.diag<- mvBM(pruned_tree, PCA_cer$x[,1:2], model="BM1",param=list(decomp="diagonal"))
model_fit.GPA.OU1.cer.diag<- mvOU(pruned_tree, PCA_cer$x[,1:2], model="OU1",param=list(decomp="diagonal"))

PCA_oper<-gm.prcomp(oper.lm)
model_fit.GPA.BM1.oper.diag<- mvBM(pruned_tree, PCA_oper$x[,1:3], model="BM1",param=list(decomp="diagonal"))
model_fit.GPA.OU1.oper.diag<- mvOU(pruned_tree, PCA_oper$x[,1:3], model="OU1",param=list(decomp="diagonal"))

PCA_palat<-gm.prcomp(palat.lm)
model_fit.GPA.BM1.palat.diag<- mvBM(pruned_tree, PCA_neuro$x[,1:2], model="BM1",param=list(decomp="diagonal"))
model_fit.GPA.OU1.palat.diag<- mvOU(pruned_tree, PCA_neuro$x[,1:2], model="OU1",param=list(decomp="diagonal"))

PCA_uppjaws<-gm.prcomp(uppjaws.lm)
model_fit.GPA.BM1.uppjaws.diag<- mvBM(pruned_tree, PCA_uppjaws$x[,1:2], model="BM1",param=list(decomp="diagonal"))
model_fit.GPA.OU1.uppjaws.diag<- mvOU(pruned_tree, PCA_uppjaws$x[,1:2], model="OU1",param=list(decomp="diagonal"))

save(model_fit.GPA.BM1.diag, model_fit.GPA.OU1.diag, model_fit.GPA.BM1.premax.diag, model_fit.GPA.OU1.premax.diag, model_fit.GPA.BM1.ang.diag, model_fit.GPA.OU1.ang.diag, model_fit.GPA.BM1.dent.diag, model_fit.GPA.OU1.dent.diag, model_fit.GPA.BM1.nas.diag, model_fit.GPA.OU1.nas.diag, model_fit.GPA.BM1.neuro.diag, model_fit.GPA.OU1.neuro.diag, model_fit.GPA.BM1.pjaws.diag,  model_fit.GPA.OU1.pjaws.diag, model_fit.GPA.BM1.max.diag, model_fit.GPA.OU1.max.diag, model_fit.GPA.BM1.hyo.diag,  model_fit.GPA.OU1.hyo.diag,  model_fit.GPA.BM1.uro.diag, model_fit.GPA.OU1.uro.diag, model_fit.GPA.BM1.cer.diag,  model_fit.GPA.OU1.cer.diag, model_fit.GPA.BM1.oper.diag, model_fit.GPA.OU1.oper.diag, model_fit.GPA.BM1.palat.diag, model_fit.GPA.OU1.palat.diag,  model_fit.GPA.BM1.uppjaws.diag, model_fit.GPA.OU1.uppjaws.diag,  file="final_mvMORPH_GPA_diag.Rdata")


### Comparing the results of model-fitting and extracting the model parameters ###

library(mvMORPH)

#Loading the data from model-fitting
#Full model
load("final_mvMORPH_GPA.Rdata")

#Diagonal only model
load("final_mvMORPH_GPA_diag.Rdata")


### Extracting the model parameters

#Without diagonal
load("final_mvMORPH_GPA.Rdata")
load("final_mvMORPH_locsup.Rdata")

#With diagonal
load("final_mvMORPH_GPA_diag.Rdata")
load("final_mvMORPH_locsup_diag.Rdata")


#a) BM models

#Extracting the list of rate matrices
BM_results_sigmat<-list(model_fit.GPA.BM1$sigma, model_fit.locsup.BM1$sigma, model_fit.GPA.BM1.premax$sigma, model_fit.GPA.BM1.ang$sigma, model_fit.GPA.BM1.dent$sigma, model_fit.GPA.BM1.nas$sigma, model_fit.GPA.BM1.neuro$sigma, model_fit.GPA.BM1.pjaws$sigma, model_fit.GPA.BM1.max$sigma, model_fit.GPA.BM1.hyo$sigma, model_fit.GPA.BM1.uro$sigma, model_fit.GPA.BM1.cer$sigma, model_fit.GPA.BM1.oper$sigma, model_fit.GPA.BM1.palat$sigma, model_fit.GPA.BM1.uppjaws$sigma, model_fit.locsup.BM1.premax$sigma, model_fit.locsup.BM1.ang$sigma, model_fit.locsup.BM1.dent$sigma, model_fit.locsup.BM1.nas$sigma, model_fit.locsup.BM1.neuro$sigma, model_fit.locsup.BM1.pjaws$sigma, model_fit.locsup.BM1.max$sigma, model_fit.locsup.BM1.hyo$sigma, model_fit.locsup.BM1.uro$sigma, model_fit.locsup.BM1.cer$sigma, model_fit.locsup.BM1.oper$sigma, model_fit.locsup.BM1.palat$sigma, model_fit.locsup.BM1.uppjaws$sigma)

#Creating a function that calculates the sum of the diagonal elements to get the multivariate sigma
sum.diag <- function(sigma_mat) {
  sigma_mult <- sum(diag(sigma_mat))
  return(sigma_mult)
}

#Running the function throughout the list of sigma matrices
sig_vals<-sapply(BM_results_sigmat, FUN = sum.diag)
names(sig_vals)<-c("GPA.BM1", "locsup.BM1", "GPA.BM1.premax", "GPA.BM1.ang", "GPA.BM1.dent", "GPA.BM1.nas", "GPA.BM1.neuro", "GPA.BM1.pjaws", "GPA.BM1.max", "GPA.BM1.hyo", "GPA.BM1.uro", "GPA.BM1.cer", "GPA.BM1.oper", "GPA.BM1.palat", "GPA.BM1.uppjaws", "locsup.BM1.premax", "locsup.BM1.ang", "locsup.BM1.dent", "locsup.BM1.nas", "locsup.BM1.neuro", "locsup.BM1.pjaws", "locsup.BM1.max", "locsup.BM1.hyo", "locsup.BM1.uro", "locsup.BM1.cer", "locsup.BM1.oper", "locsup.BM1.palat", "locsup.BM1.uppjaws")
sig_vals

#b) OU models 

#Creating a list of the rate matrices
OU_results_sigmat<-list(model_fit.GPA.OU1$sigma, model_fit.locsup.OU1$sigma, model_fit.GPA.OU1.premax$sigma, model_fit.GPA.OU1.ang$sigma, model_fit.GPA.OU1.dent$sigma, model_fit.GPA.OU1.nas$sigma, model_fit.GPA.OU1.neuro$sigma, model_fit.GPA.OU1.pjaws$sigma, model_fit.GPA.OU1.max$sigma, model_fit.GPA.OU1.hyo$sigma, model_fit.GPA.OU1.uro$sigma, model_fit.GPA.OU1.cer$sigma, model_fit.GPA.OU1.oper$sigma, model_fit.GPA.OU1.palat$sigma, model_fit.GPA.OU1.uppjaws$sigma, model_fit.locsup.OU1.premax$sigma, model_fit.locsup.OU1.ang$sigma, model_fit.locsup.OU1.dent$sigma, model_fit.locsup.OU1.nas$sigma, model_fit.locsup.OU1.neuro$sigma, model_fit.locsup.OU1.pjaws$sigma, model_fit.locsup.OU1.max$sigma, model_fit.locsup.OU1.hyo$sigma, model_fit.locsup.OU1.uro$sigma, model_fit.locsup.OU1.cer$sigma, model_fit.locsup.OU1.oper$sigma, model_fit.locsup.OU1.palat$sigma, model_fit.locsup.OU1.uppjaws$sigma)

#Running the sum of diagonals function throughout the list of sigma matrices
sig_vals<-sapply(OU_results_sigmat, FUN = sum.diag)
names(sig_vals)<-c("GPA.BM1", "locsup.BM1", "GPA.BM1.premax", "GPA.BM1.ang", "GPA.BM1.dent", "GPA.BM1.nas", "GPA.BM1.neuro", "GPA.BM1.pjaws", "GPA.BM1.max", "GPA.BM1.hyo", "GPA.BM1.uro", "GPA.BM1.cer", "GPA.BM1.oper", "GPA.BM1.palat", "GPA.BM1.uppjaws", "locsup.BM1.premax", "locsup.BM1.ang", "locsup.BM1.dent", "locsup.BM1.nas", "locsup.BM1.neuro", "locsup.BM1.pjaws", "locsup.BM1.max", "locsup.BM1.hyo", "locsup.BM1.uro", "locsup.BM1.cer", "locsup.BM1.oper", "locsup.BM1.palat", "locsup.BM1.uppjaws")
sig_vals

#Extracting the list of stationary variance matrices
OU_results_statvarmats<-list(stationary(model_fit.GPA.OU1), stationary(model_fit.locsup.OU1), stationary(model_fit.GPA.OU1.premax), stationary(model_fit.GPA.OU1.ang), stationary(model_fit.GPA.OU1.dent), stationary(model_fit.GPA.OU1.nas), stationary(model_fit.GPA.OU1.neuro), stationary(model_fit.GPA.OU1.pjaws), stationary(model_fit.GPA.OU1.max), stationary(model_fit.GPA.OU1.hyo), stationary(model_fit.GPA.OU1.uro), stationary(model_fit.GPA.OU1.cer), stationary(model_fit.GPA.OU1.oper), stationary(model_fit.GPA.OU1.palat), stationary(model_fit.GPA.OU1.uppjaws), stationary(model_fit.locsup.OU1.premax), stationary(model_fit.locsup.OU1.ang), stationary(model_fit.locsup.OU1.dent), stationary(model_fit.locsup.OU1.nas), stationary(model_fit.locsup.OU1.neuro), stationary(model_fit.locsup.OU1.pjaws), stationary(model_fit.locsup.OU1.max), stationary(model_fit.locsup.OU1.hyo), stationary(model_fit.locsup.OU1.uro), stationary(model_fit.locsup.OU1.cer), stationary(model_fit.locsup.OU1.oper), stationary(model_fit.locsup.OU1.palat), stationary(model_fit.locsup.OU1.uppjaws))

#Creating a function that calculates the sum of the diagonal elements for the stationary variance
sum.diag2 <- function(statvarmats) {
  stat_var <- sum(diag(statvarmats))
  return(stat_var)
}

#Computing stationary variances to solve for alpha
stat_vars<-sapply(OU_results_statvarmats, FUN = sum.diag2)
names(stat_vars)<-c("GPA.OU1", "locsup.OU1", "GPA.OU1.premax", "GPA.OU1.ang", "GPA.OU1.dent", "GPA.OU1.nas", "GPA.OU1.neuro", "GPA.OU1.pjaws", "GPA.OU1.max", "GPA.OU1.hyo", "GPA.OU1.uro", "GPA.OU1.cer", "GPA.OU1.oper", "GPA.OU1.palat", "GPA.OU1.uppjaws", "locsup.OU1.premax", "locsup.OU1.ang", "locsup.OU1.dent", "locsup.OU1.nas", "locsup.OU1.neuro", "locsup.OU1.pjaws", "locsup.OU1.max", "locsup.OU1.hyo", "locsup.OU1.uro", "locsup.OU1.cer", "locsup.OU1.oper", "locsup.OU1.palat", "locsup.OU1.uppjaws")
stat_vars

#Solving for alpha
alpha<-(sig_vals/stat_vars)*0.5
