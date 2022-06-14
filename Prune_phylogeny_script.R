require(ape)
require(phytools)
require(geiger)


###################################################################
#Read In Tree in this case I'm using the Rabosky 2018 tree
tree<-read.tree("actinopt_12k_treePL.tre")
#plot(tree,cex=0.5)

####Read In data to match the tips of the tree witht the data points
data <- read.csv("caran_coords.csv",header=T, row.names=1)
TreeOnly <- setdiff(tree$tip.label,rownames(data))
TreeOnly # Enter the name of the object we just created to see what's in it.
DataOnly <- setdiff(rownames(data), tree$tip.label)
DataOnly # Enter to see what species are in the data set but not the tree.
# In our case, we have overlap issues in both directions. Because we have data for fewer taxa than we have in our phylogeny, let's first prune our tree to just those species in the tree that were also measured before proceeding further.
# We'll prune the tree using drop.tip. We need to give it our tree, and a list of species to prune. We'll use the TreeOnly list of species names we just made to prune these species from the tree.

pruned_tree <- drop.tip(tree,TreeOnly)

##############################################################
#Ladderize Tree########################
##############################
#Load Tree
##############################
phyloTime <- pruned_tree # Load a ultrametric tree
phyloTimeLadderized <- (ladderize(phyloTime))  # Ladderization
phyloTimeLadderized <- rescale(phyloTimeLadderized, "depth", 1) #This rescaling will make subsequent plotting functions somewhat easier. Even more importantly, it will often improve the performance of likelihood functions
plot(phyloTimeLadderized, cex=0.5, no.margin = T) #Plot ladderized and rescaled tree
add.scale.bar() # Add a simple scale bar indicating the scale for the branches in your tree

write.tree(phyloTimeLadderized, "phyloTimeLadderized_1.nwk")
phyloTimeLadderized <- read.tree("phyloTimeLadderized_1.nwk")


tree<-phyloTimeLadderized