
#########################################################
###   Ancestral State Reconstruction brief tutorial   ###
###             by Olivier Larouche                   ### 
#########################################################

# These are the three packages that contain many of the most useful functions for phylogenetic comparative methods:

library(ape)
library(phytools)
library(geiger)


# We need to create a fictitious dataset to play around with the functions for ancestral state reconstruction. First a phylogeny as a backbone for the analysis. For the purpose of the exercise, I'm creating a tree with a relatively simple process here, i.e. without any extinctions. 
test.tree<-pbtree(b=1, d=0, n=100)
plot(test.tree)

# Simulating categorical data on the tree. I am running two models here, just for fun and to see what happens, an equal rates and an all rates different model. In both cases, I am feeding the rates of transition between states (the Q matrix that we talked about) through the rate parameter. It will be interesting to see if the reconstruction methods estimate these rates accurately. For simplification purpose, I am storing the results in an object called "traits" in both cases. So use one or the other of these two lines when running the rest of the script.

trait<-rTraitDisc(test.tree, model="ER", rate=0.2, k=3, ancestor=T)
trait<-rTraitDisc(test.tree, model="ARD", rate=c(0.5, 0.25, 0.1, 0.3, 0.2, 0.15), k=3, ancestor=T)

# We can plot the data on the tree. Here I am including the ancestral states so we can compare, at least qualitatively, the simulated states at the nodes to the estimated ones using the reconstruction methods.
plot(test.tree, show.tip.label = FALSE)
co <- c("blue", "yellow", "green")
Y <- trait[1:100]
A <- trait[-(1:100)]
nodelabels(A, bg = co[A])
tiplabels(Y, bg = co[Y])

# Finally, with an empirical dataset, you would most likely only have tip data, and no info at the nodes. So let's extract this tip data from our simulations:
tip.trait<-trait[1:100]

# Now that we have out data and our tree, we need to make sure that we're using the correct model of evolution for our ancestral state reconstructions. This can be done with the fitDiscrete function in the package geiger to distinguish between ER (equal rates), SYM (symmetrical rates) or ARD (all rates different) based on AICc.
result1<-fitDiscrete(test.tree, tip.trait, model="ER", transform="none")
result2<-fitDiscrete(test.tree, tip.trait, model="SYM", transform="none")
result3<-fitDiscrete(test.tree, tip.trait, model="ARD", transform="none")


# First, let's try ancestral state reconstruction using the maximum likelihood method. This method uses a continuous chain Markov model to infer evolutionary histories (as does the other method we will discuss next). The likelihood-based method  of ACR reconstructs the states (more specifically posterior probabilities for these states) at the nodes, given the distribution of the tip data. The max likelihood method will yield a single "most parsimonious" solution. Remember to change the model according to what is best supported by fitDiscrete.

trait_recon<-ace(tip.trait,test.tree, type="discrete", model="ER")

# The "trait_recon" object where we stored the results will show you the Q matrix estimates, as well as the likelihood for the root states.
trait_recon

#You can also access the posterior probabilities for each state along the nodes:
trait_recon$lik.anc

#The following block of code will use these posterior probabilities in order to generate pie charts for each nodes of the phylogeny. You can compare this to our original distribution of the trait using the simulation data. Hopefully it looks pretty similar. This plot can help determine if there are nodes that are more ambiguous in the reconstruction of their evolutionary history. And of course, if there are many of these, it can be concerning. 

plotTree(test.tree,fsize=0.8,ftype="i")

cols<-setNames(palette()[1:length(unique(tip.trait))],sort(unique(tip.trait)))
nodelabels(node=1:test.tree$Nnode+Ntip(test.tree),
           pie=trait_recon$lik.anc,piecol=cols,cex=0.5)
tiplabels(pie=to.matrix(tip.trait,sort(unique(tip.trait))),piecol=cols,cex=0.3)


# The other method for ancestral state reconstruction of discrete traits is stochastic character mapping. This also uses a continuous chain Markov model to reconstruct evolutionary histories. Among the differences, changes can occur along a branch and not only at nodes, and moreover changes are more likely to happen along longer branches. Using stochastic character mapping may be a more honest way of acknowledging that there is uncertainty in our data and in our tree. However, it does increase the analytical pipeline quite a bit. Indeed, you would normally generate 500 or 1000 of these character maps, and then loop your analyses through all of these evolutionary histories sequentially and average the results.

# You can use the make.simmap function to generate the character maps. Again, choose the model that is returned as best fitting by fitDiscrete. Also for simplification purposes, here I am only generating 100 of these evolutionary histories.

simmap_tree<-make.simmap(test.tree, tip.trait, pi="estimated" , model="ER", Q="empirical", nsim=100)

#You can plot individual character maps by changing the number in the double brackets. Does the reconstruction of evolutionary histories look like our simulated data? It is also always good, and reassuring, to check if the maps look pretty consistent.
plot(simmap_tree[[25]], ftype="off")

#We can use the describe.simmap function to summarize the data across our stochastic character maps. This will give you information such as the average number of transitions (as well as the distribution of these transitions i.e. changes of each type should sum to the the number of average changes between states). Finally you also get the absolute and relative time spent in each state.
test.stats<-describe.simmap(simmap_tree)
test.stats

#As with the max-likelihood method, you can overlay pie charts on the tree to show where the character maps differ in their reconstruction at nodes. This can be done by plotting the object generated by the describe.simmmap function. 
plot(test.stats)

# Well, that's it for now. I think this tutorial should be sufficient to get you started. We'll fill in whatever blanks are left as we see what else you will need for your analyses.

# For some additional options, you can take a look at this tutorial by Liam Revell about ACR for both continuous and discrete traits. Liam is a great methods developer of phylogenetic comparative methods in R and he has a mostly useful blog about his package phytools. I say "mostly useful" because although it contains a wealth of resources for coding analyses in R, once on the domain of his blog, it's not so easy to find what you may be looking for. Usually I am luckier finding useful tutorials from his blog via targeted google searches.

#http://www.phytools.org/eqg2015/asr.html


