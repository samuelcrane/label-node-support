# labelNodeSupport: map node support values from parsimony, likelihood, and Bayesian phylogenies onto a single target tree.
# This script was put together by Rich Glor and first published (as far as I can tell) in 2008 on the Dechronization blog.
# http://treethinkers.blogspot.com/2008/10/labeling-trees-posterior-probability.html
# I've updated the code (some of the syntax was depreciated), made some corrections, and am in the process of expanding the function by
# (a) adding in the ability to map a third set of support values and
# (b) deal with the verbose (as well as the simple) output of MrBayes sumt consensus command [SNC]

#Fist we need to open some necessary libraries
library(ape)
library(geiger)
library(phyloch) # Can be ignored if using the conformat=simple option in MrBayes sumt command

# The getAllSubTrees function below is a necessary subfunction that atomizes a tree into each 
# individual subclade and was provided compliments of Luke Harmon. [RG]
getAllSubtrees <- function(phy, minSize=2) 
{
    res <- list() 
    count = 1 
    ntip <- length(phy$tip.label) 
    for(i in 1:phy$Nnode) 
    { 
        l <- tips(phy, ntip+i) 
        bt <- match(phy$tip.label, l) 
        if(sum(is.na(bt)) == 0) 
        {
            st <- phy 
        } 
        else st <- drop.tip(phy, phy$tip.label[is.na(bt)]) 
        if(length(st$tip.label)>=minSize) 
        { 
            res[[count]] <- st 
            count <- count+1 
        }
    } 
    res
}

# The plotBayesBoot function below plots both posterior probability and bootstrap values on each 
# node of the consensus tree obtained from your Bayesian analysis. Bootstrap values will appear in 
# bold text immediately below and to the left of the node they support, whereas Bayesian posterior 
# probabilies will appear in regular face above and to the left of the node. [RG]
plotBayesBoot <- function(bayesTree,bootTree) 
{
    getAllSubtrees(bayesTree) -> bayesSub
    getAllSubtrees(bootTree) -> bootSub
    bootList <- matrix("<50", Nnode(bayesTree), 1)

    #The commands below compare all the subclades in the Bayes tree to all the subclades in the bootstrap tree, and vice versa, and identifies all those clades that are identical.
    for(i in 1:Nnode(bayesTree)) 
    {
        for(j in 1:Nnode(bootTree)) 
        {
            match(bayesSub[[i]]$tip.label[order(bayesSub[[i]]$tip.label)], bootSub[[j]]$tip.label[order(bootSub[[j]]$tip.label)]) -> shared
            match(bootSub[[j]]$tip.label[order(bootSub[[j]]$tip.label)], bayesSub[[i]]$tip.label[order(bayesSub[[i]]$tip.label)]) -> shared2
            if(sum(is.na(c(shared,shared2)))==0) 
            {
                bootTree$node.label[j] -> bootList[i]
            }
        }
    }
    plot(ladderize(bayesTree), cex=0.5, lwd=0.5, use.edge.length=FALSE, direction='r') #Plots your Bayesian consensus tree
    nodelabels(bayesTree$node.label, adj=c(1.2, -0.3), frame="n", cex=0.5, font=1) #Adds posterior probability values to the tree. Change the 'cex' value to make the fond smaller or larger. A value of 1 will give you a readable result in the R quartz window, but a value closer to 0.25 might be better for publication)
    nodelabels(bootList, adj=c(1.4, 1.3), frame="n", cex=0.5, font=2) #Adds bootstrap values.
}

# Use plotBayesBootTwo to map three sets of suport values onto the target tree. This works the same
# as plotBayesBoot but iterates across a third tree to compare with the first tree. The values are 
# drawn in a different format as well. 
plotBayesBootTwo <- function(bayesTree, bootTree1, bootTree2) 
{
    getAllSubtrees(bayesTree) -> bayesSub
    getAllSubtrees(bootTree1) -> bootSub1
    getAllSubtrees(bootTree2) -> bootSub2
    bootList1 <- matrix("<50", Nnode(bayesTree), 1)
    bootList2 <- matrix("<50", Nnode(bayesTree), 1) #Watch this <50 value if using a scrict consensus tree
    
    #The commands below compare all the subclades in the Bayes tree to all the subclades in the bootstrap tree, and vice versa, and identifies all those clades that are identical.
    for(i in 1:Nnode(bayesTree)) 
    {
        for(j in 1:Nnode(bootTree1)) 
        {
            match(bayesSub[[i]]$tip.label[order(bayesSub[[i]]$tip.label)], bootSub1[[j]]$tip.label[order(bootSub1[[j]]$tip.label)]) -> shared
            match(bootSub1[[j]]$tip.label[order(bootSub1[[j]]$tip.label)], bayesSub[[i]]$tip.label[order(bayesSub[[i]]$tip.label)]) -> shared2
            if(sum(is.na(c(shared,shared2)))==0) 
            {
                bootTree1$node.label[j] -> bootList1[i]
            }
        }
        for(j in 1:Nnode(bootTree2)) 
        {
            match(bayesSub[[i]]$tip.label[order(bayesSub[[i]]$tip.label)], bootSub2[[j]]$tip.label[order(bootSub2[[j]]$tip.label)]) -> shared
            match(bootSub2[[j]]$tip.label[order(bootSub2[[j]]$tip.label)], bayesSub[[i]]$tip.label[order(bayesSub[[i]]$tip.label)]) -> shared2
            if(sum(is.na(c(shared,shared2)))==0) 
            {
                bootTree2$node.label[j] -> bootList2[i]
            }
        }
    }
    plot(ladderize(bayesTree), cex=0.5, lwd=0.5, direction='r', use.edge.length=FALSE, label.offset=1, x.lim=c(0.0000, 300), y.lim=c(0, 110)) #Plots your Bayesian consensus tree
    nodelabels(bayesTree$node.label, adj=c(-0.1, -0.1), frame="n", cex=0.5, font=2) #Adds posterior probability values to the tree. Change the 'cex' value to make the fond smaller or larger. A value of 1 will give you a readable result in the R quartz window, but a value closer to 0.25 might be better for publication)
    nodelabels(bootList1, adj=c(1.2, -0.5), frame="n", cex=0.5, font=1) #Adds bootstrap values.
    nodelabels(bootList2, adj=c(1.2, 1.5), frame="n", cex=0.5, font=3) #Adds bootstrap values.
}


# Read in two trees
read.nexus("yourBayesTree.con")->bayesTree #Reads in the .con file that results from analyses in MrBayes.
bayesTree[[1]]->bayesTree #Extracts one of the two trees in the .con file.
read.nexus("yourBootTree.nex")->bootTree #Reads in the consensus tree from a bootstrap analysis in PAUP.
plotBayesBoot(bayesTree, bootTree) # For two trees

# Read in three trees
read.nexus("yourBayesTree.con")->bayesTree #Reads in the .con file that results from analyses in MrBayes.
bayesTree[[1]]->bayesTree #Extracts one of the two trees in the .con file.
read.nexus("yourBootTree.nex")->bootTree1
read.nexus("yourOtherBootTree.nex")->bootTree2
plotBayesBootTwo(bayesTree, bootTree1, bootTree2) # For three trees (e.g. Bayesian, ML, and MP)
