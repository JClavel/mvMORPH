################################################################################
##                                                                            ##
##                               mvMORPH: util.r                              ##
##                                                                            ##
##   Internal functions for the mvMORPH package                               ##
##                                                                            ##
##  Created by Julien Clavel - 26-05-2016                                     ##
##  (julien.clavel@hotmail.fr/ julien.clavel@biologie.ens.fr)                 ##
##   require: phytools, ape, corpcor, subplex, spam                           ##
##                                                                            ##
################################################################################

# Function to extract the aic-weights
aicw <- function(x,...){
    
    args <- list(...)
    if(is.null(args[["aicc"]])) args$aicc <- FALSE
    
    if(inherits(x, "list")){
        if(inherits(x[[1]],"mvmorph")){
            
            if(args$aicc==TRUE){
                aic_model <- sapply(1:length(x),function(i) x[[i]]$AICc)
            }else{
                aic_model <- sapply(x,AIC)
            }
        
                models_names <- sapply(1:length(x),function(i){
                    if(!is.null(x[[i]]$param[["constraint"]])){
                        paste(x[[i]]$param$model[length(x[[i]]$param$model)],x[[i]]$param$constraint,i)
                    }else{
                        paste(x[[i]]$param$model[length(x[[i]]$param$model)],i)}
                })
        }else{
                aic_model <- unlist(x)
                models_names <- as.character(1:length(aic_model))
        }
        
        aics <- data.frame(models=models_names, AIC=aic_model, diff=aic_model)
        row.names(aics) <- as.character(models_names)
        
    }else{
        if(is.null(names(x))){
            models_names <- as.character(1:length(x))
        }else{
            models_names <- names(x)
        }
        
        aics <- data.frame(models=models_names, AIC=x, diff=x)
        row.names(aics) <- as.character(models_names)
    }
    
    aics <- aics[order(-aics$AIC),]
    for(i in 1:length(x)){aics$diff[i] <- aics$AIC[i]-min(aics$AIC)}
    aics$wi <- exp(-0.5*aics$diff)
    aics$aicweights <- aics$wi/sum(aics$wi)
    aics <- aics[models_names,] # reorder the results to the original order

    
    class(aics) <- c("mvmorph.aicw")
   return(aics)
}

# ------------------------------------------------------------------------- #
# Function to transform ancestral states reconstruction to a tree of class  #
# "simmap"  (push to mvMORPH? current version is a bit quick & dirty)       #
#                                                      J. CLAVEL - 2020     #
# ------------------------------------------------------------------------- #

# tree = the phylogenetic tree
# ancestral = either a "describe.simmap" object from phytools or an "ace" object from ape.
# tips = states at the tips. Check that the tree and data are in same order!

mapping.asr <- function(tree, ancestral, tips){
  
  if(inherits(ancestral, "describe.simmap")){
    names_grps <- colnames(ancestral$ace)
    statesNodes <- names_grps[apply(ancestral$ace, 1, which.max)]
  }else{
    names_grps <- colnames(ancestral$lik.anc)
    statesNodes <- names_grps[apply(ancestral$lik.anc, 1, which.max)]
  }
  
  combined = as.character(c(tips, statesNodes))
  treebis=tree
  for(i in sort(tree$edge[,2])){
    treebis <- paintBranches(treebis, i, combined[i], anc.state=combined[Ntip(tree)+1])
  }
  return(treebis)
}
