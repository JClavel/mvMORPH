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
    
    if(class(x)=="list"){
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
        models_names <- 1:length(aic_model)
        }
        aics <- data.frame(models=models_names, AIC=aic_model, diff=aic_model)
        row.names(aics) <- as.character(models_names)
        
    }else{
        if(is.null(names(x))){
            models_names <- names(x)
        }else{
            models_names <- 1:length(x)
        }
        aics <- data.frame(models=models_names, AIC=x, diff=x)
        row.names(aics) <- as.character(models_names)
    }
    
    aics <- aics[order(-aics$AIC),]
    for(i in 1:length(x)){aics$diff[i] <- aics$AIC[i]-min(aics$AIC)}
    aics$wi <- exp(-0.5*aics$diff)
    aics$aicweights <- aics$wi/sum(aics$wi)
    aics <- aics[sort(row.names(aics), decreasing=FALSE),]
    
    class(aics) <- c("mvmorph.aicw")
   return(aics)
}


