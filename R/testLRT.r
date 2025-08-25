################################################################################
##                                                                            ##
##                               mvMORPH: test.LRT                            ##
##                                                                            ##
##  Created by Julien Clavel - 02-02-2015                                     ##
##  (julien.clavel@hotmail.fr/ clavel@biologie.ens.fr)                        ##
##   require: phytools, ape, corpcor, spam                                    ##
##                                                                            ##
################################################################################

## ------------------------ S3 Method for switching between LRT options ------ ##
LRT <- function(model1, model2, echo=TRUE, plot=TRUE, ...) UseMethod("LRT")

## LRT for class mvMORPH
LRT.mvmorph<-function(model1, model2, echo=TRUE, plot=TRUE,...){

## Options (TODO)
args <- list(...)
if(is.null(args[["simulations"]])){ simulations <- TRUE }else{ simulations <- FALSE}


##-------------------LRT comparison of the models-----------------------------##

if(any(class(model1)[2]!=class(model2)[2])){warning("You are using an LR test which is probably not at its nominal level!! You should use simulations-based LRT instead")}

if(model1$param$nparam<=model2$param$nparam){
    mod1<-model1;mod2<-model2
    L2<-model1$LogLik
    L1<-model2$LogLik
    model1<-mod2;model2<-mod1
    warning("The order of the models has been changed for the comparison")
}else{
L1<-model1$LogLik
L2<-model2$LogLik
}
# Set names
if(class(model1)[2]=="mvmorph.bm"){
    model1$param$model<-if(model1$param$constraint=="equal"){
        paste(model1$param$model," equal variance/rates")
    }else if(model1$param$constraint=="shared"){
        paste(model1$param$model," shared eigenvectors")
    }else if(model1$param$constraint=="proportional"){
        paste(model1$param$model," proportional")
    }else if(model1$param$constraint=="correlation"){
        paste(model1$param$model," shared correlation")
    }else if(model1$param$constraint=="diagonal"){
        paste(model1$param$model," diagonal")
    }else if(model1$param$constraint=="variance"){
        paste(model1$param$model," shared variance")
    }else if(model1$param$constraint=="equaldiagonal"){
        paste(model1$param$model," equal diagonal")
    }else{model1$param$model}
}

if(class(model2)[2]=="mvmorph.bm"){
    model2$param$model<-if(model2$param$constraint=="equal"){
        paste(model2$param$model," equal variance/rates")
    }else if(model2$param$constraint=="shared"){
        paste(model2$param$model," shared eigenvectors")
    }else if(model2$param$constraint=="proportional"){
        paste(model2$param$model," proportional")
    }else if(model2$param$constraint=="correlation"){
        paste(model2$param$model," shared correlation")
    }else if(model2$param$constraint=="diagonal"){
        paste(model2$param$model," diagonal")
    }else if(model2$param$constraint=="variance"){
        paste(model2$param$model," shared variance")
    }else if(model2$param$constraint=="equaldiagonal"){
        paste(model2$param$model," equal diagonal")
    }else{model2$param$model}
}

if(class(model1)[2]=="mvmorph.ou"){
    model1$param$model<-if(model1$param$decomp=="diagonal"){
        paste(model1$param$model," diagonal")
    }else if(model1$param$decomp=="equal"){
        paste(model1$param$model," equal")
    }else if(model1$param$decomp=="eigen"){
        paste(model1$param$model," symmetric")
    }else if(model1$param$decomp=="eigen+" | model1$param$decomp=="cholesky" | model1$param$decomp=="spherical"){
        paste(model1$param$model," symmetric positive")
    }else if(model1$param$decomp=="lower"){
        paste(model1$param$model," lower triangular")
    }else if(model1$param$decomp=="upper"){
        paste(model1$param$model," upper triangular")
    }else if(model1$param$decomp=="qr" | model1$param$decomp=="svd" | model1$param$decomp=="schur"){
        paste(model1$param$model," non-symmetric")
    }else if(model1$param$decomp=="qr+" | model1$param$decomp=="svd+" | model1$param$decomp=="schur+"){
        paste(model1$param$model," non-symmetric positive")
    }else{model1$param$model}
}

if(class(model2)[2]=="mvmorph.ou"){
    model2$param$model<-if(model2$param$decomp=="diagonal"){
        paste(model2$param$model," diagonal")
    }else if(model2$param$decomp=="equal"){
        paste(model2$param$model," equal")
    }else if(model1$param$decomp=="eigen"){
        paste(model1$param$model," symmetric")
    }else if(model2$param$decomp=="eigen+" | model2$param$decomp=="cholesky" | model2$param$decomp=="spherical"){
        paste(model2$param$model," symmetric positive")
    }else if(model2$param$decomp=="lower"){
        paste(model2$param$model," lower triangular")
    }else if(model2$param$decomp=="upper"){
        paste(model2$param$model," upper triangular")
    }else if(model2$param$decomp=="qr" | model2$param$decomp=="svd" | model2$param$decomp=="schur"){
        paste(model2$param$model," non-symmetric")
    }else if(model2$param$decomp=="qr+" | model2$param$decomp=="svd+" | model2$param$decomp=="schur+"){
        paste(model2$param$model," non-symmetric positive")
    }else{model2$param$model}
    
}

#
LRT<-(2*((L1-L2)))
#difference in degrees of freedom
ddf<-model1$param$nparam-model2$param$nparam

LRT.prob<-pchisq(LRT,ddf,lower.tail=FALSE)

if(echo==TRUE){
    if(LRT.prob<0.001){signif<-c("***")}else if(LRT.prob<0.01){
        signif<-c("**") }else if(LRT.prob<0.05){signif<-c("*")}else if(LRT.prob<0.1){signif<-c(".")}else{signif<-""}
    cat("-- Log-likelihood Ratio Test --","\n")
    cat("Model",model1$param$model," versus ",model2$param$model,"\n")
    cat("Number of degrees of freedom :",ddf,"\n")
    cat("LRT statistic:",LRT," p-value:",LRT.prob,signif,"\n")
    cat("---","\n")
    cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1")
}

results<-list(pval=LRT.prob, ratio=LRT, ddf=ddf, model1=model1$param$model, model2=model2$param$model)
class(results)<-c("mvmorph.lrt")
invisible(results)
#End
}


# ------------------------------------------------------------------------- #
# LRTsim                                                                    #
# options: model1, model2, nsim=100, plot=TRUE, nbcores=1, ...              #
# Compute the log-likelihood ratio test (using asymptotic                   #
# solution (nested models) and simulation)                                  #
# ------------------------------------------------------------------------- #


LRT.mvgls <- function(model1, model2, echo=TRUE, plot=TRUE, ...){
  
  args <- list(...)
  if(is.null(args[["nsim"]])) nsim = 100 else nsim = args$nsim
  if(is.null(args[["nbcores"]])) nbcores = 1L else nbcores = args$nbcores
  if(is.null(args[["alternative"]])) alternative = FALSE else alternative = args$alternative
  if(is.null(args[["parametric"]])) parametric = TRUE else parametric = args$parametric
  if(is.null(args[["REML"]])) REML = FALSE else REML = args$REML
  
  # # to remove
  #require(pbmcapply)
  
  # # generate the test statistic:
  if((model1$REML | model2$REML) & REML==FALSE){
    flag = TRUE
    ll1 = .reml_to_ml(model1)
    ll2 = .reml_to_ml(model2)
    lrt <- -2*(ll1 - ll2)
  }else{
    flag = FALSE
    ll1 = model1$logLik
    ll2 = model2$logLik
    lrt <- -2*(ll1 - ll2)
  }
  
  # just for coercion
  lrt <- as.numeric(lrt)
  
  # generate random samples from model 1
  if(parametric) new_data <- simulate(model1, nsim=nsim) else new_data <- sbootstrap(model1, nboot=nsim)
 
  names_data = model1$corrSt$phy$tip.label

  # prepar models
  modelNull <- model1$call
  modelAlt <- model2$call
  modelNull$start <- model1$opt$par
  modelAlt$start <- model2$opt$par
  modelNull$grid.search <- modelAlt$grid.search <- quote(FALSE)
  if(modelNull$method == "EmpBayes"){
      modelNull$MMSE <- modelAlt$MMSE <- quote(FALSE)
      modelNull$FCI <- modelAlt$FCI <- quote(FALSE)
      modelNull$Hessian <- modelAlt$Hessian <- quote(FALSE)
  }
  
  # loop over new datasets
  stat_dist <- pbmcmapply(function(x){
    rownames(new_data[[x]]) = names_data
    modelNull$response <- modelAlt$response <- quote(new_data[[x]]);
    estimModelNull <- eval(modelNull);
    estimModelalt <- eval(modelAlt);
    if(flag){
      -2*(.reml_to_ml(estimModelNull) - .reml_to_ml(estimModelalt))
    }else{
      -2*(estimModelNull$logLik -estimModelalt$logLik)
    }
  },1:nsim, mc.cores= getOption("mc.cores", nbcores))
  
  # compute the p-value
  lrtpval <- mean(as.numeric(lrt)<=stat_dist)
  
  ## Check if comparison to the alternative is wanted?
  if(alternative){
    if(parametric) new_data2 <- simulate(model2, nsim=nsim) else new_data2 <- sbootstrap(model2, nboot=nsim)

    # loop over new datasets
    stat_dist2 <- pbmcmapply(function(x){
      rownames(new_data2[[x]]) = names_data
      modelNull$response <- modelAlt$response <- quote(new_data2[[x]]);
      estimModelNull <- eval(modelNull);
      estimModelalt <- eval(modelAlt);
      if(flag){
        -2*(.reml_to_ml(estimModelNull) - .reml_to_ml(estimModelalt))
      }else{
        -2*(estimModelNull$logLik -estimModelalt$logLik)
      }
    },1:nsim, mc.cores= getOption("mc.cores", nbcores))
    
  }
  
  # plot
  if(plot & alternative==FALSE){
      limits_x = c(min(stat_dist,lrt),max(stat_dist,lrt))
    hist(stat_dist, freq = FALSE, breaks=50, las=1, main=paste("LRT:",round(lrt, digits=3), "p-value", round(lrtpval, digits=5)),
         xlab="Null distribution", xlim=limits_x, ...); abline(v=lrt)
  }else if(plot & alternative){
      limits_x = c(min(stat_dist,stat_dist2,lrt),max(stat_dist,stat_dist2,lrt))
    hist(stat_dist, freq = FALSE, breaks=50, las=1, main=paste("LRT:",round(lrt, digits=3), "p-value", round(lrtpval, digits=5)),
                xlab="Likelihood ratio", xlim = limits_x, ...)
    hist(stat_dist2, freq = FALSE, breaks=50, las=1, add=TRUE, col="red"); abline(v=lrt, lty=2)
  }
  
  # TODO :> use the print options from LRT (define it as class(results)<-c("mvmorph.lrt"))
  # print
  cat("LRT test (non-parametric)", lrt," p-value:",lrtpval, "log-lik model 1:",ll1, "log-lik model 2:", ll2)
  
  # results
  if(alternative){
      results = list(ratio=lrt, dist=stat_dist, pval=lrtpval, dist_alt=stat_dist2)
  } else{
      results = list(ratio=lrt, dist=stat_dist, pval=lrtpval)
  }
  
  invisible(results)
  
}

## Function to perform semi-parametric bootstrap from model fit by mvgls
sbootstrap <- function(object, nboot, ...){
  
  # compute the square root matrix for the model fit
  phyloSqrt <- pruning(object$corrSt$phy, trans=FALSE, inv=FALSE)$sqrtM
  
  # retrieve the expected values (ancestral states)
  expect = fitted(object)
  
  # the normalized/decorrelated residuals
  resid <- residuals(object, type="normalized")
  
 
     # generate datasets
     sim <- lapply(1:nboot, function(j){
       bootstrap_residuals <- phyloSqrt %*% resid[c(sample(Ntip(object$corrSt$phy)-1, replace = TRUE),Ntip(object$corrSt$phy)),] # remove intercept from boostrapping if the pruning algorithm is used
       data = expect + bootstrap_residuals
       rownames(data) = object$corrSt$phy$tip.label
       data
       })
  
  return(sim)
}
