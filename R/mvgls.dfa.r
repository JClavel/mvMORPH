################################################################################
##                                                                            ##
##                       mvMORPH: mvgls.dfa.r                                 ##
##                                                                            ##
##   Discriminant Function Analysis based on GLS model fit                    ##
##                                                                            ##
##  Created by Julien Clavel - 01-11-2020                                     ##
##  (julien.clavel@hotmail.fr/ julien.clavel@univ-lyon.fr)                    ##
##   require: phytools, ape, corpcor, subplex, spam, glassoFast, stats        ##
##                                                                            ##
################################################################################


# ------------------------------------------------------------------------- #
# mvgls.dfa                                                                 #
# options: object, ...                                                      #
#                                                                           #
# ------------------------------------------------------------------------- #

# return the priors, the discriminant axes (coeffs), the scores, the % explained by all factors, all original fit needed parameter or the original object fit?
mvgls.dfa <- function(object, ...){
  
  # options
  args <- list(...)
  if(is.null(args[["term"]])){
    term <- object$dims$assign[which(attr(object$terms,"dataClasses")=="factor")[1]]
    if(is.na(term)) stop("there's no factor terms in the model")
  } else term <- args$term
  if(is.null(args[["type"]])) type <- "I" else type <- args$type
  if(is.null(args[["normalized"]])) normalized <- FALSE else normalized <- args$normalized
  
  # Number of degrees of freedom
  if(object$REML) ndimCov = object$dims$n - object$dims$m else ndimCov = object$dims$n
  asgn <- object$dims$assign
  
  # variables and parameters
  Y <- object$corrSt$Y
  X <- object$corrSt$X
  N <- nrow(Y)
  
  # error matrix
  E <- object$sigma$Pinv*ndimCov
  
  # Projection matrix.
  H <- .hypothesis(object, Y, X, type=type, term=term)
  
  # decomposition
  decompE <- eigen(E)
  U <- decompE$vectors %*% diag(1/sqrt(decompE$values)) # Rencher 2002, p. 279 => U^-1
  UHU <- t(U) %*% H %*% U # Transform the hypothesis matrix with the inverse square root matrix (~cholesky factor)
  avalues <- eigen(UHU, symmetric=TRUE)
  pct <- 100 * avalues$values / sum(avalues$values)
  if(is.null(args[["tol"]])) tol <- max(dim(UHU))*max(avalues$values)*.Machine$double.eps else tol <- args$tol
  rank <- min(qr(H)$rank, sum(avalues$values>tol))
  
  # canonical discriminants:
  coeffs.raw <- U %*% avalues$vectors * sqrt(object$dims$n - object$dims$m) 
  coeffs.std <- diag(sqrt(diag(object$sigma$Pinv))) %*% coeffs.raw #8.17 Rencher 2002
  discrim <- paste("Discr.",1:ncol(coeffs.raw))
  colnames(coeffs.raw) = colnames(coeffs.std) = discrim
  rownames(coeffs.raw) = rownames(coeffs.std) = colnames(object$variable$Y)
  
  # scores (where Y are the residuals to the grand mean)
  if(any(asgn<=0)){
    var_reduc=which(asgn<=0) # zero for now FIXME for more general model
    Ystand <- object$variable$Y - object$variable$X[,var_reduc, drop=FALSE]%*%pseudoinverse(X[,var_reduc, drop=FALSE])%*%Y
  }else{
    # FIXME:  is producing two discriminants...
    warning("The function is not working yet for models without intercepts")
    Xdesign <- matrix(1,ncol=1, nrow=N)
    if(!is.null(object$corrSt$diagWeight)) Xmean <- crossprod(pruning(object$corrSt$phy, trans=FALSE)$sqrtMat, Xdesign*(1/object$corrSt$diagWeight))
    else Xmean <- crossprod(pruning(object$corrSt$phy, trans=FALSE)$sqrtMat, Xdesign)
    Ystand <- object$variable$Y - Xdesign%*%pseudoinverse(Xmean)%*%Y
  }

  # compute the scores
  scores <- Ystand %*% coeffs.raw
  
  # define priors
  nclass <- sum(asgn==term)
  if(any(asgn<=0)) nclass = nclass + 1
  if(is.null(args[["prior"]])){
    prior <-  rep(1/nclass, nclass)
  }else{
    prior <- args$prior
  }
  
  # classes ID for covariate models
  if(any(asgn<=0)) classid <- asgn%in%c(0,term) else classid <- asgn==term
  
  # checks
  if(length(prior)!=nclass) warning("The number of classes and  priors doesn't match")
  
  # results
  results <- list(coeffs=coeffs.raw, coeffs.std=coeffs.std
                  , scores=scores, H=H, E=E, residuals=Ystand
                  ,rank=rank, pct=pct, prior=prior, nclass=nclass, classid=classid, term=term, fit=object)
  
  class(results) <- "mvgls.dfa"
  return(results)
  
}

# ------------------------------------------------------------------------- #
# .hypothesis                                                               #
# options: object, Y, X, type="I", term=1, ...                              #
#                                                                           #
# ------------------------------------------------------------------------- #
.hypothesis <- function(object, Y, X, type="I", term=1, ...){
  
  # Number of degrees of freedom
  asgn <- object$dims$assign
  # FIXME handle singular design
  
  switch (type,
          "I" = {
            # QR decomposition for the Hypothesis matrix
            Q_r <- qr(X)
            Q <- qr.Q(Q_r)
            
            variables=which(asgn==term) # set to term 1 for now: FIXME: do a function
            QQl <- Q[, variables] %*% t(Q[, variables])
            H <- t(Y) %*% QQl %*% Y
          },
          "II" = {
            
            asgn <- asgn[asgn!=0]
            intercept = TRUE # (TODO: handle cases where a model without intercept is fitted => type III)
            model_terms <- terms(object$formula)
            facTerms <- crossprod(attr(model_terms, "factors"))
            facTerms <- facTerms[,asgn,drop=FALSE]
            
            # Indexing of the design matrix
            var_reduc <- var_full <- facTerms[term,]<unique(facTerms[term,which(asgn==term)])
            var_full[which(asgn==term)] <- TRUE
            
            # full model
            Proj_full <- X[,c(intercept,var_full)] %*% pseudoinverse(X[,c(intercept,var_full), drop=FALSE])
            
            # reduced model
            Proj_reduc <- X[,c(intercept,var_reduc)] %*% pseudoinverse(X[,c(intercept,var_reduc), drop=FALSE])
            
            # Hypothesis SSCP matrix
            H <- crossprod(Y, (Proj_full - Proj_reduc) %*% Y)
          },
          "III"={
            intercept = NULL  
            
            # Indexing of the design matrix
            var_reduc=which(asgn!=(term-1))
            
            # reduced model
            Proj_reduc <- X[,c(intercept,var_reduc)] %*% pseudoinverse(X[,c(intercept,var_reduc), drop=FALSE])
            
            # Hypothesis SSCP matrix
            H <- crossprod(Y, (Proj_full - Proj_reduc) %*% Y)
          })
  
  return(H)
}


# ------------------------------------------------------------------------- #
# plot.mvgls.dfa                                                            #
# options: x, ..., dims=c(1,2), type=c("raw","std"), biplot=FALSE           #
#                                                                           #
# ------------------------------------------------------------------------- #
plot.mvgls.dfa <- function(x, ..., dims=c(1,2), type=c("raw","std"), biplot=FALSE){
  type = match.arg(type[1], c("raw","std"))
  if(x$rank==1) dims = 1
  if(length(dims)>2) stop("The maximum number of discriminant axes to plot is 2")
  names_variables <- attr(x$fit$variables$X,"dimnames")[[2]][x$classid]
  grp <- as.factor(x$fit$variables$X[,x$classid]%*%(1:x$nclass))
  names(grp)= names_variables
  if(length(dims)==1){
    boxplot(x$scores[,dims]~grp, xlab="Classes", ylab="Discr. 1", names=names_variables,...) #FIXME
  }else{
    if(type=="raw") scores = x$scores[,dims] 
    else scores = (x$residuals%*%x$coeffs.std)[,dims]
    
    # plot the scores
    plot(scores, xlab=paste("Discr.",dims[1]), ylab=paste("Discr.",dims[2]), ...)
    
    # plot the mean for each classes. FIXME
    # coeff(fit)%*%x$coeffs[,dims] # if the design matrix is not a treatment contrast?
    mu <- pseudoinverse(model.matrix(~grp+0))%*%scores
    points(mu, pch=8, col="red")
    
    if(biplot){
        if(type=="raw") coefficients = x$coeffs else if(type=="std") coefficients = x$coeffs.std
        scaling <- min(apply(scores,2,range)/apply(coefficients[,dims],2,range))
        if(scaling<1) scaling <- scaling - 0.05*scaling else scaling <- 1
        arrows(x0=0, y0=0, x1=scaling*coefficients[,dims[1]], y1=scaling*coefficients[,dims[2]])
        text(scaling*coefficients[,dims], labels=rownames(coefficients))
        abline(h=0,v=0,lty=2)
    }
  }
}

# ------------------------------------------------------------------------- #
# print.mvgls.dfa                                                           #
# options: x, ...                                                           #
#                                                                           #
# ------------------------------------------------------------------------- #
print.mvgls.dfa <- function(x, ...){
  cat("\n")
  message("-- Canonical/Linear Discriminant Analysis --","\n")
  cat("Prior probabilities for grouping factors:","\n")
  cat(x$prior,"\n")
  cat("\n")
  cat("Number of discriminant axes:",x$rank,"\n")
  discrim <- paste("Discr.",1:x$rank,":")
  pct = signif(x$pct[1:x$rank], digits = 2)
  names(pct) = discrim
  cat(paste(discrim, pct,"%"),"\n")
  cat("\n")
  if(nrow(x$coeffs)>10){
    cat("Coefficients of linear discriminants [truncated]:","\n")
    print(x$coeffs[1:10,1:x$rank])
  }else{
    cat("Coefficients of linear discriminants:","\n")
    print(x$coeffs[,1:x$rank])
  }
  cat("\n")
  if(nrow(x$coeffs.std)>10){
    cat("Standardized coefficients of linear discriminants [truncated]:","\n")
    print(x$coeffs.std[1:10,1:x$rank])
  }else{
    cat("Standardized coefficients of linear discriminants:","\n")
    print(x$coeffs.std[,1:x$rank])
  }
  cat("\n")
}

# ------------------------------------------------------------------------- #
# predict.mvgls.dfa                                                         #
# options: object, newdata, prior, ...                                      #
#                                                                           #
# ------------------------------------------------------------------------- #
predict.mvgls.dfa <- function(object, newdata, prior = object$prior, ...){
    
    args <- list(...)
    if(is.null(args[["tree"]])) tree <- NULL else tree <- args$tree
    # Checks that the model is based on dummy contrasts
    contrasts_types <- object$fit$contrasts[attr(object$fit$terms, "term.labels")]
    if(contrasts_types[object$term]!="contr.treatment") warning("Model assumed that the contrasts for the term of interest is of type \"contr.treatment\"")
    # Check if there's extra terms in the model
    if(sum(object$classid)!=ncol(object$fit$variables$X)) stop("The DFA function is not currently working for multi-predictors models") # FIXME
    
    # Retrieve coefficients
    B <- coef(object$fit)
    index_B <- (1:nrow(B))[object$classid]
    
    # FIXME: works for treatment contrasts. To be generalized
    if(any(object$fit$dims$assign==0)){
        B[object$fit$dims$assign==object$term,] <-  sweep(B[object$fit$dims$assign==object$term,,drop=FALSE], 2, B[object$fit$dims$assign==0,], "+")
    }
    
    # estimated (inverse) covariance matrix
    Rinv <- object$fit$sigma$P
    
    # scale the covariance by the determinant of the evolutionary model
    # this make the vcv to the same scale as the traits.
    # should make the distance measure consistant with the Bayes rule between both OLS and GLS approaches.
    if(!all(prior==prior[1])){
        if(object$fit$REML){ #TODO handle "const" in REML determinant for OUM
            Rinv <- Rinv * (1/exp( (object$fit$corrSt$det - determinant(crossprod(object$fit$corrSt$X))$modulus) * (1/object$fit$dims$n)))
        }else{
            Rinv <- Rinv * (1/exp(object$fit$corrSt$det * (1/object$fit$dims$n)))
        }
    }
    
    # log-sum-exp trick to avoid over/under flow
    logsumexp <- function(v) max(v) + log(sum(exp(v - max(v))))
    
    if(is.null(tree)){
        #if no datasets provided, classification rule is applied to the training dataset > see 9.12 in Rencher 2002
        # Bayes classifier for group "g" is d(x) = t(x)S^-1u - 0.5t(u)S^-1u + log(P(g)) => can be used to precompute quantities for multiple values to predict
        if(missing(newdata)){
            prediction <- sapply(index_B, function(i){
                SB <- Rinv%*%B[i,]
                const_prior <- -0.5*t(B[i,])%*%SB + log(prior[i])
                sapply(1:nrow(object$fit$variable$Y), function(x){
                    t(object$fit$variable$Y[x,])%*%SB + const_prior
                })
            })
            
        }else{
            prediction <- sapply(index_B, function(i){
                SB <- Rinv%*%B[i,]
                const_prior <- -0.5*t(B[i,])%*%SB + log(prior[i])
                sapply(1:nrow(newdata), function(x){
                    t(newdata[x,])%*%SB + const_prior
                })
            })
            if(!is.matrix(prediction))  prediction <- matrix(prediction, nrow=nrow(newdata))
        }
        
        # Compute the posterior and classification
        
        posterior <- exp(sweep(prediction, 1, apply(prediction,1,logsumexp), "-")) # FIXME to check
        classif <- apply(posterior, 1, which.max)
        
        # rename
        names_variables <- attr(object$fit$variables$X,"dimnames")[[2]][object$classid]
        colnames(posterior) <- names_variables
        if(missing(newdata)) rownames(posterior) = rownames(object$fit$variable$Y) else rownames(posterior) = rownames(newdata)
        if(missing(newdata)) names(classif) = rownames(object$fit$variable$Y) else names(classif) = rownames(newdata)
        names(prior) <- names_variables
        
    }else{
        
        # checks
        if(!inherits(tree, "phylo")) stop("You must provide a tree object of class \"phylo\"")
        if(!is.data.frame(newdata) & !is.matrix(newdata)) stop("the \"newdata\" should be a data.frame object with column names matching predictors names, and row names matching names in the tree ")
        # prep.
        predicted_names <- rownames(newdata)
        rcov <- .resid_cov_phylo(tree, object$fit, predicted_names)
        X1 <- matrix(1,ncol=1,nrow=nrow(newdata))
        rownames(X1) <- predicted_names
        
        # modify the design matrix FIXME, is it necessary?
        if(any(object$fit$dims$assign==0)){
            grp <- as.factor(object$fit$variables$X[,object$classid]%*%(1:object$nclass))
            Xn <- model.matrix(~grp+0)
            resid <- object$fit$variable$Y - Xn%*%B
        }else{
            resid <- residuals(object$fit)
        }
        
        # we compute the residuals across all subjects
        # the "bias" term can be computed before hand?
        bias <- rcov$w%*%solve(rcov$Vt)%*%resid[rcov$train,,drop=FALSE]
        
        # compute prediction scores
        # prediction <- sapply(1:nrow(newdata), function(x){
        prediction <- sapply(predicted_names, function(x){
            sapply(1:nrow(B), function(i){
                predicted <- X1%*%B[i,] + bias
                -0.5*( t(as.numeric(newdata[x,] - predicted[x,]))%*%Rinv%*%as.numeric(newdata[x,] - predicted[x,])) + log(prior[i]) #eg eq 26 in Hastie et al. 1994 - PDA
            })
        })
        
        # Compute the posterior and classification
        posterior <- t(exp(sweep(prediction, 2, apply(prediction, 2, logsumexp), "-"))) # FIXME to check
        classif <- apply(posterior, 1, which.max)
        
        # rename
        #colnames(posterior) <- attr(object$fit$variables$X,"dimnames")[[2]][object$classid]
        rownames(posterior) <- predicted_names
        # names(prior) <- attr(object$fit$variables$X,"dimnames")[[2]][object$classid]
        names(classif) <- predicted_names
    }
    
    
    # results
    if(missing(newdata)){
        confusion <- table(classif, object$fit$variables$X[,object$classid]%*%(1:object$nclass) ) # FIXME provide the categorical variable instead?
        rownames(confusion) = colnames(confusion) = names_variables
        results <- list(class=classif, posterior=posterior, prior=prior, confusion=confusion, training=TRUE)
    }else{
        results <- list(class=classif, posterior=posterior, prior=prior, training=FALSE)
    }
    class(results) <- "mvgls.dfa.predict"
    return(results)
}

# ------------------------------------------------------------------------- #
# print.mvgls.dfa.predict                                                   #
# options: x, ...                                                           #
#                                                                           #
# ------------------------------------------------------------------------- #
print.mvgls.dfa.predict <- function(x, ...){
    cat("\n")
    message("-- DFA predictions --","\n")
    cat("Prior probabilities of groups:","\n")
    cat(signif(x$prior, digits=4),"\n")
    cat("\n")
    if(x$training){
        error = signif(sum(diag(x$confusion))/length(x$class), digits = 3)
        cat("Training misclassification error rate:",1-error,"\n")
    }
    cat("\n")
}

