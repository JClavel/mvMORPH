# ------------------------------------------------------------------------- #
# manova.gls                                                                #
# options: object, test, type, permutations, L...                           #
#                                                                           #
# ------------------------------------------------------------------------- #

manova.gls <- function(object, test=c("Pillai", "Wilks", "Hotelling-Lawley", "Roy"), type=c("I","II","III"), nperm=999L, L=NULL, ...){
  
  # options
  type <- match.arg(type)[1]
  test <- match.arg(test)[1]
  args <- list(...)
  if(is.null(args[["nbcores"]])) nbcores <- 1L else nbcores <- args$nbcores
  if(is.null(args[["parametric"]])) param <- TRUE else param <- args$parametric
  if(is.null(args[["permutation"]])) penalized <- NULL else penalized <- args$permutation
  if(is.null(args[["rhs"]])) rhs <- NULL else rhs <- args$rhs
  if(is.null(args[["verbose"]])) verbose <- FALSE else verbose <- args$verbose
  
  # Performs the tests
  if(!inherits(object, "mvgls")) stop("Please provide an object of class \"mvgls\", see ?mvgls ")
    
    # TEMPORARY?
    if(object$penalty!="LL" & object$penalty!="RidgeArch") stop("sorry, currently only the ML method or the \"RidgeArch\" penalized method is allowed")
    if(!is.null(L)){
        type <- "glh"
        if(!is.matrix(L)) warning("\n","The supplied contrasts vector L has been formatted to a matrix")
        L <- matrix(L, ncol=nrow(object$coefficients))
       
    } 
       
    # if ML we can use the parametric tests or permutations
    if(object$method=="LL" & param==TRUE){
        
      if(type=="I" | type==1) paramTest <- .aov.mvgls.I(object, test)
      if(type=="II" | type==2) paramTest <- .aov.mvgls.mar(object, test, type=type)
      if(type=="III" | type==3) paramTest <- .aov.mvgls.mar(object, test, type=type)
      if(type=="glh") paramTest <- .linearhypothesis.gls(object, test, L, rhs=rhs, nperm=nperm, nbcores=nbcores, parametric=TRUE, penalized=FALSE)
        
      # terms labels
      if(type=="III" | type=="3") terms <- c("(Intercept)",attr(terms(object$formula),"term.labels")) else terms <- attr(terms(object$formula),"term.labels")
      
      summary_tests <- list(test=test, type=type, stat=paramTest[2,], approxF=paramTest[3,],
                            Df=paramTest[1,], NumDf=paramTest[4,], DenDf=paramTest[5,], pvalue=paramTest[6,],  param=param, terms=terms, dims=object$dims)
    }else{
      # we use permutation methods
      if(object$method!="LL" & is.null(penalized)) penalized <- "approx" else if(object$method=="LL") penalized <- "none"

        
      param = FALSE # we use permutation rather than parametric test
      
      if(type=="I" | type==1) permTests <- .aov.mvgls.perm.I(object, test, nperm=nperm, nbcores=nbcores, penalized=penalized, verbose=verbose)
      if(type=="II" | type==2) permTests <- .aov.mvgls.perm.mar(object, test, nperm=nperm, nbcores=nbcores, type=type, penalized=penalized, verbose=verbose)
      if(type=="III" | type==3) permTests <- .aov.mvgls.perm.mar(object, test, nperm=nperm, nbcores=nbcores, type=type, penalized=penalized, verbose=verbose)
      if(type=="glh") permTests <- .linearhypothesis.gls(object, test, L, rhs=rhs, nperm=nperm, nbcores=nbcores, parametric=param, penalized=penalized, verbose=verbose)
         
      # compute the p-values for the tests
      if(test=="Wilks"){ # Wilks's lambda is an inverse test (reject H0 for small values of lambda)
        p_val <- sapply(1:ncol(permTests$simulated), function(i) sum(permTests$observed[i]>=c(permTests$simulated[,i],permTests$observed[i]))/(nperm+1) )
      }else{
        p_val <- sapply(1:ncol(permTests$simulated), function(i) sum(permTests$observed[i]<=c(permTests$simulated[,i],permTests$observed[i]))/(nperm+1) )
      }
      
      # terms labels
      if(type=="III" | type=="3") terms <- c("(Intercept)",attr(terms(object$formula),"term.labels")) else terms <- attr(terms(object$formula),"term.labels")
      
      # retrieve the statistic
      summary_tests <- list(test=test, type=type, stat=permTests$observed, pvalue=p_val, param=param, terms=terms, nperm=nperm, nullstat=permTests$simulated, dims=object$dims)
        
    }
  
  class(summary_tests) = "manova.mvgls"
  invisible(summary_tests)
}


# ------------------------------------------------------------------------- #
# .multivTests                                                              #
# options: eig, q, df.res, test="Pillai" - modified from the Stats package  #
#                                                                           #
# ------------------------------------------------------------------------- #

.multivTests <- function(eig, q=NULL, df.res=NULL, test="Pillai", stat=FALSE){
  
  switch(test,
         "Pillai"={
           
           test <- sum(eig/(1 + eig))
           if(!stat){
             p <- length(eig)
             s <- min(p, q)
             n <- 0.5 * (df.res - p - 1)
             m <- 0.5 * (abs(p - q) - 1)
             tmp1 <- 2 * m + s + 1
             tmp2 <- 2 * n + s + 1
             stats <- c(test, (tmp2/tmp1 * test)/(s - test), s * tmp1, s *
                          tmp2)
           }else{
             stats = test 
           }
         },
         "Wilks"={
           test <- prod(1/(1 + eig))
           if(!stat){
             p <- length(eig)
             tmp1 <- df.res - 0.5 * (p - q + 1)
             tmp2 <- (p * q - 2)/4
             tmp3 <- p^2 + q^2 - 5
             tmp3 <- if (tmp3 > 0)
               sqrt(((p * q)^2 - 4)/tmp3)
             else 1
             stats <- c(test, ((test^(-1/tmp3) - 1) * (tmp1 * tmp3 - 2 * tmp2))/p/q,
                        p * q, tmp1 * tmp3 - 2 * tmp2)
           }else{
             stats = test
           }
         },
         "Hotelling-Lawley"={
           test <- sum(eig)
           if(!stat){
             p <- length(eig)
             m <- 0.5 * (abs(p - q) - 1)
             n <- 0.5 * (df.res - p - 1)
             s <- min(p, q)
             tmp1 <- 2 * m + s + 1
             tmp2 <- 2 * (s * n + 1)
             stats <- c(test, (tmp2 * test)/s/s/tmp1, s * tmp1, tmp2)
           }else{
             stats=test
           }
         },
         "Roy"={
           p <- length(eig)
           test <- max(eig)
           if(!stat){
             tmp1 <- max(p, q)
             tmp2 <- df.res - tmp1 + q
             stats <- c(test, (tmp2 * test)/tmp1, tmp1, tmp2)
           }else{
             stats = test
           }
         })
  
  return(stats)
}
                        
                        
# ------------------------------------------------------------------------- #
# aov.mvgls.I                                                               #
# options: object, test                                                     #
#                                                                           #
# ------------------------------------------------------------------------- #

.aov.mvgls.I <- function(object, test){
  
  # variables and parameters
  Y <- object$corrSt$Y
  X <- object$corrSt$X
  N = nrow(Y)
  
  # QR decomposition of the Hypothesis matrix
  Q_r <- qr(X)
  Q <- qr.Q(Q_r)
  
  # Hypothesis (projection matrix)
  Pf  <- X %*% pseudoinverse(X)
  Id  <- diag(N)
  WW  <- solve(t(Y) %*% (Id - Pf) %*% Y)
  
  # Number of degrees of freedom
  nb.resid <- N - Q_r$rank
  asgn <- object$dims$assign
  nterms <- sum(unique(asgn)!=0)
  
  # Tests on each variables (type I SS type)
  tests_stats <- sapply(1:nterms,
                        FUN = function(i) {
                          # Projection matrix.
                          variables=which(asgn==i)
                          QQl <- Q[, variables] %*% t(Q[, variables])
                          S <- t(Y) %*% QQl %*% Y
                          # Compute the test statistic. 
                          HE=S%*%WW
                          eig=eigen(HE, only.values = TRUE)
                          Stats <- .multivTests(Re(eig$values), length(variables), nb.resid, test=test)
                          Pval<-pf(Stats[2],Stats[3],Stats[4],lower.tail=FALSE)
                          results <- c(length(variables), Stats[1], Stats[2],
                                       Stats[3], Stats[4], Pval)
                          results
                        })
  
  return(tests_stats)
  
}


# ------------------------------------------------------------------------- #
# aov.mvgls.perm.I                                                          #
# options: object, test, nperm, nbcores, type, penalized,...                #
#                                                                           #
# ------------------------------------------------------------------------- #

.aov.mvgls.perm.I <- function(object, test, nperm=100, nbcores=1L, type="I", penalized=TRUE, ...){
  
  # options
  args <- list(...)
  if(is.null(args[["verbose"]])) verbose <- TRUE else verbose <- args$verbose
      
  # variables and parameters
  Y <- object$corrSt$Y
  X <- object$corrSt$X
  B0 <- object$coefficients
  N = nrow(Y)
  p = object$dims$p
  if(object$REML) ndimCov = object$dims$n - object$dims$m else ndimCov = object$dims$n
  tuning <- object$tuning
  target <- object$target
  penalty <- object$penalty
  
  
  if(penalized=="full"){
    Dsqrt <- pruning(object$corrSt$phy, trans=FALSE, inv=FALSE)$sqrtM
    modelPerm <- object$call
    modelPerm$grid.search <- quote(FALSE)
    modelPerm$start <- quote(object$opt$par)
  } 
  
  if(penalty=="RidgeArch") upPerm <-  1 else upPerm <- Inf 
  
  # QR decomposition of the Hypothesis matrix
  Q_r <- qr(X)
  Q <- qr.Q(Q_r)
  
  # Hypothesis (projection matrix)
  Pf  <- X %*% pseudoinverse(X)
  In  <- diag(N)
  
  # Error SSCP matrix (i.e. inverse of the unscaled covariance)
  WW  <- object$sigma$P/ndimCov
  
  # Number of degrees of freedom
  asgn <- object$dims$assign
  uasgn <- unique(asgn)
  nterms <- sum(uasgn!=0)
  nb.resid <- N - Q_r$rank
  
  # Tests on each variables (type I SS type)
  tests_stats <- sapply(1:nterms,
                        FUN = function(i) {
                          # Projection matrix.
                          variables=which(asgn==i)
                          QQl <- Q[, variables] %*% t(Q[, variables])
                          S <- t(Y) %*% QQl %*% Y
                          # Compute the test statistic. 
                          HE=S%*%WW
                          eig=eigen(HE, only.values = TRUE)
                          Stats <- .multivTests(Re(eig$values), length(variables), nb.resid, test=test)
                          Stats[1]
                        })
  
  # Loop across the terms of the model
  permutations <- sapply(1:nterms,
                         FUN = function(k){
                           
                           # Indexing of the design matrix
                           var_reduc=which(asgn<=(k-1))
                           var_full=which(asgn<=k)
                           
                           # reduced model (test for intercept) = the null matrix is zero we are comparing the mean to zero
                           if(!any(uasgn==0) & k==1) Proj_reduc <- matrix(0, nrow=N, ncol=N) else Proj_reduc <- X[,var_reduc] %*% pseudoinverse(X[,var_reduc, drop=FALSE])
                           Proj_full <-  X[,var_full] %*% pseudoinverse(X[,var_full, drop=FALSE])
                           
                           # compute the residuals under the reduced model
                           if(!any(uasgn==0) & k==1) resnull <- (In - Proj_full)%*%Y  else resnull <- (In - Proj_reduc)%*%Y
                           
                           # Fitted values under the reduced model
                           if(penalized=="full") MeanNull <- object$variables$X[,var_reduc] %*% pseudoinverse(X[,var_reduc, drop=FALSE]) %*% Y else MeanNull <- Proj_reduc%*%Y
                           #MeanNull <- object$variables$X[,var_reduc] %*% pseudoinverse(X[,var_reduc, drop=FALSE]) %*% Y
                           
                           # Permutations
                           Null <- .parallel_mapply(function(i){
                             
                             if(penalized=="approx"){ 
                               
                               # randomize the residuals of the reduced model
                               Yp <- MeanNull + resnull[c(sample(N-1),N),]
                               
                               # residuals with permuted data for the error matrix
                               XB <- Pf %*% Yp 
                               residuals <- (Yp - XB)
                               
                               # Preconditioning
                               ZZD   <- Evalues <- matrix(nrow=p, ncol=(N-1))
                               for(j in 1:(N-1)){
                                 Bi <- pseudoinverse(X[-j,,drop=FALSE])%*%Yp[-j,,drop=FALSE]
                                 resid_j <- Yp-X%*%Bi
                                 S_j <- crossprod(resid_j[-j,])/(ndimCov-1)
                                 # Rotate the results
                                 eig <- eigen(S_j)
                                 Pj <- eig$vectors
                                 Zi <- (residuals[j,]%*%Pj)^2
                                 Evalues[,j]  <- eig$values
                                 ZZD[,j] <- Zi
                               }
                               
                               # target (general but not optimized)
                               targetb <- diag(.targetM(crossprod(residuals)/ ndimCov, target, penalty=penalty))
                               #targetb <- mean(diag(crossprod(residuals)/ ndimCov ))
                               
                               ## Estim the regularization parameter (we reuse the "residuals" vector created here)
                               estimModelNull <- optim(par = tuning, fn = .loocvVect, gr=.derivllik, method="L-BFGS-B",
                                                       upper=upPerm, lower=1e-8, Evalues=Evalues, ZZD=ZZD, target=targetb, const=ndimCov/(N-1), p=p)
                               
                               # param if(penalty=="RidgeArch")
                               tuningNull <- estimModelNull$par[1] 
                               
                               # Hypothesis SSCP matrix
                               Hp <- crossprod(Yp, (Proj_full - Proj_reduc) %*% Yp)
                               
                               # SSCP matrix
                               SSCP <- crossprod(residuals)
                                 
                               # Error matrix - Shrinkage estimator
                               targetE <- .targetM(SSCP, target, penalty)
                               Ep <- (1-tuningNull)*SSCP + tuningNull*targetE
                               
                               # HE matrix
                               HE <- Hp%*%solve(Ep)
                               
                             }else if(penalized=="full"){
                               
                               # randomize the residuals of the reduced model and transform it with the appropriate structure
                               Yp <- MeanNull + Dsqrt%*%(resnull[c(sample(N-1),N),]) 
                               rownames(Yp) <- rownames(object$variables$Y)
                               
                               # retrieve the model and refit with successive eval... time consuming but it's the more general and convenient way...
                               modelPerm$response <- quote(Yp)
                               estimModelNull <- eval(modelPerm)
                               residuals <- .mvGLS(estimModelNull$corrSt)$residuals
                               
                               # Hypothesis SSCP
                               # we have to recompute the projection matrices with new X matrix
                               if(!any(uasgn==0) & k==1) Proj_reduc <- matrix(0, nrow=N, ncol=N) else Proj_reduc <- estimModelNull$corrSt$X[,var_reduc] %*% pseudoinverse(estimModelNull$corrSt$X[,var_reduc, drop=FALSE])
                               Proj_full <-  estimModelNull$corrSt$X[,var_full] %*% pseudoinverse(estimModelNull$corrSt$X[,var_full, drop=FALSE])
                               
                               Hp <- crossprod(estimModelNull$corrSt$Y, (Proj_full - Proj_reduc) %*% estimModelNull$corrSt$Y)
                               
                               # Error SSCP matrix
                               SSCP <- crossprod(residuals)
                               # Error matrix - Shrinkage estimator
                               targetE <- .targetM(SSCP, target, penalty)
                               
                               switch(penalty,
                                      "RidgeAlt"={Ep <- .makePenaltyQuad(SSCP,estimModelNull$tuning,targetE,target)$P},
                                      "RidgeArch"={Ep <- solve((1-estimModelNull$tuning)*SSCP + estimModelNull$tuning*targetE)},
                                      "LASSO"={ Ep <- glassoFast(SSCP,tuning)$wi})
                               
                               # HE matrix
                               HE <- Hp%*%Ep
                               
                             }else if(penalized=="none"){
                               
                               # randomize the residuals of the reduced model
                               Yp <- MeanNull + resnull[sample(N),] # -1 because we use the constrasts
                               
                               # Hypothesis SSCP
                               Hp <- crossprod(Yp, (Proj_full - Proj_reduc) %*% Yp)
                               
                               # Error SSCP
                               Ep <- crossprod(Yp, (In - Pf)%*%Yp)
                               # HE matrix
                               HE <- Hp%*%solve(Ep)
                             }
                             
                             # compute the statistic
                             eig <- eigen(HE, only.values = TRUE)$values
                             .multivTests(Re(eig), test=test, stat=TRUE)[1]
                             
                           }, 1:nperm, mc.cores = getOption("mc.cores", nbcores), verbose=verbose) # END loop across permutations
                         }) # END loop across terms
  
  
  # Return the simulations
  results <- list(observed=tests_stats, simulated=permutations)
  
  return(results)
}



# ------------------------------------------------------------------------- #
# aov.mvgls.mar | Type II & III SSCP                                        #
# options: object, test, type                                               #
#                                                                           #
# ------------------------------------------------------------------------- #

.aov.mvgls.mar <- function(object, test, type="III"){
  
  # variables and parameters
  Y <- object$corrSt$Y
  X <- object$corrSt$X
  N = nrow(Y)
  
  # QR decomposition
  Q_r <- qr(X)
  
  # Hypothesis (projection matrix)
  Id  <- diag(N)
  Proj_full  <- X %*% pseudoinverse(X)
  WW  <- solve(t(Y) %*% (Id - Proj_full) %*% Y)
  
  # Number of degrees of freedom
  asgn <- object$dims$assign
  if(type=="II"){
    asgn <- asgn[asgn!=0]
    intercept = TRUE # (TODO: handle cases where a model without intercept is fitted => type III)
    model_terms <- terms(object$formula)
    facTerms <- crossprod(attr(model_terms, "factors"))
    facTerms <- facTerms[,asgn,drop=FALSE]
  }else{
    intercept = NULL  
  } 
  nb.resid <- N - Q_r$rank
  nterms <- length(unique(asgn))
  
  # Tests on each variables (type III SS type)
  tests_stats <- sapply(1:nterms,
                        FUN = function(k) {
                            
                          if(type=="III"){
                                # Indexing of the design matrix
                                var_reduc=which(asgn!=(k-1))
                                
                            }else{
                                # Indexing of the design matrix
                                var_reduc <- var_full <- facTerms[k,]<unique(facTerms[k,which(asgn==k)])
                                var_full[which(asgn==k)] <- TRUE
                                # full model
                                Proj_full <- X[,c(intercept,var_full)] %*% pseudoinverse(X[,c(intercept,var_full), drop=FALSE])
                            }
                            
                          # reduced model
                          Proj_reduc <- X[,c(intercept,var_reduc)] %*% pseudoinverse(X[,c(intercept,var_reduc), drop=FALSE])
                          # Hypothesis SSCP matrix
                          S <- crossprod(Y, (Proj_full - Proj_reduc) %*% Y)
                          # Compute the test statistic. 
                          HE=S%*%WW
                          eig=eigen(HE, only.values = TRUE)
                          Stats <- .multivTests(Re(eig$values), length(which(asgn==asgn[k])), nb.resid, test=test)
                          Pval<-pf(Stats[2],Stats[3],Stats[4],lower.tail=FALSE)
                          results <- c(length(which(asgn==asgn[k])), Stats[1], Stats[2],
                                       Stats[3], Stats[4], Pval)
                          results
                        })
  
  return(tests_stats)
}

                        
# ------------------------------------------------------------------------- #
# aov.mvgls.perm.mar | type II & III SSCP                                   #
# options: object, test, nperm, nbcores, type, penalized, ...               #
#                                                                           #
# ------------------------------------------------------------------------- #
                        
.aov.mvgls.perm.mar <- function(object, test, nperm=100, nbcores=1L, type="III", penalized=TRUE, ...){
  
    
  # options
  args <- list(...)
  if(is.null(args[["verbose"]])) verbose <- TRUE else verbose <- args$verbose
      
  # variables and parameters
  Y <- object$corrSt$Y
  X <- object$corrSt$X
  B0 <- object$coefficients
  N = nrow(Y)
  p = object$dims$p
  if(object$REML) ndimCov = object$dims$n - object$dims$m else ndimCov = object$dims$n
  tuning <- object$tuning
  target <- object$target
  penalty <- object$penalty
  
  
  if(penalized=="full"){
    Dsqrt <- pruning(object$corrSt$phy, trans=FALSE, inv=FALSE)$sqrtM
    modelPerm <- object$call
    modelPerm$grid.search <- quote(FALSE)
    modelPerm$start <- quote(object$opt$par)
  } 
  
  if(penalty=="RidgeArch") upPerm <-  1 else upPerm <- Inf 
  
  # QR decomposition
  Q_r <- qr(X)
  
  # Hypothesis (projection matrix)
  Proj_full  <- X %*% pseudoinverse(X)
  In  <- diag(N)
  
  # Error SSCP matrix (i.e. inverse of the unscaled covariance)
  WW  <- object$sigma$P/ndimCov
  
  # Number of degrees of freedom
  asgn <- object$dims$assign
  if(type=="II"){
    asgn <- asgn[asgn!=0]
    intercept = TRUE # we removed the intercept (TODO: handle cases where a model without intercept is fitted => type III)
    model_terms <- terms(object$formula)
    facTerms <- crossprod(attr(model_terms, "factors"))
    facTerms <- facTerms[,asgn,drop=FALSE]
  }else{
    intercept = NULL  
  } 
  uasgn <- unique(asgn)
  nterms <- length(uasgn)
  nb.resid <- N - Q_r$rank
  
  # Tests on each terms (type II & III SS)
  tests_stats <- sapply(1:nterms,
                        FUN = function(k) {
                        
                            if(type=="III"){
                                # Indexing of the design matrix
                                var_reduc=which(asgn!=(k-1))
                                
                            }else{
                                # Indexing of the design matrix
                                var_reduc <- var_full <- facTerms[k,]<unique(facTerms[k,which(asgn==k)])
                                var_full[which(asgn==k)] <- TRUE
                                # full model
                                Proj_full <- X[,c(intercept,var_full)] %*% pseudoinverse(X[,c(intercept,var_full), drop=FALSE])
                            }
                          
                          # reduced model
                          Proj_reduc <- X[,c(intercept,var_reduc)] %*% pseudoinverse(X[,c(intercept,var_reduc), drop=FALSE])
                          # Hypothesis SSCP matrix
                          S <- crossprod(Y, (Proj_full - Proj_reduc) %*% Y)
                          # Compute the test statistic. 
                          HE=S%*%WW
                          eig=eigen(HE, only.values = TRUE)$values
                          Stats <- .multivTests(Re(eig), test=test, stat=TRUE)
                          Stats[1]
                        })
  
  # Loop across the terms of the model
  permutations <- sapply(1:nterms,
                         FUN = function(k){
                           
                            if(type=="III"){
                                # Indexing of the design matrix
                                var_reduc=which(asgn!=(k-1))
                                
                            }else{
                                # Indexing of the design matrix
                                var_reduc <- var_full <- facTerms[k,]<unique(facTerms[k,which(asgn==k)])
                                var_full[which(asgn==k)] <- TRUE
                                
                                # full model
                                Proj_full <- X[,c(intercept,var_full)] %*% pseudoinverse(X[,c(intercept,var_full), drop=FALSE])
                            }
                           
                           # reduced model
                           Proj_reduc <- X[,c(intercept,var_reduc)] %*% pseudoinverse(X[,c(intercept,var_reduc), drop=FALSE])
                           
                           # compute the residuals under the reduced model
                           resnull <- (In - Proj_reduc)%*%Y
                           
                           # Fitted values under the reduced model (probably not useful)
                           if(penalized=="full") MeanNull <- object$variables$X[,c(intercept,var_reduc)] %*% pseudoinverse(X[,c(intercept,var_reduc), drop=FALSE]) %*% Y else MeanNull <- Proj_reduc%*%Y
                           
                           # Permutations
                           Null <- .parallel_mapply(function(i){
                             
                             if(penalized=="approx"){ 
                               
                               # randomize the residuals of the reduced model
                               Yp <- MeanNull + resnull[c(sample(N-1),N),]
                               
                               # residuals with permuted data for the error matrix
                               XB <- Proj_full %*% Yp 
                               residuals <- (Yp - XB)
                               
                               # Preconditioning
                               ZZD   <- Evalues <- matrix(nrow=p, ncol=(N-1))
                               for(j in 1:(N-1)){
                                 Bi <- pseudoinverse(X[-j,,drop=FALSE])%*%Yp[-j,,drop=FALSE]
                                 resid_j <- Yp-X%*%Bi
                                 S_j <- crossprod(resid_j[-j,])/(ndimCov-1)
                                 # Rotate the results
                                 eig <- eigen(S_j)
                                 Pj <- eig$vectors
                                 Zi <- (residuals[j,]%*%Pj)^2
                                 Evalues[,j]  <- eig$values
                                 ZZD[,j] <- Zi
                               }
                               
                               # target (general but not optimized)
                               targetb <- diag(.targetM(crossprod(residuals)/ ndimCov, target, penalty=penalty))
                               #targetb <- mean(diag(crossprod(residuals)/ ndimCov ))
                               
                               ## Estim the regularization parameter (we reuse the "residuals" vector created here)
                               estimModelNull <- optim(par = tuning, fn = .loocvVect, gr=.derivllik, method="L-BFGS-B",
                                                       upper=upPerm, lower=1e-8, Evalues=Evalues, ZZD=ZZD, target=targetb, const=ndimCov/(N-1), p=p)
                               
                               # param if(penalty=="RidgeArch")
                               tuningNull <- estimModelNull$par[1] 
                               
                               # Hypothesis SSCP matrix
                               Hp <- crossprod(Yp, (Proj_full - Proj_reduc) %*% Yp)
                               
                               # SSCP matrix
                               SSCP <- crossprod(residuals)
                               
                               # Error matrix - Shrinkage estimator
                               targetE <- .targetM(SSCP, target, penalty)
                               Ep <- (1-tuningNull)*SSCP + tuningNull*targetE
                               
                               # HE matrix
                               HE <- Hp%*%solve(Ep)
                               
                             }else if(penalized=="full"){
                               
                               # randomize the residuals of the reduced model and transform it with the appropriate structure
                               Yp <- MeanNull + Dsqrt%*%(resnull[c(sample(N-1),N),]) 
                               rownames(Yp) <- rownames(object$variables$Y)
                               
                               # retrieve the model and refit with successive eval... time consuming but it's the more general and convenient way...
                               modelPerm$response <- quote(Yp)
                               estimModelNull <- eval(modelPerm)
                               residuals <- .mvGLS(estimModelNull$corrSt)$residuals
                               
                               # Hypothesis SSCP
                               # we have to recompute the projection matrices with new X matrix
                               Proj_reduc <- estimModelNull$corrSt$X[,c(intercept,var_reduc)] %*% pseudoinverse(estimModelNull$corrSt$X[,c(intercept,var_reduc), drop=FALSE])
                               Hp <- crossprod(estimModelNull$corrSt$Y, (Proj_full - Proj_reduc) %*% estimModelNull$corrSt$Y)
                               
                               # Error SSCP matrix
                               SSCP <- crossprod(residuals)
                               # Error matrix - Shrinkage estimator
                               targetE <- .targetM(SSCP, target, penalty)
                               
                               switch(penalty,
                                      "RidgeAlt"={Ep <- .makePenaltyQuad(SSCP,estimModelNull$tuning,targetE,target)$P}, # to standardize to the SSCP for futur developments
                                      "RidgeArch"={Ep <- solve((1-estimModelNull$tuning)*SSCP + estimModelNull$tuning*targetE)},
                                      "LASSO"={ Ep <- glassoFast(SSCP,tuning)$wi})
                               
                               # HE matrix
                               HE <- Hp%*%Ep
                               
                             }else if(penalized=="none"){
                               
                               # randomize the residuals of the reduced model
                               Yp <- MeanNull + resnull[c(sample(N-1),N),] # -1 because we use the constrasts
                               
                               # Hypothesis SSCP
                               Hp <- crossprod(Yp, (Proj_full - Proj_reduc) %*% Yp)
                               
                               # Error SSCP
                               Ep <- crossprod(Yp, (In - Proj_full)%*%Yp)
                               # HE matrix
                               HE <- Hp%*%solve(Ep)
                             }
                             
                             # compute the statistic
                             eig <- eigen(HE, only.values = TRUE)$values
                             .multivTests(Re(eig), test=test, stat=TRUE)[1]
                             
                           }, 1:nperm, mc.cores = getOption("mc.cores", nbcores), verbose = verbose) # END loop across permutations
                         }) # END loop across terms
  
  
  # Return the simulations
  results <- list(observed=tests_stats, simulated=permutations)
  
  return(results)
}
                        
# ------------------------------------------------------------------------- #
# .linearhypothesis.gls                                                     #
# options: object, test, L, rhs, nperm, nbcores, parametric, penalized,...  #
#                                                                           #
# ------------------------------------------------------------------------- #
                     
.linearhypothesis.gls <- function(object, test, L, rhs=NULL, nperm=100, nbcores=1L, parametric=FALSE, penalized=TRUE, ...){
  
  # options
  args <- list(...)
  if(is.null(args[["verbose"]])) verbose <- TRUE else verbose <- args$verbose
      
  # variables and parameters
  Y <- object$corrSt$Y
  X <- object$corrSt$X
  B0 <- object$coefficients
  N = nrow(Y)
  p = object$dims$p
  if(object$REML) ndimCov = object$dims$n - object$dims$m else ndimCov = object$dims$n
  tuning <- object$tuning
  target <- object$target
  penalty <- object$penalty
  
  if(penalty=="RidgeArch") upPerm <-  1 else upPerm <- Inf 
  if(is.null(rhs)){
      rhs <- matrix(0, ncol=p, nrow=nrow(L))
  }else if(length(rhs)==1){
      rhs <- matrix(rhs, ncol=p, nrow=nrow(L))
  }
      
  if(penalized==TRUE) penalized <- "approx"
  # QR decomposition
  Q_r <- qr(X)
  
  # Hypothesis
  Xc <- X%*%pseudoinverse(t(X)%*%X)%*%t(L)
  XCXC <- pseudoinverse(t(Xc)%*%Xc)
  H <- t(L%*%B0 - rhs)%*%XCXC%*%(L%*%B0 - rhs) # The Hypothesis matrix under LB=rhs
  
  # Error SSCP matrix (i.e. inverse of the unscaled covariance)
  WW  <- object$sigma$P/ndimCov
  
  # Compute the statistic
  HE <- H%*%WW
  eig=eigen(HE, only.values = TRUE)$values
  nb.df <- qr(L)$rank
  nb.resid <- N - Q_r$rank
  Stats <- .multivTests(Re(eig), nb.df, nb.resid, test=test) 
    
  if(parametric){
      # Perform parametric test
      Pval <- pf(Stats[2], Stats[3], Stats[4], lower.tail=FALSE)
      results <- as.matrix(c(nb.df, Stats[1], Stats[2], Stats[3], Stats[4], Pval))
      
    }else{
      
      # Perform simulations for both regularized and conventional
      
      # Define the reduced model design
      Z = X - X%*%t(L)%*%t(pseudoinverse(L))
      # Define the hat matrix for the reduced model and for the full model
      Pz = Z%*%pseudoinverse(Z)
      XtX1Xt = pseudoinverse(X)
      # Construct the residuals making matrix
      Rz = (diag(N) - Pz)
  
        # Perform the permutations and show progress bar
        permutations <- .parallel_mapply(function(i){
                             
                             if(penalized=="approx"){ # use the approximation = assumes C is fixed
                               
                               # randomize the residuals of the reduced model
                               Yp <- (Rz[c(sample(N-1),N),] + Pz)%*%Y
                                
                               # residuals with permuted data for the error matrix
                               Bp <- XtX1Xt %*% Yp 
                               residuals <- (Yp - X%*%Bp)
                               
                               # Preconditioning
                               ZZD   <- Evalues <- matrix(nrow=p, ncol=(N-1))
                               for(j in 1:(N-1)){
                                 Bi <- pseudoinverse(X[-j,,drop=FALSE])%*%Yp[-j,,drop=FALSE]
                                 resid_j <- Yp-X%*%Bi
                                 S_j <- crossprod(resid_j[-j,])/(ndimCov-1)
                                 # Rotate the results
                                 eig <- eigen(S_j)
                                 Pj <- eig$vectors
                                 Zi <- (residuals[j,]%*%Pj)^2
                                 Evalues[,j]  <- eig$values
                                 ZZD[,j] <- Zi
                               }
                               
                               # target (general but not optimized)
                               targetb <- diag(.targetM(crossprod(residuals)/ ndimCov, target, penalty=penalty))
                               #targetb <- mean(diag(crossprod(residuals)/ ndimCov ))
                               
                               ## Estim the regularization parameter (we reuse the "residuals" vector created here)
                               estimModelNull <- optim(par = tuning, fn = .loocvVect, gr=.derivllik, method="L-BFGS-B",
                                                       upper=upPerm, lower=1e-8, Evalues=Evalues, ZZD=ZZD, target=targetb, const=ndimCov/(N-1), p=p)
                               
                               # param if(penalty=="RidgeArch")
                               tuningNull <- estimModelNull$par[1] 
                               
                               # Hypothesis SSCP matrix
                               Hp <- t(L%*%Bp)%*%XCXC%*%(L%*%Bp)
                               
                               # SSCP matrix
                               SSCP <- crossprod(residuals)
                               
                               # Error matrix - Shrinkage estimator
                               targetE <- .targetM(SSCP, target, penalty)
                               Ep <- (1-tuningNull)*SSCP + tuningNull*targetE
                               
                               # HE matrix
                               HE <- Hp%*%solve(Ep)
                               
                             }else if(penalized=="none"){
                               
                               # randomize the residuals of the reduced model
                               Yp <- Rz[c(sample(N-1),N), ]%*%Y # -1 because we use the constrasts; we also remove the effect it's not necessary
                               
                               # Hypothesis SSCP
                               Bp <- XtX1Xt %*% Yp
                               Hp <- t(L%*%Bp)%*%XCXC%*%(L%*%Bp)
                               
                               # Error SSCP
                               Ep <- crossprod(Yp - X%*%Bp)
                               
                               # HE matrix
                               HE <- Hp%*%solve(Ep)
                             }
                             
                             # compute the statistic
                             eig <- eigen(HE, only.values = TRUE)$values
                             .multivTests(Re(eig), test=test, stat=TRUE)[1]
                             
                           }, 1:nperm, mc.cores = getOption("mc.cores", nbcores), verbose = verbose) # END loop across permutations
                  
        # Return the simulations
        results <- list(observed=Stats[1], simulated=as.matrix(permutations))
    }
    
    
    return(results)
}                        

# ------------------------------------------------------------------------- #
# .loocvVect (vectorized log-likelihood) with the RidgeArch penalty         #
# options: alpha, Evalues, ZZD, target, nobs, p                             #
#                                                                           #
# ------------------------------------------------------------------------- #

.loocvVect <- function(alpha, Evalues, ZZD, target, const=1, p){
  res <- const*log(((1-alpha)*Evalues + alpha*target)) + (ZZD)/((1-alpha)*Evalues + alpha*target)
  LL <- 0.5 * sum(res)
  return(LL)
}

# ------------------------------------------------------------------------- #
# .derivllik (vectorized derivatives) with the RidgeArch penalty            #
# options: alpha, Evalues, ZZD, target, nobs                                #
#                                                                           #
# ------------------------------------------------------------------------- #

.derivllik <- function(alpha, Evalues, ZZD, target, const=1, p=NULL){
  res <- ((target - Evalues)*(target*alpha*const + Evalues*const - Evalues*alpha*const - ZZD))/(Evalues - Evalues*alpha + alpha*target)^2
  return(0.5*sum(res))
}

# ------------------------------------------------------------------------- #
# effectsize - multivariate measures of association                         #
# options: x, ...                                                           #
#                                                                           #
# ------------------------------------------------------------------------- #

effectsize <- function(x, ...){
    # Multivariate measures of association
    args <- list(...)
    if(is.null(args[["normalized"]])) normalized <- TRUE else normalized <- args$normalized
    if(is.null(args[["adjusted"]])) adjusted <- FALSE else adjusted <- args$adjusted
    if(x$param){
        s = x$Df #min(vh,p)
        switch(x$test,
        "Wilks"={
            # generalized eta^2 => take as normalized the conservative estimate based on the geometric mean of the canonical correlations
            if(normalized) mult <- 1 - x$stat^(1/s) else mult <- 1 - x$stat
            if(adjusted){
                # Tatsuoka adjustment
                N <- x$dims$n
                mult <- abs( 1 - N*x$stat/((N - s - 1) + x$stat))
            }
        },
        "Pillai"={
            mult <- x$stat/s
            
            if(adjusted){
                # Serlin adjustment
                N <- x$dims$n
                mult <- abs(1 - ((N-1)/(N - max(x$Df,x$dims$p) - 1))*(1 -mult))
            }
        },
        "Hotelling-Lawley"={
            mult <- (x$stat/s)/(1 + x$stat/s)
        },
        "Roy"={
            mult <- x$stat/(1+x$stat)
        })
        mult <- matrix(mult,nrow=1)
        colnames(mult) = x$terms
        rownames(mult) = paste("A(",x$test,")",sep = "")
    }else{
        # Provide a standardized size effect instead if permutation were used
        if(normalized){
            if(x$test=="Wilks"){
                # first compute SES (use a z-transform for Wilks' lambda?)
                ses <- (atanh(x$stat) - colMeans(atanh(x$nullstat))) / apply(atanh(x$nullstat), 2, sd)
            }else{
                ses <- (log(x$stat) - colMeans(log(x$nullstat))) / apply(log(x$nullstat), 2, sd)
            }
        }else{
            ses <-  (x$stat - colMeans(x$nullstat)) / apply(x$nullstat, 2, sd)
        }
        # can return ses/(1+ses) for a readable measure; for Wilks we can derive a suitable measure
        if(x$test=="Wilks") multRel <- abs(1 - x$stat/tanh(colMeans(atanh(x$nullstat)))) else  multRel <- abs(ses)/(1+abs(ses))
        mult <- rbind(abs(ses), multRel)
        colnames(mult) <- x$terms
        if(x$test=="Wilks") rownames(mult) <- c("SES","A(perm.)") else rownames(mult) <- c("SES","Rel.")
        
    }
    message("Multivariate measure(s) of association","\n")
    return(mult)
}
