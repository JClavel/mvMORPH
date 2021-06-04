# ------------------------------------------------------------------------- #
# manova.gls                                                                #
# options: object, test, type, permutations, L...                           #
#                                                                           #
# ------------------------------------------------------------------------- #

manova.gls <- function(object, test=c("Pillai", "Wilks", "Hotelling-Lawley", "Roy"), type=c("I","II","III"), nperm=1000L, L=NULL, ...){
  
  # options
  type <- match.arg(type)[1]
  test <- match.arg(test)[1]
  args <- list(...)
  if(is.null(args[["nbcores"]])) nbcores <- 1L else nbcores <- args$nbcores
  if(is.null(args[["parametric"]])) param <- TRUE else param <- args$parametric
  if(is.null(args[["permutation"]])) penalized <- NULL else penalized <- args$permutation
  if(is.null(args[["rhs"]])) rhs <- NULL else rhs <- args$rhs
  if(is.null(args[["verbose"]])) verbose <- FALSE else verbose <- args$verbose
  if(is.null(args[["P"]])) P <- NULL else P <- args$P
  
  # Performs the tests
  if(!inherits(object, "mvgls")) stop("Please provide an object of class \"mvgls\", see ?mvgls ")
    
    # TEMPORARY?
    if(object$penalty!="LL" & object$penalty!="RidgeArch") stop("sorry, currently only the ML method or the \"RidgeArch\" penalized method is allowed")
    if(!is.null(L)){
        if(!is.null(P)) type <- "glhrm" else type <- "glh"
        if(!is.matrix(L)) warning("\n","The supplied contrasts vector L has been formatted to a matrix")
        L <- matrix(L, ncol=nrow(object$coefficients))
    } 
    # check if P is provided but not L?
    if(!is.null(P) & is.null(L)){
      type <- "glhrm"
      warning("\n","No contrasts vector L supplied for the Hypothesis. Assumes comparison to the intercept")
      if(!any(object$dims$assign==0)) stop("You're model doesn't includes an intercept term. Please provide a contrast matrix for L")
      intercept_X <- rep(0, length(object$dims$assign))
      intercept_X[object$dims$assign==0] = 1
      L <- matrix(intercept_X, ncol=nrow(object$coefficients))
    }
    if(!is.null(P) & !is.matrix(P)){
        warning("\n","The supplied contrasts vector P has been formatted to a matrix")
        P <- matrix(P, nrow=ncol(object$coefficients))
    }
       
    # if ML we can use the parametric tests or permutations
    if(object$method=="LL" & param==TRUE){
        
      if(type=="I" | type==1) paramTest <- .aov.mvgls.I(object, test)
      if(type=="II" | type==2) paramTest <- .aov.mvgls.mar(object, test, type=type)
      if(type=="III" | type==3) paramTest <- .aov.mvgls.mar(object, test, type=type)
      if(type=="glh" | type=="glhrm") paramTest <- .linearhypothesis.gls(object, test, L, rhs=rhs, nperm=nperm, nbcores=nbcores, parametric=TRUE, penalized=FALSE, P=P)
        
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
      if(type=="glh" | type=="glhrm") permTests <- .linearhypothesis.gls(object, test, L, rhs=rhs, nperm=nperm, nbcores=nbcores, parametric=param, penalized=penalized, verbose=verbose, P=P)
         
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
  return(summary_tests)
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
  #pivots <- Q_r$pivot[1L:Q_r$rank] # retrieve full rank pivots
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
  #pivots <- Q_r$pivot[1L:Q_r$rank] # retrieve full rank pivots
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
  
  # QR decomposition # FIXME: perform the QR ahead in the mvgls call
  Q_r <- qr(X)
  
  # Hypothesis (projection matrix)
  Id  <- diag(N)
  Proj_full  <- X %*% pseudoinverse(X)
  WW  <- solve(t(Y) %*% (Id - Proj_full) %*% Y)
  
  # Number of degrees of freedom
  #pivots <- Q_r$pivot[1L:Q_r$rank]
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
                          df <- ifelse(type=="III", length(which(asgn==(k-1))), length(which(asgn==k)))
                          Stats <- .multivTests(Re(eig$values), df, nb.resid, test=test)
                          Pval<-pf(Stats[2],Stats[3],Stats[4],lower.tail=FALSE)
                          results <- c(df, Stats[1], Stats[2],
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
  #pivots <- Q_r$pivot[1L:Q_r$rank]
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
  if(is.null(args[["P"]])) P <- NULL else P <- args$P
  
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
    if(!is.null(P)) rhs <- matrix(0, ncol=ncol(P), nrow=nrow(L))
  }else if(length(rhs)==1){
    rhs <- matrix(rhs, ncol=p, nrow=nrow(L))
    if(!is.null(P))  rhs <- matrix(rhs, ncol=ncol(P), nrow=nrow(L))
  }
      
  if(penalized==TRUE) penalized <- "approx"
  # QR decomposition
  Q_r <- qr(X)
  
  # Hypothesis contrasts matrix
  Xc <- X%*%pseudoinverse(t(X)%*%X)%*%t(L)
  XCXC <- pseudoinverse(t(Xc)%*%Xc)
  
  if(!is.null(P)){
    # Hypothesis
    H <- t(L%*%B0%*%P - rhs)%*%XCXC%*%(L%*%B0%*%P - rhs) # The Hypothesis matrix under LB=rhs
    
    # Error SSCP matrix (i.e. inverse of the unscaled covariance)
    WW  <- solve(t(P)%*%(object$sigma$Pinv*ndimCov)%*%P)
    # FIXME permutations for the intercept-only model
  }else{
    # Hypothesis
    H <- t(L%*%B0 - rhs)%*%XCXC%*%(L%*%B0 - rhs) # The Hypothesis matrix under LB=rhs
  
    # Error SSCP matrix (i.e. inverse of the unscaled covariance)
    WW  <- object$sigma$P/ndimCov
  }
  
  # Compute the statistic
  HE <- H%*%WW
  eig=eigen(HE, only.values = TRUE)$values
  if(!is.null(P)) nb.df <- min(qr(L)$rank, qr(P)$rank) else nb.df <- qr(L)$rank
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
                              
                               # SSCP matrix
                               SSCP <- crossprod(residuals)
                               
                               # Error matrix - Shrinkage estimator
                               targetE <- .targetM(SSCP, target, penalty)
                               Ep <- (1-tuningNull)*SSCP + tuningNull*targetE
                               
                              if(is.null(P)){
                                  # Hypothesis SSCP matrix
                                  LB <- L%*%Bp
                                  Hp <- t(LB)%*%XCXC%*%(LB)
                               
                                  # HE matrix
                                  HE <- Hp%*%solve(Ep)
                               }else{
                                  # Hypothesis SSCP matrix under RM deisgn
                                  LBP <- L%*%Bp%*%P
                                  Hp <- t(LBP)%*%XCXC%*%(LBP)
                                  
                                  # HE matrix
                                  HE <- Hp%*%solve(t(P)%*%Ep%*%P)
                               } 

                               
                             }else if(penalized=="none"){
                               
                               # randomize the residuals of the reduced model
                               Yp <- Rz[c(sample(N-1),N), ]%*%Y # -1 because we use the contrasts; we also remove the effect it's not necessary
                               
                               # coefficients under randomized sets
                               Bp <- XtX1Xt %*% Yp
                               
                               # Error SSCP
                               Ep <- crossprod(Yp - X%*%Bp)
          
                               if(is.null(P)){
                                 # Hypothesis SSCP matrix
                                 LB <- L%*%Bp
                                 Hp <- t(LB)%*%XCXC%*%(LB)
                                 
                                 # HE matrix
                                 HE <- Hp%*%solve(Ep)
                               }else{
                                 # Hypothesis SSCP matrix under RM deisgn
                                 LBP <- L%*%Bp%*%P
                                 Hp <- t(LBP)%*%XCXC%*%(LBP)
                                 
                                 # HE matrix
                                 HE <- Hp%*%solve(t(P)%*%Ep%*%P)
                               } 
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
    if(is.null(args[["tatsuoka"]])) tatsuoka <- FALSE else tatsuoka <- args$tatsuoka
    
    # check that an object of class manova.mvgls is provided
    if(!inherits(x,c("manova.mvgls","pairs.mvgls"))) stop("effectsize() can be used only with objects of class \"manova.mvgls\" or \"pairs.mvgls\"")
    
    if(x$param){
        s = x$Df #min(vh,p)
        switch(x$test,
        "Wilks"={
            N <- x$dims$n
            if(tatsuoka){
                # Tatsuoka 1973 w^2
                mult <- abs( 1 - N*x$stat/((N - s - 1) + x$stat))
                if(adjusted) mult <- (mult - ((x$Df^2 + x$dims$p^2)/(3*N))*(1-mult))
                row_names <- paste("\U03C9^2 [",x$test,"]",sep = "")
            }else{
                # generalized eta^2 - Cramer-Nicewander 1979
                mult <- 1 - x$stat^(1/s)
                # Serlin (1982) adjustment
                if(adjusted) mult <- (1 - ((N-1)/(N - max(x$Df,x$dims$p) - 1))*(1 - mult))
                row_names <- paste("\U03C4^2 [",x$test,"]",sep = "")
            }
        },
        "Pillai"={
            mult <- x$stat/s
            row_names <- paste("\U03BE^2 [",x$test,"]",sep = "")
            if(adjusted){
                # Serlin (1982) adjustment
                N <- x$dims$n
                mult <- (1 - ((N-1)/(N - max(x$Df,x$dims$p) - 1))*(1 - mult))
            }
        },
        "Hotelling-Lawley"={
            mult <- (x$stat/s)/(1 + x$stat/s)
            row_names <- paste("\U03B6^2 [",x$test,"]",sep = "")
            if(adjusted){
                # Serlin adjustment (see Kim & Olejnik 2005, Huberty & Olejnik 2006)
                N <- x$dims$n
                mult <- (1 - ((N-1)/(N - max(x$Df,x$dims$p) - 1))*(1 -mult))
            }
        },
        "Roy"={
            mult <- x$stat/(1+x$stat)
            row_names <- paste("\U03B7^2 [",x$test,"]",sep = "")
            if(adjusted){
                # Serlin adjustment for Roy > correspond to H&L with s=1
                N <- x$dims$n
                mult <- (1 - ((N-1)/(N - max(x$Df,x$dims$p) - 1))*(1 -mult))
            }
        })
        mult <- matrix(mult,nrow=1)
        if(x$type=="glh") colnames(mult) = "contrast" else colnames(mult) = x$terms
        rownames(mult) = row_names # paste("A(",x$test,")",sep = "")
    }else{
        
        # retrieve expectations and theoretical bounds
        Anull <- colMeans(x$nullstat)
        if(x$type=="III") s <- table(x$dims$assign) else s <- table(x$dims$assign)[-1L]
        
        switch(x$test,
        "Wilks"={
            
            if(normalized){ # Kramer-Nicewander 1979 > default as for the parametric case?
                mult <- (1 - (x$stat/tanh(colMeans(atanh(x$nullstat))))^(1/s))
                row_names <- paste("\U03C4^2 [",x$test,"]",sep = "")
            }else{
                mult <- (1 - x$stat/tanh(colMeans(atanh(x$nullstat))))
                row_names <- paste("\U03B7^2 [",x$test,"]",sep = "")
            }
            
        },
        "Pillai"={
            mult <- ((x$stat - Anull)/(s - Anull))
            row_names <- paste("\U03BE^2 [",x$test,"]",sep = "")
        },
        "Hotelling-Lawley"={
            rootb <- (x$stat - Anull)
            mult <- rootb/(s + rootb)
            row_names <- paste("\U03B6^2 [",x$test,"]",sep = "")
        },
        "Roy"={
            rootb <- (x$stat - Anull)
            mult <- rootb/(1 + rootb)
            row_names <- paste("\U03B7^2 [",x$test,"]",sep = "")
        })
        
        # We return the chance-corrected metrics. Should we return 0 when mult <0?
        mult <- matrix(mult, nrow=1)
        if(x$type=="glh") colnames(mult) = "contrast" else colnames(mult) = x$terms
        rownames(mult) = row_names #paste("Aperm.(",x$test,")",sep = "")
    }
    if(x$test=="Wilks") attr(adjusted, "tatsuoka") <- tatsuoka else attr(adjusted, "tatsuoka") <- FALSE
    results = list(effect=mult, adjusted=adjusted, parametric=x$param)
    class(results) = "effects.mvgls"
    return(results)
}

# ------------------------------------------------------------------------- #
# pairs.contrasts                                                           #
# options: object, term, ...                                                #
#                                                                           #
# ------------------------------------------------------------------------- #
pairs.contrasts <- function(object, term=1, ...){
    if(object$contrasts[term]!="contr.treatment") stop("object fit must use dummy coding - see ?contr.treatment")
    names_variables <- paste(names(object$xlevels[term]), object$xlevels[[term]], sep="")
    indice_predictor <- attr(object$variables$X,"dimnames")[[2]]%in%c("(Intercept)",names_variables)
    names_pred <- attr(object$variables$X,"dimnames")[[2]][indice_predictor]
    combinations <- combn(unique(names_pred),2)
    combinations_names <- combn(unique(object$xlevels[[term]]),2)
    nb_comb <- ncol(combinations)
    
    # prepare the contrasts matrix
    matrices_residuals <- matrix(0, nrow=nb_comb, ncol=object$dims$m)
    colnames(matrices_residuals) <- attr(object$variables$X,"dimnames")[[2]]
    names_contrasts <- vector()
    for(i in 1:nb_comb) names_contrasts[i] <- paste(combinations_names[1,i], combinations_names[2,i], sep=" - ")
    rownames(matrices_residuals) <- names_contrasts
    for(i in 1:nb_comb){
        if(combinations[1,i]=="(Intercept)"){
            matrices_residuals[i,combinations[2,i]] <- 1
        }else{
            matrices_residuals[i,combinations[1,i]] <- 1
            matrices_residuals[i,combinations[2,i]] <- -1
        }
    }
    return(matrices_residuals)
}


# ------------------------------------------------------------------------- #
# pairwise.glh                                                              #
# options: object, term, test, adjust, nperm, ...                           #
#                                                                           #
# ------------------------------------------------------------------------- #
pairwise.glh <- function(object, term=1, test=c("Pillai", "Wilks", "Hotelling-Lawley", "Roy"), adjust="holm", nperm=1000L, ...){
    
    # options
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
    
    # build a contrast matrix
    L <- pairs.contrasts(object, term=term)
    nb_contrasts <- nrow(L)
    
    # if ML we can use the parametric tests or permutations
    if(object$method=="LL" & param==TRUE){
        
        permTests <- sapply(1:nb_contrasts, function(contx){
            .linearhypothesis.gls(object, test, L=L[contx,,drop=FALSE], rhs=rhs, nperm=nperm, nbcores=nbcores, parametric=TRUE, penalized=FALSE)
        })
        
        terms <- attr(terms(object$formula),"term.labels")
        
        summary_tests <- list(test=test, stat=permTests[2,], approxF=permTests[3,],
        Df=permTests[1,], NumDf=permTests[4,], DenDf=permTests[5,], pvalue=permTests[6,],  param=param, terms=terms,
        dims=object$dims, adjust=p.adjust(permTests[6,], method = adjust), L=L)
        
    }else{
        param = FALSE # we use permutation rather than parametric test
        
        permTests <- lapply(1:nb_contrasts, function(contx){
            .linearhypothesis.gls(object, test, L=L[contx,,drop=FALSE], rhs=rhs, nperm=nperm, nbcores=nbcores, parametric=param, penalized=TRUE, verbose=verbose)
        })
        
        # compute the p-values for the tests
        if(test=="Wilks"){ # Wilks's lambda is an inverse test (reject H0 for small values of lambda)
            p_val <- sapply(1:nb_contrasts, function(i) sum(permTests[[i]]$observed>=c(permTests[[i]]$simulated,permTests[[i]]$observed))/(nperm+1) )
        }else{
            p_val <- sapply(1:nb_contrasts, function(i) sum(permTests[[i]]$observed<=c(permTests[[i]]$simulated,permTests[[i]]$observed))/(nperm+1) )
        }
        
        # terms labels
        terms <- attr(terms(object$formula),"term.labels")
        
        # statistic
        stats <- sapply(1:nb_contrasts, function(i) permTests[[i]]$observed)
        
        # retrieve the statistic
        summary_tests <- list(test=test, stat=stats, pvalue=p_val, param=param, terms=terms, nperm=nperm, nullstat=permTests, dims=object$dims,
        adjust=p.adjust(p_val, method = adjust), L=L)
    }
    
    # retrieve results
    class(summary_tests) <- "pairs.mvgls"
    return(summary_tests)
}

# ------------------------------------------------------------------------- #
# print option for MANOVA tests  (output borrowed from "car" package)       #
# options: x, digits, ...                                                   #
#                                                                           #
# ------------------------------------------------------------------------- #


print.pairs.mvgls <- function(x, digits = max(3L, getOption("digits") - 3L), ...){
    
    # select the appropriate output
    if(x$param){
        cat("General Linear Hypothesis Test:",x$test,"test statistic","\n")
        signif <- sapply(x$adjust, function(i) if(i<0.001){"***"}else if(i<0.01){
        "**"}else if(i<0.05){"*"}else if(i<0.1){"."}else{""})
        
        table_results <- data.frame(Df=x$Df, stat=x$stat, approxF=x$approxF, numDf=x$NumDf, denDf=x$DenDf, pval=x$pvalue, p.adj=x$adjust, signif=signif)
        if(is.null(rownames(x$L))) rownames(table_results) <- "Contrasts L" else rownames(table_results) <- rownames(x$L)
        colnames(table_results) <- c("Df", "test stat", "approx F", "num Df", "den Df", "Pr(>F)","adjusted", "")
        print(table_results, digits = digits, ...)
        cat("---","\n")
        cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1","\n")
        
        
        
    }else{ # permutation methods
        
        cat("General Linear Hypothesis Test with",x$nperm,"permutations:",x$test,"test statistic","\n")
        signif <- sapply(x$adjust, function(i) if(i<0.001){"***"}else if(i<0.01){
        "**"}else if(i<0.05){"*"}else if(i<0.1){"."}else{""})
        
        table_results <- data.frame(stat=x$stat, pval=x$pvalue, p.adj=x$adjust, signif=signif)
        if(is.null(rownames(x$L))) rownames(table_results) <- "Contrasts L" else rownames(table_results) <- rownames(x$L)
        colnames(table_results) <- c("Test stat", "Pr(>Stat)", "adjusted", "")
        print(table_results, digits = digits, ...)
        cat("---","\n")
        cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1","\n")
        
    }
    
}
