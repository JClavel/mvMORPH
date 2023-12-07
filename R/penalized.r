################################################################################
##                                                                            ##
##                       mvMORPH: penalized.r                                 ##
##                                                                            ##
##   Internal functions for penalized methods in the mvMORPH package          ##
##                                                                            ##
##  Created by Julien Clavel - 31-07-2018                                     ##
##  (julien.clavel@hotmail.fr/ julien.clavel@biologie.ens.fr)                 ##
##   require: phytools, ape, corpcor, subplex, spam, glassoFast, stats        ##
##                                                                            ##
################################################################################

# ------------------------------------------------------------------------- #
# .loocvPhylo                                                               #
# options: par, cvmethod, targM, corrStr, penalty, error, nobs              #
#                                                                           #
# ------------------------------------------------------------------------- #

.loocvPhylo <- function(par, cvmethod, targM, corrStr, penalty, error, nobs){
    
    if(corrStr$REML) n <- nobs-ncol(corrStr$X) else n <- nobs
    p = corrStr$p
    # parameters
    alpha = corrStr$bounds$trTun(par)
    if(!is.null(error)) corrStr$mserr = corrStr$bounds$trSE(par)
    
    # model
    mod_par = .corrStr(corrStr$bounds$trPar(par), corrStr);
    
    # GLS estimates
    XtX <- pseudoinverse(mod_par$X)
    B <- XtX%*%mod_par$Y
    residuals <- mod_par$Y - mod_par$X%*%B
    Ccov <- mod_par$det
    
    # Covariance matrix
    Sk <- crossprod(residuals)/n
    
    # Switch between LOOCV approaches
    switch(cvmethod,
    "H&L"={
        
        # target matrix
        target <- .targetM(Sk, targM, penalty="RidgeArch")
        
        # Hoffbeck & Landgrebe (1996) efficient LOOCV
        beta <- (1 - alpha)/(n - 1)
        G <- n*beta*Sk + alpha * target
        Gi <- try(chol(G), silent=TRUE)
        if(inherits(Gi, 'try-error')) return(1e6)
        
        llik <- sapply(1:(nobs-1), function(x){
            rk <- sum(backsolve(Gi, residuals[x,], transpose = TRUE)^2)
            (n/(nobs-1))*log(1 - beta*rk) + (rk/(1 - beta*rk))
        })
        
        ll <- 0.5 * (n*p*log(2*pi) + p*Ccov +
        n*sum(2*log(diag(Gi))) + sum(llik))
    },
    "Mahalanobis"={
        
        # target matrix
        target <- .targetM(Sk, targM, penalty="RidgeArch")
        
        # Mahalanobis approximation of the LOOCV
        beta <- (1 - alpha)/(n - 1)
        G <- n*beta*Sk + alpha * target
        Gi <- try(chol(G), silent=TRUE)
        if(inherits(Gi, 'try-error')) return(1e6)
        r0 <- sum(backsolve(Gi, t(residuals), transpose = TRUE)^2)/n
        
        ll <- 0.5 * (n*p*log(2*pi) + p*Ccov + n*sum(2*log(diag(Gi))) +
        n*log(1 - beta*r0) + n*(r0/(1 - beta*r0)))
    },
    "LOOCV"={
        
        # target matrix
        target <- .targetM(Sk, targM, penalty)
        
        # hat matrix
        h <- diag(mod_par$X%*%XtX)
        
        # check for hat score of 1 (e.g. MANOVA design)
        nloo <- corrStr$nloo[!h+1e-8>=1]
        const <- n/length(nloo)
        
        llik <- sapply(nloo, function(x){
            Bx <- B - tcrossprod(XtX[,x], residuals[x,])/(1-h[x]) # rank-1 update
            # update the residuals
            residuals2 <- mod_par$Y - mod_par$X%*%Bx
            Skpartial <- crossprod(residuals2[-x,])/(n-1)
            .regularizedLik(Skpartial, residuals[x,], alpha, targM, target, penalty, const) # try instead with current residual observation?
        })
        
        ll <- 0.5 * (n*p*log(2*pi) + p*Ccov + sum(llik))
        
    },
    "LL"={
        
        # Maximum Likelihood
        Gi <- try(chol(Sk), silent=TRUE)
        if(inherits(Gi, 'try-error')) return(1e6)
        detValue <- sum(2*log(diag(Gi)))
        quadprod <- sum(backsolve(Gi, t(residuals), transpose = TRUE)^2)
        ll <- 0.5 * (n*p*log(2*pi) + p*Ccov + n*detValue + quadprod)
        
    },
    stop("You must specify \"LOOCV\", \"H&L\" or \"Mahalanobis\" method for computing the LOOCV score and \"LL\" for the log-likelihood")
    )
    
    if (!is.finite(ll)) return(1e6)
    return(ll)
}

# ------------------------------------------------------------------------- #
# .mvGLS                                                                    #
# options: corrstruct object                                                #
#                                                                           #
# ------------------------------------------------------------------------- #
.mvGLS <- function(corrstruct){
    
    # GLS Estimate
    B <- pseudoinverse(corrstruct$X)%*%corrstruct$Y
    residuals <- corrstruct$Y - corrstruct$X%*%B
    
    return(list(residuals=residuals, B=B))
}

# ------------------------------------------------------------------------- #
# .scaleStuct                                                               #
# options: structure object                                                 #
#                                                                           #
# ------------------------------------------------------------------------- #
.scaleStruct <- function(structure){
    # wrapper for future developments
    if(inherits(structure, "phylo")){
        structure$edge.length <- structure$edge.length/max(node.depth.edgelength(structure))
    }
    
    return(structure)
}

# ------------------------------------------------------------------------- #
# .targetM                                                                  #
# options: S, targM, penalty, I                                             #
#                                                                           #
# ------------------------------------------------------------------------- #
.targetM <- function(S, targM, penalty="RidgeArch", I = NULL, ...){
    
    p <- dim(S)[1]
    args <- list(...)
    if(is.null(args[["userMatrix"]])) userMatrix <- NULL else userMatrix <- args$userMatrix
    # If the identity is not provided
    if(is.null(I)) I = diag(p)
    target = NULL
    
    if(penalty=="RidgeArch"){
        switch(targM,
        "Variance" = {target <- diag(diag(S))},
        "unitVariance" = {target <- I*mean(diag(S))},
        "null" = {
            warning("The \"null\" target cannot be used with the \"RidgeArch\" method. The \"unitVariance\" target is used instead.")
            target <- I*mean(diag(S))
        },
        "user" = { target <- userMatrix} # TODO
        )
    }else if(penalty=="RidgeAlt"){
        switch(targM,
        "Variance" = {target <- diag(1/diag(S))},
        "unitVariance" = {target <- I*(1/mean(diag(S)))},
        "null" = {target <- matrix(0,p,p)},
        "user" = { target <- solve(userMatrix)} # TODO
        )
    }
    
    return(target)
}

# ------------------------------------------------------------------------- #
# .corrStr (wrapper to covariance structure)                                #
# options: par, timeObject                                                  #
#                                                                           #
# ------------------------------------------------------------------------- #
.corrStr <- function(par, timeObject){
    
    if(timeObject$model%in%c("EB", "BM", "lambda", "OU", "OUvcv", "BMM", "OUM", "OU1")){
        # Tree transformation
        struct = .transformTree(timeObject$structure, par, model=timeObject$model, mserr=timeObject$mserr,
                                Y=timeObject$Y, X=timeObject$X, REML=timeObject$REML, precalc=timeObject$precalc)
    }else{
        stop("Currently works for phylogenetic models \"BM\", \"EB\", \"OU\", \"BMM\", \"OUM\" and \"lambda\"  only...")
    }
    return(struct)
}

# ------------------------------------------------------------------------- #
# .regularizedLik return the log-lik with the regularized estimate          #
# options: S, residuals, lambda, targM, target, penalty, const              #
#                                                                           #
# ------------------------------------------------------------------------- #
.regularizedLik <- function(S, residuals, lambda, targM, target, penalty, const=1){
    
    switch(penalty,
    "RidgeArch"={
        G <- (1-lambda)*S + lambda*target
        Gi <- try(chol(G), silent=TRUE)
        if(inherits(Gi, 'try-error')) return(1e6)
        rk <- sum(backsolve(Gi, residuals, transpose = TRUE)^2)
        llik <- const*sum(2*log(diag(Gi))) + rk
    },
    "RidgeAlt"={
        quad <- .makePenaltyQuad(S, lambda, target, targM)
        Gi <- quad$P
        detG <- sum(log(quad$ev))
        Swk <- tcrossprod(residuals)
        rk <- sum(Swk*Gi)
        llik <- const*detG + rk
    },
    "LASSO"={
        LASSO <- glassoFast(S, lambda, maxIt=500)
        G <- LASSO$w;
        Gi <- LASSO$wi;
        Swk <- tcrossprod(residuals);
        rk <- sum(Swk*Gi);
        llik <- const*as.numeric(determinant(G)$modulus) + rk
    })
    
    return(llik)
}


# ------------------------------------------------------------------------- #
# .makePenaltyQuad   (for quadratic ridge)                                  #
# options: S,lambda,target,targM                                            #
#                                                                           #
# ------------------------------------------------------------------------- #
.makePenaltyQuad <- function(S,lambda,target,targM){
    
    switch(targM,
    "Variance"={
        D <- (S - lambda * target)
        D2 <- D %*% D
        sqrtM <- .sqM(D2/4 + lambda * diag(nrow(S)))
        Alt <- D/2 + sqrtM
        AltInv <- (1/lambda)*(Alt - D)
        evalues <- eigen(Alt, symmetric=TRUE, only.values = TRUE)$values
    },
    "unitVariance"={
        eig  <- eigen(S, symmetric = TRUE)
        Q <- eig$vectors
        d <- eig$values - lambda*target[1]
        evalues <- sqrt(lambda + d^2/4) + d/2
        D1 <- evalues
        D2 <- 1/evalues # Inverse
        Alt <- Q %*% (D1 * t(Q))
        AltInv <- Q %*% (D2 * t(Q))
    },
    "null"={
        eig  <- eigen(S, symmetric = TRUE)
        Q <- eig$vectors
        d <- eig$values
        evalues <- sqrt(lambda + d^2/4) + d/2
        D1 <- evalues
        D2 <- 1/evalues
        Alt <- Q %*% (D1 * t(Q))
        AltInv <- Q %*% (D2 * t(Q))
    }
    )
    pen <- list(S=Alt, P=AltInv, ev=evalues)
    return(pen)
}


# ------------------------------------------------------------------------- #
# .sqM                                                                      #
# Matrix square root using eigen-decomposition                              #
#                                                                           #
# ------------------------------------------------------------------------- #
.sqM <- function(x){
    if(!all(is.finite(x))) return(Inf)
    eig <- eigen(x, symmetric = TRUE)
    sqrtM <- eig$vectors %*% (sqrt(eig$values) * t(eig$vectors))
    return(sqrtM)
}

# Build the matrix square root inverse
.sqM1 <- function(x){
    if(inherits(x, "phylo")) x <- vcv.phylo(x)
    if(!all(is.finite(x))) return(Inf)
    eig <- eigen(x, symmetric = TRUE)
    # check for singular dimensions =>  hack from corpcor package. Just retain the dimensions with positive eigenvalues
    tol = max(dim(x))*max(eig$values)*.Machine$double.eps
    Positive = eig$values > tol
    if(sum(Positive)<length(eig$values)) warning("The phylogenetic covariance matrix was singular. Check the results carefully and consider using 'eigSqm=FALSE' option and 'error=TRUE'")
    sqrtM <- eig$vectors[,Positive,drop=FALSE] %*% ((1/sqrt(eig$values[Positive])) * t(eig$vectors[,Positive,drop=FALSE]))
    return(sqrtM)
}


# ------------------------------------------------------------------------- #
# .covPenalized                                                             #
# options: S, penalty, targM="null", tuning=0                               #
#                                                                           #
# ------------------------------------------------------------------------- #
.penalizedCov <- function(S, penalty, targM="null", tuning=0){
    
    # dim of S
    p = ncol(S)
    
    # target matrix
    Target <- .targetM(S, targM, penalty)
    
    # Construct the penalty term
    switch(penalty,
    "RidgeAlt"={
        pen <- .makePenaltyQuad(S,tuning,Target,targM)
        Pi <- pen$S
        P <- pen$P
    },
    "RidgeArch"={
        Pi <- (1-tuning)*S + tuning*Target
        eig <- eigen(Pi)
        V <- eig$vectors
        d <- eig$values
        P <- V%*%((1/d) * t(V))
    },
    "LASSO"={
        LASSO <- glassoFast(S,tuning)
        Pi <- LASSO$w
        P <- LASSO$wi
    },
    "LL"={
        Pi <- S
        eig <- eigen(Pi)
        V <- eig$vectors
        d <- eig$values
        P <- V%*%((1/d) * t(V))
    })
    
    estimate <- list(Pinv=Pi, P=P, S=S)
    return(estimate)
}

# ------------------------------------------------------------------------- #
# .transformTree                                                            #
# options: phy, param, model, mserr=NULL, Y=NULL, X=NULL, REML=TRUE,        #
#      precalc=NULL                                                         #
# ------------------------------------------------------------------------- #

.transformTree <- function(phy, param, model=c("EB", "BM", "lambda", "OU", "BMM", "OUM"), mserr=NULL, Y=NULL, X=NULL, REML=TRUE, precalc=NULL){
    
    # pre-compute and checks | TODO reduce computational burden by avoiding reordering and recomputing distances, ages...etc
    n <- Ntip(phy)
    parent <- phy$edge[,1]
    descendent <- phy$edge[,2]
    extern <- (descendent <= n)
    N <- 2*n-2
    diagWeight <- NULL
    const <- 0
    flag <- FALSE
    
    # Model
    switch(model,
    "OU"={
        D = numeric(n)
        
        # check first for ultrametric tree (see Ho & Ane 2014 - Systematic Biology; R code based on "phylolm" package implementation. Courtesy of L. Ho and C. Ane)
        if(!is.ultrametric(phy)){
            dis = node.depth.edgelength(phy) # has all nodes
            D = max(dis[1:n]) - dis[1:n]
            D = D - mean(D)
            phy$edge.length[extern] <- phy$edge.length[extern] + D[descendent[extern]]
            flag <- TRUE
        }
        
        # Branching times (now the tree is ultrametric)
        times <- branching.times(phy)
        Tmax <- max(times)
        # compute the branch lengths
        distRoot <-  exp(-2*param*times)*(1 - exp(-2*param*(Tmax-times)))
        d1 = distRoot[parent-n]
        d2 = numeric(N)
        d2[extern] = exp(-2*param*D[descendent[extern]]) * (1-exp(-2*param*(Tmax-D[descendent[extern]])))
        d2[!extern] = distRoot[descendent[!extern]-n]
        
        # weights for a 3 points structured matrix
        diagWeight = exp(param*D)
        phy$edge.length = (d2 - d1)/(2*param) # scale the tree for the stationary variance
        names(diagWeight) = phy$tip.label
        
        # transform the variables
        w <- 1/diagWeight
        Y <- matrix(w*Y, nrow=n)
        X <- matrix(w*X, nrow=n)
        
        # Adjust errors
        if(!is.null(mserr)) mserr = mserr*exp(-2*param*D[descendent[extern]])
    },
    "OUM"={
        # Weight matrix OUM
        W <- .Call(mvmorph_weights, nterm=as.integer(n), epochs=precalc$epochs, lambda=param, S=1, S1=1, beta=precalc$listReg, root=as.integer(precalc$root_std))
        
        # transform the tree
        D = numeric(n)
        
        # check first for ultrametric tree (see Ho & Ane 2014 - Systematic Biology; R code based on "phylolm" package implementation. Courtesy of L. Ho and C. Ane)
        if(!is.ultrametric(phy)){
            dis = node.depth.edgelength(phy) # has all nodes
            D = max(dis[1:n]) - dis[1:n]
            D = D - mean(D)
            phy$edge.length[extern] <- phy$edge.length[extern] + D[descendent[extern]]
            flag <- TRUE
        }
        
        # Branching times (now the tree is ultrametric)
        times <- branching.times(phy)
        Tmax <- max(times)
        # compute the branch lengths
        if(precalc$randomRoot){
            distRoot <-  exp(-2*param*times)
            d1 = distRoot[parent-n]
            d2 = numeric(N)
            d2[extern] = exp(-2*param*D[descendent[extern]])
            d2[!extern] = distRoot[descendent[!extern]-n]
        }else{
            distRoot <-  exp(-2*param*times)*(1 - exp(-2*param*(Tmax-times)))
            d1 = distRoot[parent-n]
            d2 = numeric(N)
            d2[extern] = exp(-2*param*D[descendent[extern]]) * (1-exp(-2*param*(Tmax-D[descendent[extern]])))
            d2[!extern] = distRoot[descendent[!extern]-n]
        }
        
        # weights for a "3 points" structured matrix
        diagWeight = exp(param*D)
        phy$edge.length = (d2 - d1)/(2*param) # scale the tree for the stationary variance
        names(diagWeight) = phy$tip.label
        
        # transform the variables
        w <- 1/diagWeight
        Y <- matrix(w*Y, nrow=n)
        X <- matrix(w*W, nrow=n) # Here X is replaced by the weighted matrix
        
        # REML "constant"
        if(REML) const <- determinant(crossprod(W))$modulus # TODO: check for n-ultrametric trees
        
        # Adjust errors
        if(!is.null(mserr)) mserr = mserr*exp(-2*param*D[descendent[extern]])
        
    },
    "OU1"={
        # Weight matrix OU1
        W <- .Call(mvmorph_weights, nterm=as.integer(n), epochs=precalc$epochs, lambda=param, S=1, S1=1, beta=precalc$listReg, root=as.integer(precalc$root_std))
        
        # transform the tree
        D = numeric(n)
        
        # check first for ultrametric tree (see Ho & Ane 2014 - Systematic Biology; R code based on "phylolm" package implementation. Courtesy of L. Ho and C. Ane)
        if(!is.ultrametric(phy)){
            dis = node.depth.edgelength(phy) # has all nodes
            D = max(dis[1:n]) - dis[1:n]
            D = D - mean(D)
            phy$edge.length[extern] <- phy$edge.length[extern] + D[descendent[extern]]
            flag <- TRUE
        }
        
        # Branching times (now the tree is ultrametric)
        times <- branching.times(phy)
        Tmax <- max(times)
        # compute the branch lengths
        if(precalc$randomRoot){
            distRoot <-  exp(-2*param*times)
            d1 = distRoot[parent-n]
            d2 = numeric(N)
            d2[extern] = exp(-2*param*D[descendent[extern]])
            d2[!extern] = distRoot[descendent[!extern]-n]
        }else{
            distRoot <-  exp(-2*param*times)*(1 - exp(-2*param*(Tmax-times)))
            d1 = distRoot[parent-n]
            d2 = numeric(N)
            d2[extern] = exp(-2*param*D[descendent[extern]]) * (1-exp(-2*param*(Tmax-D[descendent[extern]])))
            d2[!extern] = distRoot[descendent[!extern]-n]
        }
        
        # weights for a "3 points" structured matrix
        diagWeight = exp(param*D)
        phy$edge.length = (d2 - d1)/(2*param) # scale the tree for the stationary variance
        names(diagWeight) = phy$tip.label
        
        # transform the variables
        w <- 1/diagWeight
        Y <- matrix(w*Y, nrow=n)
        X <- matrix(w*W, nrow=n) # Here X is replaced by the weighted matrix
        
        # REML "constant"
        if(REML) const <- determinant(crossprod(W))$modulus
        
        # Adjust errors
        if(!is.null(mserr)) mserr = mserr*exp(-2*param*D[descendent[extern]])
        
    },
    "EB"={
        if (param!=0){
            distFromRoot <- node.depth.edgelength(phy)
            phy$edge.length = (exp(param*distFromRoot[descendent])-exp(param*distFromRoot[parent]))/param
        }
    },
    "lambda"={
        # Pagel's lambda tree transformation
        if(param!=1) {
            root2tipDist <- node.depth.edgelength(phy)[1:n] # for non-ultrametric trees. The 'up' limit should be exactly 1 to avoid singularity issues
            phy$edge.length <- phy$edge.length * param
            phy$edge.length[extern] <- phy$edge.length[extern] + (root2tipDist * (1-param))
        }
    },
    "OUvcv"={
        V<-.Call("mvmorph_covar_ou_fixed", A=vcv.phylo(phy), alpha=param, sigma=1, PACKAGE="mvMORPH")
        C<-list(sqrtM=t(chol(solve(V))), det=determinant(V)$modulus)
    },
    "OUTS"={
       stop("Not yet implemented. The time-series models are coming soon, please be patient")
    },
    "RWTS"={
       stop("Not yet implemented. The time-series models are coming soon, please be patient")
    },
    "BMM"={
        # multirates model - proportional scaling or explicit estimation?
        #phy$edge.length <- phy$mapped.edge %*% param
        phy$edge.length <- phy$mapped.edge %*% c(1,param)
    })
    
    # Add measurment error
    if(is.numeric(mserr)) phy$edge.length[extern] = phy$edge.length[extern] + mserr
    
    # Compute the independent contrasts scores
    if(inherits(phy, "phylOLS")){
        if((sum(phy$edge.length) - n)<=.Machine$double.eps){
            # Return the determinant
            deterM <- 0
        }else{
            sqrtM <- 1/sqrt(phy$edge.length[extern])
            X <- X*sqrtM
            Y <- Y*sqrtM
            # Return the determinant => variance terms  of the 'star' tree
            deterM <- sum(log(phy$edge.length[extern]))
        }
        
    }else{
        if(model!="OUvcv") C <- pruning(phy, trans=FALSE) # FIXME -> to remove the call to OUvcv?
        #if(any(phy$edge.length<=.Machine$double.eps)) C<-list(sqrtM=t(.sqM1(phy)), det=determinant(vcv(phy))$modulus) # FIXME => remove problems with the pruning algorithms on zero branch lengths?
        X <- crossprod(C$sqrtM, X)
        Y <- crossprod(C$sqrtM, Y)
        
        # Return the determinant
        deterM <- C$det
    }
    
    # Adjust the determinant for non-ultrametric OU (see Ho & Ane 2014 - Syst. Bio., p. 401)
    if(flag) deterM <- deterM + 2*sum(log(diagWeight))
    if(REML) deterM <- deterM + determinant(crossprod(X))$modulus - const
    
    # Return the score, variances, and scaled tree
    return(list(phy=phy, diagWeight=diagWeight, X=X, Y=Y, det=deterM, const=const))
}

# Vec operator
.vec <- function(x) as.numeric(x)


# ------------------------------------------------------------------------- #
# .setBounds                                                                #
# options: penalty, model, lower, upper, tol, mserr=NULL, penalized         #
#                                                                           #
# ------------------------------------------------------------------------- #

.setBounds <- function(penalty, model, lower, upper, tol=1e-10, mserr=NULL, penalized=TRUE, corrModel=NULL, k=NULL){
    
    if(is.null(upper)){
        switch(model,
        "EB"={up <- 0},
        "OU"={up <- 30/max(node.depth.edgelength(corrModel$structure))}, # ~ 30 half-lifes upper limit for phylogenetic trees. Should change the format for general models
        "OU1"={up <- 30/max(node.depth.edgelength(corrModel$structure))},
        "OUM"={up <- 30/max(node.depth.edgelength(corrModel$structure))},
        "lambda"={up <- 1},
        "BM"={up <- Inf},
        "BMM"={up <- rep(Inf,k-1)},
        up <- Inf)
    }else{
        up <- upper
    }
    
    if(is.null(lower)){
        switch(model,
        "EB"={low <- -10},
        "OU"={low <- 1e-10},
        "OU1"={low <- 1e-10},
        "OUM"={low <- 1e-10},
        "lambda"={low <- 1e-8},
        "BM"={low <- -Inf},
        "BMM"={low <- rep(-Inf,k-1)},
        low <- -Inf)
    }else{
        low <- lower
    }
    
    # Default tolerance for the parameter search
    if(is.null(tol)){
        if(penalty=="RidgeArch"){
            tol = 1e-8
        }else{
            tol = 0
        }
    }
    
    # Set the default bounds
    if(penalized){
        if(penalty%in%c("RidgeAlt","LASSO")){
            upperBound <- c(log(10e6),up)
            lowerBound <- c(log(tol),low)
        }else if(penalty=="RidgeArch"){
            upperBound <- c(1,up)
            lowerBound <- c(tol,low)
        }
        
        # parameters
        if(model=="BMM"){
            #id1 <- 1; id2 <- 2:(k+1); id3 <- k+2
            id1 <- 1; id2 <- 2:k; id3 <- k+1
        }else if(model=="BM"){
            id1 <- id2 <- 1; id3 <- 2
        }else{
            id1 <- 1; id2 <- 2; id3 <- 3
        }
        
    }else{
        upperBound <- up
        lowerBound <- low
        
        # parameters
        if(model=="BMM"){
            #id1 <- 1; id2 <- 1:k; id3 <- k+1
            id1 <- 1; id2 <- 1:(k-1); id3 <- k
        }else if(model=="BM"){
            id1 <- id2 <- id3 <- 1
        }else{
            id1 <- 1; id2 <- 1; id3 <- 2
        }
    }
    
    # Bounds for error
    if(!is.null(mserr)){
        lowerBound = c(lowerBound,0)
        upperBound = c(upperBound,Inf)
    }
    
    # regularization parameter
    switch(penalty,
    "RidgeArch"={ transformTun <- function(x) (x[id1])},
    "RidgeAlt" ={ transformTun <- function(x) exp(x[id1])},
    "LASSO" ={ transformTun <- function(x) exp(x[id1])},
    transformTun <- function(x) (x[id1])
    )
    
    # model parameter
    switch(model,
    "OU"={ transformPar <- function(x) (x[id2])},
    "BM" ={ transformPar <- function(x) (x[id2])},
    "EB" ={ transformPar <- function(x) (x[id2])},
    "lambda" ={ transformPar <- function(x) (x[id2])},
    "BMM"={transformPar <- function(x) (x[id2]*x[id2])},
    transformPar <- function(x) (x[id2])
    )
    
    # mserr parameter
    transformSE <- function(x) x[id3]*x[id3]
    
    
    bounds <- list(upper=upperBound, lower=lowerBound, trTun=transformTun, trPar= transformPar, trSE=transformSE)
    return(bounds)
}

# ------------------------------------------------------------------------- #
# .startGuess                                                               #
# options: corrModel, cvmethod, mserr, target, penalty, echo, penalized     #
#                                                                           #
# ------------------------------------------------------------------------- #

.startGuess <- function(corrModel, cvmethod, mserr=NULL, target, penalty, echo=TRUE, penalized=TRUE,...){
    if(echo==TRUE) message("Initialization via grid search. Please wait...")
    # Penalization parameters guesses
    if(penalized){
        if(penalty=="RidgeArch"){
            range_val <- c(1e-6, 0.01, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5, 0.7, 0.9)
        }else if(penalty=="RidgeAlt"){
            range_val <- log(c(1e-12, 1e-9, 1e-6, 0.01, 0.1, 1, 10, 100, 1000, 10000))
        }else{
            range_val <- log(c(1e-6, 0.01, 0.1, 1, 10, 100, 1000))
        }
    }else{
        range_val <- NULL
    }
    
    # prepare the list
    list_param <- list()
    list_param[[1]] <- range_val
    list_param[[2]] <- 1 # dummy starting value for BM...
    
    # Models starting guesses
    switch(corrModel$model,
        "OU"={
            mod_val <- log(2)/(max(node.depth.edgelength(corrModel$structure))/c(0.1,0.5,1.5,3,8))
            list_param[[2]] <- mod_val
            index_err <- 3
        },
        "OU1"={
            mod_val <- log(2)/(max(node.depth.edgelength(corrModel$structure))/c(0.1,0.5,1.5,3,8))
            list_param[[2]] <- mod_val
            index_err <- 3
        },
        "OUM"={
            mod_val <- log(2)/(max(node.depth.edgelength(corrModel$structure))/c(0.1,0.5,1.5,3,8))
            list_param[[2]] <- mod_val
            index_err <- 3
        },
        "OUvcv"={
            mod_val <- log(2)/(max(node.depth.edgelength(corrModel$structure))/c(0.1,0.5,1.5,3,8))
            list_param[[2]] <- mod_val
            index_err <- 3
        },
        "lambda"={
            mod_val <- c(0.2,0.5,0.8)
            list_param[[2]] <- mod_val
            index_err <- 3
        },
        "EB"={
            mod_val <- -log(2)/(max(node.depth.edgelength(corrModel$structure))/c(0.1,0.5,1.5,3,8))
            list_param[[2]] <- mod_val
            index_err <- 3
        },
        "BMM"={
            
            # guess starting values
            start_values <- function(tree, data, predictors){
                tip_values <- 1:Ntip(tree)
                index_tips <- tree$edge[,2]%in%tip_values
                maps <- sapply(tree$maps[index_tips], function(x) names(x[length(x)]))
                # check if any tips are missing?
                k = ncol(tree$mapped.edge)
                if(length(unique(maps))<k) {
                    mod_val <- mean(diag(rate_pic(tree, data)))
                    guesses <- as.list(rep(sqrt(mod_val), k))
                }else{
                    guesses <- lapply(colnames(tree$mapped.edge), function(map_names) {
                        dat_red <- which(maps==map_names)
                        tree_red=drop.tip(tree, tree$tip.label[!tree$tip.label%in%tree$tip.label[dat_red]] )
                        if(Ntip(tree_red)<=1){
                                 sqrt(mean(diag(.rate_guess(tree, data[tree$tip.label,], predictors[tree$tip.label,])))) # simple estimate on the whole tree
                             }else{
                                 sqrt(mean(diag(.rate_guess(tree_red, data[tree_red$tip.label,], predictors[tree_red$tip.label,]))))
                             }
                    })
                }
                
                return(guesses)
            }
            
            mod_val <- start_values(corrModel$structure, corrModel$Y, corrModel$X)[-1]
            list_param <- c(list(range_val), mod_val)
            index_err <- length(list_param) + 1
        },
        index_err <- 3
        )
    
    # Prepare the grid search
    if(!is.null(corrModel$mserr)){
        if(corrModel$model=="BMM") list_param[[index_err]] <- sqrt(c(0.001,0.01,0.1,1,10)*mean(unlist(mod_val)^2)) else list_param[[index_err]] <- c(0.001,0.01,0.1,1,10)
        list_param[sapply(list_param, is.null)] <- NULL
        brute_force <- expand.grid(list_param)
    }else{
        list_param[sapply(list_param, is.null)] <- NULL
        brute_force <- expand.grid(list_param)
    }
    
    start <- brute_force[which.min(apply(brute_force, 1, .loocvPhylo,
        cvmethod=cvmethod, # options
        targM=target,
        corrStr=corrModel,
        penalty=penalty,
        error=mserr,
        nobs=corrModel$nobs )),]
    
    if(echo==TRUE & penalized==TRUE)  cat("Best starting for the tuning: ",as.numeric(corrModel$bounds$trTun(start[1])))
    return(start)
}


# ------------------------------------------------------------------------- #
# .check_par_results  (TODO)                                                #
# options: model, par, penalized                                            #
#                                                                           #
# ------------------------------------------------------------------------- #

.check_par_results <- function(corrModel, par, penalized=TRUE){
    
    if(penalized) indice = 2 else indice = 1
    switch(corrModel$model,
        "OU"={
           if(par==corrModel$bounds$upper[indice]) warning("Parameter search reached the upper bound. You should consider increasing the \"upper\" argument value ")
        },
        "OU1"={
           if(par==corrModel$bounds$upper[indice]) warning("Parameter search reached the upper bound. You should consider increasing the \"upper\" argument value ")
        },
        "OUM"={
           if(par==corrModel$bounds$upper[indice]) warning("Parameter search reached the upper bound. You should consider increasing the \"upper\" argument value ")
        },
    )
}

# ------------------------------------------------------------------------- #
# .rate_guess                                                               #
# options: phylo, Y, X                                                      #
#                                                                           #
# ------------------------------------------------------------------------- #

.rate_guess <- function(phylo, Y, X){
    # model
    C <- pruning(phylo, trans=FALSE)
    X <- crossprod(C$sqrtM, X)
    Y <- crossprod(C$sqrtM, Y)
    
    # GLS estimates
    XtX <- pseudoinverse(X)
    B <- XtX%*%Y
    residuals <- Y - X%*%B
    
    # Covariance matrix
    return(crossprod(residuals)/Ntip(phylo))
    }
