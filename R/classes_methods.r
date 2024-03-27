################################################################################
##                                                                            ##
##                       mvMORPH: classes_methods.r                           ##
##                                                                            ##
##   Internal functions for the mvMORPH package                               ##
##                                                                            ##
##  Created by Julien Clavel - 31-07-2018                                     ##
##  (julien.clavel@hotmail.fr/ julien.clavel@biologie.ens.fr)                 ##
##   require: phytools, ape, corpcor, subplex, spam, glassoFast, stats        ##
##                                                                            ##
################################################################################


# ------------ S3 Methods ------------------------------------------------- #
GIC <- function(object, ...) UseMethod("GIC")
EIC <- function(object, nboot=100L, nbcores=1L, ...) UseMethod("EIC")

# ------------------------------------------------------------------------- #
# BIC.mvgls                                                                 #
# options: object, ..., k = 2                                               #
# S3 method - Bayesian Information Criterion - generic from stats           #
# ------------------------------------------------------------------------- #

BIC.mvgls <- function(object, ...){
    # TODO generalize to mvXX functions
    if(object$method=="LL"){
        p <- object$dims$p
        n <- object$dims$n # should take n or n*p?

        LL = object$logLik
        nparam = if(object$model=="BM") (length(object$start_values)-1) + length(object$coefficients) + p*(p + 1)/2 else length(object$start_values) + length(object$coefficients) + p*(p + 1)/2
        # BIC
        BIC = -2*LL+ log(n)*nparam
    }else{
        stop("BIC works only for models fit by Maximum Likelihood (method=\"LL\")")
    }
    
    # return the results
    results <- list(LogLikelihood=LL, BIC=BIC, nparam=nparam, k=log(n))
    class(results) <- c("bic.mvgls","bic")
    return(results)
}

# ------------------------------------------------------------------------- #
# AIC.mvgls                                                                 #
# options: object, ..., k = 2                                               #
# S3 method - Akaike Information Criterion - generic from stats             #
# ------------------------------------------------------------------------- #

AIC.mvgls <- function(object, ..., k = 2){
    
    if(object$method=="LL"){
        p <- object$dims$p
        LL = object$logLik
        nparam = if(object$model=="BM") (length(object$start_values)-1) + length(object$coefficients) + p*(p + 1)/2 else length(object$start_values) + length(object$coefficients) + p*(p + 1)/2
        # AIC
        AIC = -2*LL+k*nparam
    }else{
        stop("AIC works only for models fit by Maximum Likelihood (method=\"LL\")")
    }
    
    # return the results
    results <- list(LogLikelihood=LL, AIC=AIC, nparam=nparam, k=k)
    class(results) <- c("aic.mvgls","aic")
    return(results)
}

# ------------------------------------------------------------------------- #
# GIC.mvgls                                                                 #
# options: model,...                                                        #
# S3 method from "RPANDA"  package                                          #
# ------------------------------------------------------------------------- #
GIC.mvgls <- function(object, ...){
    
    # retrieve arguments
    args <- list(...)
    if(is.null(args[["eigSqm"]])) eigSqm <- TRUE else eigSqm <- args$eigSqm
    if(is.null(args[["REML"]])) args$forceREML <- TRUE else args$forceREML <- args$REML # force to true to follow what has been done in the first paper?
     
    method <- object$method
    penalty <- object$penalty
    target <- object$target
    n <- object$dims$n
    p <- object$dims$p
    m <- object$dims$m
    tuning <- object$tuning
    P <- object$sigma$P # The precision matrix
    Pi <- object$sigma$Pinv # the covariance matrix
    S <- object$sigma$S # the sample estimate
    Target <- .targetM(S=S, targM=target, penalty=penalty)
    beta <- object$coefficients
    
    if(eigSqm){ # to follow the scheme in RPANDA
        sqM1 <- .sqM1(object$corrSt$phy)
        if(!is.null(object$corrSt$diagWeight)){
            w <- 1/object$corrSt$diagWeight
            Y <- crossprod(sqM1, matrix(w*object$variables$Y, nrow=n))
            X <- crossprod(sqM1, matrix(w*object$variables$X, nrow=n))
        }else{
            X <- crossprod(sqM1, object$variables$X)
            Y <- crossprod(sqM1, object$variables$Y)
        }
        residuals <- Y - X%*%beta
    }else{
        residuals <- residuals(object, type="normalized")
        X <- object$corrSt$X
        Y <- object$corrSt$Y
    }
    
    if(object$model=="BM"){
       mod.par=0
    }else if(object$model=="BMM"){
       mod.par=(ncol(object$corrSt$phy$mapped.edge)-1) # should we consider k parameters or k-1 (i.e. relative scaling to the first group)
    }else{
       mod.par=1
    }
    if(is.numeric(object$mserr)) mod.par = mod.par + 1 # already included in the covariance matrix structure?
    if(object$REML & args$forceREML==TRUE) ndimCov = n - m else ndimCov = n
    # Nominal loocv
    XtX <- pseudoinverse(crossprod(X))
    # hat matrix
    h <- diag(X%*%pseudoinverse(X))
    
    # check for hat score of 1 (e.g. MANOVA design)
    nloo <- 1:n
    nloo <- nloo[!h+1e-8>=1]
    nC = length(nloo)
    
    if(penalty=="RidgeArch"){
        
        # First and second derivative of the functional (we can use patterned matrix to target some matrix elements)
        # We use the Kronecker-vec identity to speed up the computations
        T1 <- sapply(nloo, function(i){
            Sk <- tcrossprod(residuals[i,]) ;
            VSV <- 0.5*(Pi - (1-tuning)*Sk - tuning*Target);
            VSV2 <- 0.5*(Pi - Sk);
            sum(VSV * 2*(P%*%VSV2%*%P))
        })
        
        df = sum(T1)/nC
        sigma_df <- df
        
    }else if(penalty=="LASSO" | penalty=="LL"){
        
        # LASSO or ML
        Tf2 <- function(S, P) {
            I <- ifelse(P==0,0,1) ;
            t(.vec(S*I))%*%.vec(P%*%(S*I)%*%P)
        }
        
        sigma_df <- (1/(2*nC))*sum(sapply(nloo, function(i){ Tf2(tcrossprod(residuals[i,]) , P)})) - (1/2)*Tf2(S,P)
        
    }else if(penalty=="RidgeAlt"){
        # Alternative Ridge
        eig <- eigen(Pi)
        V <- eig$vectors
        d <- eig$values
        H <- (1/(0.5*(kronecker(d,d)+tuning)))
        
        # 2) First derivative of the functional
        T1 <- sapply(nloo, function(i){
            Sk <- tcrossprod(residuals[i,]) ;
            VSV <- .vec(crossprod(V, (0.5*(Pi - (Sk - tuning*Target) - tuning*P))%*%V));
            VSV2 <- .vec(crossprod(V, (0.5*(Pi - Sk))%*%V));
            sum(VSV * (H*VSV2))
        })
        
        df = sum(T1)/nC
        sigma_df <- df
    }
    
    # Number of parameters for the root state:
    # The Information matrix from the Hessian and gradients scores
    T2 <- sapply(nloo, function(i){
        gradient <- (X[i,])%*%t(P%*%t(Y[i,]-X[i,]%*%beta))
        sum(gradient * (XtX%*%gradient%*%Pi))
    })
    beta_df <- sum(T2)
    
    if( min(m, sum(object$dims$assign!=0))>1 ) warning("GIC criterion with multiple predictors has not been fully tested. Please use it with care and consider EIC or simulations instead")
    
    # LogLikelihood (minus)
    DP <- as.numeric(determinant(Pi)$modulus)
    if(object$REML==TRUE & args$forceREML==FALSE) Ccov <- as.numeric(object$corrSt$det - determinant(crossprod(object$corrSt$X))$modulus + object$corrSt$const) else Ccov <- as.numeric(object$corrSt$det)
    llik <- 0.5 * (ndimCov*p*log(2*pi) + p*Ccov + ndimCov*DP + ndimCov*sum(S*P))
    GIC <- 2*llik + 2*(sigma_df+beta_df+mod.par)
    
    # return the results
    results <- list(LogLikelihood=-llik, GIC=GIC, p=p, n=n, bias=sigma_df+beta_df+mod.par, bias_cov=sigma_df, args=args)
    class(results) <- c("gic.mvgls","gic")
    return(results)
}

# ------------------------------------------------------------------------- #
# EIC.mvgls                                                                 #
# options: object, nboot, nbcores, ...                                      #
# S3 method - Extended/Efron Information Criterion                          #
# ------------------------------------------------------------------------- #

EIC.mvgls <- function(object, nboot=100L, nbcores=1L, ...){
    
    # retrieve arguments
    args <- list(...)
    if(is.null(args[["eigSqm"]])) eigSqm <- TRUE else eigSqm <- args$eigSqm
    if(is.null(args[["restricted"]])) restricted <- FALSE else restricted <- args$restricted
    if(is.null(args[["REML"]])) args$forceREML <- FALSE else args$forceREML <- args$REML
    
    # retrieve data to simulate bootstrap samples
    beta <- object$coefficients
    if(eigSqm){ # to follow the scheme in RPANDA
        sqM1 <- .sqM1(object$corrSt$phy)
        if(!is.null(object$corrSt$diagWeight)){
            w <- 1/object$corrSt$diagWeight
            Y <- crossprod(sqM1, matrix(w*object$variables$Y, nrow=object$dims$n))
            X <- crossprod(sqM1, matrix(w*object$variables$X, nrow=object$dims$n))
        }else{
            X <- crossprod(sqM1, object$variables$X)
            Y <- crossprod(sqM1, object$variables$Y)
        }
        residuals <- Y - X%*%beta
    }else{
        residuals <- residuals(object, type="normalized")
        X <- object$corrSt$X
        Y <- object$corrSt$Y
    }
    
    N = nrow(Y)
    p = object$dims$p
    if(object$REML & args$forceREML==TRUE) ndimCov = object$dims$n - object$dims$m else ndimCov = object$dims$n
    tuning <- object$tuning
    target <- object$target
    penalty <- object$penalty
    if(is.null(object$corrSt$diagWeight)){
        diagWeight <- 1; is_weight = FALSE
    }else{
        diagWeight <- object$corrSt$diagWeight; is_weight = TRUE
        diagWeightInv <- 1/diagWeight
    }
    Dsqrt <- .pruning_general(object$corrSt$phy, trans=FALSE, inv=FALSE)$sqrtM # return warning message if n-ultrametric tree is used with OU?
    # TODO (change to allow n-ultrametric and OU) > just need to standardize the data by the weights
    # if(object$model=="OU" & !is.ultrametric(object$variables$tree)) stop("The EIC method does not handle yet non-ultrametric trees with OU processes")
    
    DsqrtInv <- .pruning_general(object$corrSt$phy, trans=FALSE, inv=TRUE)$sqrtM
    modelPerm <- object$call
    modelPerm$grid.search <- quote(FALSE)
    modelPerm$start <- quote(object$opt$par)
    
    # Mean and residuals for the model
    MeanNull <- object$variables$X%*%beta
    
    
    # Estimate the bias term
    D1 <- function(objectBoot, objectFit, ndimCov, p, sqM, Ccov2){ # LL(Y*|param*) - LL(Y*| param)
        
        # Y*|param*
        residualsBoot <- residuals(objectBoot, type="normalized")
        
        # For boot "i" LL1(Y*|param*)
        if(objectFit$REML==TRUE & args$forceREML==FALSE) Ccov1 <- as.numeric(objectBoot$corrSt$det - determinant(crossprod(objectBoot$corrSt$X))$modulus + objectBoot$corrSt$const) else Ccov1 <- as.numeric(objectBoot$corrSt$det)
        Gi1 <- try(chol(objectBoot$sigma$Pinv), silent=TRUE)
        if(inherits(Gi1, 'try-error')) return("error")
        quadprod <- sum(backsolve(Gi1, t(residualsBoot), transpose = TRUE)^2)
        detValue <- sum(2*log(diag(Gi1)))
        llik1 <- -0.5 * (ndimCov*p*log(2*pi) + p*Ccov1 + ndimCov*detValue + quadprod)
        
        # Y*|param
        #if(!restricted) residualsBoot <- objectBoot$corrSt$Y - objectBoot$corrSt$X%*%objectFit$coefficients # does not account for the phylo model of the original fit
        if(!restricted){
            if(is_weight){
                residualsBoot <- crossprod(sqM, (objectBoot$variables$Y - objectBoot$variables$X%*%objectFit$coefficients)*diagWeightInv)
            }else{
                residualsBoot <- crossprod(sqM, objectBoot$variables$Y - objectBoot$variables$X%*%objectFit$coefficients)
                
            }
        }
        
        # For boot "i" LL2(Y*|param)
        # if(objectFit$REML==TRUE & args$forceREML==FALSE) Ccov2 <- as.numeric(objectFit$corrSt$det - determinant(crossprod(objectFit$corrSt$X))$modulus + objectFit$corrSt$const) else Ccov2 <- as.numeric(objectFit$corrSt$det)
        Gi2 <- try(chol(objectFit$sigma$Pinv), silent=TRUE)
        if(inherits(Gi2, 'try-error')) return("error")
        quadprod <- sum(backsolve(Gi2, t(residualsBoot), transpose = TRUE)^2)
        detValue <- sum(2*log(diag(Gi2)))
        llik2 <- -0.5 * (ndimCov*p*log(2*pi) + p*Ccov2 + ndimCov*detValue + quadprod)
        
        # Return the difference in LL for D1
        return(llik1 - llik2)
    }
    
    D3 <- function(objectBoot, objectFit, loglik, ndimCov, p){ # LL(Y|param) - LL(Y| param*)
        
        # Y|param*
        if(!restricted) {
            sqM_temp <- .pruning_general(objectBoot$corrSt$phy, trans=FALSE, inv=TRUE)$sqrtM
            if(is_weight){
                residualsBoot <- try(crossprod(sqM_temp, (objectFit$variables$Y - objectFit$variables$X%*%objectBoot$coefficients)/objectBoot$corrSt$diagWeight), silent=TRUE)
            } else {
                residualsBoot <- try(crossprod(sqM_temp, objectFit$variables$Y - objectFit$variables$X%*%objectBoot$coefficients), silent=TRUE)
                
            }
        }else{ residualsBoot <- objectFit$corrSt$Y - objectFit$corrSt$X%*%objectFit$coefficients}
        
        #if(!restricted) residualsBoot <- objectFit$corrSt$Y - objectFit$corrSt$X%*%objectBoot$coefficients
        #else residualsBoot <- objectFit$corrSt$Y - objectFit$corrSt$X%*%objectFit$coefficients
        
        # For boot "i" LL2(Y|param*)
        if(objectFit$REML==TRUE & args$forceREML==FALSE) Ccov1 <- as.numeric(objectBoot$corrSt$det - determinant(crossprod(objectBoot$corrSt$X))$modulus + objectBoot$corrSt$const) else Ccov1 <- as.numeric(objectBoot$corrSt$det)
        Gi1 <- try(chol(objectBoot$sigma$Pinv), silent=TRUE)
        if(inherits(Gi1, 'try-error')) return("error")
        quadprod <- sum(backsolve(Gi1, t(residualsBoot), transpose = TRUE)^2)
        detValue <- sum(2*log(diag(Gi1)))
        llik2 <- -0.5 * (ndimCov*p*log(2*pi) + p*Ccov1 + ndimCov*detValue + quadprod)
        
        # Return the difference in LL for D1
        return(loglik - llik2)
    }
    
    # Estimate EIC: LL+bias
    
    # Maximum Likelihood
    if(object$REML==TRUE & args$forceREML==FALSE) Ccov <- as.numeric(object$corrSt$det - determinant(crossprod(object$corrSt$X))$modulus + object$corrSt$const) else Ccov <- as.numeric(object$corrSt$det)
    Gi <- try(chol(object$sigma$Pinv), silent=TRUE)
    if(inherits(Gi, 'try-error')) return("error")
    quadprod <- sum(backsolve(Gi, t(residuals), transpose = TRUE)^2)
    detValue <- sum(2*log(diag(Gi)))
    llik <- -0.5 * (ndimCov*p*log(2*pi) + p*Ccov + ndimCov*detValue + quadprod)
    
    # Estimate parameters on bootstrap samples
    bias <- pbmcmapply(function(i){
        
        # generate bootstrap sample
        Yp <- MeanNull + Dsqrt%*%(residuals[sample(N, replace=TRUE),])*diagWeight # sampling with replacement for bootstrap
        rownames(Yp) <- rownames(object$variables$Y)
        
        modelPerm$response <- quote(Yp);
        estimModelNull <- eval(modelPerm);
        d1res <- D1(objectBoot=estimModelNull, objectFit=object, ndimCov=ndimCov, p=p, sqM=DsqrtInv, Ccov2=Ccov)
        d3res <- D3(objectBoot=estimModelNull, objectFit=object, loglik=llik, ndimCov=ndimCov, p=p)
        d1res+d3res
        
    }, 1:nboot, mc.cores = getOption("mc.cores", nbcores))
    
    # check for errors first?
    bias <- .check_samples(bias)
    nboot_eff <- length(bias)
    
    # compute the EIC
    pboot <- mean(bias)
    EIC <- -2*llik + 2*pboot
    
    # standard-error
    se <- sd(bias)/sqrt(nboot_eff)
    
    # concatenate the results
    results <- list(EIC=EIC, bias=bias, LogLikelihood=llik, se=se, p=p, n=N)
    class(results) <- c("eic.mvgls","eic")
    
    return(results)
}

# ------------------------------------------------------------------------- #
# fitted.values.mvgls  / fitted.mvgls                                       #
# options: object, ...                                                      #
# S3 method from "stats" package                                            #
# ------------------------------------------------------------------------- #
fitted.mvgls <- function(object, ...){
    return(object$fitted)
}

# ------------------------------------------------------------------------- #
# residuals.mvgls                                                           #
# options: type = c("response","normalized")                                #
# S3 method "mvgls" class                                                   #
# ------------------------------------------------------------------------- #
residuals.mvgls <- function(object, type=c("response","normalized"), ...){
    type <- match.arg(type)[1]
    if(type=="response"){
        residuals <- object$residuals
    }else{
        residuals <- .mvGLS(object$corrSt)$residuals
    }
    return(residuals)
}

# ------------------------------------------------------------------------- #
# vcov.mvgls                                                                #
# options: object, ...                                                      #
# S3 method from "stats" package                                            #
# ------------------------------------------------------------------------- #
vcov.mvgls <- function(object, ...){
    args <- list(...)
    if(is.null(args[["type"]])) type <- "coef" else type <- args$type
    
    switch(type,
    "covariance"={return(object$sigma$Pinv)}, # inverse of the precision matrix
    "precision"={return(object$sigma$P)},     # precision matrix
    "coef"={
        XtX <- solve(crossprod(object$corrSt$X))
        covBeta <- kronecker(object$sigma$Pinv, XtX)
        rownames(covBeta) <- colnames(covBeta) <- rep(attr(object$variables$X,"dimnames")[[2]], object$dims$p)
       
        return(covBeta)})
}

# ------------------------------------------------------------------------- #
# coef.mvgls     / coefficients.mvgls                                       #
# options: object, ...                                                      #
# S3 method from "stats" package                                            #
# ------------------------------------------------------------------------- #
coef.mvgls <- function(object, ...){
    
    coeffs <- object$coefficients
    rownames(coeffs) <- attr(object$variables$X,"dimnames")[[2]]
    
    return(coeffs)
}

# ------------------------------------------------------------------------- #
# logLik.mvgls                                                              #
# options: object, ...                                                      #
# S3 method from "stats" package                                            #
# ------------------------------------------------------------------------- #
logLik.mvgls<-function(object,...){
    
    if(object$method=="LL"){
        LL = -object$logLik
    }else{
        # param
        n <- object$dims$n
        p <- object$dims$p
        m <- object$dims$rank
        if(object$REML) ndimCov = n - m else ndimCov = n
        DP <- as.numeric(determinant(object$sigma$Pi)$modulus)
        Ccov <- object$corrSt$det
        LL <- -0.5 * (ndimCov*p*log(2*pi) + p*Ccov + ndimCov*DP + ndimCov*sum(object$sigma$S*object$sigma$P))
    }
    return(LL)
}

# ------------ S3 Printing Methods ----------------------------------------- #

# Generic S3 print for linear models in R stats library (R core team).
print.mvgls <- function(x, digits = max(3L, getOption("digits") - 3L), ...){
    # model call
    cat("\nCall:\n",
    paste(deparse(x$call), sep = "", collapse = "\n"), "\n\n", sep = "")
    
    # loocv or LL
    meth <- ifelse(x$REML, "REML", "ML")
    if(x$method=="LL"){
        if(inherits(x, "mvols")) cat("\nOrdinary least squares fit by",meth,"\n") else cat("\nGeneralized least squares fit by",meth,"\n")
        if(x$REML) cat("Log-restricted-likelihood:",round(x$logLik, digits=digits), "\n\n") else cat("Log-likelihood:",round(x$logLik, digits=digits), "\n\n")
    }else{
        if(inherits(x, "mvols")) cat("\nOrdinary least squares fit by penalized",meth,"\n") else cat("\nGeneralized least squares fit by penalized",meth,"\n")
        if(x$REML){
            cat("LOOCV of the log-restricted-likelihood:",round(x$logLik, digits=digits), "\n\n")
        }else{
            cat("LOOCV of the log-likelihood:",round(x$logLik, digits=digits), "\n\n")
        }
    }
    
    
    # Model parameters
    cat("\nParameter estimate(s):\n")
    if(!any(is.na(x$param))){
        switch(x$model,
        "OU"={ cat("alpha:",round(x$param, digits=digits),"\n\n")},
        "EB"={ cat("r:",round(x$param, digits=digits),"\n\n")},
        "lambda"={cat("lambda:",round(x$param, digits=digits),"\n\n")},
        "BMM"={print(round(x$param, digits=digits)); cat("\n")},
        cat("parameter(s):",round(x$param, digits=digits),"\n\n")
        )
    }
    
    # Regularization parameter
    if(!is.na(x$tuning)){
        cat("Regularization parameter (gamma):", round(x$tuning, digits=digits), "\n\n")
    }
    
    # size of the evolutionary covariance matrix
    cat("\nCovariance matrix of size:",x$dims$p,"by",x$dims$p,"\n")
    cat("for",x$dims$n,"observations","\n\n")
    
    # coefficients of the linear model
    if(ncol(coef(x))<10) {
        cat("Coefficients:\n")
        print.default(format(coef(x), digits = digits),
        print.gap = 2L, quote = FALSE)
    } else {
        cat("Coefficients (truncated):\n")
        coefHead<-coef(x)[,1:10,drop=FALSE]
        print(coefHead, digits = digits, quote = FALSE)
        cat("Use \"coef\" to display all the coefficients\n")}
    cat("\n")
    invisible(x)
}


# Generic S3 print for linear models in R stats library (R core team).
print.summary.mvgls <- function(x, digits = max(3, getOption("digits") - 3), ...){
   
    # model call
    cat("\nCall:\n",
    paste(deparse(x$call), sep = "", collapse = "\n"), "\n\n", sep = "")
    
    # loocv or LL
    meth <- ifelse(x$REML, "REML", "ML")
    
    if(x$method=="LL"){
        if(x$GLS) cat("\nGeneralized least squares fit by",meth,"\n") else cat("\nOrdinary least squares fit by",meth,"\n")
        print(x$results.fit,  quote = FALSE )
    }else{
        
        if(x$GLS) cat("\nGeneralized least squares fit by penalized",meth,"\n") else cat("\nOrdinary least squares fit by penalized",meth,"\n")
        print(x$results.fit,  quote = FALSE )
    }
    
    
    # Model parameters
    cat("\nParameter estimate(s):\n")
    if(!any(is.na(x$param))){
        switch(x$model,
        "OU"={ cat("alpha:",round(x$param, digits=digits),"\n\n")},
        "EB"={ cat("r:",round(x$param, digits=digits),"\n\n")},
        "lambda"={cat("lambda:",round(x$param, digits=digits),"\n\n")},
        "BMM"={print(round(x$param, digits=digits)); cat("\n")},
        cat("parameter(s):",round(x$param, digits=digits),"\n\n")
        )
    }
    
    # Regularization parameter
    if(!is.na(x$tuning)){
        cat("Regularization parameter (gamma):", round(x$tuning, digits=digits), "\n\n")
    }
    
    # Error parameter
    if(is.numeric(x$mserr)){
        cat("Nuisance parameter (error variance):", round(x$mserr, digits=digits), "\n\n")
    }
    
    # size of the evolutionary covariance matrix
    cat("\nCovariance matrix of size:",x$dims$p,"by",x$dims$p,"\n")
    cat("for",x$dims$n,"observations","\n\n")
    
    # coefficients of the linear model
    if(ncol(coef(x))<10) {
        cat("Coefficients:\n")
        print.default(format(coef(x), digits = digits),
        print.gap = 2L, quote = FALSE)
    } else {
        cat("Coefficients (truncated):\n")
        coefHead<-coef(x)[,1:10,drop=FALSE]
        print(coefHead, digits = digits, quote = FALSE)
        cat("Use \"coef\" to display all the coefficients\n")}
    cat("\n")
    invisible(x)
}

summary.mvgls <- function(object, ...){
    
    # param
    n <- object$dims$n
    p <- object$dims$p
    m <- object$dims$rank
    if(object$REML) ndimCov = n - m else ndimCov = n
    
    # loocv or LL
    meth <- ifelse(object$REML, "REML", "ML")
    
    if(object$method=="LL"){
        LL = object$logLik
        nparam = if(object$model=="BM") (length(object$start_values)-1) + length(object$coefficients) + p*(p + 1)/2 else length(object$start_values) + length(object$coefficients) + p*(p + 1)/2 
        # AIC
        AIC = -2*LL+2*nparam
        # GIC
        GIC = GIC(object)$GIC
        
        results.fit <- data.frame("AIC"=AIC, "GIC"=GIC, "logLik"=LL, row.names = " ")
        
    }else{
        # LogLikelihood (minus)
        DP <- as.numeric(determinant(object$sigma$Pi)$modulus)
        Ccov <- object$corrSt$det
        LL <- -0.5 * (ndimCov*p*log(2*pi) + p*Ccov + ndimCov*DP + ndimCov*sum(object$sigma$S*object$sigma$P))
        # GIC
        GIC = GIC(object)$GIC
        results.fit <- data.frame("GIC"=GIC, "logLik"=LL, row.names = " ")
    }
    
    
    object$results.fit <- results.fit
    object$GLS <- if(inherits(object,"mvols")) FALSE else TRUE
    class(object) <- c("summary.mvgls","mvgls")
    object
}

# AIC printing options
print.aic.mvgls<-function(x,...){
    cat("\n")
    message("-- Akaike Information Criterion --","\n")
    cat("AIC:",x$AIC,"| Log-likelihood",x$LogLikelihood,"\n")
    cat("\n")
}

# GIC printing options
print.gic.mvgls<-function(x,...){
    cat("\n")
    message("-- Generalized Information Criterion --","\n")
    cat("GIC:",x$GIC,"| Log-likelihood",x$LogLikelihood,"\n")
    cat("\n")
}

# EIC printing options
print.eic.mvgls<-function(x,...){
    cat("\n")
    message("-- Extended Information Criterion --","\n")
    cat("EIC:",x$EIC,"| +/-",3.92*x$se,"| Log-likelihood",x$LogLikelihood,"\n")
    cat("\n")
}

# ------------------------------------------------------------------------- #
# .parallel_mapply wrapper switch options for parallel calculation          #
# options: ...                                                              #
#                                                                           #
# ------------------------------------------------------------------------- #

.parallel_mapply <- function(FUN,..., MoreArgs = NULL, mc.style = "ETA", mc.substyle = NA,
                       mc.cores = getOption("mc.cores", 2L),
                       ignore.interactive = getOption("ignore.interactive", F),
                       mc.preschedule = TRUE, mc.set.seed = TRUE, mc.cleanup = TRUE, verbose=TRUE){
    
    if(verbose){
        return(pbmcmapply(FUN, ..., MoreArgs = MoreArgs, mc.style = mc.style, mc.substyle = mc.substyle, mc.cores = mc.cores,
                          ignore.interactive = ignore.interactive, mc.preschedule = mc.preschedule, mc.set.seed = mc.set.seed, mc.cleanup = mc.cleanup))
    }else{
        return(mcmapply(FUN, ..., MoreArgs = MoreArgs, mc.cores = mc.cores,
                          mc.preschedule = mc.preschedule, mc.set.seed = mc.set.seed, mc.cleanup = mc.cleanup))
    }
}

# ------------------------------------------------------------------------- #
# .pruning_general wrapper switch options for OLS vs GLS and various models #
# options: tree, inv, scaled, trans, check                                  #
#                                                                           #
# ------------------------------------------------------------------------- #

.pruning_general <- function(tree, inv=TRUE, scaled=TRUE, trans=TRUE, check=TRUE){
    if(inherits(tree, "phylOLS")){
        n <- Ntip(tree)
        if((sum(tree$edge.length) - n)<=.Machine$double.eps){
            # Return the determinant
            det <- 0
            sqrtMat <- diag(n)
        }else{
            descendent <- tree$edge[,2]
            extern <- (descendent <= n)
            if(inv) sqrt_phy <- 1/sqrt(tree$edge.length[extern]) else sqrt_phy <- sqrt(tree$edge.length[extern])
            sqrtMat <- diag(sqrt_phy)
            # Return the determinant => variance terms of the 'star' tree
            det <- sum(2*log(sqrt_phy))
        }
        return(list(sqrtMat=sqrtMat, det=det))
    }else{
        return(pruning(tree, inv=inv, scaled=scaled, trans=trans, check=check))
    }
}
                        
# ------------------------------------------------------------------------- #
# print option for MANOVA tests  (output borrowed from "car" package)       #
# options: x, digits, ...                                                   #
#                                                                           #
# ------------------------------------------------------------------------- #


print.manova.mvgls <- function(x, digits = max(3L, getOption("digits") - 3L), ...){
  
  # select the appropriate output
  if(x$param){
    
    if(x$type=="I") cat("Sequential MANOVA Tests:",x$test,"test statistic","\n")
    if(x$type=="II") cat("Type II MANOVA Tests:",x$test,"test statistic","\n")
    if(x$type=="III") cat("Type III MANOVA Tests:",x$test,"test statistic","\n")
    if(x$type=="glh") cat("General Linear Hypothesis Test:",x$test,"test statistic","\n")
    if(x$type=="glhrm") cat("General Linear Hypothesis Test (repeated measures design):",x$test,"test statistic","\n")
    
    signif <- sapply(x$pvalue, function(i) if(i<0.001){"***"}else if(i<0.01){
      "**"}else if(i<0.05){"*"}else if(i<0.1){"."}else{""})
    
    table_results <- data.frame(Df=x$Df, stat=x$stat, approxF=x$approxF, numDf=x$NumDf, denDf=x$DenDf, pval=x$pvalue, signif=signif)
    if(x$type!="glh" & x$type!="glhrm"){
        if(x$type=="III") rownames(table_results) <- x$terms[unique(x$dims$assign)+1] else rownames(table_results) <- x$terms[unique(x$dims$assign)]
    }else{
        rownames(table_results) <- "Contrasts L"
    }
    colnames(table_results) <- c("Df", "test stat", "approx F", "num Df", "den Df", "Pr(>F)", "")
    print(table_results, digits = digits, ...)
    cat("---","\n")
    cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1","\n")
    
    
    
  }else{ # permutation methods
    
    if(x$type=="I") cat("Sequential MANOVA Tests with",x$nperm,"permutations:",x$test,"test statistic","\n")
    if(x$type=="II") cat("Type II MANOVA Tests with",x$nperm,"permutations:",x$test,"test statistic","\n")
    if(x$type=="III") cat("Type III MANOVA Tests with",x$nperm,"permutations:",x$test,"test statistic","\n")
    if(x$type=="glh")  cat("General Linear Hypothesis Test with",x$nperm,"permutations:",x$test,"test statistic","\n")
    if(x$type=="glhrm")  cat("General Linear Hypothesis Test (repeated measures design) with",x$nperm,"permutations:",x$test,"test statistic","\n")
    signif <- sapply(x$pvalue, function(i) if(i<0.001){"***"}else if(i<0.01){
      "**"}else if(i<0.05){"*"}else if(i<0.1){"."}else{""})
    
    table_results <- data.frame(stat=x$stat, pval=x$pvalue, signif=signif)
    if(x$type!="glh" & x$type!="glhrm"){
        if(x$type=="III") rownames(table_results) <- x$terms[unique(x$dims$assign)+1] else rownames(table_results) <- x$terms[unique(x$dims$assign)]
    }else{
        rownames(table_results) <- "Contrasts L"
    }
    colnames(table_results) <- c("Test stat", "Pr(>Stat)", "")
    print(table_results, digits = digits, ...)
    cat("---","\n")
    cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1","\n")
    
  }
  
}
                     
                    
# ------------------------------------------------------------------------- #
# plot option for MANOVA tests  (distribution of the test statistic)        #
# options: x, ...                                                           #
#                                                                           #
# ------------------------------------------------------------------------- #

plot.manova.mvgls <- function(x,...){
    
    args <- list(...)
    if(is.null(args[["density"]])) density = FALSE else density = args$density
    if(is.null(args[["breaks"]])) breaks = 50 else breaks = args$breaks
    
    nterms <- length(x$terms)
    
    if(x$param==TRUE){
        
        for(i in 1:nterms){
            df_mod <- x
            d1=df_mod$NumDf[i]
            d2=df_mod$DenDf[i]
            
        curve(df(x, df1=d1, df2=d2), 0, qf(0.9999, d1, d2), las=1, 
           main=paste("F test:",x$terms[i]), 
           xlab=paste("F value","(",round(x$approxF[i],3),")","p-value :", round(x$pvalue[i],3)), ylab="density" );
            abline(v=x$approxF[i], col="red")
        }
        
    }else{
    
    
    # plot histogram with permuted statistics
    for(i in 1:nterms){
      
      if(density){
          plot(density(x$nullstat[,i]), main=paste("Statistic distribution:",x$terms[i]),xlab=paste(x$test,"(",round(x$stat[i],3),")","p-value :",
                round(x$pvalue[i],3)), las=1, xlim=range(c(x$nullstat[,i],x$stat[i])))
      }else{
          hist(x$nullstat[,i], main=paste("Statistic distribution:",x$terms[i]),
        xlab=paste(x$test,"(",round(x$stat[i],3),")","p-value :",
                round(x$pvalue[i],3)), las=1, breaks=breaks, border=NA, col="lightgrey", xlim=range(c(x$nullstat[,i],x$stat[i]))); 
      }
        abline(v=x$stat[i], col="red", lwd=2)
        }
    }

}

# ------------------------------------------------------------------------- #
# plot option for Pairwise tests  (distribution of the test statistic)      #
# options: x, ...                                                           #
#                                                                           #
# ------------------------------------------------------------------------- #

plot.pairs.mvgls <- function(x,...){
  
  args <- list(...)
  if(is.null(args[["density"]])) density = FALSE else density = args$density
  if(is.null(args[["breaks"]])) breaks = 50 else breaks = args$breaks
  
  nterms <- nrow(x$L)
  if(is.null(rownames(x$L))) namesContrasts <- "contrasts L" else namesContrasts <- rownames(x$L)
  
  if(x$param==TRUE){
    
    for(i in 1:nterms){
      df_mod <- x
      d1=df_mod$NumDf[i]
      d2=df_mod$DenDf[i]
      
      curve(df(x, df1=d1, df2=d2), 0, qf(0.9999, d1, d2), las=1,
            main=paste("F test:", namesContrasts[i]),
            xlab=paste("F value","(",round(x$approxF[i],3),")","p-value :", round(x$pvalue[i],3)), ylab="density" );
      abline(v=x$approxF[i], col="red")
    }
    
  }else{
    
    
    # plot histogram with permuted statistics
    for(i in 1:nterms){
      
      if(density){
        plot(density(x$nullstat[,i]), main=paste("Statistic distribution:",x$terms[i]),xlab=paste(x$test,"(",round(x$stat[i],3),")","p-value :",
                                                                                                  round(x$pvalue[i],3)), las=1, xlim=range(c(x$nullstat[,i],x$stat[i])))
      }else{
        hist(x$nullstat[,i], main=paste("Statistic distribution:",namesContrasts[i]),
             xlab=paste(x$test,"(",round(x$stat[i],3),")","p-value :",
                        round(x$pvalue[i],3)), las=1, breaks=breaks, border=NA, col="lightgrey", xlim=range(c(x$nullstat[,i],x$stat[i])));
      }
      abline(v=x$stat[i], col="red", lwd=2)
    }
  }
  
}



# ------------------------------------------------------------------------- #
# plot.mvgls                                                                #
# options: x, term, ..., fitted=TRUE                                        #
#                                                                           #
# ------------------------------------------------------------------------- #
plot.mvgls <- function(x, term, ..., fitted=FALSE, residuals=FALSE){
    
    if(missing(term)){
        term <- which(attr(x$variables$X,"dimnames")[[2]]!="(Intercept)")[1]
        term <- attr(x$variables$X,"dimnames")[[2]][term]
    }
    
    if(!is.numeric(term) & !term%in%attr(x$variables$X,"dimnames")[[2]]) stop("Unknown predictor name.","\n")
    
    # based on Drake & Klingenberg 2008 shape score
    betas <- coefficients(x)[term,,drop=TRUE]
    standBeta <- betas %*% sqrt(solve(crossprod(betas)))
    if(residuals) scoreVar <- (x$residuals)%*% standBeta else scoreVar <- (x$variables$Y)%*% standBeta
    
    # plot
    plot(scoreVar ~ x$variables$X[,term], xlab=term, ylab="mvScore", ...)
    
    # plot predictions on the same space?
    if(fitted){
        scoreVar2 <- (x$fitted ) %*% standBeta
        points(scoreVar2 ~ x$variables$X[,term], col="red", pch=16)
    }
    
    # loess on the residuals?
    if(residuals){
        abline(h=0, lty=2)
        scores_residuals <- data.frame(score=scoreVar, xvar=x$variables$X[,term])
        loess_fit <- loess(score ~ xvar, data=scores_residuals)
        xseq <- seq(from=min(scores_residuals$xvar), to=max(scores_residuals$xvar), length=80)
        pred <- predict(loess_fit, newdata=data.frame(xvar=xseq))
        lines(pred~xseq, col="red", xpd=FALSE)
    }
    
    
    results <- list(scores = scoreVar, standBeta=standBeta, betas=betas, term=term)
    invisible(results)
}

# ------------------------------------------------------------------------- #
# predict.mvgls                                                             #
# options: object, newdata, ...                                             #
#                                                                           #
# ------------------------------------------------------------------------- #
predict.mvgls <- function(object, newdata, ...){
    
    args <- list(...)
    # if "tree" is provided
    if(!is.null(args[["tree"]])){
        if(!inherits(args$tree, "phylo")) stop("the provided tree is not of class \"phylo\" ") else tree <- args$tree
        if(!is.data.frame(newdata)) stop("the \"newdata\" should be a data.frame object with column names matching predictors names, and row names matching names in the tree ")
    } else tree <- NULL
    if(is.null(args[["na.action"]])) na.action <- na.pass else na.action <- args$na.action
    
    # check if newdata is provided
    if(missing(newdata) || is.null(newdata)) {
        X <- object$variables$X # simply return fitted values when newdata is empty
    }else{
        
        Terms <- delete.response(object$terms)
        # as in "stats v3.3.0"
        m <- model.frame(Terms, newdata, xlev = object$xlevels, na.action = na.action)
        
        # check the arguments
        if(!is.null(cl <- attr(Terms, "dataClasses"))) .checkMFClasses(cl, m)
        X <- model.matrix(Terms, m, contrasts.arg = object$contrasts)
        
        # FIXME allow lists
        predictors_names <- rownames(newdata)
    }
    
    
    # GLS/OLS prediction
    if(is.null(tree)){
        predicted <- X%*%object$coefficients # simply return fitted values when newdata is empty
    }else{
        rcov <- .resid_cov_phylo(tree, object, predictors_names)
        predicted <- X%*%object$coefficients + rcov$w%*%solve(rcov$Vt)%*%object$residuals[rcov$train,,drop=FALSE] # FIXME account for (multivariate) variance scaling? Rao & Toutenberg => no just the correlation structure
    }
    
    return(predicted)
}

# ------------------------------------------------------------------------- #
# .resid_cov_phylo                                                          #
# options: tree, object, sp_name, ...                                       #
#                                                                           #
# ------------------------------------------------------------------------- #
.resid_cov_phylo <- function(tree, object, sp_name, ...){
    
    if(is.null(sp_name)) stop("You must provide species names to \"newdata\"")
    if(any(!sp_name%in%tree$tip.label)) stop("the \"newdata\" names does not matches names in the tree ")
    train_sample <- tree$tip.label[!tree$tip.label%in%sp_name]
    
    # check first that species in the training sample are the same as in the model fit object
    if(any(!train_sample%in%object$corrSt$phy$tip.label)) train_sample <- object$corrSt$phy$tip.label
    
    # helper to obtain the covariances between data used in a model and newdata
    switch(object$model,
    "BM"={ V <- vcv.phylo(tree)},
    "OU"={
        V <- .Call("mvmorph_covar_ou_fixed", A=vcv.phylo(tree), alpha=as.double(object$param), sigma=1, PACKAGE="mvMORPH")
        rownames(V) <- colnames(V) <- tree$tip.label
    },
    "EB"={ V <- vcv.phylo(.transformPhylo(tree, model="EB", param=object$param)) },
    "lambda"={ V <- vcv.phylo(.transformPhylo(tree, model="lambda", param=object$param)) },
    #FIXME -- add BMM
    "BMM"={stop("BMM model is not handled yet. Please contact the author for further assistance.")},
    )
    
    # If error=TRUE, we add it to the covariance matrix here
    #if(!is.na(object$mserr)) diag(V) = diag(V) + object$mserr
    
    # Build the covariance matrices
    w <- V[sp_name, train_sample, drop=FALSE]
    Vt <- V[train_sample, train_sample, drop=FALSE]
    
    # return the covariances
    results <- list(w=w, Vt=Vt, train=train_sample)
    return(results)
}

# ------------------------------------------------------------------------- #
# .transformPhylo                                                           #
# options: phy, model, param, ...                                           #
#                                                                           #
# ------------------------------------------------------------------------- #
.transformPhylo <- function(phy, model, param, ...){
    
    # precomputations
    n <- Ntip(phy)
    parent <- phy$edge[,1]
    descendent <- phy$edge[,2]
    extern <- (descendent <= n)
    
    switch(model,
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
    },)
    
    return(phy)
}

# ------------------------------------------------------------------------- #
# print option for effect/association of multivariate tests                 #
# options: x, digits, ...                                                   #
#                                                                           #
# ------------------------------------------------------------------------- #

print.effects.mvgls <- function(x, digits = max(3L, getOption("digits") - 3L), ...){
    
    if(x$adjusted & x$parametric==TRUE){
        tatsuoka <- attr(x$adjusted, "tatsuoka")
        cat("## Multivariate measure(s) of association ##","\n")
        if(tatsuoka) cat("## Tatsuoka' bias adjustment              ##","\n") else cat("## Serlin' bias adjustment                ##","\n")
        print(x$effect, digits=digits)
    }else{
        cat("## Multivariate measure(s) of association ##","\n")
        print(x$effect, digits=digits)
        if(x$adjusted) cat("## Note: bias is empirically adjusted     ##","\n")
    }
    if(any(x$effect<0)) message("## Values < 0 represent no association    ##","\n")
}

# ------------------------------------------------------------------------- #
# ancestral.mvgls                                                           #
# options: object, ...                                                      #
#                                                                           #
# ------------------------------------------------------------------------- #


## S3 Method for ancestral states estimation
ancestral <- function(object, ...) UseMethod("ancestral")

# core function
ancestral.mvgls <- function(object, ...){
    
    # arguments
    args <- list(...)
    
    # extract objects
    if(!inherits(object,"mvgls")){
        
        # wrapper to "estim"
        if(is.null(args[["data"]]) | is.null(args[["tree"]])) stop("Need a \"tree\" object and a new \"data\" matrix to predict ancestral states. See ?estim")
        estim(tree = args$tree, data = args$data, object = object, asr=TRUE)
        
    }else if(inherits(object,"mvgls")){
        
        # If regression, must provide a regressor for the ancestral (nodes) states
        # check if newdata is provided
        if(any(object$dims$assign>0)){
            
            if(is.null(args[["na.action"]])) na.action <- na.pass else na.action <- args$na.action
            if(is.null(args[["newdata"]])) stop("Regression model. You must provide a new \"dataset\" of predictors for each nodes. See also ?predict")
            
            Terms <- delete.response(object$terms)
            # as in "stats v3.3.0"
            m <- model.frame(Terms, args$newdata, xlev = object$xlevels, na.action = na.action)
            
            # check the arguments
            if(!is.null(cl <- attr(Terms, "dataClasses"))) .checkMFClasses(cl, m)
            X <- model.matrix(Terms, m, contrasts.arg = object$contrasts)
            predicted_fit <- X %*% object$coefficients
            
        }else{
            # Just use the grand mean - i.e. the ancestral states at the root
            predicted_fit <- object$variables$X %*% object$coefficients
            predicted_fit <- predicted_fit[-1,,drop=TRUE]
        }
        
        # start estimating ancestral states using GLS
        n <- object$dims$n
        p <- object$dims$p
        
        # covariance for the nodes
        if(!is.null(object$corrSt$diagWeight)){
            V<-.Call("mvmorph_covar_ou_fixed", A=.vcvPhyloInternal(object$variables$tree), alpha=as.double(object$param), sigma=1, PACKAGE="mvMORPH")
        }else{
            V <- .vcvPhyloInternal(object$corrSt$phy)
        }
        indice <- (1:n)
        AY <- V[-indice,indice]
        vY <- V[indice,indice]
        
        # states at the nodes
        residuals_fit <- object$residuals
        recons_t <- (AY%*%pseudoinverse(vY)%*%residuals_fit)+predicted_fit
        colnames(recons_t) = colnames(object$variables$Y)
        rownames(recons_t) = paste("node_",n+1:Nnode(object$variables$tree), sep="")
        #class(recons_t) = "anc.mvgls"
        return(recons_t)
        
    }else{
        stop("only works with \"mvgls\" class objects. See ?mvgls, or use instead \"estim\" function")
    }
    
}

