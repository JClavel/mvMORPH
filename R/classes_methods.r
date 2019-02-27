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
# GIC.mvgls                                                                 #
# options: model,...                                                        #
# S3 method from "RPANDA"  package                                          #
# ------------------------------------------------------------------------- #
GIC.mvgls <- function(object, ...){
    
    # retrieve arguments
    args <- list(...)
    if(is.null(args[["eigSqm"]])) eigSqm <- TRUE else eigSqm <- args$eigSqm
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
      mod.par=ncol(object$corrSt$phy$mapped.edge)
    }else{
      mod.par=1
    }
    if(is.numeric(object$mserr)) mod.par = mod.par + 1 # already included in the covariance matrix structure?
    if(object$REML) ndimCov = n - m else ndimCov = n
    # Nominal loocv
    XtX <- solve(crossprod(X))
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
    XtX <- solve(t(X)%*%(X))
    T2 <- sapply(nloo, function(i){
        gradient <- (X[i,])%*%t(P%*%t(Y[i,]-X[i,]%*%beta))
        sum(gradient * (XtX%*%gradient%*%Pi))
    })
    beta_df <- sum(T2)
    
    if(m>1) warning("GIC criterion with multiple predictors has not been fully tested. Please use it with cautions and consider simulations instead")
    
    # LogLikelihood (minus)
    DP <- as.numeric(determinant(Pi)$modulus)
    Ccov <- object$corrSt$det
    llik <- 0.5 * (ndimCov*p*log(2*pi) + p*Ccov + ndimCov*DP + ndimCov*sum(S*P))
    GIC <- 2*llik + 2*(sigma_df+beta_df+mod.par)
    
    # return the results
    results <- list(LogLikelihood=-llik, GIC=GIC, p=p, n=n, bias=sigma_df+beta_df+mod.par, bias_cov=sigma_df)
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
    
    # retrieve data to simulate bootstrap samples
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
    
    N = nrow(Y)
    p = object$dims$p
    if(object$REML) ndimCov = object$dims$n - object$dims$m else ndimCov = object$dims$n
    tuning <- object$tuning
    target <- object$target
    penalty <- object$penalty
    Dsqrt <- pruning(object$corrSt$phy, trans=FALSE, inv=FALSE)$sqrtM # return warning message if n-ultrametric tree is used with OU?
    DsqrtInv <- pruning(object$corrSt$phy, trans=FALSE, inv=TRUE)$sqrtM
    modelPerm <- object$call
    modelPerm$grid.search <- quote(FALSE)
    modelPerm$start <- quote(object$opt$par)
    
    # Mean and residuals for the model
    MeanNull <- object$variables$X%*%beta
    
    # Generate bootstrap samples
    boots <- mclapply(1:nboot, function(i){
        #Yp <- MeanNull + Dsqrt%*%(residuals[c(sample(N-1, replace=TRUE),N),]) # sampling with replacement for bootstrap
        Yp <- MeanNull + Dsqrt%*%(residuals[sample(N, replace=TRUE),]) # sampling with replacement for bootstrap
        rownames(Yp) <- rownames(object$variables$Y)
        Yp
    }, mc.cores = getOption("mc.cores", nbcores))
    
    # Estimate parameters on bootstrap samples
    estimPar <- mclapply(1:nboot, function(i){
        temp_data <- boots[[i]];
        modelPerm$response <- quote(temp_data);
        estimModelNull <- eval(modelPerm);
        estimModelNull
        
    }, mc.cores = getOption("mc.cores", nbcores))
    
    # Estimate the bias term
    D1 <- function(objectBoot, objectFit, ndimCov, p, sqM){ # LL(Y*|param*) - LL(Y*| param)
        
        # Y*|param*
        residualsBoot <- residuals(objectBoot, type="normalized")
        
        # For boot "i" LL1(Y*|param*)
        Ccov1 <- as.numeric(objectBoot$corrSt$det)
        Gi1 <- try(chol(objectBoot$sigma$Pinv), silent=TRUE)
        if(inherits(Gi1, 'try-error')) return("error")
        quadprod <- sum(backsolve(Gi1, t(residualsBoot), transpose = TRUE)^2)
        detValue <- sum(2*log(diag(Gi1)))
        llik1 <- -0.5 * (ndimCov*p*log(2*pi) + p*Ccov1 + ndimCov*detValue + quadprod)
        
        # Y*|param
        #if(!restricted) residualsBoot <- objectBoot$corrSt$Y - objectBoot$corrSt$X%*%objectFit$coefficients
        if(!restricted) residualsBoot <- crossprod(sqM, objectBoot$variables$Y - objectBoot$variables$X%*%objectFit$coefficients)
        
        # For boot "i" LL2(Y*|param)
        Ccov2 <- as.numeric(objectFit$corrSt$det)
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
            sqM_temp <- pruning(objectBoot$corrSt$phy, trans=FALSE, inv=TRUE)$sqrtM
            residualsBoot <- try(crossprod(sqM_temp, objectFit$variables$Y - objectFit$variables$X%*%objectBoot$coefficients), silent=TRUE)
             
        }else{ residualsBoot <- objectFit$corrSt$Y - objectFit$corrSt$X%*%objectFit$coefficients}
        
        #if(!restricted) residualsBoot <- objectFit$corrSt$Y - objectFit$corrSt$X%*%objectBoot$coefficients
        #else residualsBoot <- objectFit$corrSt$Y - objectFit$corrSt$X%*%objectFit$coefficients
        
        # For boot "i" LL2(Y|param*)
        Ccov1 <- as.numeric(objectBoot$corrSt$det)
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
    Ccov <- as.numeric(object$corrSt$det)
    Gi <- try(chol(object$sigma$Pinv), silent=TRUE)
    if(inherits(Gi, 'try-error')) return("error")
    quadprod <- sum(backsolve(Gi, t(residuals), transpose = TRUE)^2)
    detValue <- sum(2*log(diag(Gi)))
    llik <- -0.5 * (ndimCov*p*log(2*pi) + p*Ccov + ndimCov*detValue + quadprod)
    
    # Estimate the bias term (variance reduction method)
    D1_simul <- sapply(estimPar, function(x) D1(objectBoot=x, objectFit=object, ndimCov=ndimCov, p=p, sqM=DsqrtInv) )
    D3_simul <- sapply(estimPar, function(x) D3(objectBoot=x, objectFit=object, loglik=llik, ndimCov=ndimCov, p=p) )
    
    # combine the values for D1 and D3
    bias <- D1_simul+D3_simul
    
    # compute the EIC
    EIC <- -2*llik + 2*mean(bias)
    
    # concatenate the results
    results <- list(EIC=EIC, bias=bias, LogLikelihood=llik, boot=boots, estimPar=estimPar, se=sd(bias)/sqrt(nboot), p=p, n=N)
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
        m <- object$dims$m
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
        cat("\nGeneralized least squares fit by",meth,"\n")
        if(x$REML) cat("Log-restricted-likelihood:",round(x$logLik, digits=digits), "\n\n") else cat("Log-likelihood:",round(x$logLik, digits=digits), "\n\n")
    }else{
        cat("\nGeneralized least squares fit by penalized",meth,"\n")
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
        #"BMM"={cat(names(x$param),"\n",round(x$param, digits=digits),"\n\n")},
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
        cat("\nGeneralized least squares fit by",meth,"\n")
        print(x$results.fit,  quote = FALSE )
    }else{
        
        cat("\nGeneralized least squares fit by penalized",meth,"\n")
        print(x$results.fit,  quote = FALSE )
    }
    
    
    # Model parameters
    cat("\nParameter estimate(s):\n")
    if(!any(is.na(x$param))){
        switch(x$model,
        "OU"={ cat("alpha:",round(x$param, digits=digits),"\n\n")},
        "EB"={ cat("r:",round(x$param, digits=digits),"\n\n")},
        "lambda"={cat("lambda:",round(x$param, digits=digits),"\n\n")},
        #"BMM"={cat(names(x$param),"\n",round(x$param, digits=digits),"\n\n")},
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
    m <- object$dims$m
    if(object$REML) ndimCov = n - m else ndimCov = n
    
    # loocv or LL
    meth <- ifelse(object$REML, "REML", "ML")
    
    if(object$method=="LL"){
        LL = object$logLik
        nparam = sum(!is.na(object$param)) + p + p*(p + 1)/2
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
    class(object) <- c("summary.mvgls","mvgls")
    object
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
    cat("EIC:",x$EIC,"| SE",x$se, "| Log-likelihood",x$LogLikelihood,"\n")
    cat("\n")
}

