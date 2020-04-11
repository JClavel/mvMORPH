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
    
    if(object$model=="BM") mod.par=0 else mod.par=1
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
    
    if(m>1) warning("GIC criterion with multiple predictors has not been fully tested. Please use it with cautions and consider EIC or simulations instead")
    
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
    Dsqrt <- pruning(object$corrSt$phy, trans=FALSE, inv=FALSE)$sqrtM # return warning message if n-ultrametric tree is used with OU?
    # TODO (change to allow n-ultrametric and OU
    if(object$model=="OU" & !is.ultrametric(object$variables$tree)) stop("The EIC method does not handle yet non-ultrametric trees with OU processes")
    
    DsqrtInv <- pruning(object$corrSt$phy, trans=FALSE, inv=TRUE)$sqrtM
    modelPerm <- object$call
    modelPerm$grid.search <- quote(FALSE)
    modelPerm$start <- quote(object$opt$par)
    
    # Mean and residuals for the model
    MeanNull <- object$variables$X%*%beta
    
    
    # Estimate the bias term
    D1 <- function(objectBoot, objectFit, ndimCov, p, sqM){ # LL(Y*|param*) - LL(Y*| param)
        
        # Y*|param*
        residualsBoot <- residuals(objectBoot, type="normalized")
        
        # For boot "i" LL1(Y*|param*)
        if(objectFit$REML==TRUE & args$forceREML==FALSE) Ccov1 <- as.numeric(objectBoot$corrSt$det - determinant(crossprod(objectBoot$corrSt$X))$modulus) else Ccov1 <- as.numeric(objectBoot$corrSt$det)
        Gi1 <- try(chol(objectBoot$sigma$Pinv), silent=TRUE)
        if(inherits(Gi1, 'try-error')) return("error")
        quadprod <- sum(backsolve(Gi1, t(residualsBoot), transpose = TRUE)^2)
        detValue <- sum(2*log(diag(Gi1)))
        llik1 <- -0.5 * (ndimCov*p*log(2*pi) + p*Ccov1 + ndimCov*detValue + quadprod)
        
        # Y*|param
        #if(!restricted) residualsBoot <- objectBoot$corrSt$Y - objectBoot$corrSt$X%*%objectFit$coefficients # does not account for the phylo model of the original fit
        if(!restricted) residualsBoot <- crossprod(sqM, objectBoot$variables$Y - objectBoot$variables$X%*%objectFit$coefficients)
        
        # For boot "i" LL2(Y*|param)
        if(objectFit$REML==TRUE & args$forceREML==FALSE) Ccov2 <- as.numeric(objectFit$corrSt$det - determinant(crossprod(objectFit$corrSt$X))$modulus) else Ccov2 <- as.numeric(objectFit$corrSt$det)
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
        if(objectFit$REML==TRUE & args$forceREML==FALSE) Ccov1 <- as.numeric(objectBoot$corrSt$det - determinant(crossprod(objectBoot$corrSt$X))$modulus) else Ccov1 <- as.numeric(objectBoot$corrSt$det)
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
    if(object$REML==TRUE & args$forceREML==FALSE) Ccov <- as.numeric(object$corrSt$det - determinant(crossprod(object$corrSt$X))$modulus) else Ccov <- as.numeric(object$corrSt$det)
    Gi <- try(chol(object$sigma$Pinv), silent=TRUE)
    if(inherits(Gi, 'try-error')) return("error")
    quadprod <- sum(backsolve(Gi, t(residuals), transpose = TRUE)^2)
    detValue <- sum(2*log(diag(Gi)))
    llik <- -0.5 * (ndimCov*p*log(2*pi) + p*Ccov + ndimCov*detValue + quadprod)
    
    # Estimate parameters on bootstrap samples
    bias <- pbmcmapply(function(i){
        
        # generate bootstrap sample
        Yp <- MeanNull + Dsqrt%*%(residuals[sample(N, replace=TRUE),]) # sampling with replacement for bootstrap
        rownames(Yp) <- rownames(object$variables$Y)
        
        modelPerm$response <- quote(Yp);
        estimModelNull <- eval(modelPerm);
        d1res <- D1(objectBoot=estimModelNull, objectFit=object, ndimCov=ndimCov, p=p, sqM=DsqrtInv)
        d3res <- D3(objectBoot=estimModelNull, objectFit=object, loglik=llik, ndimCov=ndimCov, p=p)
        d1res+d3res
    }, 1:nboot, mc.cores = getOption("mc.cores", nbcores))
    
    # compute the EIC
    pboot <- mean(bias)
    EIC <- -2*llik + 2*pboot
    
    # standard-error
    se <- sd(bias)/sqrt(nboot)
    
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
    if(!is.na(x$param)){
        switch(x$model,
        "OU"={ cat("alpha:",round(x$param, digits=digits),"\n\n")},
        "EB"={ cat("r:",round(x$param, digits=digits),"\n\n")},
        "lambda"={cat("lambda:",round(x$param, digits=digits),"\n\n")},
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
    if(!is.na(x$param)){
        switch(x$model,
        "OU"={ cat("alpha:",round(x$param, digits=digits),"\n\n")},
        "EB"={ cat("r:",round(x$param, digits=digits),"\n\n")},
        "lambda"={cat("lambda:",round(x$param, digits=digits),"\n\n")},
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
        nparam = length(object$start_values) + p + p*(p + 1)/2
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
    
    signif <- sapply(x$pvalue, function(i) if(i<0.001){"***"}else if(i<0.01){
      "**"}else if(i<0.05){"*"}else if(i<0.1){"."}else{""})
    
    table_results <- data.frame(Df=x$Df, stat=x$stat, approxF=x$approxF, numDf=x$NumDf, denDf=x$DenDf, pval=x$pvalue, signif=signif)
    if(x$type!="glh") rownames(table_results) <- x$terms else rownames(table_results) <- "Contrasts L"
    colnames(table_results) <- c("Df", "test stat", "approx F", "num Df", "den Df", "Pr(>F)", "")
    print(table_results, digits = digits, ...)
    cat("---","\n")
    cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1")
    
    
    
  }else{ # permutation methods
    
    if(x$type=="I") cat("Sequential MANOVA Tests with",x$nperm,"permutations:",x$test,"test statistic","\n")
    if(x$type=="II") cat("Type II MANOVA Tests with",x$nperm,"permutations:",x$test,"test statistic","\n")
    if(x$type=="III") cat("Type III MANOVA Tests with",x$nperm,"permutations:",x$test,"test statistic","\n")
    if(x$type=="glh")  cat("General Linear Hypothesis Test with",x$nperm,"permutations:",x$test,"test statistic","\n")
 
    signif <- sapply(x$pvalue, function(i) if(i<0.001){"***"}else if(i<0.01){
      "**"}else if(i<0.05){"*"}else if(i<0.1){"."}else{""})
    
    table_results <- data.frame(stat=x$stat, pval=x$pvalue, signif=signif)
    if(x$type!="glh") rownames(table_results) <- x$terms else rownames(table_results) <- "Contrasts L"
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

