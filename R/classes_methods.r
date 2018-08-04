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
        isIntercept <- attr(terms(formula(object)),"intercept")
        covBeta <- kronecker(object$sigma$Pinv, XtX)
        if(isIntercept){
            rownames(covBeta) <- colnames(covBeta) <- rep(c("(Intercept)",attr(terms(formula(object)),"term.labels")), object$dims$p)
        }else{
            rownames(covBeta) <- colnames(covBeta) <- rep(attr(terms(formula(object)),"term.labels"), object$dims$p)
        }
        return(covBeta)})
}

# ------------------------------------------------------------------------- #
# coef.mvgls     / coefficients.mvgls                                       #
# options: object, ...                                                      #
# S3 method from "stats" package                                            #
# ------------------------------------------------------------------------- #
coef.mvgls <- function(object, ...){
    
    coeffs <- object$coefficients
    isIntercept <- attr(terms(formula(object)),"intercept")
    
    if(isIntercept){
        rownames(coeffs) <- c("(Intercept)",attr(terms(formula(object)),"term.labels"))
    }else{
        rownames(coeffs) <- attr(terms(formula(object)),"term.labels")
    }
    
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
        if(x$REML) cat("Log-restricted-likelihood:",round(-x$logLik, digits=digits), "\n\n") else cat("Log-likelihood:",round(-x$logLik, digits=digits), "\n\n")
    }else{
        cat("\nGeneralized least squares fit by penalized",meth,"\n")
        if(x$REML){
            cat("LOOCV of the log-restricted-likelihood:",round(-x$logLik, digits=digits), "\n\n")
        }else{
            cat("LOOCV of the log-likelihood:",round(-x$logLik, digits=digits), "\n\n")
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
        LL = -object$logLik
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


