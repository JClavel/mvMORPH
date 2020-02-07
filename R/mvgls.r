################################################################################
##                                                                            ##
##                       mvMORPH: mvgls.r                                     ##
##                                                                            ##
##   Multivariate Generalized Least Squares Linear Models by ML and PL        ##
##                                                                            ##
##  Created by Julien Clavel - 31-07-2018                                     ##
##  (julien.clavel@hotmail.fr/ julien.clavel@biologie.ens.fr)                 ##
##   require: phytools, ape, corpcor, subplex, spam, glassoFast, stats        ##
##                                                                            ##
################################################################################


mvgls <- function(formula, data=list(), tree, model, method=c("PL-LOOCV","LL"), REML=TRUE, ...){
    
    # Recover options
    args <- list(...)
    if(is.null(args[["scale.height"]])) scale.height <- FALSE else scale.height <- args$scale.height
    if(is.null(args[["echo"]])) echo <- FALSE else echo <- args$echo
    if(is.null(args[["grid.search"]])) grid_search <- TRUE else grid_search <- args$grid.search
    if(is.null(args[["target"]])) target <- "unitVariance" else target <- args$target
    if(is.null(args[["error"]])) mserr <- NULL else mserr <- args$error
    if(is.null(args[["penalty"]])) penalty <- "RidgeArch" else penalty <- args$penalty
    if(is.null(args[["optimization"]])) optimization <- "L-BFGS-B" else optimization <- args$optimization
    if(is.null(args[["ncores"]])) ncores <- 1L else ncores <- args$ncores
    if(is.null(args[["upper"]])) up <- NULL else up <- args$upper
    if(is.null(args[["lower"]])) low <- NULL else low <- args$lower
    if(is.null(args[["tol"]])) tol <- NULL else tol <- args$tol
    if(is.null(args[["start"]])) start <- NULL else start <- args$start
    if(is.null(args[["contrasts"]])) contrasts.def <- NULL else contrasts.def <- args$contrasts
    
    # check for coercion issues
    data_format = sapply(data, function(x) inherits(x,"phylo"))
    if(any(data_format)){
        index <- which(data_format==TRUE)
        data[[index]] <- NULL
    }
    
    # retrieve data and formula as in lm
    model_fr = model.frame(formula=formula, data=data)
    X = model.matrix(attr(model_fr, "terms"), data=model_fr, contrasts.arg=contrasts.def)
    Y = model.response(model_fr)
    assign <- attr(X, "assign")
    
    # Option for bootstrap and permutation method
    if(!is.null(args[["response"]])) Y <- args$response
    
    # Warnings & checks
    method = match.arg(method[1], c("PL-LOOCV","LOOCV","LL","H&L","Mahalanobis"))
    if(method=="PL-LOOCV") method = "LOOCV" # to keep the explicit name with 'PL'
    if(missing(tree)) stop("Please provide a phylogenetic tree of class \"phylo\" ")
    if(any(is.na(Y))) stop("Sorry, the PL approach do not handle yet missing cases.")
    if(missing(model)) stop("Please provide a model (e.g., \"BM\", \"OU\", \"EB\", or \"lambda\" ")
    if(ncol(as.matrix(Y))==1) stop("mvgls can be used only with multivariate datasets. See \"gls\" function in \"nlme\" or \"phylolm\" package instead.")
    if(!penalty%in%c("RidgeArch","RidgeAlt","LASSO")) stop("The penalization method must be \"RidgeArch\", \"RidgeAlt\" or \"LASSO\"");
    if(!target%in%c("unitVariance","Variance","null")) warning("Default target are \"unitVariance\", \"null\" or \"Variance\". Check the target matrix provided");
    if(nrow(model_fr)!=length(tree$tip.label)) stop("number of rows in the data does not match the number of tips in the tree.")
    if (all(rownames(model_fr) %in% tree$tip.label)){ # to be changed for TS
        Y <- Y[tree$tip.label,,drop=FALSE]
        X <- X[tree$tip.label,,drop=FALSE]
    }else if(is.null(args[["response"]])){
        warning("the data has no names, order assumed to be the same as tip labels in the tree.\n")
    }
    if(!inherits(tree, "phylo")) stop("object \"tree\" is needed if no custom correlation structure provided.")
    if(method%in%c("H&L","Mahalanobis") & penalty%in%c("RidgeAlt","LASSO")) stop("\"H&L\" and \"Mahalanobis\" works only with \"RidgeArch\" penalization")
    if(!is.ultrametric(tree) & model=="OU" & !method%in%c("LOOCV","LL")) warning("The nominal LOOCV method should be preferred with OU on non-ultrametric trees.\n")
    
    # Miscellanous - pre-calculations
    precalc = .prepModel(tree, model)
    
    # dimensions
    n = nobs = nrow(Y)
    p = ncol(Y)
    m = ncol(X) # dim of predictors
    nloo = 1:n
    if(REML) ndimCov = n - m else ndimCov = n
    if(method=="LL") penalized=FALSE else penalized=TRUE
    if(n<p & method=="LL") stop("There are more variables than observations. Please try instead the penalized methods \"RidgeArch\", \"RidgeAlt\" or \"LASSO\"")
    
    # Set bounds for parameter search
    bounds <- .setBounds(penalty=penalty, model=model, lower=low, upper=up, tol=tol, mserr=mserr, penalized=penalized)
    
    # CorrStruct object (include data, model, covariance...)
    if(scale.height) tree <- .scaleStruct(tree)
    corrModel <- list(Y=Y, X=X, REML=REML, mserr=mserr,
                    model=model, structure=tree, p=p, nobs=nobs,
                    nloo=nloo, bounds=bounds, precalc=precalc)
    
    # Starting values & parameters ID
    if(grid_search & is.null(start)){
      
        start <- .startGuess(corrModel, cvmethod=method, mserr=mserr, target=target, penalty=penalty, echo=echo, penalized)
        
    }else if(is.null(start)){
        if(method=="LL" | model=="BM") start <- 0.5 else start <- c(0.5,0.5)
        if(!is.null(mserr)) start <- c(start,1e-4)
    }
    
    # Optimization
    if(echo==TRUE) message("Start optimization. Please wait...")
    estimModel <- optim(start,
                        fn = .loocvPhylo,
                        method=optimization,
                        upper=bounds$upper,
                        lower=bounds$lower,
                        cvmethod=method,
                        targM=target,
                        corrStr=corrModel,
                        penalty=penalty,
                        error=mserr,
                        nobs=nobs)
    
    # Estimates
    tuning <- bounds$trTun(estimModel$par)
    mod_par <- bounds$trPar(estimModel$par)
    if(!is.null(mserr)) corrModel$mserr <- mserr_par <- bounds$trSE(estimModel$par) else mserr_par <- NA
    ll_value <- -estimModel$value # either the loocv or the regular likelihood (minus because we minimize)
    
    
    # List of results to return
    corrSt = .corrStr(mod_par, corrModel);
    par_estimates <- .mvGLS(corrSt)
    residuals <- par_estimates$residuals # normalized residuals
    coefficients <- par_estimates$B
    fitted.values <- X%*%coefficients # raw coefficients
    residuals_raw <- Y - fitted.values
    call <- formula
    model.frame <- model_fr
    glsStruct <- corrSt
    method <- method
    numIter <- estimModel$count[1]
    S <- crossprod(residuals)/ndimCov
    R <- .penalizedCov(S, penalty=ifelse(method=="LL", method, penalty), targM=target, tuning=tuning)
    ndims <- list(n=n, p=p, m=m, assign=assign)
    
    # Exceptions [to improve]
    X <- .make.x(tree, mod_par, X, model)
    
    # End
    if(echo==TRUE) message("Done in ", numIter," iterations.")
    
    # Return the results
    results = list(formula=formula,
        call = match.call(),
        coefficients=coefficients,
        variables=list(Y=Y, X=X, tree=tree),
        dims=ndims,
        fitted=fitted.values,
        logLik=ll_value,
        method=method,
        model=model,
        numIter=numIter,
        residuals=residuals_raw,
        sigma=R,
        tuning=if(method=="LL") NA else tuning,
        param=if(model=="BM") NA else mod_par,
        mserr=mserr_par,
        start_values=start,
        corrSt=corrSt,
        penalty=if(method=="LL") "LL" else penalty,
        target=if(method=="LL") "LL" else target,
        REML=REML,
        opt=estimModel)
    
    class(results) <- "mvgls"
    return(results)
}
