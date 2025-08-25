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
    if(is.null(args[["randomRoot"]])) randomRoot <- TRUE else randomRoot <- args$randomRoot
    if(is.null(args[["root"]])) root <- "stationary" else root <- args$root
    if(root=="stationary") root_std <- 1L else root_std <- 0L
    if(is.null(args[["hessian"]])) hessian <- "FALSE" else hessian <- args$hessian
    if(is.null(args[["MMSE"]])) MMSE <- "TRUE" else MMSE <- args$MMSE
    if(is.null(args[["FCI"]])) FCI <- "FALSE" else FCI <- args$FCI
    if(is.null(args[["comp_ll"]])) comp_ll <- "TRUE" else comp_ll <- args$comp_ll
    
    # check for coercion issues
    data_format = sapply(data, function(x) inherits(x,"phylo"))
    if(any(data_format)){
        index <- which(data_format==TRUE)
        data[[index]] <- NULL
    }
    
    # retrieve data and formula as in lm
    model_fr = model.frame(formula=formula, data=data)
    X = model.matrix(attr(model_fr, "terms"), data=model_fr, contrasts.arg=contrasts.def)
    X.formula = X
    Y = model.response(model_fr)
    assign <- attr(X, "assign")
    terms <- attr(model_fr, "terms")
    xlevels <- .getXlevels(terms, data)
    contrasts <- attr(X, "contrasts")
    
    # Option for bootstrap and permutation method
    if(!is.null(args[["response"]])) Y <- args$response
    
    # Warnings & checks
    method = match.arg(method[1], c("PL-LOOCV","LOOCV","LL","H&L","Mahalanobis","EmpBayes"))
    if(method=="PL-LOOCV") method = "LOOCV" # to keep the explicit name with 'PL'
    if(method=="EmpBayes" | penalty=="EmpBayes") penalty = "EmpBayes"; # To use the Empirical Bayes approach
    if(penalty=="EmpBayes") method = "EmpBayes";
    if(missing(tree)) stop("Please provide a phylogenetic tree of class \"phylo\" ")
    if(!inherits(tree, "simmap") & (model=="BMM" | model=="OUM")) stop("Please provide a phylogenetic tree of class \"simmap\" for the \"BMM\" and \"OUM\" models")
   # if(any(tree$edge.length<=.Machine$double.eps)) warning("There are zero branch lengths in the supplied tree. This may cause numerical issues")
    if(any(is.na(Y))) stop("Sorry, the PL approach do not handle yet missing cases.")
    if(missing(model)) stop("Please provide a model (e.g., \"BM\", \"OU\", \"EB\", \"BMM\", \"OUM\" or \"lambda\" ")
    if(ncol(as.matrix(Y))==1) stop("mvgls can be used only with multivariate datasets. See \"gls\" function in \"nlme\" or \"phylolm\" package instead.")
    if(!penalty%in%c("RidgeArch","RidgeAlt","LASSO","EmpBayes")) stop("The penalization method must be \"RidgeArch\", \"RidgeAlt\", \"EmpBayes\", or \"LASSO\"");
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
   # if(!is.ultrametric(tree) & model=="OU" & !method%in%c("LOOCV","LL")) warning("The nominal LOOCV method should be preferred with OU on non-ultrametric trees.\n")
    if(isTRUE(mserr) & model=="lambda") warning("Pagel's lambda and measurement error cannot be distinguished.\n")
    if(isFALSE(mserr)) mserr = NULL
    if(isTRUE(FCI)) hessian = TRUE
   
    # further checks
    qrx <- qr(X)
    fullrank = ifelse(ncol(X)==qrx$rank, TRUE, FALSE)
    if(!fullrank) warning("The design matrix is not of full rank. The dimensionality has been reduced by ",ncol(X)-qrx$rank,". Check your results carefully. \n")
    if(!fullrank) assign <- assign[qrx$pivot[1L:qrx$rank]]
    X <- X[,qrx$pivot[1L:qrx$rank], drop=FALSE] # FIXME: be less strict and handle case specific issues?
    
    # dimensions
    n = nobs = nrow(Y)
    p = ncol(Y)
    m = qrx$rank # dim of predictors
    nloo = 1:n
    
    # Miscellanous - pre-calculations | TODO handle time series
    precalc = .prepModel(tree, model, root)
    precalc$randomRoot = randomRoot
    precalc$root_std = root_std

    if(inherits(tree, "simmap")){
        if(model=="BMM") k <- ncol(tree$mapped.edge)
        if(model=="OUM") m <- ncol(tree$mapped.edge)
        # TODO handle cases with covariate for OUM
    }else k <- NULL
    if(method=="LL") penalized=FALSE else penalized=TRUE
    if(n<p & method=="LL") stop("There are more variables than observations. Please try instead the penalized methods \"RidgeArch\", \"RidgeAlt\" or \"LASSO\"")
    
    if(root_std==0 & model=="OUM") m = m + 1
    if(isTRUE(REML)) ndimCov = n - m else ndimCov = n
    
    # CorrStruct object (include data, model, covariance...)
    if(scale.height) tree <- .scaleStruct(tree)
    corrModel <- list(Y=Y, X=X, REML=REML, mserr=mserr,
                    model=model, structure=tree, p=p, nobs=nobs, m=m,
                    nloo=nloo, precalc=precalc)
    
    # Set bounds for parameter search
    bounds <- corrModel$bounds <- .setBounds(penalty=penalty, model=model, lower=low, upper=up, tol=tol, mserr=mserr, penalized=penalized, corrModel=corrModel, k=k)
    
    # Starting values & parameters ID
    if(grid_search & is.null(start)){
      
        start <- .startGuess(corrModel, cvmethod=method, mserr=mserr, target=target, penalty=penalty, echo=echo, penalized, tol=tol)
        
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
                        nobs=nobs,
                        hessian=hessian)
    
    # Estimates
    tuning <- bounds$trTun(estimModel$par)
    mod_par <- bounds$trPar(estimModel$par)
    
    # convergence & bounds checks?
    .check_par_results(corrModel, mod_par, penalized);
    
    if(!is.null(mserr)) corrModel$mserr <- mserr_par <- bounds$trSE(estimModel$par) else mserr_par <- NA
    ll_value <- -estimModel$value # either the loocv or the regular likelihood (minus because we minimize)
    const_mtdist = NA # Not efficient
    if(method=="EmpBayes" & comp_ll==TRUE){
        v = p + 1 # This needs to be modified if it is allowed to select different values for v
        const_mtdist = -(lmvgamma((v+ndimCov+p-1)/2, p) - lmvgamma((v+p-1)/2, p) - 0.5*(ndimCov*p*log(pi)))
        ll_value = -(estimModel$value + const_mtdist)
    }
    
    # Exceptions [to improve]
    X <- .make.x(tree, mod_par, X, model, root, root_std)
    
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
    
    # add an option to avoid the computation of the covariance matrix with EIC - EmpBayes method
    if(MMSE==FALSE & method=="EmpBayes"){
        R <- NULL
    } else{
        S <- crossprod(residuals)/ndimCov
        R <- .penalizedCov(S, penalty=ifelse(method=="LL", method, penalty), targM=target, tuning=tuning, n=ndimCov) # change n for reml in EmpBayes
    }
    
    # Computing confidence intervals using the Fisher information matrix. Only available for method='EmpBayes'
    if(FCI){
        if(method!="EmpBayes") warning("The CI are available only for the \"EmpBayes\" method") else if(method=="EmpBayes"){
            if(model=="BM") fci = NA else{
                fisher_information<-solve(estimModel$hessian) # it's already the negative of the Hessian since the negative ll is minimized
                sigma_for_ci<-sqrt(diag(fisher_information)) # the first entry is for the "regularization" term. For ML optimization it will be the first term
                
                # Does this work for BMM model or should it be constrain?
                upper_ci<-mod_par+1.96*sigma_for_ci[2]
                lower_ci<-mod_par-1.96*sigma_for_ci[2]
                fci = c('lw'=lower_ci, 'up'=upper_ci)
            }
        }
    }
    
    # Multiple rates BMM - we scale the average rate (mean of the diagonal of the covariance matrix)
    if(inherits(tree, "simmap") && model=="BMM"){
        if(method=="EmpBayes"){
            empvar <- colSums(residuals^2)/ndimCov
            # in the PL we use the regularized estimate because it shrinks (for most of them) toward the mean or the diagonal of the empirical. The regularization of the MMSE Emp. Bayes will not.
            avg_rate <- mean(empvar) # mean((empvar*ndimCov + tuning*mean(empvar))/(v+ndimCov-2))
        }else{
            avg_rate <- mean(diag(R$Pinv)) # we can make it more efficient for the other approaches, but should we? It's better to use the EmpBayes on very high-dimensionnal datasets
        }
        mod_par <- c(avg_rate, avg_rate*mod_par)
        names(mod_par) <- attr(tree$mapped.edge,"dimnames")[[2]] # set the names of the groups for BMM. we remove the first one which is used as reference
    }
    
    # number of dimensions
    ndims <- list(n=n, p=p, m=m, assign=assign, rank=qrx$rank, pivot=qrx$pivot, fullrank=fullrank)
    if(model=='OUM') variables <- list(Y=Y, X=X, tree=tree, regimes=X.formula) else variables <- list(Y=Y, X=X, tree=tree)
    # End
    if(echo==TRUE) message("Done in ", numIter," iterations.")
    
    # Return the results
    results = list(formula=formula,
        call = match.call(),
        coefficients=coefficients,
        terms=terms,
        xlevels=xlevels,
        contrasts=contrasts,
        variables=variables,
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
        FCI=if(isTRUE(FCI)) fci else NA,
        const_mtd = const_mtdist,
        opt=estimModel)
    
    class(results) <- "mvgls"
    return(results)
}
