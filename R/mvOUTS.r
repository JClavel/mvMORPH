################################################################################
##                                                                            ##
##                               mvMORPH: mvTIME/mvTS                         ##
##                                                                            ##
##  Created by Julien Clavel - 04-12-2015                                     ##
##  (julien.clavel@hotmail.fr/ clavel@biologie.ens.fr)                        ##
##   require: phytools, ape, corpcor, spam                                    ##
##                                                                            ##
################################################################################


mvOUTS <- function(times, data, error=NULL, param=list(sigma=NULL,alpha=NULL, vcv="randomRoot", decomp=c("cholesky","spherical","eigen","qr","diagonal","upper","lower")),method=c("rpf","inverse","pseudoinverse","univarpf"),scale.height=FALSE, optimization=c("L-BFGS-B","Nelder-Mead","subplex"),control=list(maxit=20000),precalc=NULL, diagnostic=TRUE, echo=TRUE){
    
    if(missing(times)) stop("The time-series object is missing!")
    if(missing(data)) stop("You must provide a dataset along with your time-series!")
    
    # Param
    n <- length(times)
    k <- 1 # un seul groupe pour le moment
    
    # number of traits
    if(is.vector(data)){
        p<-1 }else{ p<-ncol(data)}
    
    # Set option for maximizing the likelihood in the optimization method
    # control$fnscale=-1
    
    # bind data to a vector
    if(is.matrix(data)){
        dat<-as.vector(data) }else{ dat<-as.vector(as.matrix(data))}
    
    # Match method
    # method <- match.arg(method)
    method <- method[1]
    optimization <- optimization[1]
    
    # Check if there is missing cases
    NA_val<-FALSE
    Indice_NA<-NULL
    if(any(is.na(data))){
        if(method!="pic" & method!="sparse"){
            NA_val<-TRUE
        }else{
            stop("NA values are allowed only with the \"rpf\",\"inverse\" or \"pseudoinverse\" methods")
        }
        Indice_NA<-which(is.na(as.vector(data)))
    }
    
    
    # Scale the time-serie between 0 and 1
    # check for external time series for the expectation?
    if(scale.height==TRUE){
        times <- times/max(times)
    }
    
    # Compute the covariance matrix for the time serie
    vcv_time <- vcv.ts(times)
    
    # decomposition of the pull matrix
    index.user<-NULL
    if(is.null(param[["decomp"]])==TRUE){
        decomp <- param$decomp <- "cholesky"
    }else if(is.matrix(param$decomp)){
        decomp <- "user"
        index.user <- param$decomp
    }else{
        decomp <- param$decomp[1]
    }
    
    # tolerance value
    if(is.null(param[["tol"]])){
        tol <- param$tol <- 0.1
    }else{
        tol <- param$tol
    }
    
    # scatter (sigma) matrix decomposition
    index.user2 <- NULL
    if(is.null(param[["decompSigma"]])){
        decompSigma<-param$decompSigma<-"cholesky"
    }else if(is.matrix(param$decompSigma)){
        decompSigma <- "user"
        index.user2 <- param$decompSigma
    }else{
        decompSigma <- param$decompSigma[1]
    }
    
    # function for alpha parameterization
    buildA <- function(x,p){matrixParam(x,p,decomp,index.user,tol)}
    
    # function for sigma parameterization
    sigmafun <- function(par) {symPar(par, decomp=decompSigma, p=p, index.user=index.user2, tol=tol)}
    
    # diag matrix
    matdiag <- diag(p)
    
    # option for computing the variance covariance matrix (sparse not yet included for mvOUTS)
    if(is.null(param[["vcv"]])==TRUE){
        if(method=="sparse"){
            vcvtype<-"sparse"
        }else{
            vcvtype<-"randomRoot"
        }
        
    }else{
        vcvtype<-param$vcv
    }
    
    # Prevent singular matrices
    if(is.null(error) & vcvtype=="fixedRoot" |  any(error[1,]==0) & vcvtype=="fixedRoot"){
        warning("No sampling error provided for the initial observation(s) while using the \"fixedRoot\" parameterization.","\n"," An arbitrary sampling error value is set to 0.01 for the initial observation(s) to avoid singularity issues during the likelihood computation. ")
        if(any(error[1,]==0)==TRUE){
            if(p==1){error[1] <- 0.01}else{error[1,] <- 0.01}
            
        }else{
            error <- matrix(0,ncol=p, nrow=n)
            error[1,] <- 0.01
        }
        error<-as.vector(error)
    }else{
        if(!is.null(error)){
         # transform to a vector
         error<-as.vector(as.matrix(error))
         error[is.na(error)] <- 0
        }
    }
    
    if(vcvtype=="fixedRoot" & p==1){ # a modifier pour integrer sparse et rpf?
        # si explicite ->
        if(method=="univarpf"){vcvtype<-"univarpfFixed"}else{ vcvtype<-"univarFixed"}
    }else if(vcvtype=="randomRoot" & p==1){
        if(method=="univarpf"){vcvtype<-"univarpfRandom"}else{ vcvtype<-"univarRandom"}
    }
    
    # starting values (need to be changed for taking into account the parameterization
    if(is.null(param[["alpha"]])){
      alpha<-param$alpha<-startParam(p,decomp,times,index.user)
    }else{
        if(is.matrix(param$alpha)){
            alpha<-param$alpha<-startParam(p,decomp,times,index.user,hlife=Re(eigen(param$alpha)$values))
        }else{
            alpha<-param$alpha
        }
    }
    
    if(is.null(param[["sigma"]])){
        empirical<-Stationary_to_scatter(buildA(alpha,p)$A, cov(as.matrix(data), use="complete.obs"))
         test <- try(chol(empirical), silent=TRUE)
         if(inherits(test ,'try-error')){
             empirical <- NULL
         }
        sigma<-param$sigma<-startParamSigma(p, decompSigma, times, data, empirical, index.user2)
        #sigma<-param$sigma<-startParamSigma(p, decompSigma, times, data)
    }else{
        if(is.matrix(param$sigma)){
            sigma<-startParamSigma(p, decompSigma, times, data, param$sigma, index.user2)
        }else{
            sigma<-param$sigma
        }
    }
    
    if(is.null(param[["theta0"]])){
        if(!is.vector(data)){
            theta0<-param$theta0<-colMeans(data, na.rm=TRUE)
        }else{
            theta0<-param$theta0<-mean(data, na.rm=TRUE)
        }
    }else{
        theta0<-param$theta0
    }
    
    if(is.null(param[["theta1"]])){
        if(!is.vector(data)){
            theta1<-param$theta1<-colMeans(data, na.rm=TRUE)
        }else{
            theta1<-param$theta1<-mean(data, na.rm=TRUE)
        }
    }else{
        theta1<-param$theta1
    }
    
    if(is.null(param[["root"]])){
        root<-param$root<-TRUE
    }else{
        root<-param$root
    }
    
    
    # number of parameters for each matrix
    nalpha<-length(alpha)
    nsigma<-length(sigma)
    ntheta0<-p
    if(root==FALSE){
    ntheta1<-0
    }else{
    ntheta1<-p
    }
    ntheta<-ntheta0+ntheta1
    
    
    # number of columns of the design matrix
    if(root==TRUE){
        sizeD<-p*2
    }else{
        sizeD<-p
        # precalculate the design matrix
        design_matrix <- multD(NULL,p,n,smean=TRUE)
    }
    
    ##-----------------------Precalc-on-------------------------------------------##
    if(is.null(precalc)==FALSE & class(precalc)=="mvmorph.precalc"){
        
        mt<-list(mDist=precalc$C1)
        root<-precalc$param$root
        if(method=="sparse"){
            # Yale sparse format
            JAr<-precalc$JAr
            IAr<-precalc$IAr
            ch<-precalc$ch
            precalcMat<-precalc$V
        }

        
    }else{
        
        # number of columns of the design matrix
        
        if(method=="sparse"){
            vcv_fixed<-vcv_time
            vcv_fixed[1,1]<-1e-5
            V<-kronecker((matrix(1,p,p)+diag(p)), vcv_fixed)
            # spam object
            precalcMat<-as.spam(V);
            # precal the cholesky
            if(is.null(param[["pivot"]])){pivot<-"MMD"}else{pivot<-param$pivot}
            ch<-chol(precalcMat,pivot=pivot)
            # Yale Sparse Format indices
            JAr<-precalcMat@colindices-1
            IAr<-precalcMat@rowpointers-1
        }else{
            ch<-NULL
            precalcMat<-NULL
        }
    }
    ##-----------------------Models matrices--------------------------------------##

    # Define the variance-covariance and weights matrix for the likelihood function
    
    switch(vcvtype,
    "sparse"={
        model_fun_matrix<-function(vcv,n,alphaA,sigmA,times,theta0=NULL,theta1=NULL,matdiag=NULL,theta_mle=TRUE){
            
             V<-.Call("mvmorph_covar_ou_sparse", A=as.double(precalcMat@entries), JA=as.integer(JAr), IA=as.integer(IAr), as.integer(n), bt=as.numeric(vcv), lambda=alphaA$values, S=alphaA$vectors, sigmasq=sigmA, S1=alphaA$invectors)
            if(theta_mle==TRUE & root==TRUE){
                W<-.Call("Weight_matrix", S1=alphaA$invectors, S=alphaA$vectors, lambda=alphaA$values, time=as.numeric(times), matdiag=as.numeric(matdiag))
                
            }else if(theta_mle==TRUE & root==FALSE){
                W<-design_matrix
            }else{
                W<-.Call("Expect_matrix",S1=alphaA$invectors, S=alphaA$vectors, lambda=alphaA$values,  time=as.numeric(times), theta0=theta0, theta1=theta1, matdiag=as.numeric(matdiag))
            }
            
            list(V=V, W=W)
        }
    },
    "randomRoot"={
        model_fun_matrix<-function(vcv,n,alphaA,sigmA,times,theta0=NULL,theta1=NULL,matdiag=NULL,theta_mle=TRUE){
            
            V<-.Call("simmap_covar", as.integer(n), bt=as.numeric(vcv), lambda=alphaA$values, S=alphaA$vectors, S1=alphaA$invectors, sigmasq=sigmA)
            
            if(theta_mle==TRUE & root==TRUE){
                W<-.Call("Weight_matrix", S1=alphaA$invectors, S=alphaA$vectors, lambda=alphaA$values, time=as.numeric(times), matdiag=as.numeric(matdiag))
                
            }else if(theta_mle==TRUE & root==FALSE){
                W<-design_matrix
            }else{
                W<-.Call("Expect_matrix",S1=alphaA$invectors, S=alphaA$vectors, lambda=alphaA$values,  time=as.numeric(times), theta0=theta0, theta1=theta1, matdiag=as.numeric(matdiag))
            }
            
            list(V=V, W=W)
        }
    },
    "fixedRoot"={
        model_fun_matrix<-function(vcv,n,alphaA,sigmA,times,theta0=NULL,theta1=NULL,matdiag=NULL,theta_mle=TRUE){
            
            V<-.Call("mvmorph_covar_mat", as.integer(n), bt=as.numeric(vcv), lambda=alphaA$values, S=alphaA$vectors, sigmasq=sigmA, S1=alphaA$invectors)
            
            if(theta_mle==TRUE & root==TRUE){
                W<-.Call("Weight_matrix", S1=alphaA$invectors, S=alphaA$vectors, lambda=alphaA$values, time=as.numeric(times), matdiag=as.numeric(matdiag))
                
            }else if(theta_mle==TRUE & root==FALSE){
                W<-design_matrix
            }else{
                W<-.Call("Expect_matrix",S1=alphaA$invectors, S=alphaA$vectors, lambda=alphaA$values,  time=as.numeric(times), theta0=theta0, theta1=theta1, matdiag=as.numeric(matdiag))
            }
            
            list(V=V, W=W)
        }
    },"univarFixed"={
        model_fun_matrix<-function(vcv,n,alphaA,sigmA,times,theta0=NULL,theta1=NULL,matdiag=NULL,theta_mle=TRUE){
            
            V<-.Call("mvmorph_covar_ou_fixed",A=vcv,alpha=alphaA$values, sigma=sigmA)
            
            if(theta_mle==TRUE & root==TRUE){
                W<-.Call("Weight_matrix", S1=alphaA$invectors, S=alphaA$vectors, lambda=alphaA$values, time=as.numeric(times), matdiag=as.numeric(matdiag))
                
            }else if(theta_mle==TRUE & root==FALSE){
                W<-design_matrix
            }else{
                W<-.Call("Expect_matrix",S1=alphaA$invectors, S=alphaA$vectors, lambda=alphaA$values,  time=as.numeric(times), theta0=theta0, theta1=theta1, matdiag=as.numeric(matdiag))
            }
            
            list(V=V, W=W)
        }
    },"univarRandom"={
        model_fun_matrix<-function(vcv,n,alphaA,sigmA,times,theta0=NULL,theta1=NULL,matdiag=NULL,theta_mle=TRUE){
            
            V<-.Call("mvmorph_covar_ou_random",A=vcv,alpha=alphaA$values, sigma=sigmA)
            
            if(theta_mle==TRUE & root==TRUE){
                W<-.Call("Weight_matrix", S1=alphaA$invectors, S=alphaA$vectors, lambda=alphaA$values, time=as.numeric(times), matdiag=as.numeric(matdiag))
                
            }else if(theta_mle==TRUE & root==FALSE){
                W<-design_matrix
            }else{
                W<-.Call("Expect_matrix",S1=alphaA$invectors, S=alphaA$vectors, lambda=alphaA$values,  time=as.numeric(times), theta0=theta0, theta1=theta1, matdiag=as.numeric(matdiag))
            }
            
            list(V=V, W=W)
        }
    },"univarpfFixed"={
        model_fun_matrix<-function(vcv,n,alphaA,sigmA,times,theta0=NULL,theta1=NULL,matdiag=NULL,theta_mle=TRUE){
            
             V<-.Call("mvmorph_covar_ou_rpf_fixed",A=vcv,alpha=alphaA$values, sigma=sigmA)
             
            if(theta_mle==TRUE & root==TRUE){
                W<-.Call("Weight_matrix", S1=alphaA$invectors, S=alphaA$vectors, lambda=alphaA$values, time=as.numeric(times), matdiag=as.numeric(matdiag))
                
            }else if(theta_mle==TRUE & root==FALSE){
                W<-design_matrix
            }else{
                W<-.Call("Expect_matrix",S1=alphaA$invectors, S=alphaA$vectors, lambda=alphaA$values,  time=as.numeric(times), theta0=theta0, theta1=theta1, matdiag=as.numeric(matdiag))
            }
            
            list(V=V, W=W)
        }
    },"univarpfRandom"={
        model_fun_matrix<-function(vcv,n,alphaA,sigmA,times,theta0=NULL,theta1=NULL,matdiag=NULL,theta_mle=TRUE){
            
            V<-.Call("mvmorph_covar_ou_rpf_random",A=vcv,alpha=alphaA$values, sigma=sigmA)
            
            if(theta_mle==TRUE & root==TRUE){
                W<-.Call("Weight_matrix", S1=alphaA$invectors, S=alphaA$vectors, lambda=alphaA$values, time=as.numeric(times), matdiag=as.numeric(matdiag))
                
            }else if(theta_mle==TRUE & root==FALSE){
                W<-design_matrix
            }else{
                W<-.Call("Expect_matrix",S1=alphaA$invectors, S=alphaA$vectors, lambda=alphaA$values,  time=as.numeric(times), theta0=theta0, theta1=theta1, matdiag=as.numeric(matdiag))
            }
            
            list(V=V, W=W)
        }
    })
    
    ##-----------------------Likelihood Calculation-------------------------------##
    
    devianc<-function(alpha,sigma,theta0=NULL,theta1=NULL,dat,error,theta_mle=TRUE){
        
        alphaA<-buildA(alpha,p)
        sigmA<-sigmafun(sigma)
        
        
        # avoid negative eigenvalues
        if(any(Re(alphaA$values)<=0)){return(1000000)}
        
        matEstim<-model_fun_matrix(vcv=vcv_time,n=n,alphaA,sigmA,times,theta0,theta1,matdiag,theta_mle)
        
        loglik<-loglik_mvmorph(dat,matEstim$V,matEstim$W,n,p,error=error,precalc=precalc,method=method,ch=ch,precalcMat=precalcMat,sizeD=sizeD,NA_val=NA_val,Indice_NA=Indice_NA,theta_mle=theta_mle,theta=as.matrix(1))
        
        
        if(is.infinite(loglik$logl)){
            return(10000000)
        }
        
        return(-loglik$logl)
        
    }
    
    # estim ancestral states from MLE
    estim.theta<-function(estimML){
        
        alphaA<-buildA(estimML[seq_len(nalpha)],p)
        sigmA<-sigmafun(estimML[nalpha+seq_len(nsigma)])
        
        matEstim<-model_fun_matrix(vcv=vcv_time,n=n,alphaA,sigmA,times,theta0=NULL,theta1=NULL,matdiag,theta_mle=TRUE)
        
        mvmorphEstim<-loglik_mvmorph(dat,matEstim$V,matEstim$W,n,p,error=error,precalc=precalc,method=method,ch=ch,precalcMat=precalcMat,sizeD=sizeD,NA_val=NA_val,Indice_NA=Indice_NA,theta_mle=TRUE)
        
        list(theta=mvmorphEstim$anc)
    }
    
    ##------------------function closure------------------------------------------##

    llfun <- function(par,...){
        
        args <- list(...)
        if(is.null(args[["root.mle"]])) args$root.mle <- TRUE
        if(is.null(args[["theta"]])) args$theta <- FALSE
        
        
        ## Arguments order: alpha, sigma, theta (if not the MLE)
        if(args$root.mle==TRUE){
           result <- -as.numeric(devianc(alpha=par[seq_len(nalpha)],sigma=par[nalpha+seq_len(nsigma)],error=error,dat=dat,theta_mle=TRUE))
           if(args$theta==TRUE){
               theta <- estim.theta(par)$theta
               result <- list(logl=result, theta=theta)
           }
        }else{
           result <- -as.numeric(devianc(alpha=par[seq_len(nalpha)],sigma=par[nalpha+seq_len(nsigma)], theta0=par[nalpha+nsigma+seq_len(ntheta0)], theta1=par[nalpha+nsigma+ntheta0+seq_len(ntheta1)],error=error,dat=dat,theta_mle=FALSE))
        }
        
        return(result)
    }
    # attributes to the loglik function
    attr(llfun, "model")<-"OU"
    attr(llfun, "alpha")<-nalpha
    attr(llfun, "sigma")<-nsigma
    attr(llfun, "theta")<-ntheta0+ntheta1
    
    # Matrix parameterization used for the A matrix
    
    decompfun <- function(par){
        buildA(par,p)
    }
    

    ## Check first if we want to return the log-likelihood function only
    if(optimization=="fixed"){
        message("No optimization performed, only the Log-likelihood function is returned.")
        param$nbspecies<-n
        param$ntraits<-p
        param$method<-method
        param$model<-"OUTS"
        param$optimization<-optimization
        param$traits<-colnames(data)
        param$vcv<-vcvtype
        param$alphafun<-decompfun
        param$sigmafun<-sigmafun
        theta.mat<-matrix(c(theta0,theta1),p)
        results<-list(llik=llfun, theta=theta.mat, alpha=buildA(alpha,p)$A, sigma=sigmafun(sigma), param=param)
        class(results)<-c("mvmorph")
        invisible(results)
        return(results)
    }
    ##----------------------Likelihood optimization-------------------------------##
    if(optimization!="subplex"){
        estim <- optim(par=c(alpha,sigma),fn = function (par) { devianc(alpha=par[seq_len(nalpha)],sigma=par[nalpha+seq_len(nsigma)], error=error,dat=dat)},gr=NULL,hessian=TRUE,method=optimization,control=control)# hack - use ntheta0 for the case there is a stationary mean
        hess<-eigen(estim$hessian)$values
    }else{# !subplex minimize a function, we must change the sign
        estim <- subplex(par=c(alpha,sigma),fn = function (par) { devianc(alpha=par[seq_len(nalpha)],sigma=par[nalpha+seq_len(nsigma)], error=error,dat=dat)},hessian=TRUE,control=control)
        hess<-eigen(estim$hessian)$values
    }
    
    ##---------------------Diagnostics--------------------------------------------##
    if(estim$convergence==0 & diagnostic==TRUE){
        cat("successful convergence of the optimizer","\n")
    }else if(estim$convergence==1 & diagnostic==TRUE){
        cat("maximum limit iteration has been reached, please consider increase maxit","\n")
    }else if(diagnostic==TRUE){
        cat("convergence of the optimizer has not been reached, try simpler model","\n")
    }
    
    if(any(hess<0)){
        hess.val<-1
        if(diagnostic==TRUE){
            cat("unreliable solution has been reached, check hessian eigenvectors or try simpler model","\n")}
    }else{
        hess.val<-0
        if(diagnostic==TRUE){
            cat("a reliable solution has been reached","\n")}
    }
    
    ##-------------------Summarize Results----------------------------------------##
    # names for the plot
    if(is.null(colnames(data))){
        names_data_matrix<-rep("",p)
        names_data<-list(names_data_matrix,names_data_matrix)
    }else{
        names_data_matrix<-colnames(data)
        names_data<-list(names_data_matrix,names_data_matrix)
    }
    
    # LogLik
    LL<- -estim$value 
    nparam=nalpha+nsigma+ntheta
    
    # maximum likelihood estimates of alpha and sigma
    estim.alpha<-estim$par[seq_len(nalpha)]
    estim.sigma<-estim$par[nalpha+seq_len(nsigma)]
    theta.mat<-estim.theta(estim$par)$theta
    #theta.mat<-estim$par[nalpha+nsigma+seq_len(ntheta)]
    
    # alpha matrix
    alpha.mat<-buildA(estim.alpha,p)$A
    # sigma matrix
    sigma.mat<-sigmafun(estim.sigma)
    dimnames(sigma.mat)<-names_data
    dimnames(alpha.mat)<-names_data
    
    # AIC
    AIC<- -2*LL+2*nparam
    # AIC corrected
    nobs <- length(which(!is.na(data)))
    AICc<-AIC+((2*nparam*(nparam+1))/(nobs-nparam-1)) # Hurvich et Tsai, 1989
    # matrix of estimated theta values
    theta.mat<-matrix(theta.mat,ncol=p)
    if(root==TRUE){
       rownames(theta.mat)<-c("theta_0","theta_1")
    }else{
       rownames(theta.mat)<-"theta"
    }
    colnames(theta.mat)<-names_data_matrix
    
    ##-------------------Print results--------------------------------------------##
    if(echo==TRUE){
        cat("\n")
        cat("-- Summary results for the Time-Series OU model --","\n")
        cat("LogLikelihood:","\t",LL,"\n")
        cat("AIC:","\t",AIC,"\n")
        cat("AICc:","\t",AICc,"\n")
        cat(nparam,"parameters","\n")
        cat("\n")
        cat("Estimated theta values","\n")
        cat("______________________","\n")
        print(theta.mat)
        cat("\n")
        cat("ML alpha values","\n")
        cat("______________________","\n")
        print(alpha.mat)
        cat("\n")
        cat("ML sigma values","\n")
        cat("______________________","\n")
        print(sigma.mat)
    }

    ##-------------------Save infos in parameters---------------------------------##
    param$nparam<-nparam
    param$nbspecies<-n
    param$ntraits<-p
    param$nregimes<-k
    param$method<-method
    param$model<-"OUTS"
    param$optimization<-optimization
    param$traits<-colnames(data)
    param$vcv<-vcvtype
    param$decomp<-decomp
    param$alphafun<-decompfun
    param$sigmafun<-sigmafun
    param$opt<-estim
    class(llfun) = c("mvmorph.llik")
    ##------------------List results----------------------------------------------##
    
    results<-list(LogLik=LL, AIC=AIC, AICc=AICc, theta=theta.mat, alpha=alpha.mat, sigma=sigma.mat, convergence=estim$convergence, hess.values=hess.val, param=param, llik=llfun)
    
    class(results)<-c("mvmorph","mvmorph.ou")
    invisible(results)
    
}
