################################################################################
##                                                                            ##
##                               mvMORPH: mvTIME/mvTS                         ##
##                                                                            ##
##  Created by Julien Clavel - 04-12-2015                                     ##
##  (julien.clavel@hotmail.fr/ clavel@biologie.ens.fr)                        ##
##   require: phytools, ape, corpcor, spam                                    ##
##                                                                            ##
################################################################################


mvRWTS <- function(times, data, error=NULL, param=list(sigma=NULL, trend=FALSE, decomp="cholesky"),method=c("rpf","inverse","pseudoinverse"),scale.height=FALSE, optimization=c("L-BFGS-B","Nelder-Mead","subplex"),control=list(maxit=20000),precalc=NULL, diagnostic=TRUE, echo=TRUE){
    
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
    
    # Scale the time-series to have the root (oldest point in time) at zero
    # check for external time series for the expectation?
    if(scale.height==TRUE){
        times <- times/max(times)
    }
    
    # Compute the covariance matrix for the time serie
    vcv_time <- vcv.ts(times)
    
    # define the default user constraint
    index.user <- NULL

    # constraints?
    if(is.null(param[["constraint"]])==FALSE & is.null(param[["decomp"]])==TRUE){
        
        if(class(param$constraint)=="matrix"){
            # user defined constraints
            index.user<-param$constraint
            if(!isSymmetric(unname(index.user)))
            stop("The user specified matrix constraint must be square and symmetric")
            constraint<-param$constraint<-"default"
            decomp<-param$decomp<-"user"
            message("Efficient parameterization is not yet implemented for user-defined constraints","\n","Please check the results carefully!!")
            # shrink the smallest eigenvalues
            if(!is.null(param[["tol"]])){ tol<-param$tol }else{ tol<-param$tol<-1 }
        }else if(param$constraint==TRUE){
            decomp <- param$decomp <- "equal"
        }else{
            decomp <- param$constraint
        }
    }
     
    # decomposition of the scatter matrix
    if(is.null(param[["decomp"]])==TRUE){
        decomp <- param$decomp <- "cholesky"
    }else{
        decomp <- param$decomp[1]
    }
    
    # model for now
    model <- "BM1"
    
    # sigma matrix prameterization
    
    buildSigma<-function(par,index.mat=NULL,sig=NULL,model="BM1",decomp){
 
            switch(decomp,
            "user"={
                if(model=="BMM"){
                    sig[]<-c(par)[index.mat]
                    sigma<-lapply(1:k,function(x){symPar(sig[x,], decomp="user", p=p, index.user=index.user, tol=tol)})
                }else if(model=="BM1"){
                     sigma<-symPar(par, decomp="user", p=p, index.user=index.user, tol=tol)
                }
            },
            "cholesky"={
                if(model=="BMM"){
                    sig[]<-c(par)[index.mat]
                    sigma<-lapply(1:k,function(x){ sym.par(sig[x,])})
                }else if(model=="BM1"){
                    sigma<-sym.par(par)
                }
            },
            "spherical"={
                if(model=="BMM"){
                    sig[]<-c(par)[index.mat]
                    sigma<-lapply(1:k,function(x){symPar(sig[x,], decomp="spherical", p=p)})
                }else if(model=="BM1"){
                    sigma<-symPar(par, decomp="spherical", p=p)
                }
            },
            "eigen"={
                if(model=="BMM"){
                    sig[]<-c(par)[index.mat]
                    sigma<-lapply(1:k,function(x){ symPar(sig[x,], decomp="eigen", p=p) })
                }else if(model=="BM1"){
                    sigma<-symPar(par, decomp="eigen", p=p)
                }
            },
            "eigen+"={
                if(model=="BMM"){
                    sig[]<-c(par)[index.mat]
                    sigma<-lapply(1:k,function(x){ symPar(sig[x,], decomp="eigen+", p=p) })
                }else if(model=="BM1"){
                    sigma<-symPar(par, decomp="eigen+", p=p)
                }
            },
            "diagonal"={
                if(model=="BMM"){
                    sig[] <- c(par)[index.mat]
                    sigma<-lapply(1:k,function(x){ diag(diag(sig[x,]%*%t(sig[x,])))})
                    
                }else if(model=="BM1"){
                    sigma<-diag(diag(par%*%t(par)))
                }
            },
            "equal"={
                if(model=="BMM"){
                    sig[] <- c(par)[index.mat]
                    sigma<-lapply(1:k,function(x){ symPar(sig[x,], decomp="equal", p=p)})
                    
                }else if(model=="BM1"){
                    sigma<-symPar(par, decomp="equal", p=p)
                }
            })
            
        
        return(sigma)
    }
    
    # starting values (need to be changed for taking into account the parameterization
    
    if(is.null(param[["sigma"]])){
        sigma<-param$sigma<-startParamSigma(p, decomp, times, data, index.user=index.user)
    }else{
        if(length(param$sigma)==p*(p+1)/2 & decomp!="user"){
            sigma<-param$sigma
        }else{
            sigma<-startParamSigma(p, decomp, times, data, guess=param$sigma, index.user=index.user)
        }
    }
    
    if(is.null(param[["theta"]])){
        if(!is.vector(data)){
            theta0<-param$theta0<-colMeans(data, na.rm=TRUE)
        }else{
            theta0<-param$theta0<-mean(data, na.rm=TRUE)
        }
    }else{
        theta0<-param$theta
    }
    
    # add a trend to the process
    if(is.null(param[["trend"]])==TRUE){
        istrend<-param$trend<-FALSE
        Vdiag<-NULL
        ntrend<-0
        trend_par<-rep(0,p)
        startvalue<-sigma
    }else if(any(param$trend==FALSE)){
        istrend<-FALSE
        Vdiag<-NULL
        ntrend<-0
        trend_par<-rep(0,p)
        startvalue<-sigma
    }else{
        Vdiag<-rep(diag(vcv_time),p)
        istrend<-TRUE
        theta_par<-colMeans(data, na.rm = TRUE)
        ntheta<-p
        if(!is.numeric(param$trend)){
            trend_par<-rep(0,p)
            index.mat<-1:p
            ntrend<-p
            
        }else if(is.numeric(param$trend) & length(param$trend)==p){
            ntrend<-length(unique(param$trend))
            trend_par<-rep(0,ntrend)
            index.mat<-param$trend
            
        }else{
            stop("the size of the integer indices vector for the \"trend\" argument is wrong")
        }
        startvalue<-c(sigma,trend_par,theta_par)
    }
    
    # Managing sampling errors and prevent singular matrices
    if(is.null(error) |  any(error[1,]==0)){
        warning("No sampling error provided for the initial observation(s).","\n"," An arbitrary sampling error value is set to 0.01 for the initial observation(s) to avoid singularity issues during the likelihood computation. ")
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

    # precalculate the design matrix for BM and Trend models
    design_matrix <- multD(NULL,p,n,smean=TRUE)
    
    # number of parameters for each matrix
    nsigma<-length(sigma)
    ntheta<-p
    
    # number of columns of the design matrix
    sizeD<-p

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
    
    # Define the variance-covariance functions
    switch(method,
    "rpf"={
        bm_fun_matrix<-function(C,sig,dat,D,precalcMat,n,p,error,method,theta_mle,theta,istrend,trend_val){
            V<-.Call("kronecker_mvmorph", R=sig, C=C, Rrows=as.integer(p),  Crows=as.integer(n), PACKAGE="mvMORPH")
            loglik<-loglik_mvmorph(dat,V,D,n,p,error=error,precalc=precalc,method=method,ch=ch,precalcMat=precalcMat,sizeD=ncol(D),NA_val=NA_val,Indice_NA=Indice_NA,theta_mle=theta_mle,theta=theta,istrend=istrend,trend=trend_val)
            return(loglik)
        }
    },
    "sparse"={
        bm_fun_matrix<-function(C,sig,dat,D,precalcMat,n,p,error,method,theta_mle,theta,istrend,trend_val){
            V<-.Call("kroneckerSumSpar", R=list(sig), C=list(C), Rrows=as.integer(p),  Crows=as.integer(n),  dimlist=as.integer(1), IA=IAr, JA=JAr, A=precalcMat@entries, PACKAGE="mvMORPH")
            loglik<-loglik_mvmorph(dat,V,D,n,p,error=error,precalc=precalc,method=method,ch=ch,precalcMat=precalcMat,sizeD=ncol(D),NA_val=FALSE,Indice_NA=NULL,theta_mle=theta_mle,theta=theta,istrend=istrend,trend=trend_val)
            return(loglik)
        }
    },
    "pseudoinverse"={
        bm_fun_matrix<-function(C,sig,dat,D,precalcMat,n,p,error,method,theta_mle,theta,istrend,trend_val){
            V<-.Call("kronecker_mvmorph", R=sig, C=C, Rrows=as.integer(p),  Crows=as.integer(n), PACKAGE="mvMORPH")
            loglik<-loglik_mvmorph(dat,V,D,n,p,error=error,precalc=precalc,method=method,ch=ch,precalcMat=precalcMat,sizeD=ncol(D),NA_val=NA_val,Indice_NA=Indice_NA,theta_mle=theta_mle,theta=theta,istrend=istrend,trend=trend_val)
            return(loglik)
        }
    },
    "inverse"={
        bm_fun_matrix<-function(C,sig,dat,D,precalcMat,n,p,error,method,theta_mle,theta,istrend,trend_val){
            V<-.Call("kronecker_mvmorph", R=sig, C=C, Rrows=as.integer(p),  Crows=as.integer(n), PACKAGE="mvMORPH")
            loglik<-loglik_mvmorph(dat,V,D,n,p,error=error,precalc=precalc,method=method,ch=ch,precalcMat=precalcMat,sizeD=ncol(D),NA_val=NA_val,Indice_NA=Indice_NA,theta_mle=theta_mle,theta=theta,istrend=istrend,trend=trend_val)
            return(loglik)
        }
    })
    
    ##---------------------Loglik BM1---------------------------------------------##
    lik.BM<-function(par,dat,C,D,error,method,precalcMat,n,p, constraint,theta_mle=TRUE,theta=NULL,istrend=FALSE){ ##
        if(istrend==TRUE){
            trend_val<-Vdiag*D%*%par[nsigma+seq_len(ntrend)][index.mat]
             if(NA_val==TRUE) trend_val<-trend_val[-Indice_NA]
            theta_mle<-FALSE
            theta<-par[nsigma+ntrend+seq_len(ntheta)]
        }
        loglik<-bm_fun_matrix(C,buildSigma(par[seq_len(nsigma)],NULL,NULL,model,constraint),dat,D,precalcMat,n,p,error,method,theta_mle,theta,istrend,trend_val)
        
        list(loglik=-loglik$logl, ancstate=loglik$anc)
    }
    
    ##-------------------Function log-likelihood----------------------------------##
    llfun <- function(par,...){
        
        args <- list(...)
        if(is.null(args[["root.mle"]])) args$root.mle <- TRUE
        if(is.null(args[["theta"]])) args$theta <- FALSE
        
        if(args$root.mle==TRUE){
            
            result <- -as.numeric(lik.BM(par=par,dat=data,C=vcv_time,D=design_matrix,error=error,method=method,precalcMat=precalcMat,n=n,p=p,constraint=decomp,theta_mle=TRUE,theta=NULL,istrend=istrend)$loglik)
            
            if(args$theta==TRUE){
              theta <- as.numeric(lik.BM(par=par,dat=data,C=vcv_time,D=design_matrix,error=error,method=method,precalcMat=precalcMat,n=n,p=p,constraint=decomp,theta_mle=TRUE,theta=NULL,istrend=istrend)$ancstate)
              result<-list(logl=result, theta=theta)
            }
            
        }else{
            
            result <- -as.numeric(lik.BM(par=par[seq_len(nsigma+ntrend)],dat=data,C=vcv_time,D=design_matrix,error=error,method=method,precalcMat=precalcMat,n=n,p=p,constraint=decomp,theta_mle=FALSE,theta=par[nsigma+ntrend+seq_len(ntheta)],istrend=istrend)$loglik)
            
        }
        return(result)
    }
    # attributes to the loglik function
    attr(llfun, "model")<-"RWTS"
    attr(llfun, "sigma")<-nsigma
    attr(llfun, "theta")<-ntheta
    attr(llfun, "trend")<-ntrend
    class(llfun) = c("mvmorph.llik")
    
    # Sigma fun function
    sigmafun <- function(par){ buildSigma(par,NULL,NULL,model,decomp)}
    
    ## Check first if we want to return the log-likelihood function only
    if(optimization=="fixed"){
        message("No optimization performed, only the Log-likelihood function is returned.")
        param$nbspecies<-n
        param$ntraits<-p
        param$method<-method
        param$model<-"RWTS"
        param$optimization<-optimization
        param$traits<-colnames(data)
        param$sigmafun<-sigmafun
        theta.mat<-matrix(theta0,nrow=1)
        results<-list(llik=llfun, theta=theta.mat, sigma=sigmafun(sigma), trend=trend_par, param=param)
        class(results)<-c("mvmorph")
        invisible(results)
        return(results)
    }
    ##----------------------Likelihood optimization-------------------------------##

    if(optimization!="subplex"){
        estim<-optim(par=startvalue,fn=function(par){lik.BM(par=par,dat=data,C=vcv_time,D=design_matrix,error=error,method=method,precalcMat=precalcMat,n=n,p=p,constraint=decomp,theta_mle=TRUE,theta=NULL,istrend=istrend)$loglik},control=control,hessian=TRUE,method=optimization)
    }else{
        estim<-subplex(par=startvalue,fn=function(par){lik.BM(par=par,dat=data,C=vcv_time,D=design_matrix,error=error,method=method,precalcMat=precalcMat,n=n,p=p,constraint=decomp,theta_mle=TRUE,theta=NULL,istrend=istrend)$loglik},control=control,hessian=TRUE)
    }
    
    ##---------------------Diagnostics--------------------------------------------##
    if(estim$convergence==0 & diagnostic==TRUE){
        cat("successful convergence of the optimizer","\n")
    }else if(estim$convergence==1 & diagnostic==TRUE){
        cat("maximum limit iteration has been reached, please consider increase maxit","\n")
    }else if(diagnostic==TRUE){
        cat("convergence of the optimizer has not been reached, try simpler model","\n")
    }
    
    # Hessian eigen decomposition to check the derivatives
    hess<-eigen(estim$hessian)$values

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
        if(p==1){
            names_data<-list("sigma^2:","")
            
        }else{
            names_data_matrix<-rep("",p)
            names_data<-list(names_data_matrix,names_data_matrix)
        }
    }else{
        names_data_matrix<-colnames(data)
        names_data<-list(names_data_matrix,names_data_matrix)
    }
    
    # LogLik
    LL<- -estim$value
    nparam=nsigma+ntheta+ntrend
    
    # maximum likelihood estimates sigma
    estim.sigma<-estim$par[seq_len(nsigma)]
    
    #ancestral states estimates
    theta.mat<-matrix(lik.BM(par=estim$par,dat=data,C=vcv_time,D=design_matrix,error=error, p=p, n=n, precalcMat=precalcMat, method=method,constraint=decomp,theta_mle=TRUE,theta=NULL,istrend=istrend)$ancstate,nrow=1)
    rownames(theta.mat)<-c("theta_0:")
    colnames(theta.mat)<-names_data_matrix
   
    # sigma matrix
    sigma.mat<-sigmafun(estim.sigma)
    dimnames(sigma.mat)<-names_data
    
    # trend matrix
    if(istrend==TRUE){
        trend.mat<-matrix(estim$par[nsigma+seq_len(ntrend)][index.mat], nrow=1)
        rownames(trend.mat)<-c("drift:")
        colnames(trend.mat)<-names_data_matrix
    }
    
    # AIC
    AIC<- -2*LL+2*nparam
    
    # AIC corrected
    nobs <- length(which(!is.na(data)))
    AICc<-AIC+((2*nparam*(nparam+1))/(nobs-nparam-1)) # Hurvich et Tsai, 1989
    
    ##-------------------Print results--------------------------------------------##
    if(echo==TRUE){
        
        cat("\n")
        cat("-- Summary results for the Time-Series BM/RW model --","\n")
        cat("LogLikelihood:","\t",LL,"\n")
        cat("AIC:","\t",AIC,"\n")
        cat("AICc:","\t",AICc,"\n")
        cat(nparam,"parameters","\n")
        cat("\n")
        cat("Estimated rate matrix","\n")
        cat("______________________","\n")
        print(sigma.mat)
        cat("\n")
        cat("Estimated theta values","\n")
        cat("______________________","\n")
        print(theta.mat)
        cat("\n")
        if(istrend==TRUE){
        cat("Estimated trend values","\n")
        cat("______________________","\n")
        print(trend.mat)
        cat("\n")
        }
    }
    
    ##-------------------Save infos in parameters---------------------------------##
    param$nparam<-nparam
    param$constraint<-FALSE
    param$nbspecies<-n
    param$ntraits<-p
    param$nregimes<-k
    param$method<-method
    param$model<-"RWTS"
    param$optimization<-optimization
    param$traits<-colnames(data)
    param$decomp<-decomp
    param$sigmafun<-sigmafun
    param$opt<-estim
    
    ##------------------List results----------------------------------------------##
    
    if(istrend==TRUE){
        
    results<-list(LogLik=LL, AIC=AIC, AICc=AICc, theta=theta.mat, sigma=sigma.mat, trend=trend.mat, convergence=estim$convergence, hess.values=hess.val, param=param, llik=llfun)
    }else{
    results<-list(LogLik=LL, AIC=AIC, AICc=AICc, theta=theta.mat, sigma=sigma.mat,  convergence=estim$convergence, hess.values=hess.val, param=param, llik=llfun)
    }
    
    class(results)<-c("mvmorph","mvmorph.bm")
    invisible(results)
    
}
