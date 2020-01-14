################################################################################
##                                                                            ##
##                               mvMORPH: mvSIM                               ##
##                                                                            ##
##  Created by Julien Clavel - 22-11-2014                                     ##
##  (julien.clavel@hotmail.fr/ clavel@biologie.ens.fr)                        ##
##   require: phytools, ape, corpcor, spam                                    ##
##                                                                            ##
################################################################################



mvSIM<-function(tree,nsim=1,error=NULL,model=c("BM1","BMM","OU1","OUM","EB"), param=list(theta=0,sigma=0.1,alpha=1,beta=0)){
    
    if(any(class(param)=="mvmorph")){
    ## Parameters
    model<-param$param$model[1]
    p<-param$param$ntraits
    n<-param$param$nbspecies
    k<-param$param$nregimes
    mu<-param$theta
    sigma<-param$sigma
    names_traits<-param$param$traits
    
    if(model=="BM1" | model=="BMM"){
        param$smean<-param$param$smean
    }else{
        param$smean<-TRUE
    }
    
    if(is.null(param$param[["root"]])){
        root<-FALSE
    }else{
        root<-param$param$root
    }
    
    if(is.null(param$param[["trend"]])){
        istrend<-FALSE
    }else if(param$param$trend==FALSE){
        istrend<-FALSE
    }else{
        istrend<-TRUE
        trend<-as.vector(param$trend)
    }
    
    if(any(class(param)=="mvmorph.shift")){
        before<-param$param$before
        after<-param$param$after
    }
    
    if(model=="ER"){
        alpha<-param$alpha
        sig<-param$sigma
    }else if(model=="RR"){
        alpha<-param$alpha
        sig<-param$sig
    }
    
    if(model=="OV"){
        alpha<-param$alpha
        beta<-matrix(param$beta,p,p)
        sig<-param$sigma
    }else if(model=="OVG"){
        alpha<-param$alpha
        sig<-param$sig
        beta<-matrix(param$beta,p,p)
    }
    
    if(model=="CV"){
        beta<-matrix(param$beta,p,p)
        sig<-param$sigma
    }else if(model=="CVG"){
        beta<-matrix(param$beta,p,p)
        sig<-param$sig
    }
    
    if(model=="OUM" | model=="OU1" | model=="OUTS"){
        alpha<-param$alpha
        param$vcv<-param$param$vcv
    }
    
    if(model=="BMM"){
        sigma<-lapply(1:k,function(x){sigma[,,x]})
    }
    
    if(any(class(param)=="mvmorph.acdc")){
        beta<-matrix(param$beta,p,p)
        model<-"EB"
    }
    
    }else{
        ## Default values for simulating data
        # Choose the model
        model<-model[1]
        # check
        if(missing(tree)) stop("The tree or the time-series object is missing")
        
        if(model=="ACDC" | model=="AC"){model<-"EB"}
        
        # traits names
        if(is.null(param[["names_traits"]])){
            names_traits<-NULL
        }else{
            names_traits<-param$names_traits
        }
        
        # number of regimes
        if(model=="BM1" | model=="OU1" | model=="EB"){
            k<-1
        }else{
             if(inherits(tree,"phylo")){
                 if(!is.null(tree[["mapped.edge"]])){
                    k<-length(colnames(tree$mapped.edge))
                 }
             }else{
                 k<-1
             }
             
        }
        
        ## Root state for the weight matrix
        if(is.null(param[["root"]])){
            root<-param$root<-FALSE
            if(model=="OUTS") root<-param$root<-TRUE
        }else{ root<-param$root}

        ## Nombre de traits / number of traits
        if(is.null(param[["ntraits"]])){
            
            if(!is.null(param[["alpha"]])){
                p<-dim(as.matrix(param$alpha))[1]
            }else if(!is.null(param[["sigma"]])){
                if(model=="BMM"){
                p<-dim(as.matrix(param$sigma[[1]]))[1]
                }else{
                p<-dim(as.matrix(param$sigma))[1]
                }
            }else{
                p<-1
            cat("default simulations for 1 trait. Use \"ntraits\" argument in the \"param\" list  ","\n")
            }
            
        }else{
            p<-param$ntraits
        }
        
        if(root==TRUE){k<-k+1}
        
        # mu set default values (theta)
        if(!is.null(param[["mu"]])){
            param$theta<-param$mu
        }
        
        # theta value
        if(is.null(param[["theta"]])==TRUE){
            if(model=="OUM"){
                mu<-rep(0,k)
                if(p!=1){
                    mu<-rep(0,k*p)
                }
            }else{
                mu<-0
                if(p!=1){
                    mu<-rep(0,p)
                }
            }
                cat("No values provided for theta, default value fixed to:",mu,"\n")
            
        }else{
            if(p==1){
                if(k==1){
                        mu<-param$theta
                }else if(length(param$theta)!=k & model=="OUM"){
                        mu<-rep(0,k)
                    cat("Problems with theta values provided, default value fixed to:",mu,"\n")

                }else{
                    mu<-param$theta
                }
            }else{
                if(k==1){
                    mu<-param$theta
                }else if(length(param$theta)!=k*p & model=="OUM"){
                    
                        mu<-rep(0,k*p)
                 
                    cat("Problems with theta values provided, default value fixed to:",mu,"\n")
                    
                }else{
                    mu<-param$theta
                }

            }
        }
        
        # sigma
        if(is.null(param[["sigma"]])==TRUE){
            if(model=="BMM"){
                    sigma<-lapply(1:k,function(x){ runif(n=1)})
  
                if(p!=1){
                    sigma<-lapply(1:k,function(x){ diag(p)})
                }
                
            }else{
            sigma<-sig<-0.1
                if(p!=1){
                    sigma<-sig<-diag(p)
                }
            }
         print("No values provided for sigma, default value fixed to:",unlist(sigma),"\n")
        }else{
            if(k==1){
                sigma<-sig<-param$sigma
            }else if(length(param$sigma)!=k & model=="BMM"){
                    sigma<-lapply(1:k,function(x){ runif(n=1)})
                if(p!=1){
                    sigma<-lapply(1:k,function(x){ diag(p)})
                }
            }else if(length(param$sigma)==k & is.list(param$sigma)==FALSE){
                    sigma<-lapply(1:k,function(x){ runif(n=1)})
                if(p!=1){
                    sigma<-lapply(1:k,function(x){ diag(p)})
                }
            }else if(model!="OUBMi" & model!="BMOUi" & model!="OUEBi" & model!="EBOUi" & model!="BMEBi" & model!="EBBMi" & model!="RR" & model!="RC"){
                sigma<-sig<-param$sigma
            }else{
                sigma<-param$sigma[[1]]
                sig<-param$sigma[[2]]
                
            }
        }
        
        ## Parameters
        
        if(model=="EC"|| model=="RC" || model=="BMOU" || model=="BMOUi"){
            before<-2
            after<-1
            if(model=="EC"||model=="BMOU"){
                model<-"ER"
            }else if(model=="RC" || model=="BMOUi"){
                model<-"RR"
            }
        }else if(model=="ER"|| model=="RR" || model=="OUBM" || model=="OUBMi"){
            before<-1
            after<-2
            if(model=="RR" || model=="OUBMi"){
                model<-"RR"
            }else if(model=="ER" || model=="OUBM"){
                model<-"ER"
            }
        }else if(model=="EBOU" || model=="EBOUi"){# OU-ACDC models
            before<-2
            after<-1
            if(model=="EBOU"){
                model<-"OV"
            }else if(model=="EBOUi"){
                model<-"OVG"
            }
        }else if(model=="OUEB" || model=="OUEBi"){
            before<-1
            after<-2
            if(model=="OUEB"){
                model<-"OV"
            }else if(model=="OUEBi"){
                model<-"OVG"
            }
        }else if(model=="EBBM" || model=="EBBMi"){ # BM-ACDC models
            before<-2
            after<-1
            if(model=="EBBM"){
                model<-"CV"
            }else if(model=="EBBMi"){
                model<-"CVG"
            }
        }else if(model=="BMEB" || model=="BMEBi"){
            before<-1
            after<-2
            if(model=="BMEB"){
                model<-"CV"
            }else if(model=="BMEBi"){
                model<-"CVG"
            }
        }
        # tree with shift
        
        if(model=="RR" | model=="ER" | model=="OVG" | model=="OV" | model=="CV" | model=="CVG"){
            #set age shift with make.era from phytools if the age is provided
            if(!is.null(param[["age"]])){
                tot<-max(nodeHeights(tree))
                limit=tot-param$age
                tree<-make.era.map(tree, c(0,limit))
            }
        }
        
        # alpha
        if(model=="OU1" | model=="OUM" | model=="RR" | model=="ER" | model=="OVG" | model=="OV" | model=="OUTS"){
            if(is.null(param[["alpha"]])==TRUE){
                alpha<-1
                if(p!=1){
                    alpha<-diag(p)
                }
                cat("No values provided for alpha, default value fixed to:",alpha,"\n")
            }else{
                alpha<-param$alpha
            }
        }
        
        # beta
        if(model=="EB" | model=="CV" | model=="CVG" | model=="OVG" | model=="OV"){
            if(is.null(param[["beta"]])==TRUE){
                beta<-1
                cat("No values provided for beta, default value fixed to:",beta,"\n")
                if(p!=1){
                    beta<-matrix(1,p,p)
                }
            }else{
                beta<-param$beta
                if(p!=1){
                    beta<-matrix(param$beta,p,p)
                }
            }
        }
        
        # trend options
        if(is.null(param[["trend"]])){
            istrend<-FALSE
        }else{
            istrend<-TRUE
            trend<-param$trend
            
            if(length(trend)<p){
                stop("The number of specified values for the slopes of the trend do not match the number of traits!","\n")
            }
        }
        
        
        # verif nombre de traits
        if(model=="BMM"){
            sigmadim<-length(sigma[[1]])
        }else{
            sigmadim<-length(sigma)
        }
        
        if(sigmadim!=p*p){
            stop("The number of specified traits do not match the parameters matrix dimensions!","\n")
        }
        
    }


# Root
if(root=="stationary"){
   mod_stand<-1
}else{
   mod_stand<-0 ## a modif
}

##Species covariance matrix
# switch methods depending on the nature of the tree object
if(inherits(tree,"phylo")){
    
    if(model=="OUTS" | model=="RWTS") stop("Error! you must provides a time-series, not a phylo object.")
    ## Nombre d'especes / number of species
    n<-length(tree$tip.label)
      if(model=="OU1" | model=="BM1" | model=="EB"){
            C<-vcv.phylo(tree)
        }else{
            C<-vcvSplit(tree)
        }
         names_rows <- tree$tip.label
    }else{
        # to change?
            tree <- tree-min(tree)
            C<-vcv.ts(tree)
            ## Nombre d'especes / number of species
            n<-length(tree)
            if(is.matrix(tree)){
                names_rows <- rownames(tree)
            }else{
                names_rows <- names(tree)
            }
            
}

# vcv mtrix for Ornstein-Uhlenbeck model
if(model=="OUM" | model=="OU1" | model=="OUTS"){
    if(is.null(param[["vcv"]])){
        vcv<-param$vcv<-"randomRoot"
    }else{
        vcv<-param$vcv
    }
    
    if(model=="OUTS"){
        if(is.null(error) & vcv=="fixedRoot" |  any(error[1,]==0) & vcv=="fixedRoot"){
            warning("No sampling error provided for the initial observation(s) while using the \"fixedRoot\" parameterization.","\n","An arbitrary sampling error value is set to 0.01 for the initial observation(s) to avoid singularity issues. ")
            if(any(error[1,]==0)==TRUE){
                if(p==1){error[1] <- 0.01}else{error[1,] <- 0.01}
            }else{
                error <- matrix(0,ncol=p, nrow=n)
                error[1,] <- 0.01
            }
        }
    }
    
}

# sampling error default for the RWTS model
if(model=="RWTS"){
    if(is.null(error) |  any(error[1,]==0)){
        warning("No sampling error provided for the initial observation(s).","\n","An arbitrary sampling error value is set to 0.01 for the initial observation(s) to avoid singularity issues. ")
        if(any(error[1,]==0)==TRUE){
            if(p==1){error[1] <- 0.01}else{error[1,] <- 0.01}
        }else{
            error <- matrix(0,ncol=p, nrow=n)
            error[1,] <- 0.01
        }
    }
}

# sampling error values:
if(!is.null(error)){
    error<-as.vector(as.matrix(error))
    error[is.na(error)] <- 0
}

# Select the method for drawing values from the multivariate normal distribution; default is "cholesky"
if(is.null(param[["method"]])){ methodSim <- "cholesky" }else{ methodSim <- param$method }

# Compute the VCV and Design matrix for each models

switch(model,
"RWTS"={
    # Compute the design matrix
    param$smean<-TRUE
    W<-multD(tree,p,n,smean=param$smean)
    V<-.Call(kronecker_mvmorph, R=sigma, C=C, Rrows=as.integer(p),  Crows=as.integer(n))
    # Add measurment error?
    if(is.null(error)==FALSE){
        diag(V)<-diag(V)+error
    }
    if(ncol(W)!=length(mu)) stop("\n","The number of parameters for theta is wrong","\n",ncol(W)," values are expected")

    if(istrend==TRUE){
        Vdiag<-rep(diag(C),p)
        W<-W%*%as.numeric(mu) + Vdiag*(W%*%trend)
        mu<-as.numeric(1)
    }
},
"BM1"={
    # Compute the design matrix
    if(is.null(param[["smean"]])==TRUE){ param$smean<-TRUE }
    
    # If there is no trends/errors/multiple means. We can simplify the computation of BM1 (for large dimensions...)
    if(istrend==FALSE & is.null(error) & p!=1 & param$smean){
        
        Croot <- chol(C); # the cholesky factor to correlate the data
        W<-matrix(1,nrow=n,ncol=1) # multiple mean? To be improved later
        
        if(nsim==1 & p!=1){
            X <- rmvnorm_simul(n=n, mean=rep(0,p), var=sigma, method=methodSim)
            deviates <- t(X%*%Croot)
            traits <- matrix(W%*%as.numeric(mu),ncol=p) + deviates
            rownames(traits)<-names_rows
            colnames(traits)<-names_traits
             return(traits) # if we return we exit the function.
        }else if(nsim>1 & p!=1){
            traits<-lapply(1:nsim,function(x){
                X <- rmvnorm_simul(n=n, mean=rep(0,p), var=sigma, method=methodSim)
                deviates<-t(X%*%Croot);
                traits <- matrix(W%*%as.numeric(mu),ncol=p) + deviates;
                rownames(traits)<-names_rows; colnames(traits)<-names_traits;
                traits})
             return(traits) # if we return we exit the function.
        }
    
    } # else we use the full kronecker matrix...
    
    W<-multD(tree,p,n,smean=param$smean)
    V<-.Call(kronecker_mvmorph, R=sigma, C=C, Rrows=as.integer(p),  Crows=as.integer(n))
    # Add measurment error?
    if(is.null(error)==FALSE){
        diag(V)<-diag(V)+error
    }
    if(ncol(W)!=length(mu)) stop("\n","The number of parameters for theta is wrong","\n",ncol(W)," values are expected")
    if(istrend==TRUE){
        Vdiag<-rep(diag(C),p)
        W<-W%*%as.numeric(mu) + Vdiag*(W%*%trend)
        mu<-as.numeric(1)
    }
},
"BMM"={
    if(is.null(param[["smean"]])==TRUE){ param$smean<-TRUE }
    W<-multD(tree,p,n,smean=param$smean)
    V<-.Call(kroneckerSum, R=sigma, C=C, Rrows=as.integer(p),  Crows=as.integer(n), dimlist=as.integer(k)) # Add measurment error?
    if(is.null(error)==FALSE){
        diag(V)<-diag(V)+error
    }
    if(ncol(W)!=length(mu)) stop("\n","The number of parameters for theta is wrong","\n",ncol(W)," values are expected")

    if(istrend==TRUE){
        Vdiag<-rep(diag(Reduce('+',C)),p)
        W<-W%*%as.numeric(mu) + Vdiag*(W%*%trend)
        mu<-as.numeric(1)
    }
},
"OUTS"={
    eig<-eigen(alpha)
    svec<-try(solve(eig$vectors), silent = TRUE)
    if(inherits(svec ,'try-error')){
        warning("An error occured with the inverse of the orthogonal matrix, the \"pseudoinverse\" has been used instead")
        svec<- pseudoinverse(eig$vectors)
    }
    matdiag<-diag(p)
    
    if(vcv=="fixedRoot" | vcv=="univarpfFixed" | vcv=="univarFixed" | vcv=="sparse"){
        V<-.Call(mvmorph_covar_mat, as.integer(n), bt=C, lambda=eig$values, S=eig$vectors, sigmasq=sigma, S1=svec)
    }else if(vcv=="randomRoot" | vcv=="univarpfRandom" | vcv=="univarRandom"){
        V<-.Call(simmap_covar, as.integer(n), bt=C, lambda=eig$values, S=eig$vectors, S1=svec, sigmasq=sigma)
    }
    if(root==TRUE){
    W<-.Call(Weight_matrix, S1=svec, S=eig$vectors, lambda=eig$values, time=as.numeric(tree), matdiag=as.numeric(matdiag))
    }else{
    W<-multD(NULL,p,n,smean=TRUE)
    }
    if(ncol(W)!=length(mu)) stop("\n","The number of parameters for theta is wrong","\n",ncol(W)," values are expected")

    # Add measurment error?
    if(is.null(error)==FALSE){
        diag(V)<-diag(V)+error
    }
},
"OU1"={
    param_ou<-mv.Precalc(tree,nb.traits=p,param=list(model="OU1", root=root))
    epochs<-param_ou$epochs
    listReg<-param_ou$listReg
    bt<-param_ou$C1
    eig<-eigen(alpha)
    svec<-try(solve(eig$vectors), silent = TRUE)
    if(inherits(svec ,'try-error')){
        warning("An error occured with the inverse of the orthogonal matrix, the \"pseudoinverse\" has been used instead")
        svec<- pseudoinverse(eig$vectors)
    }
    
    if(vcv=="fixedRoot" | vcv=="univarpfFixed" | vcv=="univarFixed" | vcv=="sparse"){
        V<-.Call(mvmorph_covar_mat, as.integer(n), bt=bt, lambda=eig$values, S=eig$vectors, sigmasq=sigma, S1=svec)
    }else if(vcv=="randomRoot" | vcv=="univarpfRandom" | vcv=="univarRandom"){
        V<-.Call(simmap_covar, as.integer(n), bt=bt, lambda=eig$values, S=eig$vectors, S1=svec, sigmasq=sigma)
    }
    
    W<-.Call(mvmorph_weights,nterm=as.integer(n), epochs=epochs,lambda=eig$values,S=eig$vectors,S1=svec,beta=listReg,root=as.integer(mod_stand))
    if(ncol(W)!=length(mu)) stop("\n","The number of parameters for theta is wrong","\n",ncol(W)," values are expected")

    # Add measurment error?
    if(is.null(error)==FALSE){
        diag(V)<-diag(V)+error
    }
},
"OUM"={
    param_ou<-mv.Precalc(tree,nb.traits=p,param=list(model="OUM", root=root))
    epochs<-param_ou$epochs
    listReg<-param_ou$listReg
    bt<-param_ou$C1
    eig<-eigen(alpha)
    svec<-try(solve(eig$vectors), silent = TRUE)
    if(inherits(svec ,'try-error')){
        warning("An error occured with the inverse of the orthogonal matrix, the \"pseudoinverse\" has been used instead")
        svec<- pseudoinverse(eig$vectors)
    }
    
    if(vcv=="fixedRoot" | vcv=="univarpfFixed" | vcv=="univarFixed" | vcv=="sparse"){
        V<-.Call(mvmorph_covar_mat, as.integer(n), bt=bt, lambda=eig$values, S=eig$vectors, sigmasq=sigma, S1=svec)
    }else if(vcv=="randomRoot" | vcv=="univarpfRandom" | vcv=="univarRandom"){
        V<-.Call(simmap_covar, as.integer(n), bt=bt, lambda=eig$values, S=eig$vectors, S1=svec, sigmasq=sigma)
    }
    
    W<-.Call(mvmorph_weights,nterm=as.integer(n), epochs=epochs,lambda=eig$values,S=eig$vectors,S1=svec,beta=listReg,root=as.integer(mod_stand))
    if(ncol(W)!=length(mu)) stop("\n","The number of parameters for theta is wrong","\n",ncol(W)," values are expected")

    # Add measurment error?
    if(is.null(error)==FALSE){
        diag(V)<-diag(V)+error
    }
},
"EB"={
    # Compute the design matrix
    if(is.null(param[["smean"]])==TRUE){ param$smean<-TRUE }
    W<-multD(tree,p,n,smean=param$smean)
    V<-.Call(kroneckerEB,R=sigma,C=C, beta=beta, Rrows=as.integer(p),  Crows=as.integer(n))
    if(ncol(W)!=length(mu)) stop("\n","The number of parameters for theta is wrong","\n",ncol(W)," values are expected")

    # Add measurment error?
    if(is.null(error)==FALSE){
        diag(V)<-diag(V)+error
    }
},

"RR"={
    # Brownian and OU models with different rates
    if(p==1){
        Vou<-.Call(mvmorph_covar_ou_fixed,A=C[[before]],alpha=alpha, sigma=sigma)
    }else{
        eig<-eigen(alpha)
        Vou<-.Call(mvmorph_covar_mat,nterm=as.integer(n),bt=C[[before]],lambda=eig$values,S=eig$vectors,sigmasq=sigma, S1=solve(eig$vectors))
    }
    V<-.Call(kronecker_shift, R=sig, C=C[[after]], Rrows=as.integer(p), Crows=as.integer(n), V=Vou)
    # Compute the design matrix
    if(is.null(param[["smean"]])==TRUE){ param$smean<-TRUE }
    W<-multD(tree,p,n,smean=param$smean)
    if(ncol(W)!=length(mu)) stop("\n","The number of parameters for theta is wrong","\n",ncol(W)," values are expected")

    # Add measurment error?
    if(is.null(error)==FALSE){
        diag(V)<-diag(V)+error
    }
},
"ER"={
    # Brownian and OU models with the same rates
    if(p==1){
        Vou<-.Call(mvmorph_covar_ou_fixed,A=C[[before]],alpha=alpha, sigma=sigma)
    }else{
        eig<-eigen(alpha)
        Vou<-.Call(mvmorph_covar_mat,nterm=as.integer(n),bt=C[[before]],lambda=eig$values,S=eig$vectors,sigmasq=sigma, S1=solve(eig$vectors))
    }
    V<-.Call(kronecker_shift, R=sig, C=C[[after]], Rrows=as.integer(p), Crows=as.integer(n), V=Vou)
    # Compute the design matrix
    if(is.null(param[["smean"]])==TRUE){ param$smean<-TRUE }
    W<-multD(tree,p,n,smean=param$smean)
    if(ncol(W)!=length(mu)) stop("\n","The number of parameters for theta is wrong","\n",ncol(W)," values are expected")

    # Add measurment error?
    if(is.null(error)==FALSE){
        diag(V)<-diag(V)+error
    }
},
"CV"={
    # Brownian & ACDC models with the same rates
    V<-.Call(kronecker_shiftEB_BM, R1=sig, R2=sigma, C1=C[[before]], C2=C[[after]], beta=alpha, Rrows=as.integer(p),  Crows=as.integer(n))
    # Compute the design matrix
    if(is.null(param[["smean"]])==TRUE){ param$smean<-TRUE }
    W<-multD(tree,p,n,smean=param$smean)
    if(ncol(W)!=length(mu)) stop("\n","The number of parameters for theta is wrong","\n",ncol(W)," values are expected")

    # Add measurment error?
    if(is.null(error)==FALSE){
        diag(V)<-diag(V)+error
    }
},
"CVG"={
    # Brownian & ACDC models with different rates
    V<-.Call(kronecker_shiftEB_BM, R1=sig, R2=sigma, C1=C[[before]], C2=C[[after]], beta=alpha, Rrows=as.integer(p),  Crows=as.integer(n))
    # Compute the design matrix
    if(is.null(param[["smean"]])==TRUE){ param$smean<-TRUE }
    W<-multD(tree,p,n,smean=param$smean)
    if(ncol(W)!=length(mu)) stop("\n","The number of parameters for theta is wrong","\n",ncol(W)," values are expected")

    # Add measurment error?
    if(is.null(error)==FALSE){
        diag(V)<-diag(V)+error
    }
},
"OV"={
    # OU & ACDC models with the same rates
    if(p==1){
        Vou<-.Call(mvmorph_covar_ou_fixed,A=C[[before]],alpha=alpha, sigma=sigma)
    }else{
        eig<-eigen(alpha)
        Vou<-.Call(mvmorph_covar_mat,nterm=as.integer(n),bt=C[[before]],lambda=eig$values,S=eig$vectors,sigmasq=sigma,S1=solve(eig$vectors))
    }
    V<-.Call(kronecker_shiftEB_OU, R=sig, C=C[[after]], beta=beta, Rrows=as.integer(p),  Crows=as.integer(n), V=Vou)
    # Compute the design matrix
    if(is.null(param[["smean"]])==TRUE){ param$smean<-TRUE }
    W<-multD(tree,p,n,smean=param$smean)
    if(ncol(W)!=length(mu)) stop("\n","The number of parameters for theta is wrong","\n",ncol(W)," values are expected")

    # Add measurment error?
    if(is.null(error)==FALSE){
        diag(V)<-diag(V)+error
    }
},
"OVG"={
    # OU & ACDC models with independent rates
    if(p==1){
        Vou<-.Call(mvmorph_covar_ou_fixed,A=C[[before]],alpha=alpha, sigma=sigma)
    }else{
        eig<-eigen(alpha)
        Vou<-.Call(mvmorph_covar_mat,nterm=as.integer(n),bt=C[[before]],lambda=eig$values,S=eig$vectors,sigmasq=sigma,S1=solve(eig$vectors))
    }
    V<-.Call(kronecker_shiftEB_OU, R=sig, C=C[[after]], beta=beta, Rrows=as.integer(p),  Crows=as.integer(n), V=Vou)
    # Compute the design matrix
    if(is.null(param[["smean"]])==TRUE){ param$smean<-TRUE }
    W<-multD(tree,p,n,smean=param$smean)
    if(ncol(W)!=length(mu)) stop("\n","The number of parameters for theta is wrong","\n",ncol(W)," values are expected")

    # Add measurment error?
    if(is.null(error)==FALSE){
        diag(V)<-diag(V)+error
    }
})


## Simulate the traits from a multivariate normal distribution

    if(nsim==1){
        traits<-matrix(rmvnorm_simul(n=1,as.numeric(W%*%as.numeric(mu)),V,methodSim),ncol=p)
        rownames(traits)<-names_rows
        colnames(traits)<-names_traits
    }else if(nsim>1 & p!=1){
        traits<-lapply(1:nsim,function(x){traits<-matrix(rmvnorm_simul(n=1,as.numeric(W%*%as.numeric(mu)),V,methodSim),ncol=p);rownames(traits)<-names_rows; colnames(traits)<-names_traits;traits})
    }else{
        traits<-matrix(rmvnorm_simul(n=nsim,as.numeric(W%*%as.numeric(mu)),V,methodSim),ncol=nsim)
        rownames(traits)<-names_rows
        #colnames(traits)<-names_traits
    }
    return(traits)
}

