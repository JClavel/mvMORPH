################################################################################
##                                                                            ##
##                       mvMORPH: estim missing traits                        ##
##                                                                            ##
##  Created by Julien Clavel - 01-01-2016                                     ##
##  (julien.clavel@hotmail.fr/ clavel@biologie.ens.fr)                        ##
##   require: phytools, ape, corpcor, spam                                    ##
##                                                                            ##
################################################################################


estim<-function(tree, data, object, error=NULL, asr=FALSE){

if(missing(tree)) stop("The tree or time-series object is missing!")
if(missing(data)) stop("You must provide a dataset along with your tree/time-series!")

#set data as a matrix if a vector is provided instead
if(!is.matrix(data)){data<-as.matrix(data)}
traits_names<-colnames(data)

# number of traits
p <- ncol(data)

# Check if there is missing cases
Indice_NA<-NULL
if(any(is.na(data))){
    Indice_NA<-which(is.na(as.vector(data)))
   }else if(asr!=TRUE){
    stop("There is no missing cases in your dataset!!")
}

# bind error to a vector
if(!is.null(error)){error<-as.vector(error[-Indice_NA])}

# First impute the dataset to compute the ancestral states?
if(any(is.na(data)) & asr==TRUE){
    warning("Missing cases were first imputed before estimating the ancestral values!!")
    data <- estim(tree,data,object,error, asr=FALSE)$estimates
}

##-------------------Models parameters----------------------------------------##
if(any(class(object)=="mvmorph")){
    ## Parameters
    model<-object$param$model[1]
    p<-object$param$ntraits
    k<-object$param$nregimes
    mu<-object$theta
    sigma<-object$sigma
    names_traits<-object$param$traits
    
    ## Ancestral states estimation??
    ## I must adapt it to the simmap format
    if(asr==TRUE){
        if(!inherits(tree,"simmap") & k!=1) stop("The tree should be an object of class \"simmap\".")
        if(!inherits(tree,"phylo")) stop("The \"asr=TRUE\" argument works only with phylogenies.")
        # number of tips
        ntip <- length(tree$tip.label)
        Indice_EXT <- indiceTip(tree,p)
    
        # species names
        sp_names <- tree$tip.label
        if(any(sp_names==rownames(data))){data<-data[sp_names,]}
        
        if(k!=1){
            C <- vcvSplit(tree, internal=TRUE)
            n <- ncol(C[[1]])
        }else{
            C <- vcvPhyloInternal(tree)
            n <- ncol(C)
        }

    
    }else{
        if(inherits(tree,"phylo")){
            
            # number of tips
            n <- length(tree$tip.label)
            # species names
            sp_names <- tree$tip.label

            if(k!=1){
                C <- vcvSplit(tree)
                if(any(sp_names==rownames(data))){C<-lapply(C,function(x) x[rownames(data),rownames(data)])}
            }else{
                C <- vcv.phylo(tree)
                if(any(sp_names==rownames(data))){C<-C[rownames(data),rownames(data)]}
            }

                
        }else{
                C <- vcv.ts(tree)
                # number of specimens
                n <- length(tree)
                # species names
                sp_names <- names(tree)
                if(any(sp_names==rownames(data))){C<-C[rownames(data),rownames(data)]}
                
        }
    }
    
    if(is.null(object$param[["root"]])){
        root<-FALSE
        if(model=="OUTS") root<-TRUE
    }else{
        root<-object$param$root
    }
    
    if(is.null(object[["trend"]])){
        istrend<-FALSE
    }else{
        istrend<-TRUE
        trend<-object$trend
    }
    
    # Root
    if(root=="stationary"){
        mod_stand<-1
    }else{
        mod_stand<-0 ## a modif
    }
    
    if(any(class(object)=="mvmorph.shift")){
        before<-object$param$before
        after<-object$param$after
    }
    
    if(model=="ER"){
        alpha<-object$alpha
        sig<-object$sigma
    }else if(model=="RR"){
        alpha<-object$alpha
        sig<-object$sig
    }
    
    if(model=="OV"){
        alpha<-object$alpha
        beta<-matrix(object$beta,p,p)
        sig<-object$sigma
    }else if(model=="OVG"){
        alpha<-object$alpha
        sig<-object$sig
        beta<-matrix(object$beta,p,p)
    }
    
    if(model=="CV"){
        beta<-matrix(object$beta,p,p)
        sig<-object$sigma
    }else if(model=="CVG"){
        beta<-matrix(object$beta,p,p)
        sig<-object$sig
    }
    
    if(model=="OUM" | model=="OU1" | model=="OUTS"){
        alpha<-object$alpha
    }
    
    if(model=="BMM"){
        sigma<-lapply(1:k,function(x){sigma[,,x]})
    }
    
    if(any(class(object)=="mvmorph.acdc")){
        beta<-matrix(object$beta,p,p)
        model<-"EB"
    }
    
}else{
stop("The imput model must be of class \"mvmorph\"!!")
}

# Compute the VCV and Design matrix for each models

switch(model,
"RWTS"={
    # Compute the design matrix
    smean<-TRUE
    W<-multD(NULL,p,n,smean=smean)
    V<-.Call("kronecker_mvmorph", R=sigma, C=C, Rrows=as.integer(p),  Crows=as.integer(n))
    # Add measurment error?
    if(is.null(error)==FALSE){
        if(asr==TRUE){
            diag(V[Indice_EXT,Indice_EXT])<-diag(V[Indice_EXT,Indice_EXT])+error
        }else{
            diag(V)<-diag(V)+error
        }
    }
    
    if(istrend==TRUE){
        Vdiag<-rep(diag(C),p)
        W<-W%*%as.numeric(mu) + Vdiag*(W%*%trend)
        mu<-as.numeric(1)
    }
},
"OUTS"={
    eig<-eigen(alpha)
    svec<-solve(eig$vectors)
    matdiag<-diag(p)
    
    if(object$param$vcv=="fixedRoot" | object$param$vcv=="univarpfFixed" | object$param$vcv=="univarFixed"){
        V<-.Call("mvmorph_covar_mat", as.integer(n), bt=C, lambda=eig$values, S=eig$vectors, sigmasq=sigma, S1=svec)
    }else if(object$param$vcv=="randomRoot" | object$param$vcv=="univarpfRandom" | object$param$vcv=="univarRandom"){
        V<-.Call("simmap_covar", as.integer(n), bt=C, lambda=eig$values, S=eig$vectors, S1=svec, sigmasq=sigma)
    }
    if(root==TRUE){
    W<-.Call("Weight_matrix", S1=svec, S=eig$vectors, lambda=eig$values, time=as.numeric(tree), matdiag=as.numeric(matdiag))
    }else{
    W<-multD(NULL,p,n,smean=TRUE)
    }
    # Add measurment error?
    if(is.null(error)==FALSE){
        if(asr==TRUE){
            diag(V[Indice_EXT,Indice_EXT])<-diag(V[Indice_EXT,Indice_EXT])+error
        }else{
            diag(V)<-diag(V)+error
        }    }
},
"BM1"={
    # Compute the design matrix
    W<-multD(tree,p,n,smean=object$param$smean)
    
    # Check if there where multiple vcv (compatibility with smean=FALSE)
    if(is.list(C)) C <- Reduce("+",C)
    
    V<-.Call("kronecker_mvmorph", R=sigma, C=C, Rrows=as.integer(p),  Crows=as.integer(n))
    # Add measurment error?
    if(is.null(error)==FALSE){
        if(asr==TRUE){
            diag(V[Indice_EXT,Indice_EXT])<-diag(V[Indice_EXT,Indice_EXT])+error
        }else{
            diag(V)<-diag(V)+error
        }
    }
    
    if(istrend==TRUE){
        Vdiag<-rep(diag(C[1:n,1:n]),p)
        W<-W%*%as.numeric(mu) + Vdiag*(W%*%trend)
        mu<-as.numeric(1)
    }
},
"BMM"={
    W<-multD(tree,p,n,smean=object$param$smean)
    V<-.Call("kroneckerSum", R=sigma, C=C, Rrows=as.integer(p),  Crows=as.integer(n), dimlist=as.integer(k)) # Add measurment error?
    if(is.null(error)==FALSE){
        if(asr==TRUE){
            diag(V[Indice_EXT,Indice_EXT])<-diag(V[Indice_EXT,Indice_EXT])+error
        }else{
            diag(V)<-diag(V)+error
        }
    }
    
    if(istrend==TRUE){
        Vdiag<-rep(diag(C[1:n,1:n]),p)
        W<-W%*%as.numeric(mu) + Vdiag*(W%*%trend)
        mu<-as.numeric(1)
    }
},
"OU1"={
    param_ou<-mv.Precalc(tree,nb.traits=p,param=list(model="OU1", root=root))
    epochs<-param_ou$epochs
    listReg<-param_ou$listReg
    bt<-param_ou$C1
    eig<-eigen(alpha)
    svec<-solve(eig$vectors)
    if(object$param$vcv=="fixedRoot" | object$param$vcv=="univarpfFixed" | object$param$vcv=="univarFixed" | object$param$vcv=="sparse"){
        V<-.Call("mvmorph_covar_mat", as.integer(n), bt=bt, lambda=eig$values, S=eig$vectors, sigmasq=sigma, S1=svec)
    }else{
        V<-.Call("simmap_covar", as.integer(n), bt=bt, lambda=eig$values, S=eig$vectors, S1=svec, sigmasq=sigma)
        #V<-.Call("simmap_covar", as.integer(n), bt=bt, lambda=eig$values, S=eig$vectors, sigmasq=sigma)
    }
    W<-.Call("mvmorph_weights",nterm=as.integer(n), epochs=epochs,lambda=eig$values,S=eig$vectors,S1=svec,beta=listReg,root=as.integer(mod_stand))
    # Add measurment error?
    if(is.null(error)==FALSE){
        if(asr==TRUE){
            diag(V[Indice_EXT,Indice_EXT])<-diag(V[Indice_EXT,Indice_EXT])+error
        }else{
            diag(V)<-diag(V)+error
        }
    }
},
"OUM"={
    param_ou<-mv.Precalc(tree,nb.traits=p,param=list(model="OUM", root=root))
    epochs<-param_ou$epochs
    listReg<-param_ou$listReg
    bt<-param_ou$C1
    eig<-eigen(alpha)
    svec<-solve(eig$vectors)
    if(object$param$vcv=="fixedRoot"| object$param$vcv=="univarpfFixed" | object$param$vcv=="univarFixed" | object$param$vcv=="sparse"){
        V<-.Call("mvmorph_covar_mat", as.integer(n), bt=bt, lambda=eig$values, S=eig$vectors, sigmasq=sigma, S1=svec)
    }else{
        V<-.Call("simmap_covar", as.integer(n), bt=bt, lambda=eig$values, S=eig$vectors, sigmasq=sigma)
    }
    W<-.Call("mvmorph_weights",nterm=as.integer(n), epochs=epochs,lambda=eig$values,S=eig$vectors,S1=svec,beta=listReg,root=as.integer(mod_stand))
    # Add measurment error?
    if(is.null(error)==FALSE){
        if(asr==TRUE){
            diag(V[Indice_EXT,Indice_EXT])<-diag(V[Indice_EXT,Indice_EXT])+error
        }else{
            diag(V)<-diag(V)+error
        }
    }
},
"EB"={
    # Compute the design matrix
    W<-multD(tree,p,n,smean=TRUE)
    V<-.Call("kroneckerEB",R=sigma,C=C, beta=beta, Rrows=as.integer(p),  Crows=as.integer(n))
    # Add measurment error?
    if(is.null(error)==FALSE){
        if(asr==TRUE){
            diag(V[Indice_EXT,Indice_EXT])<-diag(V[Indice_EXT,Indice_EXT])+error
        }else{
            diag(V)<-diag(V)+error
        }
    }
},

"RR"={
    # Brownian and OU models with different rates
    if(p==1){
        Vou<-.Call("mvmorph_covar_ou",A=C[[before]],alpha=alpha, sigma=sigma)
    }else{
        eig<-eigen(alpha)
        Vou<-.Call("mvmorph_covar_mat",nterm=as.integer(n),bt=C[[before]],lambda=eig$values,S=eig$vectors,sigmasq=sigma, S1=solve(eig$vectors))
    }
    V<-.Call("kronecker_shift", R=sig, C=C[[after]], Rrows=as.integer(p), Crows=as.integer(n), V=Vou)
    # Compute the design matrix
    W<-multD(tree,p,n,smean=TRUE)
    # Add measurment error?
    if(is.null(error)==FALSE){
        if(asr==TRUE){
            diag(V[Indice_EXT,Indice_EXT])<-diag(V[Indice_EXT,Indice_EXT])+error
        }else{
            diag(V)<-diag(V)+error
        }
    }
},
"ER"={
    # Brownian and OU models with the same rates
    if(p==1){
        Vou<-.Call("mvmorph_covar_ou",A=C[[before]],alpha=alpha, sigma=sigma)
    }else{
        eig<-eigen(alpha)
        Vou<-.Call("mvmorph_covar_mat",nterm=as.integer(n),bt=C[[before]],lambda=eig$values,S=eig$vectors,sigmasq=sigma, S1=solve(eig$vectors))
    }
    V<-.Call("kronecker_shift", R=sig, C=C[[after]], Rrows=as.integer(p), Crows=as.integer(n), V=Vou)
    # Compute the design matrix
    W<-multD(tree,p,n,smean=TRUE)
    # Add measurment error?
    if(is.null(error)==FALSE){
        if(asr==TRUE){
            diag(V[Indice_EXT,Indice_EXT])<-diag(V[Indice_EXT,Indice_EXT])+error
        }else{
            diag(V)<-diag(V)+error
        }
    }
},
"CV"={
    # Brownian & ACDC models with the same rates
    V<-.Call("kronecker_shiftEB_BM", R1=sig, R2=sigma, C1=C[[before]], C2=C[[after]], beta=alpha, Rrows=as.integer(p),  Crows=as.integer(n))
    # Compute the design matrix
    W<-multD(tree,p,n,smean=TRUE)
    # Add measurment error?
    if(is.null(error)==FALSE){
        if(asr==TRUE){
            diag(V[Indice_EXT,Indice_EXT])<-diag(V[Indice_EXT,Indice_EXT])+error
        }else{
            diag(V)<-diag(V)+error
        }
    }
},
"CVG"={
    # Brownian & ACDC models with different rates
    V<-.Call("kronecker_shiftEB_BM", R1=sig, R2=sigma, C1=C[[before]], C2=C[[after]], beta=alpha, Rrows=as.integer(p),  Crows=as.integer(n))
    # Compute the design matrix
    W<-multD(tree,p,n,smean=TRUE)
    # Add measurment error?
    if(is.null(error)==FALSE){
        if(asr==TRUE){
            diag(V[Indice_EXT,Indice_EXT])<-diag(V[Indice_EXT,Indice_EXT])+error
        }else{
            diag(V)<-diag(V)+error
        }
    }
},
"OV"={
    # OU & ACDC models with the same rates
    if(p==1){
        Vou<-.Call("mvmorph_covar_ou",A=C[[before]],alpha=alpha, sigma=sigma)
    }else{
        eig<-eigen(alpha)
        Vou<-.Call("mvmorph_covar_mat",nterm=as.integer(n),bt=C[[before]],lambda=eig$values,S=eig$vectors,sigmasq=sigma,S1=solve(eig$vectors))
    }
    V<-.Call("kronecker_shiftEB_OU", R=sig, C=C[[after]], beta=beta, Rrows=as.integer(p),  Crows=as.integer(n), V=Vou)
    # Compute the design matrix
    W<-multD(tree,p,n,smean=TRUE)
    # Add measurment error?
    if(is.null(error)==FALSE){
        if(asr==TRUE){
            diag(V[Indice_EXT,Indice_EXT])<-diag(V[Indice_EXT,Indice_EXT])+error
        }else{
            diag(V)<-diag(V)+error
        }
    }
},
"OVG"={
    # OU & ACDC models with independent rates
    if(p==1){
        Vou<-.Call("mvmorph_covar_ou",A=C[[before]],alpha=alpha, sigma=sigma)
    }else{
        eig<-eigen(alpha)
        Vou<-.Call("mvmorph_covar_mat",nterm=as.integer(n),bt=C[[before]],lambda=eig$values,S=eig$vectors,sigmasq=sigma,S1=solve(eig$vectors))
    }
    V<-.Call("kronecker_shiftEB_OU", R=sig, C=C[[after]], beta=beta, Rrows=as.integer(p),  Crows=as.integer(n), V=Vou)
    # Compute the design matrix
    W<-multD(tree,p,n,smean=TRUE)
    # Add measurment error?
    if(is.null(error)==FALSE){
        if(asr==TRUE){
            diag(V[Indice_EXT,Indice_EXT])<-diag(V[Indice_EXT,Indice_EXT])+error
        }else{
            diag(V)<-diag(V)+error
        }
    }
})

##-------------------Prepare the covariance matrices--------------------------##
# TODO: need first to compute the VCV of simmap trees

    
# stack the data as a vector
data<-as.numeric(as.matrix(data))

if(asr==TRUE){
    ## Parameters estimates
    anc <- as.numeric(as.matrix(W[Indice_EXT,])%*%as.numeric(mu))
    anc2 <- as.numeric(as.matrix(W[-Indice_EXT,])%*%as.numeric(mu))
    
    # extend Cunningham et al. 1998 to multivariate
    # Models included in the objects
    ## covariance between tip species and ancestral states
    varAY <- V[-Indice_EXT, Indice_EXT]
    ## covariance between ancestral states
    varA <- V[-Indice_EXT, -Indice_EXT]
    ## covariance between tip species
    varY <- V[Indice_EXT,Indice_EXT]
    data_complete<-data
    
    
    # Invert the covariance matrix
    invY<-solve(varY)
    ## Estimating missing cases; follow Cunningham et al. 1998
    data_estim <- as.numeric(varAY%*%invY%*%(data_complete-anc))+anc2
    
}else{
    ## Parameters estimates
    anc <- as.numeric(W%*%as.numeric(mu))
    
    # Prepare the covariances
    varY<-V[-Indice_NA,-Indice_NA]
    varAY<-V[Indice_NA,-Indice_NA]
    varA<-V[Indice_NA,Indice_NA]
    # Prepare the data
    anc<-anc[Indice_NA]
    data_complete<-data[-Indice_NA]
    
    
    # Invert the covariance matrix
    invY<-solve(varY)
    ## Estimating missing cases; follow Cunningham et al. 1998
    data_estim <- as.numeric(varAY%*%invY%*%(data_complete-anc))+anc
}



##-------------------Estimate ancestral or missing states---------------------##

## Estimate the variances of the imputed data
data_var<-diag(varA - varAY %*% invY %*% t(varAY))
data_se<-sqrt(data_var)
# data_CI<-cbind(data_estim + (data_sd * -1.959964), data_estim - (data_sd * -1.959964))
if(asr==TRUE){
    ## Names the results
    results_matrix<-matrix(data_estim, ncol=p)
    colnames(results_matrix)<-traits_names
    rownames(results_matrix)<-ntip+1:Nnode(tree)
    
    ## Names the results of variance and sd
    results_var<-matrix(data_var, ncol=p)
    results_se<-matrix(data_se, ncol=p)
    rownames(results_var)<-rownames(results_se)<-ntip+1:Nnode(tree)
    colnames(results_var)<-colnames(results_se)<-traits_names
  
}else{
    ## Names the results
    results_matrix<-matrix(data, ncol=p)
    rownames(results_matrix)<-sp_names
    colnames(results_matrix)<-traits_names
    results_matrix[Indice_NA]<-data_estim

    ## Names the results of variance and sd
    results_var<-results_se<-matrix(NA, nrow=n, ncol=p)
    rownames(results_var)<-rownames(results_se)<-sp_names
    colnames(results_var)<-colnames(results_se)<-traits_names
    results_var[Indice_NA]<-data_var
    results_se[Indice_NA]<-data_se
}
##-------------------Store results--------------------------------------------##

results<-list(estimates=results_matrix, var=results_var, se=results_se, NA_index=Indice_NA)

class(results)<-c("mvmorph.estim")
invisible(results)
#End
}
