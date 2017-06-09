################################################################################
##                                                                            ##
##                               mvMORPH: mvBM                                ##
##                                                                            ##
##  Created by Julien Clavel - 16-07-2013                                     ##
##  (julien.clavel@hotmail.fr/ julien.clavel@univ-lyon1.fr)                   ##
##   require: phytools, ape, corpcor, spam                                    ##
##                                                                            ##
################################################################################


mvBM<-function(tree, data, error=NULL, model=c("BMM","BM1"),param=list(constraint=FALSE, smean=TRUE, trend=FALSE), method=c("rpf","pic","sparse","inverse","pseudoinverse"), scale.height=FALSE, optimization=c("L-BFGS-B","Nelder-Mead","subplex"), control=list(maxit=20000), precalc=NULL, diagnostic=TRUE, echo=TRUE){

# select default model
model<-model[1]
method<-method[1]
if(missing(tree)) stop("The tree object is missing!")
if(missing(data)) stop("You must provide a dataset along with your tree!")

#set data as a matrix if a vector is provided instead
if(!is.matrix(data)){data<-as.matrix(data)}

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

# number of root states
if(is.null(param[["smean"]])==TRUE){param$smean<-TRUE}

# number of selective regimes
if(model!="BM1" | param$smean==FALSE){
    k<-length(colnames(tree$mapped.edge))
}else{
    k<-1
}

# bind error to a vector
if(!is.null(error)){
    error<-as.vector(error)
    error[is.na(error)] <- 0
}

# number of species (tip)
n<-dim(data)[1]

# number of variables
p<-dim(data)[2]

# method for the optimizer & algorithm
optimization<-optimization[1]

# user defined constraints
index.user<-NULL

# Constraints parameterizations
if(is.null(param[["constraint"]])==TRUE){
    constraint<-param$constraint<-"default"
}else{
    if(class(param$constraint)=="matrix"){
        # user defined constraints
        index.user<-param$constraint
        if(!isSymmetric(unname(index.user)))
        stop("The user specified matrix constraint must be square and symmetric")
        constraint<-param$constraint<-"default"
        decomp<-param$decomp<-"user"
        message("Efficient parameterization is not yet implemented for user-defined constraints","\n","Please check the results carefully!!")
        # tolerance for min eigenvalue shrinkage
        if(!is.null(param[["tol"]])){ tol<-param$tol }else{ tol<-param$tol<-1 }
    }else if((param$constraint=="variance" & model=="BM1") | (param$constraint=="shared" & model=="BM1") | (param$constraint=="proportional" & model=="BM1") | (param$constraint=="correlation" & model=="BM1")){
        constraint<-param$constraint<-"default"
        if(echo==TRUE) cat(" \"shared\",\"variance\",\"correlation\" and \"proportional\" can be used only with BMM model","\n")
        
    }else if((param$constraint=="shared" & k==1) | (param$constraint=="correlation" & k==1) | (param$constraint=="proportional" & k==1) | (param$constraint=="variance" & k==1)){
        if(echo==TRUE) cat(" \"shared\",\"variance\",\"correlation\" and \"proportional\" can be used only with multiple dimension","\n")
        
        constraint<-param$constraint<-"default"
    }else if(param$constraint==TRUE){
        constraint<-param$constraint<-"equal"
    }else if(param$constraint==FALSE){
        constraint<-param$constraint<-"default"
    }else{
        constraint<-param$constraint
    }
}

##------------------------Precalc---------------------------------------------##
if(is.null(precalc)==FALSE & class(precalc)=="mvmorph.precalc"){
    tree<-precalc$tree
    
    if(is.null(tree[["mapped.edge"]])==TRUE & model=="BMM"){
        model<-"BM1"
       if(echo==TRUE) cat("No selective regimes mapped on the tree, only a BM1 model could be estimated","\n")
    }
    C1<-precalc$C1
    D<-precalc$D
    if((model=="BM1" & method=="pic") | (model=="BM1" & is.null(tree[["mapped.edge"]])==TRUE) ){param$smean<-TRUE}
    if(model=="BMM"){C2<-precalc$C2}
    # number of selective regimes
    if(model!="BM1"){
        k<-length(C2)
    }else{ k<-1 }
    
    if(method=="sparse"){
        # Yale sparse format
        JAr<-precalc$JAr
        IAr<-precalc$IAr
        ch<-precalc$ch
        precalcMat<-precalc$V
    }

    
}else{
##------------------------Precalc-off-----------------------------------------##
##------------------------Create VCV matrix-----------------------------------##
if(is.null(tree[["mapped.edge"]])==TRUE & model=="BMM"){
    model<-"BM1"
    if(echo==TRUE) cat("No selective regimes mapped on the tree, only a BM1 model could be estimated","\n")
}

# Scale the tree
if(scale.height==TRUE){
    maxHeight<-max(nodeHeights(tree))
    tree$edge.length<-tree$edge.length/maxHeight
    if(model=="BMM")tree$mapped.edge<-tree$mapped.edge/maxHeight
}

# Parameterization of the covariances/tree
if(method=="pic"){
    ind<-reorder.phylo(tree,"postorder", index.only=TRUE)
    tree$edge<-tree$edge[ind,]
    tree$edge.length<-tree$edge.length[ind]
    value<-list(tree$edge.length)
    
    if(model=="BMM"){
        tree$mapped.edge<-tree$mapped.edge[ind,]
    }
    # method for computing the log-likelihood
    if(model=="BMM" & p!=1){
        #mvMORPH-1.0.3
        warning("Sorry, the \"pic\" method only works with univariate data for the BMM model ","\n","the \"rpf\" method has been used instead...","\n")
        method<-"rpf"
    }
    C1<-NULL
    C2<-NULL
}

if(method!="pic"){
 
 # Compute the vcv matrix for BM1 or sparse matrix
 if(model=="BM1"){
    C1<-vcv.phylo(tree)
	if(!is.null(rownames(data))) {
	 if(any(tree$tip.label==rownames(data))){
         C1<-C1[rownames(data),rownames(data)]
     }else if(echo==TRUE){
         cat("row names of the data matrix must match tip names of your phylogeny!","\n")
     }
    }else if(echo==TRUE){
        cat("species in the matrix are assumed to be in the same order as in the phylogeny, otherwise specify rownames of 'data'","\n")
	}
 }else if(model=="BMM"){
  multi.tre<-list()
  class(multi.tre)<-"multiPhylo"
  C2<-list()
	for(i in 1:ncol(tree$mapped.edge)){
		multi.tre[[i]]<-tree
		multi.tre[[i]]$edge.length<-tree$mapped.edge[,i]
		multi.tre[[i]]$state<-colnames(tree$mapped.edge)[i]
		temp<-vcv.phylo(multi.tre[[i]])
		if(any(tree$tip.label==rownames(data))) { 
            C2[[i]]<-temp[rownames(data),rownames(data)]
		}else{
            C2[[i]]<-temp
        }
      }
    }# end BMM
}

##------------------------Parameters------------------------------------------##

# Compute the design matrix
if(model=="BM1" & method=="pic"){param$smean<-TRUE}
D<-multD(tree,p,n,smean=param$smean)

# Decomp (take care about the constraint option)
if(is.null(param[["decomp"]])==TRUE){
    decomp <- param$decomp <- "cholesky"
}else if(param$decomp=="diagonal" & param$constraint=="default"){
    constraint <- param$constraint <- "diagonal"
    decomp <- param$decomp <- "diagonal"
}else if(param$decomp=="equal" & param$constraint=="default"){
    constraint <- param$constraint <- "equal"
    decomp <- param$decomp <- "equal"
}else{
    decomp <- param$decomp
}

if(method=="sparse"){
    if(model=="BMM"){
        C1<-Reduce('+',C2)
    }
    if(constraint=="diagonal"){
        V<-kronecker(diag(p), C1)
    }else{
        V<-kronecker((matrix(1,p,p)+diag(p)), C1)
    }
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

# add a trend to the process
if(is.null(param[["trend"]])==TRUE){
    istrend<-param$trend<-FALSE
    Vdiag<-NULL
    npartrend<-0
    trend_par<-NULL
}else if(any(param$trend==FALSE)){
    istrend<-FALSE
    Vdiag<-NULL
    npartrend<-0
    trend_par<-NULL
}else{
    if(method=="pic"){
        warning("The \"pic\" method can't be used with the \"trend\" option. The \"rpf\" method is used instead")
        method<-"rpf"
    }
    if(model=="BMM"){
        C1<-Reduce('+',C2)
    }
    # trend
    Vdiag<-rep(diag(C1),p)
    
    # we don't use the MLE for theta
    if(param$smean==FALSE){
        theta_par<-rep(colMeans(data, na.rm = TRUE), k)
    }else{
        theta_par<-colMeans(data, na.rm = TRUE)
    }
    istrend<-TRUE
    if(!is.numeric(param$trend)){
        trend_par<-c(rep(0,p),theta_par)
        tr_index<-1:p
        npartrend<-p
        
    }else if(is.numeric(param$trend) & length(param$trend)==p){
        npartrend<-length(unique(param$trend))
        trend_par<-c(rep(0,npartrend),theta_par)
        tr_index<-param$trend
        
    }else if(is.numeric(param$trend) & param$smean==FALSE){
        npartrend<-length(unique(param$trend))
        trend_par<-c(rep(0,npartrend),theta_par)
        tr_index<-param$trend
        if(length(param$trend)!=(p*k)){
            stop("the size of the integer indices vector for the \"trend\" argument is wrong")
        }
    }else{
        stop("the size of the integer indices vector for the \"trend\" argument is wrong")
    }
    
}


}#end off-precalc

# Define the variance-covariance functions
if(model=="BMM"){
switch(method,
"rpf"={
    bm_fun_matrix<-function(C,sig,dat,D,precalcMat,n,p,error,method,theta_mle,theta,istrend,trend_val){
    V<-.Call(kroneckerSum, R=sig, C=C, Rrows=as.integer(p),  Crows=as.integer(n), dimlist=as.integer(k))
    loglik<-loglik_mvmorph(dat,V,D,n,p,error=error,precalc=precalc,method=method,ch=ch,precalcMat=precalcMat,sizeD=ncol(D),NA_val=NA_val,Indice_NA=Indice_NA,theta_mle=theta_mle,theta=theta,istrend=istrend,trend=trend_val)
    return(loglik)
    }
},
"sparse"={
    bm_fun_matrix<-function(C,sig,dat,D,precalcMat,n,p,error,method,theta_mle,theta,istrend,trend_val){
    V<-.Call(kroneckerSumSpar, R=sig, C=C, Rrows=as.integer(p),  Crows=as.integer(n),  dimlist=as.integer(k), IA=IAr, JA=JAr, A=precalcMat@entries)
    loglik<-loglik_mvmorph(dat,V,D,n,p,error=error,precalc=precalc,method=method,ch=ch,precalcMat=precalcMat,sizeD=ncol(D),NA_val=FALSE,Indice_NA=NULL,theta_mle=theta_mle,theta=theta,istrend=istrend,trend=trend_val)
    return(loglik)
    }
},
"pseudoinverse"={
    bm_fun_matrix<-function(C,sig,dat,D,precalcMat,n,p,error,method,theta_mle,theta,istrend,trend_val){
    V<-.Call(kroneckerSum, R=sig, C=C, Rrows=as.integer(p),  Crows=as.integer(n), dimlist=as.integer(k))
    loglik<-loglik_mvmorph(dat,V,D,n,p,error=error,precalc=precalc,method=method,ch=ch,precalcMat=precalcMat,sizeD=ncol(D),NA_val=NA_val,Indice_NA=Indice_NA,theta_mle=theta_mle,theta=theta,istrend=istrend,trend=trend_val)
    return(loglik)
    }
},
"inverse"={
    bm_fun_matrix<-function(C,sig,dat,D,precalcMat,n,p,error,method,theta_mle,theta,istrend,trend_val){
    V<-.Call(kroneckerSum, R=sig, C=C, Rrows=as.integer(p),  Crows=as.integer(n), dimlist=as.integer(k))
    loglik<-loglik_mvmorph(dat,V,D,n,p,error=error,precalc=precalc,method=method,ch=ch,precalcMat=precalcMat,sizeD=ncol(D),NA_val=NA_val,Indice_NA=Indice_NA,theta_mle=theta_mle,theta=theta,istrend=istrend,trend=trend_val)
    return(loglik)
    }
},
"pic"={
    bm_fun_matrix<-function(C,sig,dat,D,precalcMat,n,p,error,method,theta_mle,theta,istrend,trend_val){
    # Mult only available for univariate case
    p=1
    sig<-unlist(sig)
    tree$edge.length<-tree$mapped.edge%*%sig
    # method for pic
    modpic <- ifelse(theta_mle==TRUE,7,8)
    # Compute the LLik
    res<-.Call(PIC_gen, x=dat, n=as.integer(p), Nnode=as.integer(tree$Nnode), nsp=as.integer(n), edge1=as.integer(tree$edge[,1]), edge2=as.integer(tree$edge[,2]), edgelength=list(tree$edge.length), times=1, rate=rep(0,p), Tmax=1, Model=as.integer(modpic), mu=theta, sigma=1)
    logl<- -0.5 * ( n * p * log( 2 * pi) +  res[[5]] + n * res[[6]]  + res[[4]] )
    # return the theta value provided instead of the MLE...
    if(theta_mle==TRUE){ theta<-res[[7]] }
    return(list(logl=logl,ancstate=theta, sigma=res[[2]]))
        }
})

}else{

switch(method,
"rpf"={
bm_fun_matrix<-function(C,sig,dat,D,precalcMat,n,p,error,method,theta_mle,theta,istrend,trend_val){
    V<-.Call(kronecker_mvmorph, R=sig, C=C, Rrows=as.integer(p),  Crows=as.integer(n))
    loglik<-loglik_mvmorph(dat,V,D,n,p,error=error,precalc=precalc,method=method,ch=ch,precalcMat=precalcMat,sizeD=ncol(D),NA_val=NA_val,Indice_NA=Indice_NA,theta_mle=theta_mle,theta=theta,istrend=istrend,trend=trend_val)
    return(loglik)
    }
},
"sparse"={
bm_fun_matrix<-function(C,sig,dat,D,precalcMat,n,p,error,method,theta_mle,theta,istrend,trend_val){
    V<-.Call(kroneckerSumSpar, R=list(sig), C=list(C), Rrows=as.integer(p),  Crows=as.integer(n),  dimlist=as.integer(1), IA=IAr, JA=JAr, A=precalcMat@entries)
    loglik<-loglik_mvmorph(dat,V,D,n,p,error=error,precalc=precalc,method=method,ch=ch,precalcMat=precalcMat,sizeD=ncol(D),NA_val=FALSE,Indice_NA=NULL,theta_mle=theta_mle,theta=theta,istrend=istrend,trend=trend_val)
    return(loglik)
    }
},
"pseudoinverse"={
bm_fun_matrix<-function(C,sig,dat,D,precalcMat,n,p,error,method,theta_mle,theta,istrend,trend_val){
    V<-.Call(kronecker_mvmorph, R=sig, C=C, Rrows=as.integer(p),  Crows=as.integer(n))
    loglik<-loglik_mvmorph(dat,V,D,n,p,error=error,precalc=precalc,method=method,ch=ch,precalcMat=precalcMat,sizeD=ncol(D),NA_val=NA_val,Indice_NA=Indice_NA,theta_mle=theta_mle,theta=theta,istrend=istrend,trend=trend_val)
    return(loglik)
    }
},
"inverse"={
bm_fun_matrix<-function(C,sig,dat,D,precalcMat,n,p,error,method,theta_mle,theta,istrend,trend_val){
    V<-.Call(kronecker_mvmorph, R=sig, C=C, Rrows=as.integer(p),  Crows=as.integer(n))
    
    loglik<-loglik_mvmorph(dat,V,D,n,p,error=error,precalc=precalc,method=method,ch=ch,precalcMat=precalcMat,sizeD=ncol(D),NA_val=NA_val,Indice_NA=Indice_NA,theta_mle=theta_mle,theta=theta,istrend=istrend,trend=trend_val)
    return(loglik)
    }
},
"pic"={
bm_fun_matrix<-function(C,sig,dat,D,precalcMat,n,p,error,method,theta_mle,theta,istrend,trend_val){
    modpic <- ifelse(theta_mle==TRUE,7,6)
    res<-.Call(PIC_gen, x=dat, n=as.integer(p), Nnode=as.integer(tree$Nnode), nsp=as.integer(n), edge1=as.integer(tree$edge[,1]), edge2=as.integer(tree$edge[,2]), edgelength=list(tree$edge.length), times=1, rate=rep(0,k), Tmax=1, Model=as.integer(modpic), mu=theta, sigma=sig)
    logl<- -0.5 * ( n * p * log( 2 * pi) +  res[[5]] + n * res[[6]]  + res[[4]] )
    # return the theta value provided instead of the MLE...
    if(theta_mle==TRUE){ theta<-res[[7]] }
    return(list(logl=logl,ancstate=theta, sigma=res[[2]]))
    }
})
# End if
}

##------------------Sigma_matrix_parameterization-----------------------------##
# sigma matrix prameterization

buildSigma<-function(par,index.mat=NULL,sig=NULL,model,constraint){
        
        switch(constraint,
        "default"={
            if(model=="BMM"){
                sig[]<-c(par)[index.mat]
                sigma<-lapply(1:k,function(x){ symPar(sig[x,], decomp=decomp, p=p, index.user=index.user, tol=tol) })
            }else if(model=="BM1"){
                sigma<-symPar(par, decomp=decomp, p=p, index.user=index.user, tol=tol)
            }
        },
        "equaldiagonal"={
            if(model=="BMM"){
                sig[] <- c(par)[index.mat]
                sigma <- lapply(1:k,function(x){ symPar(sig[x,], decomp=decomp, p=p) })
            }else if(model=="BM1"){
                sigma <- symPar(par, decomp=decomp, p=p)
            }
        },
        "diagonal"={
            if(model=="BMM"){
                sig[] <- c(par)[index.mat]
                sigma <- lapply(1:k,function(x){ symPar(sig[x,], decomp=decomp, p=p) })
            }else if(model=="BM1"){
                sigma <- symPar(par, decomp=decomp, p=p)
            }
        },
        "equal"={
            if(model=="BMM"){
                sig[] <- c(par)[index.mat]
                sigma<-lapply(1:k,function(x){ symPar(sig[x,], decomp="equal", p=p)})
            }else if(model=="BM1"){
                sigma<-symPar(par, decomp="equal", p=p)
            }
        },
        "variance"={
            # Comparisons between groups (ie not allowed for BM1) - shared variance
            diagval<-1:p
            variance <- par[diagval]
            param<-par[-diagval]
            sig[] <- c(param)[index.mat]
            sigma<-lapply(1:k,function(x){.Call(spherical, param=sig[x,], variance=variance, dim=as.integer(p))})
        },
        "correlation"={
            # Comparisons between groups (ie not allowed for BM1) - shared correlations
            diagval<-1:npar
            angle <- par[diagval]
            param<-par[-diagval]
            sig[] <- c(param)[index.mat]
            sigma<-lapply(1:k,function(x){.Call(spherical, param=angle, variance=sig[x,], dim=as.integer(p))})
        },
        "shared"={
            # Comparisons between groups (ie not allowed for BM1) - shared eigenvectors
            diagval<-1:npar
            angle <- par[diagval]
            param<-par[-diagval]
            sig[] <- c(param)[index.mat]
            Q<-.Call(givens_ortho, Q=diag(p), angle=angle, ndim=as.integer(p))
            sigma<-lapply(1:k,function(x){ Q%*%diag(exp(sig[x,]))%*%t(Q) })
        },
        "proportional"={
            # Comparisons between groups (ie not allowed for BM1)
            propindex<-1:(k-1)
            constant<-par[propindex]
            rates_mat<-symPar(par[-propindex], decomp=decomp, p=p)
            sigma<-lapply(1:k,function(x){ if(x==1){rates_mat}else{rates_mat*as.vector((constant[x-1]%*%t(constant[x-1])))} })
        })
        
       
    return(sigma)
}


##------------------LogLikelihood function for multiple rates per traits------##
lik.Mult<-function(par,dat,C,D,index.mat,sig,error, p, k, n, precalcMat, method,constraint,theta_mle=TRUE,theta=NULL,istrend=FALSE){
    if(istrend==TRUE){
        trend_val<-Vdiag*(D%*%par[nparsig+seq_len(npartrend)][tr_index])
        if(NA_val==TRUE) trend_val<-trend_val[-Indice_NA]
        theta_mle<-FALSE
        theta<-par[nparsig+npartrend+seq_len(npartheta)]
    }
  loglik<-bm_fun_matrix(C,buildSigma(par[seq_len(nparsig)],index.mat,sig,model,constraint),dat,D,precalcMat,n,p,error,method,theta_mle,theta,istrend,trend_val)
  
  list(loglik=-loglik$logl, ancstate=loglik$anc)
  
}

##---------------------Loglik BM1---------------------------------------------##
lik.BM1<-function(par,dat,C,D,error,method,precalcMat,n,p, constraint,theta_mle=TRUE,theta=NULL,istrend=FALSE){ ##
    if(istrend==TRUE){
        trend_val<-Vdiag*(D%*%par[nparsig+seq_len(npartrend)][tr_index])
         if(NA_val==TRUE) trend_val<-trend_val[-Indice_NA]
        theta_mle<-FALSE
        theta<-par[nparsig+npartrend+seq_len(npartheta)]
    }
  loglik<-bm_fun_matrix(C,buildSigma(par[seq_len(nparsig)],NULL,NULL,model,constraint),dat,D,precalcMat,n,p,error,method,theta_mle,theta,istrend,trend_val)
  
  list(loglik=-loglik$logl, ancstate=loglik$anc)
}


if(model=="BMM"){
        if(param$constraint=="default"){
##---------------------Optimization BMM---------------------------------------##
        # number of parameters
        if(decomp=="user"){npar=length(unique(as.numeric(index.user[!is.na(index.user)])))}else{npar=(p*(p+1)/2)}
        
        # sigma matrix
        sig<-matrix(1,k,npar)

        # index matrix of rates
        index.mat<-matrix(1:length(sig),k,npar,byrow=TRUE)

        # initial values for the optimizer
        if(is.null(param[["sigma"]])==TRUE){
            sig1<-startParamSigma(p, decomp, tree, data, index.user=index.user)
            starting<-unlist(lapply(1:k,function(x){sig1}))
        }else{
            if(length(param$sigma[[1]])==npar){
                starting<-unlist(param$sigma)
            }else if(isSymmetric(param$sigma[[1]])){
                starting<-unlist(lapply(1:length(param$sigma),function(x){
                    startParamSigma(p, decomp, tree, data, guess=param$sigma[[x]], index.user=index.user)}))
            }else{
                starting<-unlist(lapply(1:length(param$sigma),function(x){sym.unpar(param$sigma[[x]])}))
            }
            if(length(starting)!=(k*npar)){stop("The number of starting values for the rate matrix do not match, see ?mvBM for providing user specified starting values")}
        }
        
        # length of sigma parameters
        nparsig<-length(starting)
        starting<-c(starting,trend_par)
        
##---------------------Optimization BMM diagonal------------------------------##
        }else if(param$constraint=="diagonal" | param$constraint=="equaldiagonal"){
            
            # default decomp
            if(decomp!="diagonal" & decomp!="equaldiagonal"){ decomp<-param$decomp<-constraint}

            # number of parameters
            if(decomp=="diagonal"){npar=p}else{npar=1}
            
            # sigma matrix
            sig<-matrix(1,k,npar)
            
            # index matrix of rates
            index.mat<-matrix(1:length(sig),k,npar,byrow=TRUE)
            
            # initial values for the optimizer
            if(is.null(param[["sigma"]])==TRUE){
                sig1<-startParamSigma(p, decomp, tree, data)
                starting<-unlist(lapply(1:k,function(x){sig1}))
            }else{
                if(length(param$sigma[[1]])==npar){
                    starting<-unlist(param$sigma)
                }else if(isSymmetric(param$sigma[[1]])){
                    starting<-unlist(lapply(1:length(param$sigma),function(x){
                        startParamSigma(p, decomp, tree, data, guess=param$sigma[[x]], index.user=index.user)}))
                }else{
                    starting<-unlist(lapply(1:length(param$sigma),function(x){diag(param$sigma[[x]])}))
                }
                if(length(starting)!=(k*npar)){stop("The number of starting values for the rate matrix do not match, see ?mvBM for providing user specified starting values")}
            }
            
            # length of sigma parameters
            nparsig<-length(starting)
            starting<-c(starting,trend_par)
##---------------------Optimization BMM shared eigenvectors-------------------##
        }else if(param$constraint=="shared"){
            
            # number of parameters
            npar=(p*(p-1)/2)
            # sigma matrix
            sig<-matrix(1,k,p)
            
            # index matrix of rates
            index.mat<-matrix(1:length(sig),k,p,byrow=TRUE)
            
            # initial values for the optimizer
            if(is.null(param[["sigma"]])==TRUE){
                sig1<-varBM(tree,data,n,p)
                eigval<-eigen(sig1)$values
                starting<-numeric(npar)
                varval<-unlist(lapply(1:k,function(x){log(eigval)}))
                starting<-c(starting,varval)
                
            }else{
                if(length(param$sigma[[1]])==npar+p){
                    starting<-sym.unpar_off(sym.par(param$sigma[[1]]))
                    varval<-unlist(lapply(1:length(param$sigma),function(x){log(eigen(sym.par(param$sigma[[1]]))$values)}))
                }else if(isSymmetric(param$sigma[[1]])){
                    starting<-sym.unpar_off(param$sigma[[1]])
                    varval<-unlist(lapply(1:length(param$sigma),function(x){log(eigen(param$sigma[[1]])$values)}))
                }else{
                    starting<-NULL
                    varval<-unlist(lapply(1:length(param$sigma),function(x){
                        startParamSigma(p, "eigen+", tree, data, guess=param$sigma[[x]], index.user=index.user)}))

                }
                starting<-c(starting,varval)
                if(length(starting)!=p*k+npar){stop("The number of starting values for the rate matrix do not match, see ?mvBM for providing user specified starting values")}
            }
            # length of sigma parameters
            nparsig<-length(starting)
            starting<-c(starting,trend_par)

##---------------------Optimization BMM shared variance----------------------##
        }else if(param$constraint=="variance"){
            
            # number of parameters
            npar=(p*(p-1)/2)
            # sigma matrix
            sig<-matrix(1,k,npar)
            
            # index matrix of rates
            index.mat<-matrix(1:length(sig),k,npar,byrow=TRUE)
            
            # initial values for the optimizer
            if(is.null(param[["sigma"]])==TRUE){
                sig1<-varBM(tree,data,n,p)
                varval<-diag(sig1)
                sig1<-sym.unpar_off(sig1)
                starting<-unlist(lapply(1:k,function(x){sig1}))
                starting<-c(varval,starting)
                
            }else{
                if(length(param$sigma[[1]])==npar+p){
                    varval<-diag(sym.par(param$sigma[[1]]))
                    starting<-lapply(1:length(param$sigma),function(x){sym.unpar_off(sym.par(param$sigma[[x]]))})
                }else if(isSymmetric(param$sigma[[1]])){
                    varval<-NULL
                    starting<-unlist(lapply(1:length(param$sigma),function(x){
                        startParamSigma(p, decomp, tree, data, guess=param$sigma[[x]], index.user=index.user)}))
                }else{
                    starting<-unlist(lapply(1:length(param$sigma),function(x){sym.unpar_off(param$sigma[[x]])}))
                    varval<-diag(param$sigma[[1]])
                }
                starting<-c(varval,starting)
                if(length(starting)!=k*npar+p){stop("The number of starting values for the rate matrix do not match, see ?mvBM for providing user specified starting values")}
            }
            # length of sigma parameters
            nparsig<-length(starting)
            starting<-c(starting,trend_par)
        
##---------------------Optimization BMM shared eigenvectors-------------------##
        }else if(param$constraint=="correlation"){
            
            # number of parameters
            npar=(p*(p-1)/2)
                
            # sigma matrix
            sig<-matrix(1,k,p)
            
            # index matrix of rates
            index.mat<-matrix(1:length(sig),k,p,byrow=TRUE)
            
            # initial values for the optimizer
             if(is.null(param[["sigma"]])==TRUE){
                 sig1<-varBM(tree,data,n,p)
                 starting<-sym.unpar_off(sig1)
                 varval<-unlist(lapply(1:k,function(x){diag(sig1)}))
                 starting<-c(starting,varval)
             
                }else{
                if(length(param$sigma[[1]])==npar+p){
                    starting<-sym.unpar_off(sym.par(param$sigma[[1]]))
                    varval<-lapply(1:length(param$sigma),function(x){diag(sym.par(param$sigma[[1]]))})
                }else if(isSymmetric(param$sigma[[1]])){
                    starting<-NULL
                    varval<-unlist(lapply(1:length(param$sigma),function(x){
                        startParamSigma(p, decomp, tree, data, guess=param$sigma[[x]], index.user=index.user)}))
                }else{
                    starting<-sym.unpar_off(param$sigma[[1]])
                    varval<-unlist(lapply(1:length(param$sigma),function(x){diag(param$sigma[[x]])}))
                }
                starting<-c(starting,varval)
                if(length(starting)!=p*k+npar){stop("The number of starting values for the rate matrix do not match, see ?mvBM for providing user specified starting values")}
                }
            
            # length of sigma parameters
            nparsig<-length(starting)
            starting<-c(starting,trend_par)

            
##---------------------Optimization BMM proportional--------------------------##
        }else if(param$constraint=="proportional"){
            
              npar=(p*(p+1)/2)
              proportionalval<-rep(1,k-1)
            # initial values for the optimizer
             if(is.null(param[["sigma"]])==TRUE){
                starting<-startParamSigma(p, decomp, tree, data)
            }else{
                if(length(param$sigma)==npar){
                    starting<-param$sigma
                
                }else if(isSymmetric(param$sigma[[1]])){
                    starting<-unlist(lapply(1:length(param$sigma),function(x){
                        startParamSigma(p, decomp, tree, data, guess=param$sigma[[x]], index.user=index.user)}))
                    proportionalval<-NULL
                }else{
                    starting<-startParamSigma(p, decomp, tree, data, param$sigma)
                }
                if(length(starting)!=(npar)){stop("The number of starting values for the rate matrix do not match, see ?mvBM for providing user specified starting values")}
            }
         
            starting<-c(proportionalval,starting)
            index.mat<-NULL
            
            # length of sigma parameters
            nparsig<-length(starting)
            starting<-c(starting,trend_par)
       
        }else if(param$constraint=="equal"){
##---------------------Optimization BMM constrained---------------------------##
        # number of parameters for the constrained model
        npar=(p*(p-1)/2)+1
        # sigma matrix
        sig<-matrix(1,k,npar)
        # index matrix
        index.mat<-matrix(1:length(sig),k,npar,byrow=TRUE)
        
        if(is.null(param[["sigma"]])==TRUE){
            # # m?me chose mais pour chaque regimes dans la matrice index
            valstart<-startParamSigma(p, "equal", tree, data)
            starting<-NULL
            for(i in 1:k){
                starting<-c(starting,valstart)
            }
        }else{ ## starting values are provided
            if(length(param$sigma[[1]])==npar){
                starting<-unlist(param$sigma)
            }else{
                starting<-unlist(lapply(1:length(param$sigma),function(x){c(param$sigma[[x]][[1]],param$sigma[[x]][lower.tri(param$sigma[[x]])])}))
            }
            if(length(starting)!=(k*npar)){stop("The number of starting values for the rate matrix do not match, see ?mvBM for providing user specified starting values")}
        }
        # length of sigma parameters
        nparsig<-length(starting)
        starting<-c(starting,trend_par)

    }
  
  # Number of parameters
  if(param$smean==TRUE){
      npartheta<-p
      nparam=nparsig+npartheta+npartrend
  }else{
      npartheta<-k*p
      nparam=nparsig+npartheta+npartrend
  }
##-------------------Function log-likelihood----------------------------------##

llfun <- function(par,...){
    
    args <- list(...)
    if(is.null(args[["root.mle"]])) args$root.mle <- TRUE
    if(is.null(args[["theta"]])) args$theta <- FALSE
    
    if(args$root.mle==TRUE){
        
            result <- -as.numeric(lik.Mult(par=par, dat=data, C=C2, D=D, index.mat=index.mat,sig=sig, error=error, p=p, k=k, n=n, precalcMat=precalcMat, method=method, constraint=constraint,theta_mle=TRUE,theta=NULL,istrend=istrend)$loglik)
            
            if(args$theta==TRUE){
                theta <- as.numeric(lik.Mult(par=par, dat=data, C=C2, D=D, index.mat=index.mat,sig=sig, error=error, p=p, k=k, n=n, precalcMat=precalcMat, method=method, constraint=constraint,theta_mle=TRUE,theta=NULL,istrend=istrend)$ancstate)
                result<-list(logl=result, theta=theta)
            }
            
                }else{
                    
            result <- -as.numeric(lik.Mult(par=par, dat=data, C=C2, D=D, index.mat=index.mat,sig=sig, error=error, p=p, k=k, n=n, precalcMat=precalcMat, method=method, constraint=constraint,theta_mle=FALSE,theta=par[nparsig+npartrend+seq_len(npartheta)],istrend=istrend)$loglik)
                    
        }
        return(result)
    }

# attributes to the loglik function
attr(llfun, "model")<-"BM"
attr(llfun, "sigma")<-nparsig
if(istrend==TRUE) attr(llfun, "trend")<-npartrend
attr(llfun, "theta")<-npartheta


sigmafun <- function(par) {buildSigma(par,index.mat=index.mat,sig=sig,model=model,constraint=constraint)}
            
## Check first if we want to return the log-likelihood function only
if(optimization=="fixed"){
    message("No optimization performed, only the Log-likelihood function is returned with default parameters.")
    param$sigmafun<-sigmafun
    param$model<-model
    param$constraint<-constraint
    param$nparam<-nparam
    param$nbspecies<-n
    param$ntraits<-p
    param$nregimes<-k
    param$method<-method
    param$optimization<-optimization
    param$traits<-colnames(data)
    theta.mat<-rep(0,p)

    matResults<-buildSigma(starting,index.mat,sig,model,constraint)
        resultList<-array(dim = c(p, p, k))
        states=vector()
        for(i in 1:k){
            resultList[,,i]<-matResults[[i]]
            states[i]<-colnames(tree$mapped.edge)[i]
        }
    dimnames(resultList)<-list(colnames(data), colnames(data), states)
    if(istrend==FALSE){
        results<-list(llik=llfun, theta=theta.mat, sigma=resultList, param=param)
    }else{
        results<-list(llik=llfun, trend=trend_par, theta=theta.mat, sigma=resultList, param=param)
    }
    class(results)<-c("mvmorph")
    invisible(results)
    return(results)
    }

# Optimizer
if(optimization!="subplex"){
estim<-optim(par=starting,fn=function (par) { lik.Mult(par=par, dat=data, C=C2, D=D, index.mat=index.mat, sig=sig, error=error, p=p, k=k, n=n, precalcMat=precalcMat, method=method, constraint=constraint, istrend=istrend)$loglik },control=control,hessian=TRUE,method=optimization)   #mettre les options maxit et method dans le menu
}else{
estim<-subplex(par=starting,fn=function (par){lik.Mult(par=par, dat=data, C=C2, D=D, index.mat=index.mat, sig=sig, error=error, p=p, k=k, n=n, precalcMat=precalcMat, method=method, constraint=constraint, istrend=istrend)$loglik},control=control,hessian=TRUE)   #mettre les options maxit et method dans le menu
}

}else if(model=="BM1"){
##---------------------Optimization BM1---------------------------------------##
    if(constraint=="default"){
        # number of parameters
        if(decomp=="user"){npar=length(unique(as.numeric(index.user[!is.na(index.user)])))}else{npar=(p*(p+1)/2)}
        
        # initial values for the optimizer
        if(is.null(param[["sigma"]])==TRUE){
            starting<-startParamSigma(p, decomp, tree, data, index.user=index.user)
        }else{
            if(length(param$sigma)==npar){
                starting<-param$sigma
            }else{
                starting<-startParamSigma(p, decomp, tree, data, guess=param$sigma, index.user=index.user)
            }
            if(length(starting)!=(npar)){stop("The number of starting values for the rate matrix do not match, see ?mvBM for providing user specified starting values")}
        }
        # length of sigma parameters
        nparsig<-length(starting)
        starting<-c(starting,trend_par)
        
    }else if(constraint=="diagonal" | constraint=="equaldiagonal"){
##---------------------Optimization BM1 diagonal------------------------------##
        if(decomp!="diagonal" & decomp!="equaldiagonal"){ decomp<-param$decomp<-constraint}
        if(decomp=="diagonal"){npar=p}else{npar=1}
        
        if(is.null(param[["sigma"]])==TRUE){
            starting<-startParamSigma(p, decomp, tree, data)
            
        }else{
            if(length(param$sigma)==npar){
                starting<-param$sigma
            }else{
                starting<-startParamSigma(p, decomp, tree, data, guess=param$sigma, index.user=index.user)
                #starting<-c(param$sigma[[1]][[1]],param$sigma[[1]][lower.tri(param$sigma[[1]])])
            }
        if(length(starting)!=(npar)){stop("The number of starting values for the rate matrix do not match, see ?mvBM for providing user specified starting values")}
    
        }
        # length of sigma parameters
        nparsig<-length(starting)
        starting<-c(starting,trend_par)

    }else if(param$constraint=="equal"){
##---------------------Optimization BM1 constrained---------------------------##
        #Same as Adams 2012; for p>2 I use the spherical parameterization in the MEE paper
        npar=(p*(p-1)/2)+1
        
        if(is.null(param[["sigma"]])==TRUE){
            starting<-startParamSigma(p, "equal", tree, data)
            
        }else{
            if(length(param$sigma)==npar){
                starting<-param$sigma
            }else{
                starting<-startParamSigma(p, decomp, tree, data, guess=param$sigma, index.user=index.user)
                #starting<-c(param$sigma[[1]][[1]],param$sigma[[1]][lower.tri(param$sigma[[1]])])
            }
            if(length(starting)!=(npar)){stop("The number of starting values for the rate matrix do not match, see ?mvBM for providing user specified starting values")}
            
        }
        # length of sigma parameters
        nparsig<-length(starting)
        starting<-c(starting,trend_par)
    }

# Number of parameters
if(param$smean==TRUE){
    npartheta<-p
    nparam=nparsig+npartheta+npartrend
}else{
    npartheta<-k*p
    nparam=nparsig+npartheta+npartrend
}

##-------------------Function log-likelihood----------------------------------##
llfun <- function(par,...){
    
    args <- list(...)
    if(is.null(args[["root.mle"]])) args$root.mle <- TRUE
    if(is.null(args[["theta"]])) args$theta <- FALSE
    
    if(args$root.mle==TRUE){
            
            result <- -as.numeric(lik.BM1(par=par,dat=data,C=C1,D=D,error=error,method=method,precalcMat=precalcMat,n=n,p=p,constraint=constraint,theta_mle=TRUE,theta=NULL,istrend=istrend)$loglik)
            
            if(args$theta==TRUE){
                theta <- as.numeric(lik.BM1(par=par,dat=data,C=C1,D=D,error=error,method=method,precalcMat=precalcMat,n=n,p=p,constraint=constraint,theta_mle=TRUE,theta=NULL, istrend=istrend)$ancstate)
                result<-list(logl=result, theta=theta)
            }
            
        }else{
            
            result <- -as.numeric(lik.BM1(par=par,dat=data,C=C1,D=D,error=error,method=method,precalcMat=precalcMat,n=n,p=p,constraint=constraint,theta_mle=FALSE,theta=par[nparsig+npartrend+seq_len(npartheta)],istrend=istrend)$loglik)
            
        }
        return(result)
    }
# attributes to the loglik function
attr(llfun, "model")<-"BM"
attr(llfun, "sigma")<-nparsig
if(istrend==TRUE) attr(llfun, "trend")<-npartrend
attr(llfun, "theta")<-npartheta


sigmafun <- function(par){ buildSigma(par,NULL,NULL,model,constraint)}
    
    
## Check first if we want to return the log-likelihood function only
if(optimization=="fixed"){
    message("No optimization performed, only the Log-likelihood function is returned with default parameters.")
    param$sigmafun<-sigmafun
    param$model<-model
    param$constraint<-constraint
    param$nparam<-nparam
    param$nbspecies<-n
    param$ntraits<-p
    param$nregimes<-k
    param$method<-method
    param$optimization<-optimization
    param$traits<-colnames(data)
    theta.mat<-rep(0,p)

    resultList<-buildSigma(starting,NULL,NULL,model,constraint)
    colnames(resultList)<-colnames(data)
    rownames(resultList)<-colnames(data)
    
    if(istrend==FALSE){
        results<-list(llik=llfun, theta=theta.mat, sigma=resultList, param=param)
    }else{
        results<-list(llik=llfun, trend=trend_par, theta=theta.mat, sigma=resultList, param=param)
    }
    
    class(results)<-c("mvmorph")
    invisible(results)
    return(results)
}

# Optimizer
if(optimization!="subplex"){
estim<-optim(par=starting,fn=function(par){lik.BM1(par=par,dat=data,C=C1,D=D,error=error,method=method,precalcMat=precalcMat,n=n,p=p,constraint=constraint,theta_mle=TRUE,theta=NULL,istrend=istrend)$loglik},control=control,hessian=TRUE,method=optimization)
}else{
estim<-subplex(par=starting,fn=function(par){lik.BM1(par=par,dat=data,C=C1,D=D,error=error,method=method,precalcMat=precalcMat,n=n,p=p,constraint=constraint,theta_mle=TRUE,theta=NULL,istrend=istrend)$loglik},control=control,hessian=TRUE)
} 
}

##-----------------Summarizing results----------------------------------------##
# names for the plot
if(is.null(colnames(data))){
    names_data_matrix<-rep("",p)
    names_data<-list(names_data_matrix,names_data_matrix)
}else{
    names_data_matrix<-colnames(data)
    names_data<-list(names_data_matrix,names_data_matrix)
}

if(model=="BMM"){
matResults<-buildSigma(estim$par[seq_len(nparsig)],index.mat,sig,model,constraint)
resultList<-array(dim = c(p, p, k))
states=vector()
 for(i in 1:k){
 resultList[,,i]<-matResults[[i]]
 states[i]<-colnames(tree$mapped.edge)[i]#multi.tre[[i]]$state
}
 
dimnames(resultList)<-list(names_data_matrix, names_data_matrix, states)

#ancestral states estimates
theta.mat<-lik.Mult(par=estim$par,dat=data,C=C2,D=D,index.mat=index.mat,sig=sig,error=error, p=p, k=k, n=n, precalcMat=precalcMat, method=method,constraint=constraint,theta_mle=TRUE,theta=NULL,istrend=istrend)$ancstate
if(param$smean==TRUE){
    theta.mat<-matrix(theta.mat,nrow=1)
    colnames(theta.mat)<-names_data_matrix
    rownames(theta.mat)<-"theta:"
}else{
    theta.mat<-matrix(theta.mat,nrow=k)
    colnames(theta.mat)<-names_data_matrix
    rownames(theta.mat)<-colnames(tree$mapped.edge)
}

}else if(model=="BM1"){
resultList<-buildSigma(estim$par[seq_len(nparsig)],NULL,NULL,model,constraint)
dimnames(resultList)<-names_data

 #ancestral states estimates
theta.mat<-matrix(lik.BM1(par=estim$par,dat=data,C=C1,D=D,error=error, p=p, n=n, precalcMat=precalcMat, method=method,constraint=constraint,theta_mle=TRUE,theta=NULL,istrend=istrend)$ancstate,nrow=1)
if(param$smean==TRUE){
    theta.mat<-matrix(theta.mat,nrow=1)
    colnames(theta.mat)<-names_data_matrix
    rownames(theta.mat)<-"theta:"
}else{
    theta.mat<-matrix(theta.mat,nrow=k)
    colnames(theta.mat)<-names_data_matrix
    rownames(theta.mat)<-colnames(tree$mapped.edge)
}
}

# trend matrix
if(istrend==TRUE){
    if(param$smean==TRUE){
        trend.mat<-matrix(estim$par[nparsig+seq_len(npartrend)][tr_index], nrow=1)
        rownames(trend.mat)<-c("drift:")
    }else{
        trend.mat<-matrix(estim$par[nparsig+seq_len(npartrend)][tr_index], nrow=k)
        rownames(trend.mat)<-colnames(tree$mapped.edge)
    }
    colnames(trend.mat)<-names_data_matrix
}

# LogLikelihood
LL<--estim$value

# AIC
AIC<--2*LL+2*nparam

# AIC corrected
nobs <- length(which(!is.na(data)))
AICc<-AIC+((2*nparam*(nparam+1))/(nobs-nparam-1)) #Hurvich et Tsai, 1989
# Maybe n need to be changed by length(data)? Moreover it can change when there is missing cases
##---------------------Diagnostics--------------------------------------------##

if(estim$convergence==0 & diagnostic==TRUE){  
cat("successful convergence of the optimizer","\n")
}else if(estim$convergence==1 & diagnostic==TRUE){  
cat("\n","maximum limit iteration has been reached, please consider increase maxit","\n")
}else if(diagnostic==TRUE){  
cat("\n","convergence of the optimizer has not been reached, try simpler model","\n") 
}
# Hessian eigen decomposition to check the derivatives
hess<-eigen(estim$hessian)$values
if(any(hess<0)){
hess.value<-1
if(diagnostic==TRUE){
cat("unreliable solution has been reached, check hessian eigenvectors or try simpler model","\n")}
}else{
hess.value<-0
if(diagnostic==TRUE){
cat("a reliable solution has been reached","\n")}
}

##-------------------Print results--------------------------------------------##
if(echo==TRUE){
cat("\n")
if(constraint==TRUE | constraint=="diagonal" | constraint=="equal"){
cat("-- Summary results for constrained rate ",model,"model --","\n")
}else if(constraint=="default" & decomp=="user"){
cat("-- Summary results for user-defined",model,"constrained model --","\n")
}else if(constraint=="correlation"){
cat("-- Summary results for common correlation ",model,"model --","\n")
}else if(constraint=="shared"){
cat("-- Summary results for shared eigenvectors ",model,"model --","\n")
}else if(constraint=="proportional"){
cat("-- Summary results for proportional rate matrices ",model,"model --","\n")
}else if(constraint=="variance"){
cat("-- Summary results for common variance",model,"model --","\n")
}else{
cat("-- Summary results for multiple rate",model,"model --","\n")
}
cat("LogLikelihood:","\t",LL,"\n")
cat("AIC:","\t",AIC,"\n")
cat("AICc:","\t",AICc,"\n")
cat(nparam,"parameters","\n")
cat("\n")
cat("Estimated rate matrix","\n")
cat("______________________","\n")
print(resultList)
cat("\n")
cat("Estimated root state","\n")
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
param$model<-model
param$constraint<-constraint
param$nparam<-nparam
param$nbspecies<-n
param$ntraits<-p
param$nregimes<-k
param$method<-method
param$optimization<-optimization
param$traits<-colnames(data)
# param$loglik<-llfun # Add it directly to the results?
param$sigmafun<-sigmafun
param$opt<-estim
class(llfun) = c("mvmorph.llik")
##-------------------Store results--------------------------------------------##
if(istrend==FALSE){
 results<-list(LogLik=LL, AIC=AIC, AICc=AICc, theta=theta.mat, sigma=resultList ,convergence=estim$convergence, hess.values=hess.value, param=param, llik=llfun)
}else{
 results<-list(LogLik=LL, AIC=AIC, AICc=AICc, theta=theta.mat, sigma=resultList , trend=trend.mat, convergence=estim$convergence, hess.values=hess.value, param=param, llik=llfun)
}

class(results)<-c("mvmorph","mvmorph.bm")
invisible(results)
#End
}
