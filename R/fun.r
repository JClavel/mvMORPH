################################################################################
##                                                                            ##
##                               mvMORPH: fun.r                               ##
##                                                                            ##
##   Internal functions for the mvMORPH package                               ##
##                                                                            ##
##  Created by Julien Clavel - 16-07-2013 - updated 10-03-2015                ##
##  (julien.clavel@hotmail.fr/ julien.clavel@biologie.ens.fr)                 ##
##   require: phytools, ape, corpcor, subplex, spam                           ##
##                                                                            ##
################################################################################

##------------------------Fonctions necessaires-------------------------------##

## vcv build for time-series
vcv.ts <- function(times){
   vcv_mat <- outer(times,times, FUN=pmin)
   # avoid problems with integers
   V <- matrix(as.numeric(vcv_mat),length(times))
   return(V)
}

## Matrices parameterizations

# Calcul d'une matrice positive semi-definie (Choleski transform) modified from OUCH by A. King
sym.par<-function (x) {
  nchar <- floor(sqrt(2*length(x)))
  if (nchar*(nchar+1)!=2*length(x)) {
    stop("a symmetric matrix is parameterized by a triangular number of parameters",call.=FALSE)
  }
  y <- matrix(0,nchar,nchar)
  y[lower.tri(y,diag=TRUE)] <- x
  tcrossprod(y)
}

# Inverse of Sym.par, return a vector (from OUCH) by A. King
sym.unpar <- function (x) {
  y <- t(chol(x))
  y[lower.tri(y,diag=TRUE)]
}

sym.unpar_off <- function (x) {
    y <- t(chol(x))
    y[lower.tri(y,diag=FALSE)]
}

# sigma matrix parameterization
symPar <-function(par, decomp="cholesky", p=NULL, index.user=NULL, tol=0.1){

    switch(decomp,
        "cholesky"={
            sigma<-sym.par(par)
        },
        "spherical"={
            dim1<-p*(p-1)/2
            sigma <- .Call(spherical, param=par[1:dim1], variance=par[(dim1+1):(dim1+p)], dim=as.integer(p))
        },
        "eigen"={
            dim1<-p*(p-1)/2
            Q<-.Call(givens_ortho, Q=diag(p), angle=par[1:dim1], ndim=as.integer(p))
            T<-par[(dim1+1):(dim1+p)]
            invQ<-t(Q)
            sigma<-Q%*%diag(T)%*%invQ
        },
        "eigen+"={
            dim1<-p*(p-1)/2
            Q<-.Call(givens_ortho, Q=diag(p), angle=par[1:dim1], ndim=as.integer(p))
            T<-exp(par[(dim1+1):(dim1+p)])
            invQ<-t(Q)
            sigma<-Q%*%diag(T)%*%invQ
        },
        "diagonal"={
            d_par<-diag(par,p)
            sigma<-d_par%*%t(d_par)
        },
        "equaldiagonal"={
            d_par<-diag(par,p)
            sigma<-d_par%*%t(d_par)
        },
        "equal"={
            sigma<-tcrossprod(build.chol(par,p))
            # avoid issues when the matrix is not spd?
            # with the Adams parameterization I always obtain problems with more than 2 traits (even with the original code)
            if(p>2){
                 dim1<-p*(p-1)/2
                 sigma <- .Call(spherical, param=par[1:dim1], variance=rep(par[(dim1+1)],p), dim=as.integer(p))
            }
            
        },
        "user"={
            # user defined sigma matrix
            sigma <- matrix(0,p,p)
            sigma[] <- c(par)[index.user]
            sigma[is.na(sigma)] <- 0
            dval <- diag(sigma)
            diag(sigma) <- diag(dval%*%t(dval))
            #eig <- min(eigen(sigma)$values)
            eig <- min(eigen(sigma)$values)
            if(eig<0){
                # arbitrary tolerance value to shift the smallest eigenvalue
                sigma <- sigma+(abs(eig)+tol)*diag(p)
            }
        },
        # set the default error message
        stop("Sorry, the \"decomp\" or \"decompSigma\" option in the \"param\" list is wrong...","\n")
        )
        
        return(sigma)
}


# default matrix parameterization for Sigma
startParamSigma <- function(p, matrix, tree, data, guess=NULL, index.user=NULL){
    
    if(inherits(tree,"phylo")){
        if(!is.null(guess)){
            n <- length(tree$tip.label)
            sigma <- guess
            
        }else{
            n <- length(tree$tip.label)
            sigma <- varBM(tree,data,n,p)
        }
    }else{
        if(!is.null(guess)){
            sigma <- guess
        }else{
            sigma <- cov(as.matrix(data),use="complete.obs")/max(tree)
        }
    }
    
    # off-diagonal dimension
    dim1<-p*(p-1)/2
    
    switch(matrix,
    "diagonal"={ # Diagonal parameterization
        value<-sqrt(diag(sigma))
    },
    "equaldiagonal"={ # Diagonal parameterization
        value<-sqrt(mean(diag(sigma)))
    },
    "equal"={ # Cholesky constraint / spherical parameterization
        if(p>2){
        value<-vector("numeric",dim1+1)
        value[1:dim1]<-rep(pi/2,dim1)
        value[dim1+1]<-sqrt(mean(diag(sigma)))
        }else{
        value<-vector("numeric",dim1+1)
        value[1]<-mean(diag(sigma))
        }
    },
    "eigen"={
        value <- vector("numeric",dim1+p)
        eigval <- eigen(sigma)$value
        value[(dim1+1):(dim1+p)]<-eigval
    },
    "eigen+"={
        value <- vector("numeric",dim1+p)
        eigval <- eigen(sigma)$value
        value[(dim1+1):(dim1+p)]<-log(eigval)
    },
    "spherical"={
        dval <- diag(sigma)
        value <- rep(pi/2,dim1+p)
        value[(dim1+1):(dim1+p)]<-sqrt(dval)
    },
    "cholesky"={
        value <- sym.unpar(sigma)
        
    },
    "user"={
        value <- user_guess(index.user, sigma)
    },
    stop("Sorry, the sigma \"decomp\" option in the \"param\" list is wrong...","\n")
    )
   
    return(value)
}


# default matrix parameterization for A
startParam <- function(p, matrix, tree, index.user=NULL, ...){
    
    args <- list(...)
    
    if(is.null(args[["hlife"]])){
        if(inherits(tree,"phylo")){
            if(is.ultrametric(tree)){
                times <- branching.times(tree)
            }else{
                times <- nodeHeights(tree)
            }
        }else{
            times <- max(tree)
        }
    # set default hlife to 1/3
        hlife <- log(2)/(max(times)/(2+seq_len(p)))
    }else{
        hlife <- args$hlife
    }
    
    # off-diagonal dimension
    dim1<-p*(p-1)/2
    
    switch(matrix,
    "svd"={ # SVD decomposition
        alpha <- vector("numeric",p*p)
        alpha[(dim1+1):(dim1+p)]<-hlife
    },
    "svd+"={ # SVD decomposition
        alpha <- vector("numeric",p*p)
        alpha[(dim1+1):(dim1+p)]<-log(hlife)
    },
    "eigen"={ # eigen decomposition
        alpha <- vector("numeric",dim1+p)
        alpha[(dim1+1):(dim1+p)]<-hlife
    },
    "qr"={ # QR decomposition with diagonal values forced to be positive for uniqueness
        alpha <- vector("numeric",p*p)
        alpha[(dim1+1):(dim1+p)]<-hlife
    },
    "eigen+"={ # eigen decomposition
        alpha <- vector("numeric",dim1+p)
        alpha[(dim1+1):(dim1+p)]<-log(hlife)
    },
    "qr+"={ # QR decomposition with diagonal values forced to be positive for uniqueness
        alpha <- vector("numeric",p*p)
        alpha[(dim1+1):(dim1+p)]<-log(hlife)
    },
    "spherical"={ # Spherical parameterization for correlation matrix + scaling
        alpha <- rep(pi/2,dim1+p)
        alpha[(dim1+1):(dim1+p)]<-sqrt(hlife)
    },
    "cholesky"={ # Cholesky decomposition
        alpha <- diag(p)
        diag(alpha)<-sqrt(hlife)
        alpha<-alpha[lower.tri(alpha,diag=TRUE)]
    },
    "schur"={ # Schur decomposition
        alpha <- vector("numeric",p*p)
        alpha[(dim1+1):(dim1+p)]<-hlife
    },
    "schur+"={ # Schur decomposition
        alpha <- vector("numeric",p*p)
        alpha[(dim1+1):(dim1+p)]<-log(hlife)
    },
    "equal"={
        alpha<-vector("numeric",dim1+1)
        alpha[1]<-mean(hlife)
    },
    "equaldiagonal"={
        alpha<-mean(hlife)
    },
    "diagonal"={
        alpha <- vector("numeric",p)
        alpha<-hlife
    },
    "diagonalPositive"={
        alpha <- vector("numeric",p)
        alpha<-log(hlife)
    },
    "lower"={
        alpha <- vector("numeric",p+dim1)
        alpha[(dim1+1):(dim1+p)]<-log(hlife)
    },
    "upper"={
        alpha <- vector("numeric",p+dim1)
        alpha[(dim1+1):(dim1+p)]<-log(hlife)
    },
    "univariate"={
        alpha <- hlife
    },
    "user"={
        dim_user <- length(unique(c(index.user[!is.na(index.user)])))
        alpha <- vector("numeric",dim_user)
        alpha[1:p] <- log(hlife)
    },
    # set the default error message
    stop("Sorry, the \"decomp\" option in the \"param\" list is wrong...","\n")
    )

return(alpha)
}

# alpha matrix parameterization
matrixParam<-function(x,p,matrix="cholesky",index.user=NULL,tol=0.000001){
    
    switch(matrix,
    "svd"={ # svd decomposition
        dim1<-p*(p-1)/2
        U<-.Call(givens_ortho, Q=diag(p), angle=as.numeric(x[1:dim1]), ndim=as.integer(p))
        T<-x[(dim1+1):(dim1+p)]+tol
        V<-.Call(givens_ortho, Q=diag(p), angle=as.numeric(x[(dim1+p+1):(p*p)]), ndim=as.integer(p))
        A<-U%*%diag(T)%*%t(V)
        eigval<-eigen(A)
        Adecomp<-list(vectors=eigval$vectors, values=eigval$values, A=A, invectors=solve(eigval$vectors))
    },
    "svd+"={ # svd decomposition
        dim1<-p*(p-1)/2
        U<-.Call(givens_ortho, Q=diag(p), angle=as.numeric(x[1:dim1]), ndim=as.integer(p))
        T<-exp(x[(dim1+1):(dim1+p)])
        V<-.Call(givens_ortho, Q=diag(p), angle=as.numeric(x[(dim1+p+1):(p*p)]), ndim=as.integer(p))
        A<-U%*%diag(T)%*%t(V)
        eigval<-eigen(A)
        Adecomp<-list(vectors=eigval$vectors, values=eigval$values, A=A, invectors=solve(eigval$vectors))
    },
    "eigen"={ # eigen decomposition with positive eigenvalues
        dim1<-p*(p-1)/2
        Q<-.Call(givens_ortho, Q=diag(p), angle=x[1:dim1], ndim=as.integer(p))
        T<-x[(dim1+1):(dim1+p)]
        invQ<-t(Q)
        A<-Q%*%diag(T)%*%invQ
        Adecomp<-list(vectors=Q, values=T, A=A, invectors=invQ)
    },
    "eigen+"={ # eigen decomposition with positive eigenvalues
        dim1<-p*(p-1)/2
        Q<-.Call(givens_ortho, Q=diag(p), angle=x[1:dim1], ndim=as.integer(p))
        T<-exp(x[(dim1+1):(dim1+p)])
        invQ<-t(Q)
        A<-Q%*%diag(T)%*%invQ
        Adecomp<-list(vectors=Q, values=T, A=A, invectors=invQ)
    },
    "qr"={ # QR decomposition (R diagonal values are forced to be positive)
        dim1<-p*(p-1)/2
        Q<-.Call(givens_ortho, Q=diag(p), angle=x[1:dim1], ndim=as.integer(p))
        R<-diag(x[(dim1+1):(dim1+p)],p)
        R[upper.tri(R)]<-x[(dim1+p+1):(p*p)]
        A<-Q%*%R
        eigval<-eigen(A)
        Adecomp<-list(vectors=eigval$vectors, values=eigval$values, A=A, invectors=solve(eigval$vectors))
    },
    "qr+"={ # QR decomposition (R diagonal values are forced to be positive)
        dim1<-p*(p-1)/2
        Q<-.Call(givens_ortho, Q=diag(p), angle=x[1:dim1], ndim=as.integer(p))
        R<-diag(exp(x[(dim1+1):(dim1+p)]),p)
        R[upper.tri(R)]<-x[(dim1+p+1):(p*p)]
        A<-Q%*%R
        eigval<-eigen(A)
        Adecomp<-list(vectors=eigval$vectors, values=eigval$values, A=A, invectors=solve(eigval$vectors))
    },
    "spherical"={ # Spherical parameterization
        dim1<-p*(p-1)/2
        A <- .Call(spherical, param=x[1:dim1], variance=x[(dim1+1):(dim1+p)], dim=as.integer(p))
        eigval<-eigen(A)
        Adecomp<-list(vectors=eigval$vectors, values=eigval$values, A=A, invectors=solve(eigval$vectors))
    },
    "cholesky"={ # Cholesky decomposition
        y <- matrix(0,p,p)
        y[lower.tri(y,diag=TRUE)] <- x
        A<-tcrossprod(y)
        eigval<-eigen(A)
        Adecomp<-list(vectors=eigval$vectors, values=eigval$values, A=A, invectors=t(eigval$vectors))
    },
    "schur"={ # Schur decomposition
        dim1<-p*(p-1)/2
        Q<-.Call(givens_ortho, Q=diag(p), angle=x[1:dim1], ndim=as.integer(p))
        T<-diag(x[(dim1+1):(dim1+p)],p)
        T[upper.tri(T)]<-x[(dim1+p+1):(p*p)]
        A<-Q%*%T%*%t(Q)
        eigval<-eigen(A)
        Adecomp<-list(vectors=eigval$vectors, values=eigval$values, A=A, invectors=solve(eigval$vectors))
    },
    "schur+"={ # Schur decomposition with positive eigenvalues
        dim1<-p*(p-1)/2
        Q<-.Call(givens_ortho, Q=diag(p), angle=x[1:dim1], ndim=as.integer(p))
        T<-diag(exp(x[(dim1+1):(dim1+p)]),p) # exp to force positive eigenvalues
        T[upper.tri(T)]<-x[(dim1+p+1):(p*p)]
        A<-Q%*%T%*%t(Q)
        eigval<-eigen(A)
        Adecomp<-list(vectors=eigval$vectors, values=eigval$values, A=A, invectors=solve(eigval$vectors))
    },
    "equaldiagonal"={
        A<-diag(x,p)
        eigval<-eigen(A)
        Adecomp<-list(vectors=eigval$vectors, values=eigval$values, A=A, invectors=t(eigval$vectors))
    },
    "diagonal"={
        A<-diag(x,p)
        eigval<-eigen(A)
        Adecomp<-list(vectors=eigval$vectors, values=eigval$values, A=A, invectors=t(eigval$vectors))
    },
    "diagonalPositive"={
        A<-diag(exp(x)+tol,p)
        eigval<-eigen(A)
        Adecomp<-list(vectors=eigval$vectors, values=eigval$values, A=A, invectors=t(eigval$vectors))
    },
    "equal"={
        sigma<-tcrossprod(build.chol(x,p))
        # avoid issues when the matrix is not spd?
        # with the Adams parameterization I always obtain problems with more than 2 traits (even with the original code)
        if(p>2){
            dim1<-p*(p-1)/2
            sigma <- .Call(spherical, param=x[1:dim1], variance=rep(x[(dim1+1)],p), dim=as.integer(p))
        }
    },
    "lower"={
        dim1<-p*(p-1)/2
        A<-matrix(0,p,p)
        A[lower.tri(A,diag=FALSE)] <- x[1:dim1]
        diag(A) <- exp(x[(dim1+1):(dim1+p)])
        eigval<-eigen(A)
        Adecomp<-list(vectors=eigval$vectors, values=eigval$values, A=A, invectors=solve(eigval$vectors))
    },
    "upper"={
        dim1<-p*(p-1)/2
        A<-matrix(0,p,p)
        A[upper.tri(A,diag=FALSE)] <- x[1:dim1]
        diag(A) <- exp(x[(dim1+1):(dim1+p)])
        eigval<-eigen(A)
        Adecomp<-list(vectors=eigval$vectors, values=eigval$values, A=A, invectors=solve(eigval$vectors))
    },
    "univariate"={
        Adecomp<-list(vectors=1, values=x, A=x, invectors=1)
    },
    "user"={
        # user defined sigma matrix
        A <- matrix(0,p,p)
        A[] <- c(x)[index.user]
        A[is.na(A)] <- 0
        dval <- diag(A)
        diag(A) <- diag(dval%*%t(dval))
        eigval <- eigen(A)
        eig <- min(Re(eigval$values))
        if(eig<0){
            # arbitrary tolerance value to shift the smallest eigenvalue
            A <- A+(abs(eig)+tol)*diag(p)
            eigval<-eigen(A)
        }
        Adecomp<-list(vectors=eigval$vectors, values=eigval$values, A=A, invectors=solve(eigval$vectors))
    },
    # set the default error message
     stop("Sorry, the \"decomp\" option in the \"param\" list is wrong...","\n","Note: the \"symmetric\", \"symmetricPositive\", \"nsymmetric\", and \"nsymPositive\" options are deprecated")
    )
    return(Adecomp)
}


# Compute factor list / Creation d'une liste de facteurs (to improve)
newList<-function(factorVal,nv){
    val <-list()
    val[[1]] <-factorVal
    
    for(i in 1:(nv-1)){
        ind <- max(val[[i]])
        val[[i+1]] <- factorVal+ind
    }
    
    factorVal<-unlist(val)
    return(factorVal)
}

# Compute the design matrix for multiple mean / Creation de la matrice de variables indicatrices pour plusieurs moyennes

multD<-function(tree,k,nbtip,smean=TRUE){
    if(smean==TRUE){
        # design matrix
        y<-kronecker(diag(k),rep(1,nbtip))
    }else{
        if(is.null(tree[["mapped.edge"]])){
          stop("The specified tree must be in SIMMAP format with different regimes")
        }
        namestip<-sapply(1:nbtip,function(i){
            if(i==Ntip(tree)+1){
                ind <- which(tree$edge[,1]==i)
                names(tree$maps[[ind[1]]][1]) # we take only one of the root?
            }else{
                ind<-which(tree$edge[,2]==i);
                names(tree$maps[[ind]][length(tree$maps[[ind]])])
            }
        })

        group<-as.numeric(as.factor(namestip))
        
        if(k==1){
            ngr<-length(unique(group))
            y=matrix(0,ncol=ngr,nrow=nbtip)
            for (i in 1:nbtip){y[i,group[i]] <- 1}
        }else{
        gr<-newList(group,k)
        ngr<-unique(gr)
        tgr<-length(gr)
        y=matrix(0,tgr,length(ngr))
        for (i in 1:tgr){y[i,gr[i]] <- 1}
        }
    }
return(y)
}


# Compute matrix time to common ancestor
mTime<-function(phy, scale.height){
    vcv<-vcv.phylo(phy)
if(is.ultrametric(phy)){
   if(scale.height==TRUE){
   mdist<-vcv/max(vcv)
   }else{
    mdist<-vcv}
    mSpdist<-rep(max(mdist),length(phy$tip.label))
  }else{
   if(scale.height==TRUE){
   vstand<-vcv/max(vcv)
   }else{vstand<-vcv}
   mSpdist<-diag(vstand)
   mdist<-vstand
  }
  # mcoph<-diag(mdist)-mdist
list(mSpDist=mSpdist, mDist=mdist)
}

# Binary regime coding
regimeList<-function(mm,k,root=TRUE){ # use root="stationary" for ouch...
    nReg=length(mm)
    regime <- matrix(0,nrow=nReg,ncol=k) # remplacer 0?
    for(i in 1:nReg){
        regime[i,mm[i]]<-1
    }
    if(root=="stationary"){
        regime[nReg,]<-0 ## under OU1 implicitely assumed stationary, just use root=TRUE or FALSE (ultrametric trees)
    }
    return(regime)
}

# Set root regime
indiceReg<-function(n,indice,facInd, root=TRUE){
    if(root==TRUE){
        for(i in 1:n){
            val=length(indice[[i]])
            indice[[i]][val+1]<-facInd[facInd=="_root_state"]
        }
    }else{
        for(i in 1:n){
            val=length(indice[[i]])
            indice[[i]][val+1]<-indice[[i]][val]
        }
    }
    return(indice)
}


# Test for polytomies
eval_polytom<-function(tree){
    nb.tip <- length(tree$tip.label)
    nb.node <- tree$Nnode
    if (nb.node != nb.tip - 1) {
        stop("You can't use \"pic\" with polytomies, try instead \"rpf\",\"inverse\",\"pseudoinverse\" or \"sparse\". Otherwise, first transform the tree using the \"multi2di\" function")
    }
}

# Estimate starting values for the variance matrix
varBM<-function(tree,data,n,k){
    if(any(is.na(data))){
        res<-diag(0.1,k)
        return(res)
    }
    nb.tip <- length(tree$tip.label)
    nb.node <- tree$Nnode
    if (nb.node != nb.tip - 1) {
    tree <- multi2di(tree)
    }
    rate<-rep(0,k)
    tree<-reorder.phylo(tree,"postorder")
    value<-list(tree$edge.length)
    res<-.Call(PIC_gen, x=as.vector(as.matrix(data)), n=as.integer(k), Nnode=as.integer(tree$Nnode), nsp=as.integer(n), edge1=as.integer(tree$edge[,1]), edge2=as.integer(tree$edge[,2]), edgelength=value, times=1, rate=rate, Tmax=1, Model=as.integer(13), mu=1, sigma=1)
    return(res[[2]])
}

# Starting values guess for user-defined constraint matrix
user_guess <- function(user_const, sigma){
    diag(sigma) <- sqrt(diag(sigma))
    indvalue <- sort(unique(c(user_const[!is.na(user_const)])))
    rcind <- sapply(indvalue,function(x) which(user_const==x, arr.ind=TRUE))
    guess <- sapply(rcind, function(x) {
        dim_val <- nrow(x)
        mean(sapply(1:dim_val, function(i) sigma[x[i,1],x[i,2]]))})
    return(guess)
}

# Generate random multivariate distributions
rmvnorm_simul<-function(n=1, mean, var, method="cholesky"){
    
    p<-length(mean)
    if (!all(dim(var)==c(p,p)))
    stop("length of ",sQuote("mean")," must equal the dimension of the square matrix ",sQuote("var"))
    
    # try with fast cholesky
    if(method=="cholesky"){
        chol_factor <- try(t(chol(var)), silent = TRUE)
        
      # else use svd
      if(inherits(chol_factor ,'try-error')){
        warning("An error occured with the Cholesky decomposition, the \"svd\" method is used instead")
        s. <- svd(var)
        if (!all(s.$d >= -sqrt(.Machine$double.eps) * abs(s.$d[1]))){
            warning("The covariance matrix is numerically not positive definite")
        }
        R <- t(s.$v %*% (t(s.$u) * sqrt(s.$d)))
        X <- matrix(mean,p,n) + t(matrix(rnorm(n * p), nrow=n )%*%  R)

      } else {
        X <- matrix(mean,p,n) + chol_factor%*%matrix(rnorm(p*n),p,n)
      }
      
    }else{
        s. <- svd(var)
        if (!all(s.$d >= -sqrt(.Machine$double.eps) * abs(s.$d[1]))){
            warning("The covariance matrix is numerically not positive definite")
        }
        R <- t(s.$v %*% (t(s.$u) * sqrt(s.$d)))
        X <- matrix(mean,p,n) + t(matrix(rnorm(n * p), nrow=n )%*%  R)
    }
    
    return(X)
}

# Compute the variance-covariance matrix of tips as well as internal nodes
vcvPhyloInternal <- function(tree){
    nbtip <- Ntip(tree)
    dis <- dist.nodes(tree)
    MRCA <- mrca(tree, full = TRUE)
    M <- dis[as.character(nbtip + 1), MRCA]
    dim(M) <- rep(sqrt(length(M)), 2)
    return(M)
}

# Generate a multi-phylo list for SIMMAP trees
vcvSplit<-function(tree, internal=FALSE){
    if(!inherits(tree,"simmap")) stop("tree should be an object of class \"simmap\".")
    multi.tre<-list()
    class(multi.tre)<-"multiPhylo"
    #Array method - better for memory?
    #C2<-array(dim=c(nrow(C1),ncol(C1),ncol(tree$mapped.edge)))
    C<-list()
    for(i in 1:ncol(tree$mapped.edge)){
        multi.tre[[i]]<-tree
        multi.tre[[i]]$edge.length<-tree$mapped.edge[,i]
        multi.tre[[i]]$state<-colnames(tree$mapped.edge)[i]
        # hack from Liam Revell
        C[[i]]<-if(internal) vcvPhyloInternal(multi.tre[[i]]) else vcv.phylo(multi.tre[[i]])
    }
return(C)
}

# Get the indice of tip species on ancestral vcv
indiceTip <- function(tree, p){
    totsp <- Ntip(tree)+Nnode(tree) # retrieve the root state
    init  <- res <- (1:Ntip(tree))
    if(p==1) return(init)
    for(i in 1:(p-1)){
        init<-init+totsp
        res <- c(res,init)
    }
    return(res)
}

# Return the path from the root to a node
pathToNode <- function(tree, node, root=NULL){
    ntip <- Ntip(tree)
    if(is.null(root)) root <- ntip+1 # ?? always true?
    if(node==root) return(node)
    vector_of_nodes = NULL # Memory allocation may be better...
    increment = 1
    ind = node
    repeat{
        ind <- tree$edge[which(tree$edge[,2]==ind),1]
        vector_of_nodes[increment] = ind
        increment = increment+1
        if(ind==root) break
    }
    return(c(rev(vector_of_nodes),node))
}


# Prepar the epochs and regime lists to compute the weight-matrix of the OU process
prepWOU <- function(tree, n, p, k, model="OUM", root=FALSE){
    
    # Compute root to tip nodes paths
    ntip = Ntip(tree)
    nnode = tree$Nnode
    root2tip <- .Call(seq_root2tipM, tree$edge, ntip, nnode)

    # parameters
    root_node = ntip+1
    ntot = ntip+nnode
    # hack
    if(n>ntip){
        for(i in root_node:ntot){
            root2tip[[i]] <- pathToNode(tree,i)
        }
    }
    
    if(model=="OUM"){
        
        # lineages maps
        valLineage<-sapply(1:n,function(z){
            if(z!=root_node){ # ntip+1 is assumed to be the root
                rev(unlist(
                sapply(1:(length(root2tip[[z]])-1),function(x){vec<-root2tip[[z]][x:(x+1)]; val<-which(tree$edge[,1]==vec[1] & tree$edge[,2]==vec[2]); tree$maps[[val]]<-tree$maps[[val]]},simplify=FALSE)))
                
            }else{
                val <- max(nodeHeights(tree))
                ind <- which(tree$edge[,1]==z)
                names(val) <- names(tree$maps[[ind[1]]][1])
                val
            }
        },simplify=FALSE)
        
        
        # set factors
        if(root==FALSE | root=="stationary"){
            facInd<-factor(colnames(tree$mapped.edge))
        }else if(root==TRUE){
            facInd<-factor(c("_root_state",colnames(tree$mapped.edge)))
        }
        
        # indice factors
        indice<-lapply(1:n,function(z){
            if(z!=root_node){
                rev(unlist(lapply(1:(length(root2tip[[z]])-1),function(x){vec<-root2tip[[z]][x:(x+1)]; val<-which(tree$edge[,1]==vec[1] & tree$edge[,2]==vec[2]); factor(names(tree$maps[[val]]),levels=facInd)})))
            }else{
              if(root==FALSE | root=="stationary"){
                ind <- which(tree$edge[,1]==root_node)
                factor(names(tree$maps[[ind[1]]][1]), levels=facInd)
              }else{
                factor("_root_state", levels=facInd)
              }
            }
        })
        
    }else{ # model is "OU1"
        
        # change the number of regimes to 1
        k <- 1
        
        # lineage maps
        valLineage<-sapply(1:n,function(z){
            if(z!=ntip+1){
                rev(unlist(
                sapply(1:(length(root2tip[[z]])-1),function(x){vec<-root2tip[[z]][x:(x+1)]; val<-which(tree$edge[,1]==vec[1] & tree$edge[,2]==vec[2]); tree$edge.length[val]<-tree$edge.length[val]},simplify=FALSE)))
            }else{
                max(nodeHeights(tree))
            }
        } ,simplify=FALSE)
        
        if(root==TRUE){
            k<-k+1
            facInd<-factor(c("_root_state","theta_1"))
            indice<-lapply(1:n,function(z){ 
              if(z!=root_node){
                as.factor(rep(facInd[facInd=="theta_1"],length(valLineage[[z]])))
              }else{
                factor("_root_state", levels=facInd)
              }
            })
        }else{
            indice<-lapply(1:n,function(z){ as.factor(rep(1,length(valLineage[[z]])))})
        }
        
    }
    
    # Liste avec dummy matrix
    if(root==FALSE){
        indiceA<-indiceReg(n, indice, facInd, FALSE)
        mod_stand<-0 # the root is not provided nor assumed to be one of the selected regimes, so we rowstandardize (could be optional)
    }else if(root==TRUE){
        indiceA<-indiceReg(n, indice, facInd, TRUE)
        mod_stand<-0
    }else if(root=="stationary"){
        indiceA<-indiceReg(n, indice, facInd, FALSE)
        mod_stand<-1
    }
    
    ## Return the epochs and regime lists
    
    # regime lists
    listReg<-sapply(1:n,function(x){sapply(1:p,function(db){regimeList(indiceA[[x]],k=k,root)},simplify=FALSE)},simplify=FALSE)
    # mapped epochs
    epochs<-sapply(1:n,function(x){lineage<-as.numeric(c(cumsum(valLineage[[x]])[length(valLineage[[x]])],(cumsum(valLineage[[x]])[length(valLineage[[x]])]-cumsum(valLineage[[x]])))); lineage[which(abs(lineage)<1e-15)]<-0; lineage },simplify=FALSE)
    
    results <- list(epochs=epochs, listReg=listReg, mod_stand=mod_stand)
    return(results)
}

# Generate box constraints for the likelihood search
ratevalue<-function(up,low,x){
    y<-x<low | x>up
    x[y]=0
    return(x)
}

# Return the expected scatter matrix given the "alpha" matrix and the empirical covariance.
Stationary_to_scatter <- function(alpha,stationary){
    return(alpha%*%stationary + stationary%*%t(alpha))
}

# Compute the stationary multivariate normal distribution for the multivariate Ornstein-Uhlenbeck (Bartoszek et al. 2012 - B.8) - to modify to add the time
StationaryVariance <- function(alpha,sigma){
    sigma <- sigma
    eig <- eigen(alpha)
    P <- eig$vectors
    invP <- solve(P)
    eigvalues <- eig$values
    p=dim(sigma)[1]
    Mat <- matrix(0,p,p)
    for(i in 1:p){
        for(j in 1:p){
            Mat[i,j] <- 1/(eigvalues[i]+eigvalues[j])
        }
    }
    StVar <- P%*%(Mat*(invP%*%sigma%*%t(invP)))%*%t(P)
    return(StVar)
}

# Constrainded Choleski decomposition (Adams, 2013: Systematic Biology)
build.chol<-function(b,p){
 c.mat<-matrix(0,nrow=p,ncol=p)
 c.mat[lower.tri(c.mat)] <- b[-1]
 c.mat[p,p]<-exp(b[1])
 c.mat[1,1]<-sqrt(sum((c.mat[p,])^2))
 if(p>2){
 for (i in 2:(p-1)){
 c.mat[i,i]<-ifelse((c.mat[1,1]^2-sum((c.mat[i,])^2))>0,sqrt(c.mat[1,1]^2-sum((c.mat[i,])^2)), 0)
 }}
 return(c.mat)
}


# Generate a weight matrix for OUM
.make.x <- function(phy, param, X, model, root, std=0){
    
    switch(model,
    "OUM"={
        if(!inherits(phy,"simmap")) stop("A tree of class \"simmap\" is required for the OUM model")
        n <- Ntip(phy)
        precalc <- mv.Precalc(phy, nb.traits=1, param=list(model="OUM", root=root))
        X <- .Call(mvmorph_weights, nterm=as.integer(n), epochs=precalc$epochs, lambda=param, S=1, S1=1, beta=precalc$listReg, root=as.integer(std))
        if(root==TRUE) colnames(X) <- c("root", colnames(phy$mapped.edge)) else colnames(X) <- colnames(phy$mapped.edge)
    },
    "OU1"={
        n <- Ntip(phy)
        precalc <- mv.Precalc(phy, nb.traits=1, param=list(model="OU1", root=root))
        X <- .Call(mvmorph_weights, nterm=as.integer(n), epochs=precalc$epochs, lambda=param, S=1, S1=1, beta=precalc$listReg, root=as.integer(0))
        if(root==TRUE) colnames(X) <- c("root","optimum") else colnames(X) <- c("optimum")
    },
    )
    
    return(X)
}

# precalculations for mvgls structures [ for future developments ]
.prepModel <- function(phy, model, root){
    
    switch(model,
    "OUM"={
        if(!inherits(phy,"simmap") && model=="OUM") stop("A tree of class \"simmap\" is required for the OUM model")
        precalc <- mv.Precalc(phy, nb.traits=1, param=list(model="OUM", root=root))
    },
    "OU1"={
        if(is.ultrametric(phy) & root==TRUE) warning("Estimation of the root and optimum in \"OU1\" is not identifiable on ultrametric trees")
        precalc <- mv.Precalc(phy, nb.traits=1, param=list(model="OU1", root=root))
    },
    precalc <- list(randomRoot=TRUE)
    )
    return(precalc)
}



# Functions to check input values in lists > e.g. in EIC
.check_samples <- function(list_values){
    check <- sapply(list_values, function(x){
        if(inherits(x, "try-error")) NA else x
    })
    
    na_rm <- check[!is.na(check)]
    if(length(na_rm)<length(check)) warning("There were multiple issues/aborted estimations in the bootstrapped samples")
    return(na_rm)
}


##----------------------mvfit_likelihood--------------------------------------##

loglik_mvmorph<-function(data,V=NULL,D=NULL,n,k,error=NULL,precalc=precalc,method, param=list(),ch=NULL,precalcMat=NULL,sizeD=NULL,NA_val=NULL,Indice_NA=NULL,theta_mle=TRUE,theta=NULL,istrend=FALSE,trend=0){

data<-as.numeric(as.matrix(data))
size<-k*n
if(NA_val==TRUE){
    V<-V[-Indice_NA,-Indice_NA]
    D<-D[-Indice_NA,]
    data<-data[-Indice_NA]
    error<-error[-Indice_NA]
    size<-length(data)
}

switch(method,

"rpf"={
    if(is.null(error)!=TRUE){
        ms<-1
    }else{ ms<-0} # length of the error vector should be the same as the data vector
   
    
    if(theta_mle==TRUE){
       cholres<-.Call(Chol_RPF,V,D,data,as.integer(sizeD),as.integer(size),mserr=error,ismserr=as.integer(ms))
       theta<-pseudoinverse(cholres[[3]])%*%cholres[[4]]
    }else{
        # We avoid the transformation of D and data by the Cholesky factor
       cholres<-.Call(Chol_RPF_only,V,as.integer(size),mserr=error,ismserr=as.integer(ms))
    }
        det<-cholres[[2]]
        
    if(istrend==FALSE){
        residus=D%*%theta-data
    }else{
        residus=(D%*%theta+trend)-data
    }
    quad<-.Call(Chol_RPF_quadprod, cholres[[1]], residus, as.integer(size))
    logl<--.5*quad-.5*as.numeric(det)-.5*(size*log(2*pi))
    results<-list(logl=logl,anc=theta)
},
"univarpf"={
    if(is.null(error)!=TRUE){
        ms<-1
    }else{ ms<-0}
    

   if(theta_mle==TRUE){
       cholres<-.Call(Chol_RPF_univ,V,D,data,as.integer(sizeD),as.integer(size),mserr=error,ismserr=as.integer(ms))
       theta<-pseudoinverse(cholres[[3]])%*%cholres[[4]]
   }else{
       cholres<-.Call(Chol_RPF_univ_only,V,as.integer(size),mserr=error,ismserr=as.integer(ms))
   }
 
    det<-cholres[[2]]
    
    if(istrend==FALSE){
        residus=D%*%theta-data
    }else{
        residus=(D%*%theta+trend)-data
    }
    
    quad<-.Call(Chol_RPF_quadprod_column, cholres[[1]], residus, as.integer(size))
    logl<--.5*quad-.5*as.numeric(det)-.5*(size*log(2*pi))
    results<-list(logl=logl,anc=theta)
},
"sparse"={
    ## On considere que les valeurs dans @entries de V de precalc ont ete modifiees / assume that the values in @entries were updated
    if(is.null(error)==FALSE){
        if(!is.null(precalc)){
        diag(precalc$V)<-diag(precalc$V)+error
        }else{
        diag(precalcMat)<-diag(precalcMat)+error
        }
    }
    if(!is.null(precalc)){
    U<-update.spam.chol.NgPeyton(precalc$ch,precalc$V)
    }else{
    U<-update.spam.chol.NgPeyton(ch,precalcMat)
    }
    if(theta_mle==TRUE){
        vec<-forwardsolve(U,data)
        xx<-forwardsolve(U,D)
        theta<-pseudoinverse(matrix(xx,ncol=sizeD))%*%vec
    }
    
    if(istrend==FALSE){
        res<-D%*%theta-data
    }else{
        res<-(D%*%theta+trend)-data
    }
    
    vec1<-forwardsolve(U,res)
    a<-sum(vec1^2)
    DET<-determinant(U)
    logl<--.5*a-.5*as.numeric(DET$modulus*2)-.5*(n*k*log(2*pi))
    results<-list(logl=logl,anc=theta)
},
"pseudoinverse"={
    if(is.null(error)==FALSE){
        diag(V)<-diag(V)+error
    }
    
    inv<-pseudoinverse(V)
    if(theta_mle==TRUE){
         theta<-pseudoinverse(t(D)%*%inv%*%D)%*%t(D)%*%inv%*%data
    }
    
    DET<-determinant(V, logarithm=TRUE)
    
    if(istrend==FALSE){
        res<-D%*%theta-data
    }else{
        res<-(D%*%theta+trend)-data
    }
   
    logl<--.5*(t(res)%*%inv%*%(res))-.5*as.numeric(DET$modulus)-.5*(size*log(2*pi))
    results<-list(logl=logl,anc=theta)
},
"inverse"={
    if(is.null(error)==FALSE){
        diag(V)<-diag(V)+error
    }

    inv<-solve(V)
    if(theta_mle==TRUE){
    theta<-solve(t(D)%*%inv%*%D)%*%t(D)%*%inv%*%data
    }
    DET<-determinant(V, logarithm=TRUE)
    if(istrend==FALSE){
        res<-D%*%theta-data
    }else{
        res<-(D%*%theta+trend)-data
    }
    logl<--.5*(t(res)%*%inv%*%(res))-.5*as.numeric(DET$modulus)-.5*(size*log(2*pi))
    results<-list(logl=logl,anc=theta)
})


return(results)


}


##----------------------Print_functions---------------------------------------##

print.mvmorph.bm<-function(x,...){
    
    cat("\n")
    if(x$param$constraint==TRUE){
        cat("-- Summary results for multiple constrained rate",x$param$model,"model --","\n")
    }else if(x$param$constraint=="default" & x$param$decomp=="user"){
        cat("-- Summary results for user-defined",x$param$model,"constrained model --","\n")
    }else if(x$param$constraint=="correlation"){
        cat("-- Summary results for common correlation ",x$param$model,"model --","\n")
    }else if(x$param$constraint=="shared"){
        cat("-- Summary results for shared eigenvectors ",x$param$model,"model --","\n")
    }else if(x$param$constraint=="proportional"){
        cat("-- Summary results for proportional rate matrices ",x$param$model,"model --","\n")
    }else if(x$param$constraint=="variance"){
        cat("-- Summary results for common variance ",x$param$model,"model --","\n")
    }else{
        cat("-- Summary results for multiple rate",x$param$model,"model --","\n")
    }
    cat("LogLikelihood:","\t",x$LogLik,"\n")
    cat("AIC:","\t",x$AIC,"\n")
    cat("AICc:","\t",x$AICc,"\n")
    cat(x$param$nparam,"parameters","\n")
    cat("\n")
    cat("Estimated rate matrix","\n")
    cat("______________________","\n")
    print(x$sigma)
    cat("\n")
    cat("Estimated root state","\n")
    cat("______________________","\n")
    print(x$theta)
    cat("\n")
    
    if(!is.null(x[["trend"]])){
        cat("Estimated trend values","\n")
        cat("______________________","\n")
        print(x$trend)
        cat("\n")
    }


}

print.mvmorph.acdc<-function(x,...){
    cat("\n")
    cat("-- Summary results for Early Burst or ACDC model --","\n")
    cat("LogLikelihood:","\t",x$LogLik,"\n")
    cat("AIC:","\t",x$AIC,"\n")
    cat("AICc:","\t",x$AICc,"\n")
    cat(x$param$nparam,"parameters","\n")
    cat("Rate change:","\n")
    cat("______________________","\n")
    print(x$beta)
    cat("\n")
    cat("Estimated rate matrix","\n")
    cat("______________________","\n")
    print(x$sigma)
    cat("\n")
    cat("Estimated root states","\n")
    cat("______________________","\n")
    print(x$theta)
    cat("\n")
}

print.mvmorph.ou<-function(x,...){
    cat("\n")
    cat("-- Summary results --","\n")
    cat("LogLikelihood:","\t",x$LogLik,"\n")
    cat("AIC:","\t",x$AIC,"\n")
    cat("AICc:","\t",x$AICc,"\n")
    cat(x$param$nparam,"parameters","\n")
    cat("\n")
    cat("Estimated theta values","\n")
    cat("______________________","\n")
    print(x$theta)
    cat("\n")
    cat("ML alpha values","\n")
    cat("______________________","\n")
    print(x$alpha)
    cat("\n")
    cat("ML sigma values","\n")
    cat("______________________","\n")
    print(x$sigma)
}

print.mvmorph.shift<-function(x,...){
    cat("\n")
    cat("-- Summary results for the",x$param$model[2]," --","\n")
    cat("LogLikelihood:","\t",x$LogLik,"\n")
    cat("AIC:","\t",x$AIC,"\n")
    cat("AICc:","\t",x$AICc,"\n")
    cat(x$param$nparam,"parameters")
    cat("\n")
    cat("Estimated theta values","\n")
    cat("______________________","\n")
    print(x$theta)
    cat("\n")
    if(x$param$model[1]=="CV" || x$param$model[1]=="CVG"){
        cat("ML beta values","\n")
        cat("______________________","\n")
        print(x$beta)
    }else{
        cat("ML alpha values","\n")
        cat("______________________","\n")
        print(x$alpha)
    }
    cat("\n")
    cat("ML sigma values","\n")
    cat("______________________","\n")
    print(x$sigma)
    if(x$param$model[1]=="RR" || x$param$model[1]=="radiate"){
        cat("\n")
        cat("ML sigma radiation values","\n")
        cat("______________________","\n")
        print(x$sig)
    }
    if(x$param$model[1]=="CVG" || x$param$model[1]=="OVG"){
        cat("\n")
        cat("ML sigma values ( recent slice:",x$param$names_regimes[2],")","\n")
        cat("______________________","\n")
        print(x$sig)
    }
    if(x$param$model[1]=="OVG" || x$param$model[1]=="OV"){
        cat("\n")
        cat("ML beta values","\n")
        cat("______________________","\n")
        print(x$beta)
    }
}

print.mvmorph.lrt<-function(x,...){
    if(x$pval<0.000001){signif<-c("***")}else if(x$pval<0.001){
        signif<-c("**") }else if(x$pval<0.01){signif<-c("*")}else if(x$pval<0.05){signif<-c(".")}else{signif<-""}
    cat("-- Log-likelihood Ratio Test --","\n")
    cat("Model",x$model1," versus ",x$model2,"\n")
    cat("Number of degrees of freedom :",x$ddf,"\n")
    cat("LRT statistic:",x$ratio," p-value:",x$pval,signif,"\n")
    cat("---","\n")
    cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1")
}

print.mvmorph.precalc<-function(x,...){
    cat("A tree with",length(x$tree$tip.label),"species used in precalc","\n")
    cat("Optimized for:","\n")
    cat("-----------------","\n")
    cat("method:",x$param$method,"\n")
    cat("model:",x$model,"\n")
    cat("number of traits:",x$param$nbtraits,"\n")
}

print.mvmorph.llik<-function(x,...){
    cat("Log-likelihood function :","\n")
    cat("llik(par, root.mle=TRUE)","\n")
    cat("\"par\" argument vector order:","\n")
    model<-attr(x,"model")
    switch(model,
    "OU"={
        cat("1) alpha (",attr(x,"alpha")," parameters)","\n")
        cat("2) sigma (",attr(x,"sigma")," parameters)","\n")
        cat("3) theta (",attr(x,"theta")," parameters if \"root.mle=FALSE\")","\n")
    },
    "RWTS"={
        if(is.null(attr(x,"trend"))){
            cat("1) sigma (",attr(x,"sigma")," parameters)","\n")
            cat("2) theta (",attr(x,"theta")," parameters if \"root.mle=FALSE\")","\n")
        }else{
            cat("1) sigma (",attr(x,"sigma")," parameters)","\n")
            cat("2) trend (",attr(x,"trend")," parameters)","\n")
            cat("3) theta (",attr(x,"theta")," parameters) Note: \"root.mle=TRUE\" is not allowed with the \"trend\" model)","\n")
        }
    },
    "BM"={
        if(is.null(attr(x,"trend"))){
        cat("1) sigma (",attr(x,"sigma")," parameters)","\n")
        cat("2) theta (",attr(x,"theta")," parameters if \"root.mle=FALSE\")","\n")
        }else{
        cat("1) sigma (",attr(x,"sigma")," parameters)","\n")
        cat("2) trend (",attr(x,"trend")," parameters)","\n")
        cat("3) theta (",attr(x,"theta")," parameters) Note: \"root.mle=TRUE\" is not allowed with the \"trend\" model)","\n")
        }
    },
    "EB"={
        cat("1) beta (",attr(x,"beta")," parameters)","\n")
        cat("2) sigma (",attr(x,"sigma")," parameters)","\n")
        cat("3) theta (",attr(x,"theta")," parameters if \"root.mle=FALSE\")","\n")
    },
    "ER"={
        cat("1) alpha (",attr(x,"alpha")," parameters)","\n")
        cat("2) sigma (",attr(x,"sigma")," parameters)","\n")
        cat("3) theta (",attr(x,"theta")," parameters if \"root.mle=FALSE\")","\n")
    },
    "RR"={
        cat("1) alpha (",attr(x,"alpha")," parameters)","\n")
        cat("2) sigma (",attr(x,"sigma")," parameters)","\n")
        cat("3) sig (",attr(x,"sig")," parameters)","\n")
        cat("4) theta (",attr(x,"theta")," parameters if \"root.mle=FALSE\")","\n")
    },
    "CV"={
        cat("1) beta (",attr(x,"beta")," parameters)","\n")
        cat("2) sigma (",attr(x,"sigma")," parameters)","\n")
        cat("3) theta (",attr(x,"theta")," parameters if \"root.mle=FALSE\")","\n")
    },
    "CVG"={
        cat("1) beta (",attr(x,"beta")," parameters)","\n")
        cat("2) sigma (",attr(x,"sigma")," parameters)","\n")
        cat("3) sig (",attr(x,"sig")," parameters)","\n")
        cat("4) theta (",attr(x,"theta")," parameters if \"root.mle=FALSE\")","\n")
    },
    "OV"={
        cat("1) alpha (",attr(x,"alpha")," parameters)","\n")
        cat("2) sigma (",attr(x,"sigma")," parameters)","\n")
        cat("3) beta (",attr(x,"beta")," parameters)","\n")
        cat("4) theta (",attr(x,"theta")," parameters if \"root.mle=FALSE\")","\n")
    },
    "OVG"={
        cat("1) alpha (",attr(x,"alpha")," parameters)","\n")
        cat("2) sigma (",attr(x,"sigma")," parameters)","\n")
        cat("3) sig (",attr(x,"sig")," parameters)","\n")
        cat("4) beta (",attr(x,"beta")," parameters)","\n")
        cat("5) theta (",attr(x,"theta")," parameters if \"root.mle=FALSE\")","\n")
    })
    cat("-----------------","\n")
    fun<-x; attributes(fun)=NULL
    print(fun)
}

summary.mvmorph<-function(object,...){
    cat("mvMORPH model :",object$param$model," summary","\n")
    cat("AIC :",object$AIC,"\n")
    cat("AICc:",object$AICc,"\n")
    cat("Log-Likelihood:",object$LogLik,"\n")
    if(object$convergence==0){cat("Succesful convergence","\n")}else{cat("Convergence has not been reached","\n")}
    if(object$hess.value==0){cat("Reliable solution","\n")}else{cat("Unreliable solution (Likelihood at a saddle point)","\n")}
}

print.mvmorph.aicw<-function(x,...){
    cat("-- Akaike weights --","\n")
    aics <- data.frame(Rank=1:length(x$models), AIC=x$AIC, diff=x$diff, wi=x$wi, AICw=x$aicweights)
    row.names(aics) <- as.character(x$models)
    aics <- aics[order(aics$wi, decreasing=TRUE),]
    aics$Rank <- 1:length(x$models)
    aics[,c(4,5)][aics[,c(4,5)]<1e-8]<-0
    print(aics, zero.print = ".", digits=3)
}

print.mvmorph.estim<-function(x,...){
    number_of_missing <- length(x$NA_index)
    if(is.null(number_of_missing)){
    cat("Ancestral states reconstructed for",nrow(x$estimates)," nodes and",ncol(x$estimates)," traits.","\n")
    }else{
    cat("Dataset with",number_of_missing," imputed missing values.","\n")
    }
}

## Return the model AIC
AIC.mvmorph<-function(object,...,k){
    return(object$AIC)
}

## Return the model AICc
AICc <- function(x) UseMethod("AICc")
AICc.mvmorph<-function(object){
    return(object$AICc)
}

## Return the model fit log-likelihood
logLik.mvmorph<-function(object,...){
    return(object$LogLik)
}

## Change to include the tree in the analysis? == problematic for large trees... need too much storage!
simulate.mvmorph<-function(object,nsim=1,seed=NULL,...){
    mvSIM(...,param=object,nsim=nsim)
}

## Return the stationary variance for the multivariate Ornstein-Uhlenbeck

stationary.mvmorph<-function(object){
    if(any(class(object)=="mvmorph.shift") | any(class(object)=="mvmorph.ou") ){
            if(is.null(object[["alpha"]])==TRUE){
                stop("The stationary distribution can be computed only for models including Ornstein-Uhlenbeck processes.")
            }
    statMat<-StationaryVariance(object$alpha,object$sigma)
    rownames(statMat)<-rownames(object$sigma)
    colnames(statMat)<-colnames(object$sigma)
    return(statMat)   
    }else{
      warning("The stationary distribution can be computed only for Ornstein-Uhlenbeck processes.","\n")
    }
}

## Compute the phylogenetic half-life

halflife.mvmorph<-function(object){
    if(inherit(object,"mvmorph.ou") | inherit(object,"mvmorph.shift")){
        if(is.null(object[["alpha"]])==TRUE){
            stop("The phylogenetic half-life can be computed only for models including Ornstein-Uhlenbeck processes.")
        }
    lambda<-eigen(object$alpha)$values
    phyhalflife<-log(2)/Re(lambda)
        return(phyhalflife)
    }else{
        warning("The phylogenetic half-life is computed only for Ornstein-Uhlenbeck models.","\n")  
    }
}
