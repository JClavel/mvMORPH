################################################################################
##                                                                            ##
##                               mvMORPH: mvEB                                ##
##                                                                            ##
##  Created by Julien Clavel - 16-07-2013                                     ##
##  (julien.clavel@hotmail.fr/ julien.clavel@univ-lyon1.fr)                   ##
##   require: phytools, ape, corpcor, geiger                                  ##
##                                                                            ##
################################################################################

mvEB<-function(tree, data, error=NULL, param=list(up=0), method=c("rpf","sparse","inverse","pseudoinverse","pic"), scale.height=FALSE, optimization=c("Nelder-Mead","L-BFGS-B","subplex"), control=list(maxit=20000), precalc=NULL, diagnostic=TRUE, echo=TRUE){

if(missing(tree)) stop("The tree object is missing!")
if(missing(data)) stop("You must provide a dataset along with your tree!")

#set data as a matrix if a vector is provided instead
if(!is.matrix(data)){data<-as.matrix(data)}
# number of species (tip)
n<-dim(data)[1]
# number of variables
k<-dim(data)[2]
# method for the optimizer & algorithm
optimization<-optimization[1]
# select default model
method<-method[1]

# Check the order of the dataset and the phylogeny
if(!is.null(rownames(data))) {
    if(any(tree$tip.label==rownames(data))){
        data<-data[tree$tip.label,]
    }else if(echo==TRUE){
        cat("row names of the data matrix must match tip names of your phylogeny!","\n")
    }
}else if(echo==TRUE){
    cat("species in the matrix are assumed to be in the same order as in the phylogeny, otherwise specify rownames of 'data'","\n")
}


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
# bind error to a vector
if(!is.null(error)){
    error<-as.vector(error)
    error[is.na(error)] <- 0
}


# scale height of the tree
maxHeight<-1
if(scale.height==TRUE){
maxHeight<-max(nodeHeights(tree))
tree$edge.length<-tree$edge.length/maxHeight
}
# number of traits
k<-ncol(data)
if(is.null(k)){
k<-1
}

# number of parameters
npar<-(k*(k+1)/2)

if(is.null(param[["low"]])==TRUE){
    if(scale.height==FALSE) maxHeight<-max(nodeHeights(tree))
    low<-param$low<- log(10^-5)/maxHeight # Slater & Pennell 2013
   if(echo==TRUE) cat("No lower bound provided. Use of default setting \" ",low,"\"")

}else{low<-param$low}
if(is.null(param[["up"]])==TRUE){
  if(echo==TRUE)  cat("No upper bound provided. Use of default setting \"0\"")
    up<-param$up<-0
}else{up<-param$up}

# Scatter (sigma) matrix decomposition
if(is.null(param[["decomp"]])){
    decomp<-param$decomp<-"cholesky"
}else{
    decomp<-param$decomp[1]
}

# function for sigma parameterization
sigmafun <- function(par) {symPar(par, decomp=decomp, p=k)}

if(is.null(precalc)==FALSE & class(precalc)=="mvmorph.precalc"){
    
    tree<-precalc$tree
    C<-precalc$C1
    D<-precalc$D
    
    if(method=="sparse"){
        # Yale sparse format
        JAr<-precalc$JAr
        IAr<-precalc$IAr
        ch<-precalc$ch
        precalcMat<-precalc$V
    }
    
    
}else{
if(method!="pic"){
# compute the vcv
C<-vcv.phylo(tree) 
# Design matrix
D<-multD(tree,k,n,smean=TRUE)
    }else{
        C<-NULL
        D<-NULL
    }

if(method=="sparse"){
    V<-kronecker((matrix(1,k,k)+diag(k)), C)
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
##--------------Likelihood_functions--------------------------------------------------##


# switch method
switch(method,
"rpf"={
    eb_fun_matrix<-function(beta,sig,n,k,precalcMat,C,D,dat,error,method,precalc,theta_mle,theta){
        V<-.Call("kroneckerEB",R=sig,C=C, beta=matrix(beta,k,k), Rrows=as.integer(k),  Crows=as.integer(n), PACKAGE="mvMORPH")
        loglik<-loglik_mvmorph(dat,V,D,n,k,error=error,precalc=precalc,method=method,ch=ch,precalcMat=precalcMat,sizeD=ncol(D),NA_val=NA_val,Indice_NA=Indice_NA,theta_mle=theta_mle,theta=theta) ######### A modif pour sparse
        return(loglik)
    }
},
"sparse"={
    eb_fun_matrix<-function(beta,sig,n,k,precalcMat,C,D,dat,error,method,precalc,theta_mle,theta){
        V<-.Call("kroneckerSparEB",R=sig,C=C, beta=matrix(beta,k,k), Rrows=as.integer(k),  Crows=as.integer(n),  IA=as.integer(precalcMat@rowpointers-1), JA=as.integer(precalcMat@colindices-1), A=precalcMat@entries, PACKAGE="mvMORPH")
        loglik<-loglik_mvmorph(dat,V,D,n,k,error=error,precalc=precalc,method=method,ch=ch,precalcMat=precalcMat,sizeD=ncol(D),NA_val=NA_val,Indice_NA=NULL,theta_mle=theta_mle,theta=theta) ######### A modif pour sparse
        return(loglik)
    }
},
"pseudoinverse"={
    eb_fun_matrix<-function(beta,sig,n,k,precalcMat,C,D,dat,error,method,precalc,theta_mle,theta){
        V<-.Call("kroneckerEB",R=sig,C=C, beta=matrix(beta,k,k), Rrows=as.integer(k),  Crows=as.integer(n), PACKAGE="mvMORPH")
        loglik<-loglik_mvmorph(dat,V,D,n,k,error=error,precalc=precalc,method=method,ch=ch,precalcMat=precalcMat,sizeD=ncol(D),NA_val=NA_val,Indice_NA=Indice_NA,theta_mle=theta_mle,theta=theta) ######### A modif pour sparse
        return(loglik)
    }
},
"inverse"={
    eb_fun_matrix<-function(beta,sig,n,k,precalcMat,C,D,dat,error,method,precalc,theta_mle,theta){
    V<-.Call("kroneckerEB",R=sig,C=C, beta=matrix(beta,k,k), Rrows=as.integer(k),  Crows=as.integer(n), PACKAGE="mvMORPH")
    loglik<-loglik_mvmorph(dat,V,D,n,k,error=error,precalc=precalc,method=method,ch=ch,precalcMat=precalcMat,sizeD=ncol(D),NA_val=NA_val,Indice_NA=Indice_NA,theta_mle=theta_mle,theta=theta) ######### A modif pour sparse
    return(loglik)
    }
},
"pic"={
    eb_fun_matrix<-function(beta,sig,n,k,precalcMat,C,D,dat,error,method,precalc,theta_mle,theta){

        res<-.Call("PIC_gen", x=dat, n=as.integer(k), Nnode=as.integer(tree$Nnode), nsp=as.integer(n), edge1=as.integer(tree$edge[,1]), edge2=as.integer(tree$edge[,2]), edgelength=list(tree$edge.length), times=tt, rate=as.double(rep(beta,k)), Tmax=1, Model=as.integer(1), mu=theta, sigma=sig, PACKAGE="mvMORPH")
        logl<- -0.5 * ( n * k * log( 2 * pi) +  res[[5]] + n * res[[6]]  + res[[4]] )
        # return the beta value provided instead of the MLE...
        if(theta_mle==TRUE){ theta<-res[[7]] }
       return(list(logl=logl,anc=theta, sigma=res[[2]]))
    }
}
)



##--------------Maximum Likelihood estimation of EB---------------------------##  

likEB <- function(dat, error, C, D, beta, sig, k, n, method, precalc, precalcMat, theta_mle=TRUE, theta=NULL) {

    loglik<-eb_fun_matrix(beta,sig,n,k,precalcMat,C,D,dat,error,method,precalc,theta_mle,theta) ######### A modif pour sparse

    list(loglik=-loglik$logl, ancstate=loglik$anc)
}

##-------------------Function log-likelihood----------------------------------##

llfun <- function(par,...){
    
    args <- list(...)
    if(is.null(args[["root.mle"]])) args$root.mle <- TRUE
    if(is.null(args[["theta"]])) args$theta <- FALSE
    
    if(args$root.mle==TRUE){
        
        result <- -as.numeric(likEB(dat=as.matrix(data), error=error, C=C, D=D, beta=ratevalue(up,low,par[1]), sig=sigmafun(par[1+seq_len(npar)]),k=k, n=n, method=method, precalc=precalc, precalcMat=precalcMat, theta_mle=TRUE, theta=NULL)$loglik)
        if(args$theta==TRUE){
            theta <- as.numeric(likEB(dat=as.matrix(data), error=error, C=C, D=D, beta=ratevalue(up,low,par[1]), sig=sigmafun(par[1+seq_len(npar)]),k=k, n=n, method=method, precalc=precalc, precalcMat=precalcMat, theta_mle=TRUE, theta=NULL)$ancstate)
            result<-list(logl=result, theta=theta)
        }
    }else{
        
        result <- -as.numeric(likEB(dat=as.matrix(data), error=error, C=C, D=D, beta=ratevalue(up,low,par[1]), sig=sigmafun(par[1+seq_len(npar)]),k=k, n=n, method=method, precalc=precalc, precalcMat=precalcMat, theta_mle=FALSE, theta=par[1+npar+seq_len(k)])$loglik)
        
    }
    return(result)
}
# attributes to the loglik function
attr(llfun, "model")<-"EB"
attr(llfun, "beta")<-1
attr(llfun, "sigma")<-npar
attr(llfun, "theta")<-k

##---------------Optimizing function------------------------------------------##
# initial value for exponent parameter
if(is.null(param[["beta"]])==TRUE){
    ebval<- -1/maxHeight
}else{
    ebval<-param$beta
}

# initial values for the optimizer
if(is.null(param[["sigma"]])==TRUE){
sig1<-varBM(tree,data,n,k)
sig1<-sym.unpar(sig1)
}else{
sig1<-param$sigma
}

## Check first if we want to return the log-likelihood function only
if(optimization=="fixed"){
    message("No optimization performed, only the Log-likelihood function is returned with default parameters.")
    param$sigmafun<-sigmafun
    param$nbspecies<-n
    param$ntraits<-k
    param$nregimes<-1
    param$method<-method
    param$optimization<-optimization
    param$traits<-colnames(data)
    param$model<-if(ebval>0){"AC"}else if(ebval<0){"EB"}else{"ACDC"}
    theta.mat<-rep(0,k)
    class(llfun) = c("mvmorph.llik")

    results<-list(llik=llfun, theta=theta.mat, sigma=sym.par(sig1), beta=ebval, param=param)
    class(results)<-c("mvmorph")
    invisible(results)
    return(results)
}


if(method=="pic"){
    warning<-eval_polytom(tree)
    tree<-reorder.phylo(tree,"postorder")
    # times from the root
    tt<-.Call("times_root", brlength=tree$edge.length, edge1=as.integer(tree$edge[,1]), edge2=as.integer(tree$edge[,2]), ntip=as.integer(n), Nnode=as.integer(tree$Nnode), PACKAGE="mvMORPH")
    
    if(optimization=="subplex"){
        estim <- subplex(par=c(ebval,sig1), fn = function (par) {likEB(dat=as.matrix(data), error=error, C=C, D=D, beta=ratevalue(up,low,par[1]), sig=sigmafun(par[1+seq_len(npar)]),k=k, n=n, method=method, precalc=precalc, precalcMat=precalcMat)$loglik }, hessian=TRUE, control=control)
    }else{
        estim <- optim(par=c(ebval,sig1), fn = function (par) {likEB(dat=as.matrix(data), error=error, C=C, D=D, beta=ratevalue(up,low,par[1]), sig=sigmafun(par[1+seq_len(npar)]),k=k, n=n, method=method, precalc=precalc, precalcMat=precalcMat)$loglik }, gr=NULL, hessian=TRUE, method = optimization, control=control)
    }
}else{


if(optimization=="subplex"){
    estim <- subplex(par=c(ebval,sig1), fn = function (par) {likEB(dat=as.matrix(data), error=error, C=C, D=D, beta=ratevalue(up,low,par[1]), sig=sigmafun(par[1+seq_len(npar)]),k=k, n=n, method=method, precalc=precalc, precalcMat=precalcMat)$loglik }, hessian=TRUE, control=control)
   }else{
    estim <- optim(par=c(ebval,sig1), fn = function (par) {likEB(dat=as.matrix(data), error=error, C=C, D=D, beta=ratevalue(up,low,par[1]), sig=sigmafun(par[1+seq_len(npar)]),k=k, n=n, method=method, precalc=precalc, precalcMat=precalcMat)$loglik }, gr=NULL, hessian=TRUE, method = optimization, control=control)
   }
}
##-----------------Summarizing results----------------------------------------##
# names for the plot
if(is.null(colnames(data))){
    names_data_matrix<-rep("",k)
    names_data<-list(names_data_matrix,names_data_matrix)
}else{
    names_data_matrix<-colnames(data)
    names_data<-list(names_data_matrix,names_data_matrix)
}

# sigma matrix
resultList<-sigmafun(estim$par[1+seq_len(npar)])
dimnames(resultList)<-names_data

#ancestral states estimates
anc<-matrix(likEB(dat=as.matrix(data), error=error, C=C, D=D, beta=ratevalue(up,low,estim$par[1]), sig=sigmafun(estim$par[1+seq_len(npar)]),k=k, n=n, method=method, precalc=precalc, precalcMat=precalcMat)$ancstate,nrow=1)
rownames(anc)<-"theta"
colnames(anc)<-names_data_matrix

# rate parameter
r=ratevalue(up,low,estim$par[1])
# LogLikelihood
LL<--estim$value
# models parameters
nparam=k+npar+1 # length(estim$par) # AIC
AIC<--2*LL+2*nparam
# AIC corrected
nobs <- length(which(!is.na(data)))
AICc<-AIC+((2*nparam*(nparam+1))/(nobs-nparam-1)) #Hurvich et Tsai, 1989
##---------------------Diagnostics--------------------------------------------##

if(estim$convergence==0 & diagnostic==TRUE){  
cat("\n","successful convergence of the optimizer","\n") 
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
cat("-- Summary results for Early Burst or ACDC model --","\n")
cat("LogLikelihood:","\t",LL,"\n")
cat("AIC:","\t",AIC,"\n")
cat("AICc:","\t",AICc,"\n")
cat(nparam,"parameters","\n")
cat("Rate change:","\n")
cat("______________________","\n")
print(r)
cat("\n")
cat("Estimated rate matrix","\n")
cat("______________________","\n")
print(resultList)
cat("\n")
cat("Estimated root states","\n")
cat("______________________","\n")
print(anc)
cat("\n")
}
##-------------------Save infos in parameters---------------------------------##
param$nparam<-nparam
param$nbspecies<-n
param$ntraits<-k
param$nregimes<-1
param$method<-method
param$optimization<-optimization
param$traits<-colnames(data)
param$model<-if(r>0){"AC"}else if(r<0){"EB"}else{"ACDC"}
param$sigmafun<-sigmafun
param$opt<-estim
class(llfun) = c("mvmorph.llik")
##-------------------Store results--------------------------------------------##
 
results<-list(LogLik=LL, AIC=AIC, AICc=AICc, theta=anc, beta=r, sigma=resultList,  convergence=estim$convergence, hess.values=hess.value, param=param, llik=llfun)

class(results)<-c("mvmorph","mvmorph.acdc")
invisible(results)

}
