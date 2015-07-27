################################################################################
##                                                                            ##
##                               mvMORPH: mvLL                                ##
##                                                                            ##
##  Created by Julien Clavel - 22-11-2014                                     ##
##  (julien.clavel@hotmail.fr/ clavel@biologie.ens.fr)                        ##
##   require: phytools, ape, corpcor, spam                                    ##
##                                                                            ##
################################################################################



mvLL<-function(tree,data,error=NULL,method=c("pic","rpf","sparse","inverse","pseudoinverse"), param=list(estim=TRUE, mu=0,sigma=0 ,D=NULL, check=TRUE), precalc=NULL){

    # select default method depending on data type for the tree object
 method<-method[1]
 dim_data<-dim(as.matrix(data))
 ntrait<-dim_data[2]
 if(class(tree)=="phylo" | class(tree[[1]])=="phylo"){
            datatype<-"tree"
            
        if(method!="pic"){
            
            if(ntrait==1){
            V<-vcv.phylo(tree)
            }else{
          stop("Provide a multivariate vcv-matrix for multivariate dataset")
            }
        }

        }else{
            datatype<-"vcv"
            V<-tree
            D<-param$D
            nspvcv<-nrow(V)
            
            if(nrow(as.matrix(data))!=nspvcv/ntrait){
          stop("Provide a multivariate vcv-matrix for multivariate dataset")
            }
          
        }
        
     
    if(datatype=="vcv" && method=="pic"){
          stop("object \"tree\" is not of class \"phylo\" try using methods \"sparse\", \"rpf\", \"inverse\" or \"pseudoinverse\", see details ?mvLL")
    }
    
    if(method!="pic"){
        # precalc options for faster computations
        if(is.null(param[["D"]]) & is.null(precalc)){
        
            if(datatype=="tree"){ntip<-length(tree$tip.label)}else{ntip=nspvcv/ntrait}
            D<-multD(tree,ntrait,ntip,smean=TRUE)
        
        }else if(class(precalc)=="precalc.mvmorph"){
        
            D<-precalc$D
        }
    }
    

    # Data format
    if(is.vector(data)==FALSE){data<-as.vector(as.matrix(data))}

    # switch methods depending on the nature of the tree object

switch(method,
    "pic"={
        
        if(is.null(param[["estim"]])){ param$estim<-FALSE }
        if(is.null(param[["check"]])){ param$check<-TRUE }
        
        # Preparing the phylo
        if(is.null(precalc)==TRUE){
            if(class(tree[[1]])=="phylo"){
                n=length(tree[[1]]$tip.label)
                k=dim(matrix(data,nrow=n))[2]
                if(length(tree)!=k){ stop("The number of trees in the list object for the \"tree\" argument must be the same as the number of traits ")}
                if(param$check==TRUE){
                    ind<-reorder(tree[[1]],"postorder", index.only=TRUE) # implique que c'est le meme arbre
                    value<-lapply(1:k,function(x){tree[[x]]$edge.length[ind]})
                    tree$edge<-tree[[1]]$edge[ind,]
                }else{
                    value<-lapply(1:k,function(x){tree[[x]]$edge.length})
                    tree$edge<-tree[[1]]$edge
                }
                    tree$Nnode<-tree[[1]]$Nnode
                # Change the way are computed the contrasts
                # Generalized Brownian Motion?
                if(param$estim==TRUE){
                        mod<-3
                        param$sigma<-1
                        param$mu<-1
                }else{
                    if(is.null(param[["mu"]])==TRUE & is.null(param[["sigma"]])==TRUE){
                        warning("No values specified for sigma and mu in the \"param\" list, analytical MLE is computed")
                        mod<-3
                        param$sigma<-1
                        param$mu<-1
                        
                    }else if(is.null(param[["mu"]])==TRUE & is.null(param[["sigma"]])==FALSE){
                        # estimate theta / user sigma
                        mod<-5
                        param$mu<-1
                    }else if(is.null(param[["mu"]])==FALSE & is.null(param[["sigma"]])==TRUE){
                        # estimate sigma / user theta
                        mod<-9
                        param$sigma<-1
                    }else if(is.null(param[["mu"]])==FALSE & is.null(param[["sigma"]])==FALSE){
                        # user sigma and theta
                        mod<-4
                    }
                }

            }else{ # one tree
                # Generalized Brownian Motion?
                if(param$estim==TRUE){
                    mod<-10 # meme topologie (sinon 3)
                    param$sigma<-1
                    param$mu<-1
                
                }else{
                    if(is.null(param[["mu"]])==TRUE & is.null(param[["sigma"]])==TRUE){
                        warning("No values specified for sigma and mu in the \"param\" list, analytical MLE is computed")
                        mod<-10
                        param$sigma<-1
                        param$mu<-1
                    }else if(is.null(param[["mu"]])==TRUE & is.null(param[["sigma"]])==FALSE){
                        # estimate theta/ user sigma
                        mod<-7
                        param$mu<-1
                    }else if(is.null(param[["mu"]])==FALSE & is.null(param[["sigma"]])==TRUE){
                        # user theta / estimate sigma
                        mod<-8
                        param$sigma<-1
                    }else if(is.null(param[["mu"]])==FALSE & is.null(param[["sigma"]])==FALSE){
                        # user sigma and theta
                        mod<-6
                    }

                }

                n=length(tree$tip.label)
                k=dim(matrix(data,nrow=n))[2]
                if(param$check==TRUE){
                 tree<-reorder(tree,"postorder")
                }
                value<-list(tree$edge.length)
            }
        
        }else{
                tree<-precalc$tree
                # Generalized Brownian Motion?
                if(param$estim==TRUE){
                    mod<-10 # meme topologie (sinon 3)
                }else{
                    if(is.null(param[["mu"]])==TRUE & is.null(param[["sigma"]])==TRUE){
                        warning("No values specified for sigma and mu in the \"param\" list, analytical MLE is computed")
                        mod<-10
                        param$sigma<-1
                        param$mu<-1
                    }else if(is.null(param[["mu"]])==TRUE & is.null(param[["sigma"]])==FALSE){
                        # estimate theta/ user sigma
                        mod<-7
                        param$mu<-1
                    }else if(is.null(param[["mu"]])==FALSE & is.null(param[["sigma"]])==TRUE){
                        # user theta / estimate sigma
                        mod<-8
                        param$sigma<-1
                    }else if(is.null(param[["mu"]])==FALSE & is.null(param[["sigma"]])==FALSE){
                        # user sigma and theta
                        mod<-6
                    }

                }#
                
                n=length(tree$tip.label)
                k=dim(matrix(data,nrow=n))[2]
                value<-list(tree$edge.length)
            }
        
        rate<-rep(0,k)
        # Compute the LLik
        res<-.Call("PIC_gen", x=data, n=as.integer(k), Nnode=as.integer(tree$Nnode), nsp=as.integer(n), edge1=as.integer(tree$edge[,1]), edge2=as.integer(tree$edge[,2]), edgelength=value, times=1, rate=rate, Tmax=1, Model=as.integer(mod), mu=param$mu, sigma=param$sigma)
        logl<- -0.5 * ( n * k * log( 2 * pi) +  res[[5]] + n * res[[6]]  + res[[4]] )
        
        if(param$estim==TRUE){
        results<-list(logl=logl,theta=res[[7]], sigma=res[[2]])
        }else{
        results<-list(logl=logl,theta=res[[7]])
        }
        },
    "rpf"={
        if(is.null(error)!=TRUE){ ms<-1 }else{ ms<-0}
        k<-ncol(D)
        if(is.null(k)){k=1}
        n<-dim_data[1]
        ntot<-n*ntrait
        cholres<-.Call("Chol_RPF",V,D,data,as.integer(k),as.integer(ntot),mserr=error,ismserr=as.integer(ms))
        beta<-pseudoinverse(cholres[[3]])%*%cholres[[4]]
        det<-cholres[[2]]
        residus=D%*%beta-data
        quad<-.Call("Chol_RPF_quadprod", cholres[[1]], residus, as.integer(n))
        logl<--.5*quad-.5*as.numeric(det)-.5*(n*k*log(2*pi))
        results<-list(logl=logl,theta=beta)
    },
    "sparse"={

        if(is.null(precalc)==TRUE){
            spambig=as.spam(V)
            if(is.null(error)==FALSE){
                diag(spambig)<-diag(spambig)+error
            }
            U<-chol(spambig)
        }else{
            U<-update(precalc$ch,precalc$V)
        }
        n<-dim_data[1]
        if(is.null(n)){n=length(D)}
        k=ncol(D)
        if(is.null(k)){k=1}
        vec<-forwardsolve(U,data)
        xx<-forwardsolve(U,D)
        beta<-pseudoinverse(matrix(xx,ncol=k))%*%vec
        res<-D%*%beta-data
        vec1<-forwardsolve(U,res)
        a<-sum(vec1^2)
        DET<-determinant(U)
        logl<--.5*(a)-.5*as.numeric(DET$modulus*2)-.5*(n*k*log(2*pi))
        results<-list(logl=logl,theta=beta)
    },
    "pseudoinverse"={
        if(is.null(error)==FALSE){
            diag(V)<-diag(V)+error
        }

        n<-dim_data[1]
        k<-ncol(D)
        if(is.null(k)){k=1}
        inv<-pseudoinverse(V)
        beta<-pseudoinverse(t(D)%*%inv%*%D)%*%t(D)%*%inv%*%data
        DET<-determinant(V, logarithm=TRUE)
        res<-D%*%beta-data
        logl<--.5*(t(res)%*%inv%*%(res))-.5*as.numeric(DET$modulus)-.5*(n*k*log(2*pi))
        results<-list(logl=logl,theta=beta)
    },
    "inverse"={
        if(is.null(error)==FALSE){
            diag(V)<-diag(V)+error
        }

        n<-dim_data[1]
        k<-ncol(D)
        if(is.null(k)){k=1}
        inv<-solve(V)
        beta<-solve(t(D)%*%inv%*%D)%*%t(D)%*%inv%*%data
        DET<-determinant(V, logarithm=TRUE)
        res<-D%*%beta-data
        logl<--.5*(t(res)%*%inv%*%(res))-.5*as.numeric(DET$modulus)-.5*(n*k*log(2*pi))
        results<-list(logl=logl,theta=beta)
    })
        class(results)<-c("mvmorph","loglik")
    
       return(results)
}
