################################################################################
##                  P3CA: Probabilistic and Phylogenetic PCA                  ##
##                                                                            ##
## Function to fit trait evolution models (Pagel's lambda/BM) in multivariate ##
## data, including high-dimensional data, with or without missing values      ##
## A projection in a lower space is returned (phylogenetic PCA)               ##
##                                                                            ##
## tree: phylogenetic tree                                                    ##
## q: number of variables in the lower space                                  ##
##                                                                            ##
################################################################################

p3ca <- function(Y, tree, q=4, model='lambda', REML=TRUE, tol=1e-4, EM=FALSE, plot=TRUE, ...){
  
  if(model!='lambda' & model!='BM') stop("Only Pagel's lambda and Brownian motion models have been implemented so far. Please use model='lambda' or model='BM'")
  
  # options
  args <- list(...)

  if(is.null(args[["cex"]])) cex <- 0.7 else cex <- args$cex
  if(is.null(args[["start"]])) start_lambda  <- 0.5 else start_lambda <- args$start # init values for the parameter search
  if(is.null(args[["axes"]])) axes <- 1:2 else axes <- args$axes
  if(is.null(args[["maxit"]])) maxit <- 5000 else maxit <- args$maxit
  if(is.null(args[["verbose"]])) verbose <- TRUE else verbose <- args$verbose
  
  # initialize the data
  n = Ntip(tree)
  if(REML) nLL = n-1 else nLL=n
  
  # Ordering the dataset according to the tree
  Y = Y[tree$tip.label,] #what are the potential issues here? - no rownames?
  
  p = ncol(Y)
  hidden <- which(is.na(Y))
  missing <- length(hidden)
  index = sapply(1:p, function(j){ which(is.na(Y[,j]))})
  miss_val = apply(Y,2, function(x) any(is.na(x)))
  
  # Matrix column of ones
  One <- matrix(1,ncol=1,nrow=Ntip(tree))
  
  if(missing){
    Ymean = unlist(apply(Y,2,function(x) rep(mean(x,na.rm=TRUE),length(which(is.na(x))))))
    Y_init = Y
    Y_init[hidden] = Ymean
    
    if(model=='lambda'){
      ## lambda
      obj_init = .PRPCA_fs(par=start_lambda,Y=Y_init,tree=tree,q=q,n=nLL,One=One,model=model,REML=REML)
      W_init = obj_init$W
      s2_init = obj_init$s2
      
      opt <- optim(start_lambda,
                   fn=.em_for_p3ca,
                   n=nLL,
                   q=q,
                   Y=Y,
                   tree=tree,
                   s2=s2_init,
                   W=W_init,
                   One=One,
                   model=model,
                   meanY=Ymean,
                   hidden=hidden,
                   missing=missing,
                   miss_index=index,
                   tol=tol,
                   REML=REML,
                   miss_var=miss_val,
                   verbose=FALSE,
                   ll=TRUE,
                   maxit=maxit,
                   lower=1e-7,
                   upper=1,
                   method="Brent",
                   control=list(abstol=1e-9),
                   hessian=TRUE)
      
      prpca <- .em_for_p3ca(par=opt$par,n=nLL,q=q,Y=Y,tree=tree,s2=s2_init,W=W_init,One=One,model=model,meanY=Ymean,
                            hidden=hidden,missing=missing,miss_index=index,verbose=verbose,tol=tol,REML=REML,miss_var=miss_val,maxit=maxit)
      
    }else{
      obj_init = .PRPCA_fs(par=start_lambda,Y=Y_init,tree=tree,q=q,n=nLL,One=One,model=model,REML=REML)
      W_init = obj_init$W
      s2_init = obj_init$s2
      
      prpca <- .em_for_p3ca(par=1,n=nLL,q=q,Y=Y,tree=tree,s2=s2_init,W=W_init,One=One,model=model,meanY=Ymean,hidden=hidden,
                            missing=missing,miss_index=index,verbose=verbose,tol=tol,REML=REML,miss_var=miss_val,maxit=maxit)
    }
    
  }else{
    if(model=='lambda'){
      
      if(EM){
        #W_init = rustiefel(p,q)
        W_init = rstief(p,q)
        s2_init = 0.05
        
        opt <- optim(start_lambda,
                     fn=function(param){
                       ll_opt = .em_for_p3ca_nomissing(par=param,n=nLL,
                       q=q,
                       Y=Y,
                       tree=tree,
                       s2=s2_init,
                       W=W_init,
                       One=One,
                       model=model,
                       tol=tol,
                       REML=REML,
                       maxit=maxit,
                       verbose=FALSE)$ll
                       return(-ll_opt)
                     },
                     lower=1e-7,
                     upper=1,
                     method="Brent",
                     control=list(abstol=1e-9),
                     hessian=TRUE)
        
        prpca <- .em_for_p3ca_nomissing(par=opt$par,n=nLL,
                                        q=q,
                                        Y=Y,
                                        tree=tree,
                                        s2=s2_init,
                                        W=W_init,
                                        One=One,
                                        model=model,
                                        tol=tol,
                                        REML=REML,
                                        maxit=maxit,
                                        verbose=verbose)
        
      }else{
        opt <- optim(start_lambda,
                     fn = .PRPCA_fs,
                     Y=Y,
                     tree=tree,
                     q=q,
                     n=nLL,
                     One=One,
                     model=model,
                     REML=REML,
                     ll=TRUE,
                     lower=1e-7,
                     upper=1,
                     method="Brent",
                     control=list(abstol=1e-9),
                     hessian=TRUE)
        
        prpca = .PRPCA_fs(par=opt$par, Y=Y, tree=tree, q=q, n=nLL, One=One, model=model, REML=REML)
      }
        
    }else{
      if(EM){
       # W_init = rustiefel(p,q) # TODO => init with an SVD of the original dataset transformed by current lambda and with missing values replaced by average, this will speed up convergence
       W_init = rstief(p,q) # simple code to generate a matrix on Siefel manifold rather than relying on extra libraries?
        s2_init = 0.05
        
        prpca <- .em_for_p3ca_nomissing(par=1,n=nLL,
                                        q=q,
                                        Y=Y,
                                        tree=tree,
                                        s2=s2_init,
                                        W=W_init,
                                        One=One,
                                        model=model,
                                        tol=tol,
                                        REML=REML,
                                        maxit=maxit,
                                        verbose=verbose)
      }else{
        prpca = .PRPCA_fs(par=1, Y=Y, tree=tree, q=q, n=nLL, One=One, model=model, REML=REML)
      }
    }
  }
  
  sco_loa = .scores_loadings_for_p3ca(n=nLL,W=prpca$W,s2=prpca$s2,Y_cent=prpca$Y_hat,tree_trans=prpca$tree_trans)
  
  #if(verbose) cat("s2 converge on",prpca$s2,"\n")
  
  #is there missing values ?
  miss = if(missing) missing/length(Y) else 0
  
  # Plot the P3CA axes
  if(plot){
    
    # switch between methods
    switch(model,
           "BM"={
             plot_title = paste("P3CA - prop. missing (",round(miss*100, digits=3),"%)")
           },
           "lambda"={
             plot_title = paste("P3CA (lambda=",round(opt$par, digits=3),")"," prop. missing (",round(miss*100, digits=3),"%)")
           })
    
    # plotting options
    plot(sco_loa$scores[,axes], xlab=paste("P3C",axes[1],"(",round(sco_loa$var.exp[axes[1]], digits = 3),"%)"), 
         ylab=paste("P3C",axes[2],"(",round(sco_loa$var.exp[axes[2]], digits = 3),"%)"),
         main=plot_title, las=1);
    abline(h=0,v=0)
    text(sco_loa$scores[,axes], tree$tip.label, pos=2, cex=cex)
  } 
  
  # list results
  
  if(missing>0|EM==TRUE){
    results = list(par= if(model=='lambda') opt$par else NA,
                       logLik= if(model=='lambda') -opt$value else prpca$ll,
                       W=prpca$W,
                       sigma2=prpca$s2,
                       L=sco_loa$loadings,
                       scores=sco_loa$scores,
                       vectors=sco_loa$rotations,
                       eigenval=sco_loa$val.ppca,
                       varExp=sco_loa$var.exp,
                       coef=prpca$coef,
                       model=model,
                       imputed=prpca$recon,
                       count=prpca$count,
                       opt=if(model=='lambda') opt else NA,
                       tol=tol,
                       maxit=maxit,
                       prop_missing = miss)
  }else{
    results = list(par= if(model=='lambda') opt$par else NA,
                       logLik= if(model=='lambda') -opt$value else prpca$ll,
                       W=prpca$W,
                       sigma2=prpca$s2,
                       L=sco_loa$loadings,
                       scores=sco_loa$scores,
                       vectors=sco_loa$rotations,
                       eigenval=sco_loa$val.ppca,
                       varExp=sco_loa$var.exp,
                       coef=prpca$coef,
                       opt=if(model=='lambda') opt else NA,
                       model=model,
                       prop_missing = miss)
  }
  
  # set a class to p3ca objects
  class(results) <- "p3ca"
  return(results)
}


###############################################################################
# Performing the P3CA using the analytical solutions                          #
# Only for complete data                                                      #
###############################################################################
.PRPCA_fs = function(par, Y, tree, q, n, One, model, REML=REML, ll=FALSE){
  
  if(model=='lambda') tree_scaled = .transformTree(tree, param=par, model=model)$phy else tree_scaled = tree
  
  S = vcv(tree_scaled)
  
  # Centering Y
  invS = solve(S)
  S.invSqrt = chol(invS)
  
  OnePhylo = S.invSqrt%*%One
  traitPhylo = S.invSqrt%*%Y
  b = pseudoinverse(OnePhylo)%*%traitPhylo
  
  # set missing cases to 0 to initiate the algorithm, or to some different value (e.g. random residuals ?)
  Y_new = (Y - One%*%b)
  B = S.invSqrt%*%Y_new
  tB = t(B)
  
  p = ncol(Y)
  svdB = svd(B)
  sigma2_b = sum(svdB$d[(q+1):min(n,p)]^2/n)/(p-q)
  deltaq_b = svdB$d[1:q]^2/n
  W_b = svdB$v[,1:q]%*%diag(deltaq_b - sigma2_b, q)^0.5
    
  svdW = svd(W_b)
  term1_f = (sum(colSums((B%*%svdW$u)^2)/(svdW$d^2+sigma2_b)) + sum((B-B%*%svdW$u%*%t(svdW$u))^2)/sigma2_b)/n
    term2_f = sum(log(svdW$d^2+sigma2_b)) + (p-q)*log(sigma2_b)
    if(REML){
      loglik = -(n/2) * (p*log(2*pi) + term2_f +
                           term1_f) - (p/2)*(determinant(S)$modulus) -
        (p/2)*determinant(t(One)%*%invS%*%One)$modulus
    } else{
      loglik = -(n/2) * (p*log(2*pi) + term2_f +
                           term1_f) - (p/2)*(determinant(S)$modulus) }
  
  if(ll) return(-loglik) else{
      result = list(ll=loglik, s2=sigma2_b, W=W_b, Y_hat=Y_new, tree_trans=tree_scaled, coef=b)
      return(result)
  }
}



###############################################################################
# Performing the P3CA using an EM algorithm                                   #
# Used when missing data are included                                         #
###############################################################################
.em_for_p3ca = function(par=NULL, n, q, Y, tree, s2, W, One, model, meanY,
                        tol=1e-5, hidden=NULL, missing=NULL,
                        verbose=TRUE, maxit=1000, miss_index=NULL, REML=TRUE, miss_var=NULL, ll=FALSE){
  
  # tree transformation
  if(model=='lambda') tree_scaled = .transformTree(tree, param=par, model=model)$phy else tree_scaled = tree
  
  S = vcv(tree_scaled)
  S.invSqrt = pruning(tree_scaled, inv=TRUE)$sqrtMat
  invS = crossprod(S.invSqrt)
  OnePhylo = S.invSqrt%*%One
  b = apply(Y,2, function(yh) .mu_obs(tree_scaled, yh, One))
  Y[hidden] = meanY
  Y_new = (Y - One%*%b)
  trace_mat <- estimated_missing <- 0
  p = ncol(Y_new)
  
  # init the exp_z
  M = t(W)%*%W + s2*diag(q)
  invM = solve(M)
  exp_z = invM%*%t(W)%*%t(Y_new)
  
  cov_miss_terms <- sapply(1:p, function(i){
    if(miss_var[i]==TRUE){
      cov_miss_mat(V=S, VInv=invS, missing=miss_index[[i]])
    }else{
      NULL
    }
  }, simplify=FALSE)
  
  #### Getting the trace from the list
  trace_mat <- sum(unlist(lapply(cov_miss_terms, function(i){
    if(!is.null(i)) i$tr })))
  
  # start the EM loop
  s2_reg = ll_reg = c()
  count = 0
  loglik_p = 0
  while(count!=maxit){
    M = t(W)%*%W + s2*diag(q)
    invM = solve(M)
    
    # reconstruct Yh conditionnal on Yo, x given the evolutionary model C, for each trait 'p' with missing values
    estimated_missing <- lapply(1:p, function(i){
      if(miss_var[i]==TRUE){
        #conditional_miss(S, Y[,i], missing=miss_index[[i]], W=W[i,], X=exp_z)$recon
        conditional_miss(cov_miss_terms[[i]]$VcovVobs, Y_new[,i], missing=miss_index[[i]], W=W[i,], X=exp_z)$recon
      }else{ NULL }})
    
    # adding the E[Y_m] to Yh
    Y_new[hidden] <- unlist(estimated_missing)
    Yres = S.invSqrt%*%Y_new
    
    # E[z]
    exp_z = invM%*%t(W)%*%t(Y_new)
    # E[z(S^-1)z']
    exp_zSzt = n*s2*invM + exp_z%*%invS%*%t(exp_z)
    
    # Maximization
    Wnew = (t(invS%*%Y_new)%*%t(exp_z))%*%solve(exp_zSzt)
    #t1 = trace(t(Y)%*%invS%*%Y)
    t1 = sum(Yres*Yres)
    #t2 = 2*trace(t(Y)%*%invS%*%t(Wnew%*%exp_z))
    t2 = 2*(trace(S.invSqrt%*%t(Wnew%*%exp_z)%*%t(Yres)))
    #t3 = trace(t(Wnew)%*%Wnew%*%exp_zSzt)
    t3 = sum(crossprod(Wnew)*exp_zSzt)
    s2new = (t1 - t2 + t3 + trace_mat*s2)/(n*p)
    #s2new = (trace(t(Y)%*%invS%*%Y) - 2*trace(t(Y)%*%invS%*%t(Wnew%*%exp_z)) +
    #  trace(t(Wnew)%*%Wnew%*%exp_zSzt) + trace_mat*s2)/(n*p)
    
    ## REML is only for observed because the latent variable are already centered
    ## i.e. we're not integrating across possible values for the ancestral states
    # Expected joint likelihood E_q(x) ~ p(y,x) = p(y|x)*p(x)
    
    loglik_em <- -0.5*(n*p*log(s2new) +
                         (1/s2new)*(t1 - t2 + t3) + trace(exp_zSzt))
    
    
    #loglik_em <- -0.5*(n*p*log(s2new) +
    #                      (1/s2new)*(trace(t(Y)%*%invS%*%Y)-2*trace(t(Y)%*%invS%*%t(exp_z)%*%t(Wnew)) +
    #                               trace(t(Wnew)%*%Wnew%*%exp_zSzt)) + trace(exp_zSzt))
    
    ## check convergence from s2, but otherwise can be from the log-lik...
    #if(abs(s2-s2new)<tol) {
    # if(verbose) cat("covergence at:",count,"\n")
    #break;
    #}
    
    if(abs(loglik_p-loglik_em)<tol) {
      if(verbose) cat("EM covergence at:",count,"\n")
      W = Wnew
      s2 = s2new
      break;
    }
    
    ll_reg = c(ll_reg,loglik_em)
    loglik_p = loglik_em
    count=count+1
    W = Wnew
    s2_reg = c(s2_reg,s2new)
    s2 = s2new
  }
  
  #if(verbose) cat("Finished round","\n")
  if(count==maxit) cat("non convergence","\n")
  
  #if(model=='BM'){
    log_detS = determinant(S)$modulus
    
    if(REML){
      
      ## REML is only for observed because the latent variable are already centered
      ## i.e. we're not integrating across possible values for the ancestral states
      # Expected joint likelihood E_q(x) ~ p(y,x) = p(y|x)*p(x)
      det_reml = determinant(t(One)%*%invS%*%One)$modulus
      loglik_comp = -0.5*((p+q)*n*log(2*pi) + q*log_detS + p*(log_detS + det_reml) + n*p*log(s2) +
                            (1/s2)*(trace(t(Y_new)%*%invS%*%Y_new)-2*trace(t(Y_new)%*%invS%*%t(exp_z)%*%t(W)) +
                                      trace(t(W)%*%W%*%exp_zSzt)) + trace(exp_zSzt) + trace_mat)
    }else{
      # Expected joint likelihood E_q(x) ~ p(y,x) = p(y|x)*p(x)
      loglik_comp <- -0.5*((p+q)*n*log(2*pi) + q*log_detS + p*log_detS + n*p*log(s2) +
                             (1/s2)*(trace(t(Y_new)%*%invS%*%Y_new)-2*trace(t(Y_new)%*%invS%*%t(exp_z)%*%t(W)) +
                                       trace(t(W)%*%W%*%exp_zSzt)) + trace(exp_zSzt) + trace_mat)
      #The last term is the expectation of the term for p(X) in the joint likelihood - note that this term is omitted in Li et al. because it is not used to find the solution of W and sigma^2
    }
    
    # entropy for q(X) - the latent variables
    # no need to account for ancestral state - i.e REML for the latent
    det_X = n*determinant(s2*invM)$modulus + q*log_detS
    H_x <- ((n*q)/2) * (1 + log(2*pi)) + 0.5*det_X
    
    # entropy for q(Yh) - the missing observations
    H_y = compute_log_det2(missing, S, invS, s2, index=miss_index, miss_var=miss_var)
    
    # lower bound: p(y_o) >= E_q(x,y_h)[p(y,x)] + H(q(x)) + H(q(y_h))
    low_bound_ll= loglik_comp + H_x + H_y
  #}else low_bound_ll = NA
  
  if(ll) return(-low_bound_ll) else{
      recon = Y_new + One%*%b
      
      return(list(ll=low_bound_ll, W=W, s2=s2, count=count, Y_hat=Y_new, recon=recon, coef=b, tree_trans=tree_scaled))
  }
}

.em_for_p3ca_nomissing = function(par=NULL, n, q, Y, tree, s2, W, One, model,
                        tol=1e-5, verbose=TRUE, maxit=1000,REML=TRUE){
  
  # tree transformation
  tree_scaled = .transformTree(tree, param=par, model=model)$phy
  
  S = vcv(tree_scaled)
  S.invSqrt = pruning(tree_scaled, inv=TRUE)$sqrtMat
  invS = crossprod(S.invSqrt)
  OnePhylo = S.invSqrt%*%One
  traitPhylo = S.invSqrt%*%Y
  b = pseudoinverse(OnePhylo)%*%traitPhylo
  Y_new = (Y - One%*%b)
  trace_mat <- estimated_missing <- 0
  p = ncol(Y_new)
  
  # init the exp_z
  M = t(W)%*%W + s2*diag(q)
  invM = solve(M)
  exp_z = invM%*%t(W)%*%t(Y_new)
  
  # start the EM loop
  s2_reg = ll_reg = c()
  count = 0
  loglik_p = 0
  while(count!=maxit){
    M = t(W)%*%W + s2*diag(q)
    invM = solve(M)
    Yres = S.invSqrt%*%Y_new
    
    # E[z]
    exp_z = invM%*%t(W)%*%t(Y_new)
    # E[z(S^-1)z']
    exp_zSzt = n*s2*invM + exp_z%*%invS%*%t(exp_z)
    
    # Maximization
    Wnew = (t(invS%*%Y_new)%*%t(exp_z))%*%solve(exp_zSzt)
    #t1 = trace(t(Y)%*%invS%*%Y)
    t1 = sum(Yres*Yres)
    #t2 = 2*trace(t(Y)%*%invS%*%t(Wnew%*%exp_z))
    t2 = 2*(trace(S.invSqrt%*%t(Wnew%*%exp_z)%*%t(Yres)))
    #t3 = trace(t(Wnew)%*%Wnew%*%exp_zSzt)
    t3 = sum(crossprod(Wnew)*exp_zSzt)
    s2new = (t1 - t2 + t3 + trace_mat*s2)/(n*p)
    #s2new = (trace(t(Y)%*%invS%*%Y) - 2*trace(t(Y)%*%invS%*%t(Wnew%*%exp_z)) +
    #  trace(t(Wnew)%*%Wnew%*%exp_zSzt) + trace_mat*s2)/(n*p)
    
    ## REML is only for observed because the latent variable are already centered
    ## i.e. we're not integrating across possible values for the ancestral states
    # Expected joint likelihood E_q(x) ~ p(y,x) = p(y|x)*p(x)
    
    loglik_em <- -0.5*(n*p*log(s2new) +
                         (1/s2new)*(t1 - t2 + t3) + trace(exp_zSzt))
    
    
    #loglik_em <- -0.5*(n*p*log(s2new) +
    #                      (1/s2new)*(trace(t(Y)%*%invS%*%Y)-2*trace(t(Y)%*%invS%*%t(exp_z)%*%t(Wnew)) +
    #                               trace(t(Wnew)%*%Wnew%*%exp_zSzt)) + trace(exp_zSzt))
    
    ## check convergence from s2, but otherwise can be from the log-lik...
    #if(abs(s2-s2new)<tol) {
    # if(verbose) cat("covergence at:",count,"\n")
    #break;
    #}
    
    if(abs(loglik_p-loglik_em)<tol) {
      if(verbose) cat("EM covergence at:",count,"\n")
      W = Wnew
      s2 = s2new
      break;
    }
    
    ll_reg = c(ll_reg,loglik_em)
    loglik_p = loglik_em
    count=count+1
    W = Wnew
    s2_reg = c(s2_reg,s2new)
    s2 = s2new
  }
  
  #if(verbose) cat("Finished round","\n")
  if(count==maxit) cat("non convergence","\n")
  
  log_detS = determinant(S)$modulus
    
    if(REML){
      
      ## REML is only for observed because the latent variable are already centered
      ## i.e. we're not integrating across possible values for the ancestral states
      # Expected joint likelihood E_q(x) ~ p(y,x) = p(y|x)*p(x)
      det_reml = determinant(t(One)%*%invS%*%One)$modulus
      loglik_comp = -0.5*((p+q)*n*log(2*pi) + q*log_detS + p*(log_detS + det_reml) + n*p*log(s2) +
                            (1/s2)*(trace(t(Y_new)%*%invS%*%Y_new)-2*trace(t(Y_new)%*%invS%*%t(exp_z)%*%t(W)) +
                                      trace(t(W)%*%W%*%exp_zSzt)) + trace(exp_zSzt) + trace_mat)
    }else{
      # Expected joint likelihood E_q(x) ~ p(y,x) = p(y|x)*p(x)
      loglik_comp <- -0.5*((p+q)*n*log(2*pi) + q*log_detS + p*log_detS + n*p*log(s2) +
                             (1/s2)*(trace(t(Y_new)%*%invS%*%Y_new)-2*trace(t(Y_new)%*%invS%*%t(exp_z)%*%t(W)) +
                                       trace(t(W)%*%W%*%exp_zSzt)) + trace(exp_zSzt) + trace_mat)
      #The last term is the expectation of the term for p(X) in the joint likelihood - note that this term is omitted in Li et al. because it is not used to find the solution of W and sigma^2
    }
    
    # entropy for q(X) - the latent variables
    # no need to account for ancestral state - i.e REML for the latent
    det_X = n*determinant(s2*invM)$modulus + q*log_detS
    H_x <- ((n*q)/2) * (1 + log(2*pi)) + 0.5*det_X
    
    # lower bound: p(y_o) >= E_q(x,y_h)[p(y,x)] + H(q(x)) + H(q(y_h))
    low_bound_ll= loglik_comp + H_x

  return(list(ll=low_bound_ll, W=W, s2=s2, count=count, Y_hat=Y_new, tree_trans=tree_scaled, recon=NULL, coef=b))
}

###############################################################################
# Estimating the mean for observed values                                     #
###############################################################################
.mu_obs = function(phy,y,X){
  pos.rem = which(is.na(y))
  
  if(length(pos.rem)==0) newphy=phy else{
    X = as.matrix(X[-pos.rem,])
    to.rem=phy$tip.label[pos.rem]
    newphy = drop.tip(phy,tip=to.rem)
  }
  
  pt <- pruning(newphy, inv=TRUE)
  Xi <- pt$sqrtMat%*%X
  Xhat <- pseudoinverse(Xi)
  Yi <- pt$sqrtMat%*%na.omit(y)
  mu <- Xhat%*%Yi
  
  return(mu)
}

###############################################################################
# Estimating the phylogenetic covariance matrices for missing and             #
# observed values                                                             #
###############################################################################
cov_miss_mat <- function(V, VInv, missing, ...){
  
  # estimate the expectation for y_h
  Vobs <- V[-missing,-missing, drop=FALSE]
  Vcov <- V[missing,-missing, drop=FALSE]
  
  VcovVobs = Vcov%*%solve(Vobs)
  
  # estimate the conditional covariance
  Vmiss <- V[missing, missing, drop=FALSE]
  covC <- Vmiss - VcovVobs%*%t(Vcov)
  
  #VmissI = solve(V)#solve(V)[missing, missing, drop=FALSE]
  #trC <- sum(diag(solve(Vmiss)%*%covC)) # won't work, because we should take into account all the species for computing the inverse
  trC <- sum(diag(VInv[missing, missing, drop=FALSE]%*%covC))
  results = list(covC=covC, tr=trC, VcovVobs=as.numeric(VcovVobs))
  return(results)
}

###############################################################################
# Estimating the expected values for missing data conditional on observed     #
# values, W, s2, and C                                                        #
###############################################################################
conditional_miss <- function(covV, Y, missing, W, X, ...){
  
  # estimate mu
  mu_h = W%*%X[,missing]
  mu_o = W%*%X[,-missing]
  
  # estimate the expectation for y_h
  #Vobs <- V[-missing,-missing, drop=FALSE] # it can be Vobs <- solve(V)[-missing,-missing, drop=FALSE]
  #Vcov <- V[missing,-missing, drop=FALSE]
  recon = mu_h + as.numeric(matrix(covV,nrow=length(missing))%*%t(Y[-missing] - mu_o))
  
  results = list(recon=recon)
  return(results)
}

###############################################################################
# Computing the entropy for missing values                                    #
###############################################################################
compute_log_det2 <- function(missing, C_lambda, C_lambda_inv, sigma2, index, miss_var) {
  
  log_dets <- sapply(1:length(miss_var), function(i){
    if(miss_var[i]==TRUE){
      as.numeric(determinant(cov_miss_mat(C_lambda*sigma2, C_lambda_inv/sigma2, missing=index[[i]])$covC)$modulus)
    }else{
      NULL
    }
  })
  
  # Entropy for the missing cases
  total_log_det <- 0.5*sum(unlist(log_dets)) + (missing/2)*(log(2*pi) + 1)
  
  #total_log_det <- 0.5*sum(unlist(log_dets)) + (missing/2)*(log(2*pi) + 1) +  ((missing)/2)*log(sigma2)
  # return the results
  return(total_log_det)
}

###############################################################################
# Computing the scores and loadings for the P3CA                              #
###############################################################################
.scores_loadings_for_p3ca = function(n,W,s2,Y_cent,tree_trans){
  
  p = dim(W)[1]
  q = dim(W)[2]
  eig_W = svd(W)
  scores <-  Y_cent%*%eig_W$u
  rownames(scores) = rownames(Y_cent)
  
  val.ppca <- c(eig_W$d^2 + s2,rep(s2,(p-q)))
  var.exp = (val.ppca*100/sum(val.ppca))[1:q] ##Only the values corresponding to q should be returned?
  
  # WtW + sigma2*I
  diag.R <- rowSums(W^2)+s2#diag(eig_W$u%*%(eig_W$d^2*t(eig_W$u))) + s2
  sqrtM <- pruning(tree_trans, inv=TRUE)$sqrtMat
  
  # compute the loadfit_bmm_diet_meansings => correlations between axes and data [check Revell 2010 for the pPCA version]
  Ccv<-(t(Y_cent)%*%crossprod(sqrtM)%*%scores)/n # compute cross covariance matrix and loadings
  L<-matrix(0,p,min(p,q),dimnames=list(colnames(Y_cent),paste("PC",1:min(n,q),sep="")))
  for(i in 1:p) for(j in 1:min(p,q)) L[i,j]<-Ccv[i,j]/(sqrt(diag.R[i])*sqrt(val.ppca[j]))
  
  return(list(scores=scores, loadings = L, val.ppca = val.ppca, var.exp = var.exp, rotations=eig_W$u))
  
}

###############################################################################
# Miscellaneous functions: trace and vector                                   #
###############################################################################
trace = function(X) sum(diag(X))
vec = function(x) as.numeric(x)

# ------------------------------------------------------------------------- #
# print option for p3ca                                                     #
# options: x, digits, ...                                                   #
#                                                                           #
# ------------------------------------------------------------------------- #

print.p3ca <- function(x, digits = max(3L, getOption("digits") - 3L), ...){
    
    cat("## Probabilistic Phylogenetic PCA                 ##","\n")
    if(x$model=="lambda") cat("## Pagel' lambda                                  ##","\n") else cat("## Brownian motion                                ##","\n")
    print(x$par, digits=digits)
    cat("## Proportion of variance explained for the q PCs ##","\n")
    print(x$varExp, digits=digits)
    
}

# ------------------------------------------------------------------------- #
# generate a matrix on Stiefel manifold                                     #
# options: n, p                                                             #
#                                                                           #
# ------------------------------------------------------------------------- #

rstief <- function(n,p){
  X <- matrix(rnorm(n*p), ncol=p, nrow=n)
  q = qr(X)
  return(qr.Q(q))
}
