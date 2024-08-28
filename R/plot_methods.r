################################################################################
##                                                                            ##
##                       mvMORPH: plot_methods.r                              ##
##                                                                            ##
##   Internal functions for the mvMORPH package                               ##
##                                                                            ##
##  Created by Julien Clavel - 01-12-2020                                     ##
##  (julien.clavel@hotmail.fr/ julien.clavel@univ-lyon1.fr)                   ##
##   require: phytools, ape, corpcor, subplex, spam, glassoFast, stats        ##
##                                                                            ##
################################################################################

# ---------------------------------------------------------------------------- #
# mvqqplot                                                                     #
# options: object, conf, ...                                                   #
#                                                                              #
# ---------------------------------------------------------------------------- #
mvqqplot <- function(object, conf=0.95, ...){
  
  # Check if weights are given (e.g. for OU)
  if(is.null(object$corrSt$diagWeight)){
    diagWeight <- 1; is_weight = FALSE
  }else{
    diagWeight <- object$corrSt$diagWeight; is_weight = TRUE
    diagWeightInv <- 1/diagWeight
  }
  
  # Mahalanobis distances - studentized residuals (e.g., Caroni 1987)
  if(object$method=="LL"){
    # FIXME - weighted transformations
    C <- vcv.phylo(object$corrSt$phy)
    Cinv <- solve(C)
    Di <- chol(Cinv)
    if(is_weight){
      resid <- Di%*%(object$variables$Y - object$variables$X%*%object$coefficients)*diagWeightInv
      X <- Di%*%(object$variables$X*diagWeightInv)
    }else{
      resid <- Di%*%(object$variables$Y - object$variables$X%*%object$coefficients)
      X <- Di%*%object$variables$X
    }
    V=diag(X%*%pseudoinverse(X))
    
    # scaled mahalanobis distance
    r2 <- sapply(1:Ntip(object$corrSt$phy), function(i) t(resid[i,])%*%object$sigma$P%*%resid[i,] / (1-V[i]))
    
    # plot the quantile-quantile
    df <- object$dims$p
    order_stat <- order(r2)
    r2 <- r2[order_stat]
    n <- length(r2)
    p <- ppoints(n)
    chi2q <- qchisq(p, df)
    plot(chi2q, r2, main="Chi-Square Q-Q Plot", xlab="Chi-Square quantiles", ylab="Squared Mahalanobis distances", las=1)
    
    # expected multivariate normal line - or take the interquartile?
    abline(0,1)
    
    # confidence interval for the order statistic (parametric if ML)
    f <- dchisq(chi2q, df)
    se <- sqrt( p * (1-p) / (n*f^2))
    ci <- qnorm(1 - (1 - conf)/2)
    lower <- chi2q - se * ci # 95% CI
    upper <- chi2q + se * ci
    g1 = g2 = 1
    
    # Check the values using a cut-off
    crit_maha <- qchisq(conf, df)
    
  }else{
    
    # hat matrix
    n <- object$dims$n
    alpha <- object$tuning
    B <- object$coefficients
    Cinv <- solve(vcv.phylo(object$corrSt$phy))
    Di <- chol(Cinv)
    if(is_weight) X <- Di%*%(object$variables$X*diagWeightInv) else X <- Di%*%object$variables$X
    if(is_weight) Y <- Di%*%(object$variables$Y*diagWeightInv) else Y <- Di%*%object$variables$Y
    XtX <- pseudoinverse(X)
    h <- diag(X%*%XtX)
    resid <- Y - X%*%B
    
    # cross-validated Mahalanobis distance
    r2 <- sapply(1:n, function(x){
      Bx <- B - tcrossprod(XtX[,x], resid[x,])/(1-h[x]) # rank-1 update
      # update the residuals
      residuals2 <- Y - X%*%Bx
      Skpartial <- crossprod(residuals2[-x,])/(n-1)
      # target matrix
      target <- .targetM(Skpartial, object$target, object$penalty)
      # Mahalanobis distance
      .regularizedMahalanobis(Skpartial, resid[x,], alpha, object$target, target, object$penalty) / (1-h[x])
    })
    
    # calculate the chi2 distribution using moment matching (cf. Nomikos & McGregor 1995)
    # replace s2 and m by robust estimates? FIXME
    s2 <- var(r2); m <- mean(r2)
    g1 <- s2/(2*m); g2 <- (2*m^2)/s2
    
    # first check the values using a cut-off
    crit_maha <- g1*qchisq(conf,g2)
    s2 <- var(r2[r2<crit_maha]); m <- mean(r2[r2<crit_maha])
    g1 <- s2/(2*m); g2 <- (2*m^2)/s2
    
    # plot the quantile-quantile
    order_stat <- order(r2)
    r2 <- r2[order_stat]
    p <- ppoints(n)
    chi2q <- g1*qchisq(p, g2)
    plot(chi2q, r2, main="Scaled Chi-Square Q-Q Plot", xlab="Chi-Square quantiles", ylab="Squared Mahalanobis distances", las=1)
    
    # expected multivariate normal line - or take the interquartile?
    abline(0,1)
    
    # confidence interval for the order statistic - take the scaled chi2q value with momment matching df
    # If X~chisq(k) then cX~gamma(k/2; 2c)
    f <- dgamma(chi2q, g2/2, scale=2*g1)
    se <- sqrt( p * (1-p) / (n*f^2))
    
    ci <- qnorm(1 - (1 - conf)/2) # 95% CI should be around ~2
    lower <- chi2q - se * ci 
    upper <- chi2q + se * ci
  }
  
  # # plot confidence lines
  lines(chi2q, lower, ...)
  lines(chi2q, upper, ...)
  
  # plot outliers??
  ident_out <- which(r2>upper | r2<lower)
  if(!length(ident_out)==0){
    names_sp <- object$corrSt$phy$tip.label[order_stat]
    names_out <- names_sp[ident_out]
    pos <- rep(2, length(ident_out))
    pos[r2[ident_out]<lower[ident_out]] <- 4
    text(chi2q[ident_out], r2[ident_out], labels=names_out, pos=pos, offset = 0.5)
  }else{
    names_out <- NULL
  }
  
  # retrieve values above Chi2 cutoff
  names_sp <- object$corrSt$phy$tip.label[order_stat]
  crit = names_sp[r2>=crit_maha]
  
  results <- list(squared_dist=r2, chi2q=chi2q, outlier=names_out, upper=upper, lower=lower, g1=g1, g2=g2, crit=crit)
  invisible(results)
}


# ------------------------------------------------------------------------- #
# .regularizedLik return the log-lik with the regularized estimate          #
# options: S, residuals, lambda, targM, target, penalty                     #
#                                                                           #
# ------------------------------------------------------------------------- #
.regularizedMahalanobis <- function(S, residuals, lambda, targM, target, penalty){
  
  switch(penalty,
         "RidgeArch"={
           G <- (1-lambda)*S + lambda*target
           Gi <- try(chol(G), silent=TRUE)
           if(inherits(Gi, 'try-error')) return(1e6)
           rk <- sum(backsolve(Gi, residuals, transpose = TRUE)^2)
         },
         "RidgeAlt"={
           quad <- .makePenaltyQuad(S, lambda, target, targM)
           Gi <- quad$P
           Swk <- tcrossprod(residuals)
           rk <- sum(Swk*Gi)
         },
         "LASSO"={
           LASSO <- glassoFast(S, lambda, maxIt=500)
           Gi <- LASSO$wi;
           Swk <- tcrossprod(residuals);
           rk <- sum(Swk*Gi);
         })
  
  return(rk)
}


# ------------------------------------------------------------------------ #
#                                                                          #
#  pcaShape: projection of shape changes along PC axes                     #
#  J. Clavel - 2019                                                        #
#                                                                          #
# ------------------------------------------------------------------------ #

# "object" is an object fit of class "mvgls"
# "axis" is the PC axis chosen
# "ndim" the number of dimension for each landmarks (2 or 3), or 1 if not a landmark...
# "landmarks" the number of landmarks
# "spp" the name of a given species, otherwise the min and max for the PC axis is used

pcaShape <- function(object, axis=1, ndim=3, spp=NULL, plot=FALSE, ...){
 
  if(!inherits(object,"mvgls")) stop("only works with \"mvgls\" class objects. See ?mvgls or ?mvols")
  par = list(...)
  phyl_pca <- mvgls.pca(object, plot=FALSE)
  pcscores <- phyl_pca$scores
  ancestral <- coef(object)
  if(is.null(par[["landmarks"]])) landmarks = ncol(pcscores)/ndim else landmarks = par$landmarks

  if(is.null(spp)){
    pc_axis_min <- min(pcscores[, axis])%*%phyl_pca$vectors[,axis] + ancestral["(Intercept)",]
    pc_axis_max <- max(pcscores[, axis])%*%phyl_pca$vectors[,axis] + ancestral["(Intercept)",]
    shape <- list()
    shape$min <- matrix(pc_axis_min, ncol=ndim, nrow=landmarks, byrow = TRUE)
    shape$max <- matrix(pc_axis_max, ncol=ndim, nrow=landmarks, byrow = TRUE)
  }else{
    shape <- lapply(spp, function(sp_name){
        pc_specimen_shape <- pcscores[sp_name, axis]%*%phyl_pca$vectors[,axis] + ancestral["(Intercept)",]
        shape <- matrix(pc_specimen_shape, ncol=ndim, nrow=landmarks, byrow=TRUE)
    })
    names(shape) = spp
  }
  
  # plot the shape changes?
  if(plot){
      if(is.null(spp)){
          plot(shape$min, pch=16, xlab="", ylab="", main=paste("Shapes changes along PC",axis))
          points(shape$max, pch=16, col="red")
          legend("bottomright", legend=c("max scores","min scores"), pch=16, col=c("red","black"))
      }else{
          plot(shape[[1]], pch=16, xlab="", ylab="", main=paste("Shapes changes along PC",axis))
          if(length(spp)>1){
              for(i in 2:length(spp)) points(shape[[i]], pch=16, col=i)
              legend("bottomright", legend=spp, pch=16, col=1:length(spp))
          }
      }
  }
  
  # return the shapes
  return(shape)
}

## ------------------------------------------------------- ##
##                                                         ##
##  dfaShape: project shape changes along PC axes          ##
##  J. Clavel - 2020 ; require mvMORPH 1.1.4               ##
##                                                         ##
## ------------------------------------------------------- ##

# "object" is an object fit of class "mvgls"
# "reference" the reference shape used to compare the deformations. Usually the mean shape
# "axis" is the PC axis chosen
# "ndim" the number of dimension for each landmarks (2 or 3), or 1 if not a landmark...
# "landmarks" the number of landmarks
# "spp" the name of a given species, otherwise the min and max for the PC axis is used
# "scaling" is an arbitrary factor used to multiply the effects (for better visualization)

dfaShape <- function(object, reference, axis=1, ndim=3, spp=NULL, scaling=1, plot=FALSE, ...){
  
  if(!inherits(object,"mvgls.dfa")) stop("only works with \"mvgls.dfa\" class objects. See ?mvgls.dfa")
  par = list(...)
  
  pcscores <- object$scores
  if(is.null(reference) | missing(reference)) reference <- coef(object$fit) # should instead fit an intercept only model
  # we standardize the discriminant coefficients by the covariance because it has been used to standardize the between-group covariance
  # see e.g. J. Claude 2008 - chapter 6
  DF <- object$fit$sigma$Pinv%*%object$coeffs
  
  if(is.null(par[["landmarks"]])) landmarks = ncol(reference)/ndim else landmarks = par$landmarks
  if(is.null(spp)){
    pc_axis_min <- scaling*min(pcscores[, axis])%*%DF[,axis] + reference["(Intercept)",] # to replace by first line or provide an option here 
    pc_axis_max <- scaling*max(pcscores[, axis])%*%DF[,axis] + reference["(Intercept)",]
    shape <- list()
    shape$min <- matrix(pc_axis_min, ncol=ndim, nrow=landmarks, byrow = TRUE)
    shape$max <- matrix(pc_axis_max, ncol=ndim, nrow=landmarks, byrow = TRUE)
  }else{
    shape <- lapply(spp, function(sp_name){
                        pc_specimen_shape <- scaling*pcscores[sp_name, axis]%*%DF[,axis] + reference["(Intercept)",]
                        matrix(pc_specimen_shape, ncol=ndim, nrow=landmarks, byrow=TRUE)
                    })
    names(shape) = spp
  }
  
  # plot the shape changes?
  if(plot){
      if(is.null(spp)){
          plot(shape$min, pch=16, xlab="", ylab="", main=paste("Shapes changes along DF",axis))
          points(shape$max, pch=16, col="red")
          legend("bottomright", legend=c("max scores","min scores"), pch=16, col=c("red","black"))
      }else{
          plot(shape[[1]], pch=16, xlab="", ylab="", main=paste("Shapes changes along DF",axis))
          if(length(spp)>1){
              for(i in 2:length(spp)) points(shape[[i]], pch=16, col=i)
              legend("bottomright", legend=spp, pch=16, col=1:length(spp))
          }
      }
  }
  
  
  # return the shapes
  return(shape)
}
