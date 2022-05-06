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
