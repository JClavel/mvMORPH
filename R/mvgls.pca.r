################################################################################
##                                                                            ##
##                       mvMORPH: mvgls.pca.r                                 ##
##                                                                            ##
##   PCA on the Multivariate Generalized Least Squares                        ##
##   covariance matrix estimator                                              ##
##                                                                            ##
##  Created by Julien Clavel - 31-07-2018                                     ##
##  (julien.clavel@hotmail.fr/ julien.clavel@biologie.ens.fr)                 ##
##   require: phytools, ape, corpcor, subplex, spam, glassoFast, stats        ##
##                                                                            ##
################################################################################

mvgls.pca <- function(object, plot=TRUE, ...){
  
    # optional arguments
    args <- list(...)
    if(is.null(args[["axes"]])) axes <- c(1,2) else axes <- args$axes
    if(is.null(args[["col"]])) col <- "black" else col <- args$col
    if(is.null(args[["pch"]])) pch <- 19 else pch <- args$pch
    if(is.null(args[["cex"]])) cex <- 0.7 else cex <- args$cex
    if(is.null(args[["las"]])) las <- 1 else las <- args$las
    if(is.null(args[["main"]])) main <- "Regularized Phylogenetic PCA" else main <- args$main
    if(is.null(args[["mode"]])) mode <- "cov" else mode <- args$mode
    
    # if correlation matrix?
    if(!inherits(object,"mvgls")) stop("only works with \"mvgls\" class objects. See ?mvgls")
    covR <- object$sigma$Pinv
    if(mode=="corr") covR <- cov2cor(covR)

    # compute the scores
    eig <- eigen(covR)
    values <- eig$values
    U <- eig$vectors
    resids <- object$residuals
    S <- resids%*%U
    
    # plot
    if(plot){
        # contribution % variance
        tot<-sum(values)
        valX<-round(values[axes[1]]*100/tot,digits=2)
        valY<-round(values[axes[2]]*100/tot, digits=2)
        xlabel <- paste("PC",axes[1]," (",valX," %)", sep="")
        ylabel <- paste("PC",axes[2]," (",valY," %)", sep="")
        plot(S[,axes], main=main, xlab=xlabel, ylab=ylabel, pch=pch, col=col, las=las)
        abline(h=0,v=0)
        text(S[,axes],object$corrSt$phy$tip.label, pos=2, cex=cex)
    }
    
    res <- list(scores=S, values=values, vectors=U, rank=qr(covR)$rank)
    class(res) <- "mvgls.pca"
    invisible(res)
}


