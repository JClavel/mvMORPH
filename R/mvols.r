################################################################################
##                                                                            ##
##                       mvMORPH: mvlm.r/mvols.r                              ##
##                                                                            ##
##      Multivariate Ordinary Least Squares Linear Models by ML and PL        ##
##                                                                            ##
##  This is a wrapper to the mvgls function when there's no trees             ##
##  Created by Julien Clavel - 31-09-2020 (terminated 6-06-2022)              ##
##  (julien.clavel@hotmail.fr/ julien.clavel@univ-lyon1.fr)                   ##
##   require: phytools, ape, corpcor, subplex, spam, glassoFast, stats        ##
##                                                                            ##
################################################################################


mvols <- function(formula, data=list(), method=c("PL-LOOCV","LL"), REML=TRUE, ...){
  
  # Recover options
  args <- list(...)
  if(is.null(args[["weights"]])) weights <- NULL else weights <- args$weights
  
  # Retrieve the dataset
  model_fr = model.frame(formula=formula, data=data)
  Y = model.response(model_fr)
  
  # Number of observations
  n = nrow(Y)
  
  # For now I'm using a simple wrapper to mvgls. The idea is to fix the tree structure to a star tree 
  # (with possibly different lengths if a weighted least squares approach is required), and use the BM model.
  phylo_struct <- generate_tree_structure(n, names_tips = rownames(Y))
  
  # Weighted least squares
  if(!is.null(weights)){
    if(!is.null(rownames(Y))){
      if(any(!phylo_struct$tip.label%in%names(weights))) stop("weights vector should be named and match the names of the observations in the response variable.")
      weights <- weights[phylo_struct$tip.label]
      phylo_struct$edge.length[sapply(1:n,function(x) which(x==phylo_struct$edge[,2]))] <- weights
    }else{
      phylo_struct$edge.length[phylo_struct$edge[,2]<=n] <- weights
    }
  }
  
  # Call mvgls => force BM on a star tree
  results <- mvgls(formula = formula, data = data, tree = phylo_struct, model = "BM", method = method, REML = REML, ...)
  
  # Define a broader class
  results$call = match.call()
  class(results) <- c("mvols","mvgls")
  return(results)
}

# ------------------------------------------------------------------------- #
# generate_tree_structure                                                   #
# options: n, names_tips, ...                                               #
# ------------------------------------------------------------------------- #
generate_tree_structure <- function(n, names_tips=NULL, ...){
  
  # bar & foo functions from 'rtree' in 'ape' package v.5.6-2
  bar <- function(n) sample.int(n - 1L, 1L, FALSE, NULL, FALSE)
  foo <- function(n, pos) {
    n1 <- bar(n)
    n2 <- n - n1
    po2 <- pos + 2L * n1 - 1L
    edge[c(pos, po2), 1L] <<- nod
    nod <<- nod + 1L
    if (n1 > 2L) {
      edge[pos, 2L] <<- nod
      foo(n1, pos + 1L)
    }
    else if (n1 == 2L) {
      edge[pos + 1:2, 1L] <<- edge[pos, 2L] <<- nod
      nod <<- nod + 1L
    }
    if (n2 > 2L) {
      edge[po2, 2L] <<- nod
      foo(n2, po2 + 1L)
    }
    else if (n2 == 2L) {
      edge[po2 + 1:2, 1L] <<- edge[po2, 2L] <<- nod
      nod <<- nod + 1L
    }
  }
  
  # generate the edge matrix
  nbr <- 2L * n - 2L 
  edge <- matrix(NA_integer_, nbr, 2L)
  
  # are the data named?
  if (is.null(names_tips)) {
    names_tips <- paste0("obs", 1:n)
  }
  
  # populate the edge matrix
  if (n == 2L) {
    edge[] <- c(3L, 3L, 1L, 2L)
  }
  else if (n == 3L) {
    edge[] <- c(4L, 5L, 5L, 4L, 5L, 1:3)
  }
  else if (n > 3L) {
    nod <- n + 1L
    foo(n, 1L)
    i <- which(is.na(edge[, 2L]))
    edge[i, 2L] <- 1:n
  }
  
  # make the "tree" object
  phy_struct <- list(edge = edge, tip.label = names_tips)
  phy_struct$Nnode <- if (n == 1L)  1L else n - 1L
  class(phy_struct) <- "phylo"
  attr(phy_struct, "order") <- "cladewise"
  
  # create branch lengths
  phy_struct$edge.length <- numeric(nbr)
  phy_struct$edge.length[edge[,2]<=n] <- 1
  
  return(phy_struct)
}
