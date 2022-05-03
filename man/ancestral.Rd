\name{ancestral}
\alias{ancestral}

\title{
Estimation of traits ancestral states.
}
\description{
Reconstruct the ancestral states at each node of a phylogenetic tree from models fit obtained using the \code{mvgls} function. For models of the class \code{mvXX} this is a wrapper to the function \code{estim}}
\usage{

ancestral(object, ...)
  
}
\arguments{
  \item{object}{
A model fit object obtained by the \code{mvgls} function.
}
 \item{...}{
Further options to be passed through. For instance, if a regression model is used, values for the predictor(s) at each node of the tree should be given in a matrix to the \code{newdata} argument. If a model of the type \code{mvXX} is used, the argument \code{tree} and \code{data} should be provided like in \code{estim} function.
}
}
 

\value{
	a matrix of reconstructed ancestral states for each node (note that the numerotation of the ancestral states starts at "N+1" [for the root], where N is the number of species in the tree)

  }

\details{
\code{ancestral} is an S3 method that reconstruct the ancestral states at each nodes of a phylogenetic tree from the models fit obtained by the \code{mvgls} function (Clavel et al. 2019). Ancestral states are estimated using generalized least squares (GLS; Martins & Hansen 1997, Cunningham et al. 1998 ). Note that when a regression model (rather than an intercept only model of the form Y~1) is provided, an argument "newdata" with a matrix of regressor values for each node should be provided (similar to what is done in the "predict" function). 

}

\note{
	The function is similar to the one used with \code{fit_t_pl} from the RPANDA package (Clavel et al. 2019).}
\references{

Clavel, J., Aristide, L., Morlon, H., 2019. A Penalized Likelihood framework for high-dimensional phylogenetic comparative methods and an application to new-world monkeys brain evolution. Syst. Biol. 68: 93-116.

Cunningham C.W., Omland K.E., Oakley T.H. 1998. Reconstructing ancestral character states: a critical reappraisal. Trends Ecol. Evol. 13:361-366.

Martins E.P., Hansen T.F. 1997. Phylogenies and the comparative method: a general approach to incorporating phylogenetic information into the analysis of interspecific data. Am. Nat. 149:646-667.
}

\author{J. Clavel}

\seealso{
\code{\link{mvgls}}, 
\code{\link{estim}},
\code{\link{predict.mvgls}}
}

\examples{
\donttest{

set.seed(1)
n <- 32 # number of species
p <- 5 # number of traits

tree <- pbtree(n=n, scale=1) # phylogenetic tree
R <- crossprod(matrix(runif(p*p), ncol=p)) # a random covariance matrix
# simulate a BM dataset
Y <- mvSIM(tree, model="BM1", nsim=1, param=list(sigma=R, theta=rep(0,p))) 
data=list(Y=Y)

fit <- mvgls(Y~1, data=data, tree, model="BM", method="LL")

# Perform the ancestral states reconstruction
anc <- ancestral(fit)

# retrieve the ancestral states
head(anc)

}
}
