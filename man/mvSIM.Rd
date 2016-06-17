\name{mvSIM}
\alias{mvSIM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulation of (multivariate) continuous traits on a phylogeny
%%  ~~function to do ... ~~
}
\description{
This function allows simulating multivariate (as well as univariate) continuous traits evolving according to a BM (Brownian motion), OU (Ornstein-Uhlenbeck), ACDC (Accelerating rates and Decelerating rates/Early bursts), and SHIFT models of phenotypic evolution.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mvSIM(tree, nsim = 1, error = NULL, model = c("BM1", "BMM", "OU1", "OUM", "EB"),
                        param = list(theta = 0, sigma = 0.1, alpha = 1, beta = 0))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree}{
Phylogenetic tree with mapped ancestral states in SIMMAP format.
    (See make.simmap function from phytools package) or a standard "phylo" object (ape). Or a time-series
%%     ~~Describe \code{tree} here~~
}
  \item{nsim}{
The number of simulated traits (or datasets for multivariate analysis).
%%     ~~Describe \code{nsim} here~~
}
  \item{error}{
Matrix or data frame with species in rows and continuous trait standard errors (squared) in columns.
%%     ~~Describe \code{error} here~~
}
  \item{model}{
The model of trait evolution for the simulations. Could be any of the models used by the mvBM, mvEB, mvOU and mvSHIFT functions.
%%     ~~Describe \code{model} here~~
}
  \item{param}{
List of parameter arguments used for the simulations. You should provide the sigma (values or matrix), alpha (for OU and SHIFT models), beta (EB and SHIFT), theta (ancestral states), ntraits (the number of traits) or others param arguments used in the models.
Alternatively you can provide a fitted object of class "mvmorph". See details.
%%     ~~Describe \code{param} here~~
}
}
\details{
This function simulates multivariate (as well as univariate) continuous traits evolving according to a BM/RW (Brownian motion/Random walk), OU (Ornstein-Uhlenbeck), ACDC (Accelerating rates and Decelerating rates/Early Bursts), and SHIFT models of phenotypic evolution. The traits are simulated by random sampling from a Multivariate Normal Distribution (Paradis, 2012).


The mvSIM function allows simulating continuous trait (univariate or multivariate) evolution along a phylogeny (or a time-series) with user specified parameters or estimated parameters from a previous fit.

The "simulate" wrapper can also be used with a fitted object of class "mvmorph":
    simulate(object, nsim=1, tree=tree). See example.

If parameters values are not provided, the default values are fixed to 1 (sigma, sig, alpha, beta) or to 0 for the mean at the root (ancestral state).

For the "BMM" model were different parts of the tree have their own rate, a list with one rate (or matrix of rates) per selective regimes must be provided.

For the "OU1" and "OUM" models, the user can define if the ancestral state (theta0) should be computed (param$root=TRUE), assumed to be at the oldest regime state (param$root=FALSE), or that there is no root and each regimes is at the stationary point (param$root="stationary"). (See also ?mvOU).

For the "BM1", "BMM", and "RWTS" models, a trend can be simulated by providing the "trend" argument and values in the "param" list.

Traits names can be provided with the "names_traits" argument in the "param" list.
For all the shift models, if the tree is not mapped the age of the shift should be directly provided (in unit of times of the tree) by the "age" argument in the "param" list.


%%  ~~ If necessary, more details than the description above ~~
}
\value{
A matrix with simulated traits (columns) for the univariate case, or a list of matrix for the multivariate case (with nsim>1).
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Paradis E. 2012. Analysis of Phylogenetics and Evolution with R. New York: Springer.

%% ~put references to the literature/web site here ~
}
\author{
Julien Clavel
%%  ~~who you are~~
}
\note{
Ancestral states for Ornstein-Uhlenbeck processes (param$root=TRUE) should be used with non-ultrametric trees.
As this method uses Multivariate Normal distribution (MVN) for simulating the traits, it is advised to avoid its use with very large datasets/trees and rely instead on recursive algorithms (see e.g., ?rTraitCont from "ape").
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{mvMORPH}}
\code{\link{mvOU}}
\code{\link{mvEB}}
\code{\link{mvBM}}
\code{\link{mvSHIFT}}
\code{\link{mvRWTS}}
\code{\link{mvOUTS}}
\code{\link{mvLL}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

## Simulated dataset
set.seed(14)
# Generating a random tree with 50 species
tree<-pbtree(n=50)

# Setting the regime states of tip species
sta<-as.vector(c(rep("Forest",20),rep("Savannah",30))); names(sta)<-tree$tip.label

# Making the simmap tree with mapped states
tree<-make.simmap(tree,sta , model="ER", nsim=1)
col<-c("blue","orange"); names(col)<-c("Forest","Savannah")

# Plot of the phylogeny for illustration
plotSimmap(tree,col,fsize=0.6,node.numbers=FALSE,lwd=3, pts=FALSE)

## Simulate trait evolution according to a bivariate "BMM" model
# Number of traits
ntraits<-2
# Number of simulated (pairs of) traits
nsim<-10
# Rates matrices for the "Forest" and the "Savannah" regimes
sigma<-list(Forest=matrix(c(2,0.5,0.5,1),2), Savannah=matrix(c(5,3,3,4),2))
# ancestral states for each traits
theta<-c(0,0)

# Simulate
simul<-mvSIM(tree,nsim=nsim, model="BMM",param=list(ntraits=ntraits,sigma=sigma,theta=theta))

# Try to fit a "BM1" model to the first simulated dataset
model_fit<-mvBM(tree,simul[[1]],model="BM1")

# Use the estimated parameters to simulate new traits!
simul2<-mvSIM(tree,nsim=nsim,param=model_fit)

# or try with generic "simulate" function
simul3<-simulate(model_fit,nsim=nsim,tree=tree)

## Just-for-fun :Comparing parameters
\donttest{
 simul4<-simulate(model_fit,nsim=100,tree=tree)
 
 results<-lapply(simul4,function(x){
    mvBM(tree,x,model="BM1",method="pic", echo=F,diagnostic=F)
    })

 sigma_simul<-sapply(results,function(x){x$sigma})

# comparison between the simulated (black) and the observed (red) multivariate rates
 layout(matrix(1:4, ncol=2))
 for(i in 1:4){
	  hist(sigma_simul[i,], main=paste("Estimated sigma on simulated traits"),
  	xlab="estimated sigma for 100 replicates");abline(v=mean(sigma_simul[i,]),lwd=2);
  	abline(v=model_fit$sigma[i],col="red",lwd=2)
 }
	}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ simulate traits }
\keyword{ mvmorph object }% __ONLY ONE__ keyword per line
