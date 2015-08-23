\name{mvOU}
\alias{mvOU}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multivariate Ornstein-Uhlenbeck model of continuous traits evolution
%%  ~~function to do ... ~~
}
\description{
This function allows the fitting of a multivariate Ornstein-Uhlenbeck (OU) model by allowing a given tree branch to be subdivided into multiple selective regimes using SIMMAP-like mapping of ancestral states.
Species measurement errors or dispersions can also be included in the model.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mvOU(tree, data, error = NULL, model = c("OUM", "OU1"), param = list(sigma = NULL,
    alpha = NULL, vcv = "mvmorph", decomp = c("symmetricPositive", "symmetric",
    "nsymPositive", "nsymmetric", "diagonal", "upper", "lower")), method = c("rpf", 
    "sparse", "inverse", "pseudoinverse", "univarpf"), scale.height = FALSE, 
    optimization = c("L-BFGS-B", "Nelder-Mead", "subplex"), control = 
    list(maxit = 20000), precalc = NULL, diagnostic = TRUE, echo = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree}{
Phylogenetic tree with mapped ancestral states in SIMMAP format.
(See make.simmap function from phytools package). A "phylo" object can be used with model "OU1".
%%     ~~Describe \code{tree} here~~
}
  \item{data}{
Matrix or data frame with species in rows and continuous traits in columns. NA values are allowed with the "rpf", "inverse", and "pseudoinverse" methods.
%%     ~~Describe \code{data} here~~
}
  \item{error}{
Matrix or data frame with species in rows and continuous trait standard errors (squared) in columns.
%%     ~~Describe \code{error} here~~
}

  \item{model}{
Choose between "OUM" for a multiple selective regime model, or "OU1" for a unique selective regime for the whole tree.
%%     ~~Describe \code{model} here~~
}
  \item{param}{
List of arguments to be passed to the function. See details.
%%     ~~Describe \code{param} here~~
}
  \item{method}{
Choose between "rpf", "sparse", "inverse", "pseudoinverse", or "univarpf" for computing the log-likelihood during the fitting process. See details.
%%     ~~Describe \code{method} here~~
}
  \item{scale.height}{
Whether the tree should be scaled to length 1 or not.
%%     ~~Describe \code{scale.height} here~~
}
  \item{optimization}{
Methods used by the optimization routines. (See ?optim and ?subplex for details).
%%     ~~Describe \code{optimization} here~~
}
  \item{control}{
Max. bound for the number of iteration of the optimizer; other options can be fixed on the list. (See ?optim or ?subplex).
%%     ~~Describe \code{control} here~~
}
  \item{precalc}{
Optional. precalculation of fixed parameters. See ?mvmorph.Precalc for details.
%%     ~~Describe \code{precalc} here~~
}
  \item{diagnostic}{
 Whether the convergence diagnostics should be returned or not.
%%     ~~Describe \code{diagnostic} here~~
}
  \item{echo}{
Whether the results must be returned or not.
%%     ~~Describe \code{echo} here~~
}
}
\details{
The mvOU function fit a multivariate model of evolution according to an Ornstein-Uhlenbeck process. The user can incorporate measurement errors and uses SIMMAP-like mapping of ancestral states. SIMMAP mapping allows one to assign parts of branchs to different selective regimes, and allows testing for change in trait variance that is not synchronous with the species divergence events.

Mapping of ancestral states can be done using the "make.simmap", "make.era.map" or "paintSubTree" functions from the "phytools" package.

The "method" argument allows the user to try different algorithms for computing the log-likelihood. The "rpf", "univarpf" (for univariate analysis) and "sparse" methods use fast GLS algorithms based on factorization for avoiding the computation of the inverse of the variance-covariance matrix and its determinant for the log-likelihood estimation. The "inverse" approach use the "stable" standard explicit computation of the inverse and determinant of the matrix and is therefore slower. The "pseudoinverse" method uses a generalized inverse that is safer for matrix near singularity but very time consuming. See ?mvLL for details.

Arguments in the "param" list are:

"sigma" or "alpha" - Starting values for the likelihood search can be specified through the "alpha" and "sigma" arguments in the param list. It is also possible to test for the significance of the off-diagonal sigma and alpha matrix in the full model by making comparison with a constrained model (using sigma="constraint", or alpha="constraint") in the "param" argument list. You can also provide starting values for the constrained model. For instance, for two traits use sigma=list("constraint", c(0.5,0.5))  (or alpha=list("constraint", c(0.5,0.5))).

"decomp" - You can further constraint the alpha matrix by specifying the decomposition of the matrix through the "decomp" argument in the "param" list. Indeed, the multivariate Ornstein-Uhlenbeck model is described by the spectral decomposition of the alpha matrix. Thus it is possible to parameterize the alpha matrix on its eigenvalues with different biological interpretations (e.g., Sy et al. 1997, Bartoszek et al. 2012). "symmetric" uses SVD parameterization of the matrix, "symmetricPositive" uses a Cholesky parameterization of the matrix, "nsymmetric" and "nsymPositive" use a Schur parameterization of the matrix for specifying the eigenvalues. It is also possible to specify "diagonal" which is similar to the use of the "constraint" argument for "alpha" argument. Finally, one can specify that the alpha matrix is "upper" or "lower" triangular (i.e., one process affect the other unilateraly). The sigma matrix is parameterized by the Cholesky method to ensure positive definiteness (Pinheiro and Bates, 1996).

"vcv" - It is possible to specify in the "param" list what kind of variance-covariance matrix to use with the "vcv" argument. The vcv="ouch" option uses a modified version of the OUCH package code (Butler and King, 2004) to compute the variance-covariance matrix. This algorithm assumes that the process is stationary. It cannot be used with the "sparse" method to speed up the computations. The vcv="mvmorph" option follows Bartoszek et al. 2012 for computing the variance covariance matrix. This approach can be used with non-ultrametric trees, with the various matrix parameterization proposed and all the GLS methods of mvMORPH. Both approaches should converge on the same results when the traits have evolved according to an OU process.

"root" - This argument allows the user to specify if the ancestral state at the root (theta 0) should be estimated (root=TRUE), assumed to be at the oldest regime state (root=FALSE), or that there is no root and each regimes is at the stationary point (root="stationary"). The first option should be used with non-ultrametric trees (i.e., with fossil species - e.g., Hansen 1997) where information on the ancestral state is directly available from the data. Indeed, estimating shifts in the ancestral state from extant species could be problematic and it seems preferable to either assume each regime to be stationary or that the root state is the state of the older regime (as currently assumed in the "ouch" package).


Note on the returned Hessian matrix in the result list:

The hessian is the matrix of second order partial derivatives of the likelihood function with respect to the maximum likelihood parameter values. This matrix provides a measure of the steepness of the likelihood surface in the region of the optimum. The eigen-decomposition of the hessian matrix returned by the optimizing function allows assessing the reliability of the fit of the model (even if the optimizer have converged).
When the optimization function does not converge on a stable result, the user may consider increasing the "maxit" argument in the "control" option, or try a simpler model with fewer parameters to estimate. Changing the starting values ("alpha" and "sigma" options in the param list) as well as the optimizing function ("optimization" options) may help sometimes (e.g., alpha=runif(3) for two traits analysis with random starting values (i.e., the lower triangular alpha matrix)). Note that the number of starting values to provide depends on the matrix decomposition chosen for the alpha parameter (p*(p+1)/2 values for symmetric alpha matrix but p*p values for non-symmetric ones - with p the number of traits).

Note that Bartoszek et al. (2012) proposed the mvSLOUCH package dedicated to multivariate Ornstein-Uhlenbeck processes, which allows fitting models with a randomly evolving predictor variables. They also provide a detailed mathematical description of the multivariate OU model which the reader should consult.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{LogLik }{The log-likelihood of the optimal model.}
\item{AIC }{Akaike Information Criterion for the optimal model.}
\item{AICc }{Sample size-corrected AIC.}
\item{theta }{Estimated ancestral states.}
\item{alpha }{Matrix of estimated alpha values (strength of selection).}
\item{sigma }{Evolutionary rates matrix (drift).}
\item{convergence }{Convergence status of the optimizing function; "0" indicates convergence. (See ?optim for details).}
\item{hessian }{Hessian matrix of second order partial derivatives at the MLE. (See details).}
\item{hess.values }{Reliability of the likelihood estimates calculated through the eigen-decomposition of the hessian matrix. "0" means that a reliable estimate has been reached. See details.}
\item{param }{List of model fit parameters (optimization, method, model, number of parameters...).}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Bartoszek K., Pienaar J., Mostad P., Andersson S., Hansen T.F. 2012. A phylogenetic comparative method for studying multivariate adaptation. J. Theor. Biol. 314:204-215.

Beaulieu J.M., Jhwueng D.-C., Boettiger C., O'Meara B.C. 2012. Modeling stabilizing selection: Expanding the Ornstein-Uhlenbeck model of adaptive evolution. Evolution. 66:2369-2389.

Butler M.A., King A.A. 2004. Phylogenetic comparative analysis: a modeling approach for adaptive evolution. Am. Nat. 164:683-695.

Hansen T.F. 1997. Stabilizing selection and the comparative analysis of adaptation. Evolution. 51:1341-1351.

Pinheiro J.C., Bates D.M. 1996. Unconstrained parameterizations for variance-covariance matrices. Stat. Comput. 6:289-296.

Sy J.P., Taylor J.M.G., Cumberland W.G. 1997. A stochastic model for the analysis of bivariate longitudinal AIDS data. Biometrics. 53:542-555.
%% ~put references to the literature/web site here ~
}
\author{
Julien Clavel
%%  ~~who you are~~
}
\note{
This function partly uses a modified version of the C code from the "OUCH" package built by Aaron King, as well as a C code which is part of the "ape" package by Emmanuel Paradis.  I kindly thank those authors for sharing their sources.
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{mvMORPH}}
\code{\link{halflife}}
\code{\link{stationary}}
\code{\link{mvBM}}
\code{\link{mvEB}}
\code{\link{mvSHIFT}}
\code{\link{mvSIM}}
\code{\link{LRT}}
\code{\link{optim}}
\code{\link{make.simmap}}
\code{\link{make.era.map}}
\code{\link{paintSubTree}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## Toy example

# Simulated dataset
set.seed(14)
# Generating a random tree
tree<-pbtree(n=50)

# Setting the regime states of tip species
sta<-as.vector(c(rep("Forest",20),rep("Savannah",30))); names(sta)<-tree$tip.label

# Making the simmap tree with mapped states
tree<-make.simmap(tree,sta , model="ER", nsim=1)
col<-c("blue","orange"); names(col)<-c("Forest","Savannah")

# Plot of the phylogeny for illustration
plotSimmap(tree,col,fsize=0.6,node.numbers=FALSE,lwd=3, pts=FALSE)

# Simulate the traits
alpha<-matrix(c(2,0.5,0.5,1),2)
sigma<-matrix(c(0.1,0.05,0.05,0.1),2)
theta<-c(2,3,1,1.3)
data<-mvSIM(tree, param=list(sigma=sigma, alpha=alpha, ntraits=2, mu=theta,
            names_traits=c("head.size","mouth.size")), model="OUM", nsim=1)

## Fitting the models
# OUM - Analysis with multiple optima
mvOU(tree, data)

# OU1 - Analysis with a unique optimum
# mvOU(tree, data, model="OU1", method="sparse")

# various options
# mvOU(tree, data, model="OUM", method="sparse", scale.height=FALSE,
#    param=list(decomp="nsymmetric", root="stationary"))
# mvOU(tree, data, model="OUM", method="sparse", scale.height=FALSE,
#    param=list(decomp="nsymmetric", root=TRUE))
# mvOU(tree, data, model="OUM", method="sparse", scale.height=FALSE,
#    param=list(decomp="symmetricPositive", root=TRUE))
# OUCH setting
# mvOU(tree, data, model="OUM", method="rpf", scale.height=FALSE,
#    param=list(decomp="symmetricPositive", root=FALSE, vcv="ouch"))

## Univariate case - FAST with RPF
# set.seed(14)
# tree<-pbtree(n=500)

# Setting the regime states of tip species
# sta<-as.vector(c(rep("Forest",200),rep("Savannah",300))); names(sta)<-tree$tip.label

# Making the simmap tree with mapped states
# tree<-make.simmap(tree,sta , model="ER", nsim=1)
# col<-c("blue","orange"); names(col)<-c("Forest","Savannah")

# Plot of the phylogeny for illustration
# plotSimmap(tree,col,fsize=0.6,node.numbers=FALSE,lwd=3, pts=FALSE)

# Parameters
# alpha<-2.5
# sigma<-0.1
# theta<-c(0,2)
# data<-mvSIM(tree, param=list(sigma=sigma, alpha=alpha, ntraits=1, mu=theta,
#             names_traits=c("body_size")), model="OUM", nsim=1)

# Fit the model
# system.time(mvOU(tree, data, model="OUM", method="univarpf",
#                param=list(root="stationary")))
# system.time(mvOU(tree, data, model="OU1", method="univarpf",
#                param=list(root="stationary")))

# Add measurement error
# error=rnorm(500,sd=0.1)
# mvOU(tree, data, error=error^2, model="OUM", method="univarpf",
#    param=list(root="stationary"))


}% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Ornstein Uhlenbeck }
\keyword{ Hessian }
\keyword{ OU }% __ONLY ONE__ keyword per line
