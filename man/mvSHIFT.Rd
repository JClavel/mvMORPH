\name{mvSHIFT}
\alias{mvSHIFT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multivariate change in mode of continuous traits evolution
%%  ~~function to do ... ~~
}
\description{
This function fits different models of evolution after a fixed point. This allows fitting models of change in mode of evolution following a given event.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mvSHIFT(tree, data, error = NULL, param = list(age = NULL, sigma = NULL,
        alpha = NULL, sig = NULL, beta = NULL), model = c("ER", "RR", "EC",
        "RC", "SR", "EBOU", "OUEB", "EBBM", "BMEB"), method = c("rpf",
        "sparse", "inverse", "pseudoinverse"), scale.height = FALSE,
        optimization = c("L-BFGS-B", "Nelder-Mead", "subplex"), control =
        list(maxit = 20000), precalc = NULL, diagnostic = TRUE, echo = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree}{
 Phylogenetic tree with shift mapped. (See "make.era.map" function from "phytools" package). A "phylo" object can be used if the "age" argument is provided in the "param" list.
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
  \item{param}{
List of arguments to be passed to the function. See details.
%%     ~~Describe \code{param} here~~
}
  \item{model}{
Choose between the different models "OUBM", "BMOU", "EBOU", "OUEB", "BMEB", "EBBM"... See details below.
%%     ~~Describe \code{model} here~~
}
  \item{method}{
Choose between "rpf", "sparse", "inverse", or "pseudoinverse" for computing the log-likelihood during the fitting process. See details.
%%     ~~Describe \code{method} here~~
}
  \item{scale.height}{
Whether the tree should be scaled to length 1 or not.
%%     ~~Describe \code{scale.height} here~~
}
  \item{optimization}{
Methods used by the optimization routines. (See ?optim and ?subplex for details). The "fixed" method return the log-likelihood function only.
%%     ~~Describe \code{optimization} here~~
}
  \item{control}{
Max. bound for the number of iteration of the optimizer; other options can be fixed on the list (See ?optim and ?subplex for details).
%%     ~~Describe \code{control} here~~
}
  \item{precalc}{
Optional. precalculation of fixed parameters. See ?mvmorph.Precalc for details.
%%     ~~Describe \code{precalc} here~~
}
  \item{diagnostic}{
Whether the diagnostics of convergence should be returned or not.
%%     ~~Describe \code{diagnostic} here~~
}
  \item{echo}{
Whether the results must be returned or not.
%%     ~~Describe \code{echo} here~~
}
}
\details{
The mvSHIFT function fits a shift in mode or rate of evolution at a fixed point in time, as previously proposed by some authors (O'Meara et al. 2006; O'Meara, 2012; Slater, 2013). Shift in mode of evolution could be mapped on a modified "phylo" object using the "make.era.map" function from the "phytools" package.
Note that only one shift is allowed by the current version. The age of the shift could be otherwise directly provided (in unit of times of the tree) in the function by the "age" argument in the "param" list.

The function allows fitting model with shift from an Orstein-Uhlenbeck to a Brownian motion process and vice-versa ("OUBM" and "BMOU"), shifts from a Brownian motion to/from an Early Burst (ACDC) model ("BMEB" and "EBBM"), or shifts from an Orstein-Uhlenbeck to/from an Early Burst (ACDC) model ("OUEB" and "EBOU"). (-- The shift models with OU process are relevant only if you use fossil species!! --)

In all these cases it is possible to allow the drift parameter to vary after the fixed point by specifying "i" (for independent) after the model name. For instance, to fit models of "ecological release" or "ecological release and radiate" following Slater (2013), we can use  "OUBM" or "OUBMi" respectively.

Alternatively it is also possible to use the shortcuts "ER" or "RR" to fit respectively models of "ecological release" and "ecological release and radiate", and "EC" for a model of "constrained ecology" (e.g., after invasion of a competitive species in a given ecosystem) where traits are constrained in an Ornstein-Uhlenbeck process after a fixed point in time ("RC" is the same model but assumes an independent rate during the early radiative phase). The "SR" model allows fitting different (Brownian) rates/drift before and after the shift point (note that this model could also be fitted using the mvBM function).

The "param" list can be used to provide lower and upper bounds for the exponential rate parameter of the Early-Burst/ACDC model. See ?mvEB for details.

The "method" argument allows the user to try different algorithms for computing the log-likelihood. The "rpf" and "sparse" methods use fast GLS algorithms based on factorization for avoiding the computation of the inverse of the variance-covariance matrix and its determinant involved in the log-likelihood estimation. The "inverse" approach uses the "stable" standard explicit computation of the inverse and determinant of the matrix and is therefore slower. The "pseudoinverse" method uses a generalized inverse that is safer for matrix near singularity but very time consuming. See ?mvLL for details.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{LogLik }{The log-likelihood of the optimal model.}
\item{AIC }{Akaike Information Criterion for the optimal model.}
\item{AICc }{Sample size-corrected AIC.}
\item{theta }{Estimated ancestral states.}
\item{alpha }{Matrix of estimated alpha values (strength of selection).}
\item{beta }{Exponent rate (of decay or increase) for the ACDC/Early-Burst model.}
\item{sigma }{Evolutionary rates matrix (drift) for the BM process before the shift.}
\item{sig   }{Evolutionary rates matrix (drift) for the BM process after the shift (only for "i" models).}
\item{convergence }{Convergence status of the optimizing function; "0" indicates convergence. (See ?optim for details).}
\item{hess.values }{Reliability of the likelihood estimates calculated through the eigen-decomposition of the hessian matrix. "0" means that a reliable estimate has been reached. (See ?mvOU for details).}
\item{param }{List of model fit parameters (optimization, method, model, number of parameters...).}
\item{llik }{The log-likelihood function evaluated in the model fit "$llik(par, root.mle=TRUE)".}

%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Clavel J., Escarguel G., Merceron G. 2015. mvMORPH: an R package for fitting multivariate evolutionary models to morphometric data. Methods in Ecology and Evolution, 6(11):1311-1319.

O'Meara B.C. 2012. Evolutionary inferences from phylogenies: a review of methods. Annu. Rev. Ecol. Evol. Syst. 43:267-285.

O'Meara B.C., Ane C., Sanderson M.J., Wainwright P.C. 2006. Testing for different rates of continuous trait evolution. Evolution. 60:922-933.

Slater G.J. 2013. Phylogenetic evidence for a shift in the mode of mammalian body size evolution at the Cretaceous-Palaeogene boundary. Methods Ecol. Evol. 4:734-744.
%% ~put references to the literature/web site here ~
}
\author{
Julien Clavel
%%  ~~who you are~~
}
\note{
Changes in rate of evolution and optima can also be fitted using the mvBM and mvOU functions using a 'make.era.map' transformed tree.
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{mvMORPH}}
\code{\link{mvOU}}
\code{\link{mvBM}}
\code{\link{mvEB}}
\code{\link{mvOUTS}}
\code{\link{mvRWTS}}
\code{\link{mvSIM}}
\code{\link{optim}}
\code{\link{subplex}}
\code{\link{paintSubTree}}
\code{\link{make.era.map}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## Toy example

# Simulated dataset
set.seed(14)
# Generating a random tree
tree<-rtree(50)

# Providing a tree whith the shift mapped on
tot<-max(nodeHeights(tree))
age=tot-3    # The shift occured 3 Ma ago
tree<-make.era.map(tree,c(0,age))

# Plot of the phylogeny for illustration
plotSimmap(tree,fsize=0.6,node.numbers=FALSE,lwd=3, pts=FALSE)

# Simulate the traits
alpha<-matrix(c(2,0.5,0.5,1),2)
sigma<-matrix(c(0.1,0.05,0.05,0.1),2)
theta<-c(2,3)
data<-mvSIM(tree, param=list(sigma=sigma, alpha=alpha, ntraits=2, theta=theta,
            names_traits=c("head.size","mouth.size")), model="ER", nsim=1)


## Fitting the models
# "Ecological release model"
mvSHIFT(tree, data, model="OUBM") # similar to mvSHIFT(tree, data, model="ER")

# "Release and radiate model"
\donttest{
 mvSHIFT(tree, data, model="RR", method="sparse")
# similar to mvSHIFT(tree, data, model="OUBMi")

# More generally...

# OU to/from BM
 mvSHIFT(tree, data, model="OUBM", method="sparse")
 mvSHIFT(tree, data, model="BMOU", method="sparse")
 mvSHIFT(tree, data, model="OUBMi", method="sparse")
 mvSHIFT(tree, data, model="BMOUi", method="sparse")

# BM to/from EB
 mvSHIFT(tree, data, model="BMEB", method="sparse")
 mvSHIFT(tree, data, model="EBBM", method="sparse")
 mvSHIFT(tree, data, model="BMEBi", method="sparse")
 mvSHIFT(tree, data, model="EBBMi", method="sparse")

# OU to/from EB
 mvSHIFT(tree, data, model="OUEB", method="sparse")
 mvSHIFT(tree, data, model="OUEBi", method="sparse")
 mvSHIFT(tree, data, model="EBOU", method="sparse")
 mvSHIFT(tree, data, model="EBOUi", method="sparse")


## Without providing mapped tree
# The shift occured 3Ma ago (param$age=3)
 set.seed(14)
 tree<-rtree(50)
 data<-mvSIM(tree, param=list(sigma=sigma, alpha=alpha, ntraits=2, theta=theta,
            names_traits=c("head.size","mouth.size"), age=3), model="ER", nsim=1)

## Fitting the models without mapped tree but by specifying the age in the param list.
 mvSHIFT(tree, data, model="OUBM", param=list(age=3))
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Ornstein Uhlenbeck }
\keyword{ Shifts }
\keyword{ Brownian Motion }
\keyword{ Evolutionary rates }
\keyword{ Early-Burst }
\keyword{ RR }
\keyword{ EC }
\keyword{ SR }
\keyword{ ER }% __ONLY ONE__ keyword per line
