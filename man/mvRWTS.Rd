\name{mvRWTS}
\alias{mvRWTS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multivariate Brownian motion / Random Walk models of continuous traits evolution on time-serie
%%  ~~function to do ... ~~
}
\description{
This function allows the fitting of multivariate rates of evolution under a Brownian motion model. This function can also fit constrained models.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mvRWTS(times, data, error = NULL, param =
    list(sigma=NULL, trend=FALSE, decomp="cholesky"), method = c("rpf",
    "inverse", "pseudoinverse"), scale.height = FALSE,
    optimization = c("L-BFGS-B", "Nelder-Mead", "subplex"),
    control = list(maxit = 20000), precalc = NULL, diagnostic = TRUE,
    echo = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{times}{
Phylogenetic tree in SIMMAP format by default. A "phylo" object can also be used with the "BM1" model.
%%     ~~Describe \code{tree} here~~
}
  \item{data}{
Matrix or data frame with species in rows and continuous traits in columns (preferentially with names and in the same order than in the tree). NA values are allowed with the "rpf", "inverse", and "pseudoinverse" methods.
%%     ~~Describe \code{data} here~~
}
  \item{error}{
Matrix or data frame with species in rows and continuous traits standard error (squared) in columns.
%%     ~~Describe \code{error} here~~
}

  \item{param}{
List of arguments to be passed to the function. See details.
%%     ~~Describe \code{param} here~~
}
  \item{method}{
Choose between "rpf", "inverse", or "pseudoinverse" for log-likelihood computation during the fitting process. See details.
%%     ~~Describe \code{method} here~~
}
  \item{scale.height}{
Whether the time-serie should be scaled to length 1 or not.
%%     ~~Describe \code{scale.height} here~~
}
  \item{optimization}{
Methods used by the optimization routines (see ?optim and ?subplex for details). The "fixed" method return the log-likelihood function only.
%%     ~~Describe \code{optimization} here~~
}
  \item{control}{
Max. bound for the number of iteration of the optimizer; other options can be fixed on the list (see ?optim or ?subplex).
%%     ~~Describe \code{control} here~~
}
  \item{precalc}{
Optional. precalculation of fixed parameters. See ?mvmorph.Precalc.
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
The mvRWTS function fits a multivariate Random Walk (RW) (the time serie counterpart of the Brownian motion) process.

The "method" argument allows the user to try different algorithms for computing the log-likelihood. The "rpf" and "sparse" methods use fast GLS algorithms based on factorization for avoiding the computation of the inverse of the variance-covariance matrix and its determinant involved in the log-likelihood estimation. The "inverse" approach uses the "stable" standard explicit computation of the inverse and determinant of the matrix and is therefore slower. The "pseudoinverse" method uses a generalized inverse that is safer for matrix near singularity but very time consuming. See ?mvLL for more details on these computational methods.

The \bold{ "param"} emph{list} arguments:

\bold{ "constraint"} - The "constraint" argument in the "param" list allows the user to compute the joint likelihood for each trait by assuming they evolved independently (\bold{ constraint="diagonal"}, or \bold{ constraint="equaldiagonal"}). If \bold{ constraint="equal"}, the model forces the sigma values for each studied trait to be the same using the constrained Cholesky decomposition proposed by Adams (2013) or a separation strategy based on spherical parameterization when p>2 (Clavel et al. 2015), which allow considering the covariations between traits.

Finally user-defined constrained models can be specified through a numeric matrix (square and symmetric) with integer values taken as indices of the parameters. For instance, for three traits: constraint=matrix(c(1,3,3,3,2,3,3,3,2),3). Covariances constrained to be zero are introduced by NA values, e.g., constraint=matrix(c(1,4,4,4,2,NA,4,NA,3),3).
Significance of nested model fit can be assessed using the "LRT" function. See example below and ?LRT.

\bold{ "decomp"} - For the general case (unconstrained models), the sigma matrix is parameterized by various methods to ensure it's positive definiteness (Pinheiro and Bates, 1996). These methods are the: "cholesky", "eigen+", and "spherical" parameterization.

\bold{ "trend"} - Default set to FALSE. If TRUE, the ancestral state is allowed to drift leading to a directional random walk. Note that it is possible to provide a vector of integer indices to constraint the the estimated trends when p>1 (see the vignette).

\bold{ "sigma"} - Starting values for the likelihood estimation. By default the traits covariances are used as starting values for the likelihood optimization. The user can specify starting values as square symmetric matrices or a simple vector of values for the upper factor of the sigma matrix. The parameterization is done using the factorization determined through the "decomp" argument (Pinheiro and Bates, 1996). Thus, you should provide p*(p+1)/2 values, with p the number of traits (e.g., random numbers or the values from the cholesky factor of a symmetric positive definite sigma matrix - see example below). If a constrained model is used, the number of starting values is (p*(p-1)/2)+1.


%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{LogLik }{The log-likelihood of the optimal model.}
\item{AIC }{Akaike Information Criterion for the optimal model.}
\item{AICc }{Sample size-corrected AIC.}
\item{theta }{Estimated ancestral states.}
\item{sigma }{Evolutionary rates matrix for each selective regime.}
\item{convergence }{Convergence status of the optimizing function; "0" indicates convergence (See ?optim for details).}
\item{hess.values }{Reliability of the likelihood estimates calculated through the eigen-decomposition of the hessian matrix. "0" means that a reliable estimate has been reached. (See ?mvOU).}
\item{param }{List of model fit parameters (optimization, method, model, number of parameters...).}
\item{llik }{The log-likelihood function evaluated in the model fit "$llik(par, root.mle=TRUE)".}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Adams D.C. 2013. Comparing evolutionary rates for different phenotypic traits on a phylogeny using likelihood. Syst. Biol. 62:181-192.

Clavel J., Escarguel G., Merceron G. 2015. mvMORPH: an R package for fitting multivariate evolutionary models to morphometric data. Methods Ecol. Evol., 6(11):1311-1319.

Hunt G. (2012). Measuring rates of phenotypic evolution and the inseparability of tempo and mode. Paleobiology, 38(3):351-373.

Revell L.J. 2012. phytools: An R package for phylogenetic comparative biology (and other things). Methods Ecol. Evol. 3:217-223.

%% ~put references to the literature/web site here ~
}
\author{
Julien Clavel
%%  ~~who you are~~
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{mvMORPH}}
\code{\link{mvOU}}
\code{\link{mvEB}}
\code{\link{mvSHIFT}}
\code{\link{mvSIM}}
\code{\link{mvOUTS}}
\code{\link{LRT}}
\code{\link{optim}}
\code{\link{brownie.lite}}
\code{\link{evol.vcv}}
\code{\link{make.simmap}}
\code{\link{make.era.map}}
\code{\link{paintSubTree}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## Toy example

set.seed(1)
# Simulate the time serie
timeserie <- 0:49

# Simulate the traits
sigma <- matrix(c(0.01,0.005,0.005,0.01),2)
theta <- c(0,1)
error <- matrix(0,ncol=2,nrow=50);error[1,]=0.001
data<-mvSIM(timeserie, error=error, 
            param=list(sigma=sigma, theta=theta), model="RWTS", nsim=1)

# plot the time serie
matplot(data, type="o", pch=1, xlab="Time (relative)")

# model fit
mvRWTS(timeserie, data, error=error, param=list(decomp="diagonal"))
mvRWTS(timeserie, data, error=error, param=list(decomp="equal"))
mvRWTS(timeserie, data, error=error, param=list(decomp="cholesky"))

# Random walk with trend
set.seed(1)
trend <- c(0.02,0.02)
data<-mvSIM(timeserie, error=error, 
            param=list(sigma=sigma, theta=theta, trend=trend), model="RWTS", nsim=1)

# plot the time serie
matplot(data, type="o", pch=1, xlab="Time (relative)")

# model fit
mvRWTS(timeserie, data, error=error, param=list(trend=TRUE))

# we can specify a vector of indices
mvRWTS(timeserie, data, error=error, param=list(trend=c(1,1)))

}% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Brownian Motion }
\keyword{ Random walk }
\keyword{ Time serie }
\keyword{ Evolutionary rates }
\keyword{ Cholesky constraint }% __ONLY ONE__ keyword per line
