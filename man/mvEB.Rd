\name{mvEB}
\alias{mvEB}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multivariate Early Burst model of continuous traits evolution
%%  ~~function to do ... ~~
}
\description{
This function fits multivariate Early Burst (EB) or ACDC models of evolution.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mvEB(tree, data, error = NULL, param = list(up = 0), method = 
	c("rpf", "sparse", "inverse", "pseudoinverse", "pic"), scale.height =
	 FALSE, optimization = c("Nelder-Mead", "L-BFGS-B", "subplex"), 
	 control = list(maxit = 20000), precalc = NULL, diagnostic = TRUE, 
	 echo = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree}{
Phylogenetic tree (phylo object).
%%     ~~Describe \code{tree} here~~
}
  \item{data}{
Matrix or data frame with species in rows and continuous traits in columns (preferentially with names and in the same order than in the tree).
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
  \item{method}{
Choose between "rpf", "sparse", "inverse", "pseudoinverse", or "pic" for computing the log-likelihood during the fitting process. See details.
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
Max. bound for the number of iteration of the optimizer; other options can be fixed on the list (See ?optim or ?subplex for details).
%%     ~~Describe \code{control} here~~
}
  \item{precalc}{
Optional. precalculation of fixed parameters. See ?mvmorph.Precalc for details.
%%     ~~Describe \code{precalc} here~~
}
  \item{diagnostic}{
Whether the diagnostics of convergence should be returned.
%%     ~~Describe \code{diagnostic} here~~
}
  \item{echo}{
Whether the results must be returned or not.
%%     ~~Describe \code{echo} here~~
}
}
\details{
The Early Burst model (Harmon et al. 2010) is a special case of the ACDC model of Blomberg et al. (2003). Using an upper bound larger than zero transform the EB model to the accelerating rates of character evolution of Blomberg et al. (2003).

The "method" argument allows the user to try different algorithms for computing the log-likelihood. The "rpf" and "sparse" methods use fast GLS algorithms based on factorization for avoiding the computation of the inverse of the variance-covariance matrix and its determinant for the log-likelihood estimation. The "inverse" approach uses the "stable" standard explicit computation of the inverse and determinant of the matrix and is therefore slower. The "pseudoinverse" method uses a generalized inverse that is safer for matrix near singularity but very time consuming. The "pic" method uses a very fast algorithm based on independent contrasts. It should be used with strictly dichotomic trees (i.e., no polytomies) and is currently not available for the multivariate "BMM" model. See ?mvLL for more details on these computational methods.


The "param" list can be used to set the lower (low) and upper (up, default value is 0 - i.e., Early Burst model) bounds for the estimation of the exponential rate (beta). The default lower bound for decelerating rates (as assumed in Early Burst) is fixed as log(min.rate) / T, where T is the depth of the tree and min.rate is the minimum rate that could be assumed for the model (following Slater and Pennell, 2014; log(10^-5)/T). Bounds may need to be adjusted by hand for specific cases.

Starting values for "sigma" and "beta" could also be provided through the "param" list.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{LogLik }{The log-likelihood of the optimal model.}
\item{AIC }{Akaike Information Criterion for the optimal model.}
\item{AICc }{Sample size-corrected AIC.}
\item{theta }{Estimated ancestral states.}
\item{beta }{Exponent rate (of decay or increase).}
\item{sigma }{Evolutionary rates matrix for each selective regimes.}
\item{convergence }{Convergence status of the optimizing function; "0" indicates convergence (See ?optim for details).}
\item{hessian }{Hessian matrix of second order partial derivatives at the MLE. (See ?mvOU for details).}
\item{hess.values }{Reliability of the likelihood estimates calculated through the eigen-decomposition of the hessian matrix. "0" means that a reliable estimate has been reached. (See ?mvOU for details).}
\item{param }{List of model fit parameters (optimization, method, model, number of parameters...).}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Blomberg S.P., Garland T.J., Ives A.R. 2003. Testing for phylogenetic signal in comparative data: behavioral traits are more labile. Evolution. 57:717-745.

Harmon L.J., Losos J.B., Davies J.T., Gillespie R.G., Gittleman J.L., Jennings B.W., Kozak K.H., McPeek M.A., Moreno-Roark F., Near T.J., Purvis A., Ricklefs R.E., Schluter D., Schulte II J.A., Seehausen O., Sidlauskas B.L., Torres-Carvajal O., Weir J.T., Mooers A.O. 2010. Early bursts of body size and shape evolution are rare in comparative data. Evolution. 64:2385-2396.

Slater G.J., Pennell M. 2014. Robust regression and posterior predictive simulation increase power to detect early bursts of trait evolution. Syst. Biol. 63: 293-308.
%% ~put references to the literature/web site here ~
}
\author{
Julien Clavel
%%  ~~who you are~~
}
\note{
The derivative-free "Nelder-Mead" optimization method is used as default setting instead of "L-BFGS-B".
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{mvMORPH}}
\code{\link{mvOU}}
\code{\link{mvBM}}
\code{\link{mvSHIFT}}
\code{\link{mvSIM}}
\code{\link{optim}}

%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## Toy example

# Simulated dataset
set.seed(14)
# Generating a random tree
tree<-pbtree(n=50, scale=10)

# Simulate the traits
sigma<-matrix(c(0.1,0.05,0.05,0.1),2)
theta<-c(0,0)
beta<- -0.34 # 5 phylogenetic half-life ( log(2)/ (10/5) )
data<-mvSIM(tree, param=list(sigma=sigma, beta=beta, ntraits=2, mu=theta,
            names_traits=c("head.size","mouth.size")), model="EB", nsim=1)

## Fitting the models
mvEB(tree, data)
mvEB(tree, data, method="pic")
mvEB(tree, data, method="pic", param=list(low=log(10^-5)/10)) # avoid internal estimation

# ACDC
beta<- 0.34
data<-mvSIM(tree, param=list(sigma=sigma, beta=beta, ntraits=2, mu=theta,
            names_traits=c("head.size","mouth.size")), model="EB", nsim=1)

fit<-mvEB(tree, data, method="pic", param=list(up=2, low=-2))

logLik(fit)
AIC(fit)
summary(fit)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Early burst }
\keyword{ Decelerating rates }
\keyword{ Accelerating rates }
% __ONLY ONE__ keyword per line
