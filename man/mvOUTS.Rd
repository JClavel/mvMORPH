\name{mvOUTS}
\alias{mvOUTS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multivariate continuous trait evolution for a stationary time series (Ornstein-Uhlenbeck model)
%%  ~~function to do ... ~~
}
\description{
This function allows the fitting of a multivariate Ornstein-Uhlenbeck (OU) model to a time series.
Species measurement errors or dispersions can also be included in the model.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mvOUTS(times, data, error = NULL, param = list(sigma = NULL, alpha = NULL,
    vcv = "randomRoot", decomp = c("cholesky","spherical","eigen","qr",
    "diagonal","upper","lower")), method = c("rpf", "inverse", "pseudoinverse",
    "univarpf"), scale.height = FALSE, optimization = c("L-BFGS-B", "Nelder-Mead", 
    "subplex"), control = list(maxit = 20000), precalc = NULL, diagnostic = TRUE, 
    echo = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{times}{
Time series - vector of sample ages.
%%     ~~Describe \code{tree} here~~
}
  \item{data}{
Matrix or data frame with species in rows and continuous traits in columns. NA values are allowed.
%%     ~~Describe \code{data} here~~
}
  \item{error}{
Matrix or data frame with species in rows and continuous trait sampling variance (squared standard errors) in columns.
%%     ~~Describe \code{error} here~~
}

  \item{param}{
List of arguments to be passed to the function. See details below.
%%     ~~Describe \code{param} here~~
}
  \item{method}{
Choose between "rpf", "inverse", "pseudoinverse", or "univarpf" for computing the log-likelihood during the fitting process. See details below.
%%     ~~Describe \code{method} here~~
}
  \item{scale.height}{
Whether the time series should be scaled to unit length or not.
%%     ~~Describe \code{scale.height} here~~
}
  \item{optimization}{
Methods used by the optimization routines (see ?optim and ?subplex for details). The "fixed" method returns the log-likelihood function only.
%%     ~~Describe \code{optimization} here~~
}
  \item{control}{
Max. bound for the number of iteration of the optimizer; other options can be fixed in the list (see ?optim or ?subplex).
%%     ~~Describe \code{control} here~~
}
  \item{precalc}{
Optional. Precalculation of fixed parameters. See ?mvmorph.Precalc for details.
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
The mvOUTS function fits a multivariate model of trait evolution on a time series according to an Ornstein-Uhlenbeck process. The user can include measurement errors to the analyzed dataset.

The "method" argument allows the user to try different algorithms for computing the log-likelihood. The \emph{"rpf"}, \emph{"univarpf"} (for univariate analysis) methods use fast GLS algorithms based on factorization for avoiding the computation of the inverse of the variance-covariance matrix and its determinant for the log-likelihood estimation. The \emph{"inverse"} approach uses the \emph{"stable"} standard explicit computation of the inverse and determinant of the matrix and is therefore slower. The \emph{"pseudoinverse"} method uses a generalized inverse that is safer for matrix near singularity but highly time consuming. See ?mvLL for details.

Arguments in the \bold{"param"} \emph{list} are:

\bold{"sigma"} or \bold{"alpha"} - Starting values for the likelihood search can be specified through the "alpha" and "sigma" arguments in the param list. It is also possible to test for the significance of the off-diagonal sigma (scatter) and alpha (drift) matrix in the full model by making comparison with a constrained model (using sigma="constraint", or alpha="constraint") in the "param" argument list. You can also provide starting values for the constrained model. For instance, for two traits use sigma=list("constraint", c(0.5,0.5))  (or alpha=list("constraint", c(0.5,0.5))).

\bold{"decomp"} - You can further constrain the alpha matrix by specifying the decomposition of the matrix through the "decomp" argument in the "param" list. Indeed, the multivariate Ornstein-Uhlenbeck model is described by the spectral decomposition of the alpha matrix. Thus it is possible to parameterize the alpha matrix to be decomposable using various parameterizations (e.g., on its eigenvalues with different biological interpretations; Sy et al. 1997, Bartoszek et al. 2012). For a symmetric matrix parameterization the user can choose the \emph{"cholesky"}, \emph{"eigen"}, or \emph{"spherical"} option. 
For general square (non-symmetric) matrices the \emph{"svd"}, \emph{"qr"} and \emph{"schur"} parameterizations can be used. The \emph{"schur"} parameterization constrains the eigenvalues of the alpha matrix to be real numbers. The \emph{"svd+"}, \emph{"qr+"} or \emph{"eigen+"} options forces the eigenvalues to be positives by taking their logarithm. It is also possible to specify \emph{"diagonal"} which is similar to the use of the "constraint" argument for the "alpha" argument, or to use \emph{"equal"} and \emph{"equaldiagonal"}. Finally, one can specify that the alpha matrix is \emph{"upper"} or \emph{"lower"} triangular (i.e., one process affect the other unilateraly). Details can be found in the package vignette: browseVignettes("mvMORPH").

\bold{"decompSigma"} - The sigma matrix is parameterized by various methods to ensure its positive definiteness (Pinheiro and Bates, 1996). These methods can be accessed through the "decompSigma" argument and are the \emph{"cholesky"}, \emph{"eigen+"}, and \emph{"spherical"} parameterization. The sigma matrix can also be forced to be diagonal using \emph{"diagonal"} or \emph{"equaldiagonal"} and forced to have the same variances using \emph{"equal"}. Details can be found in the package vignette: browseVignettes("mvMORPH").

\bold{"vcv"} - It is possible to specify in the "param" list what kind of variance-covariance matrix to use with the "vcv" argument, depending on how the root is treated. 
The \emph{vcv="randomRoot"} option assumes that the value at the root is a random variable with the stationary distribution of the process. The \emph{vcv="fixedRoot"} option assumes that the root is a fixed parameter.

\bold{"root"} - If root=TRUE, the ancestral state and the optimum (stationary mean) are estimated, otherwise (root=FALSE) the ancestral (initial) state and the optimum (long-term expectation) are assumed to be the same.

Note: for the \bold{"decomp"} and \bold{"decompSigma} arguments, an user-defined matrix with integer values taken as indices of the parameters to be estimated can be provided. See ?mvBM and ?mvRWTS.

%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{LogLik }{The log-likelihood of the optimal model.}
\item{AIC }{Akaike Information Criterion for the optimal model.}
\item{AICc }{Sample size-corrected AIC.}
\item{theta }{Estimated ancestral states.}
\item{alpha }{Matrix of estimated alpha values (strength of selection, drift matrix).}
\item{sigma }{Evolutionary rate matrix (scatter).}
\item{convergence }{Convergence status of the optimizing function; "0" indicates convergence. (See ?optim for details).}
\item{hess.values }{Reliability of the likelihood estimates calculated through the eigen-decomposition of the hessian matrix. "0" means that a reliable estimate has been reached. See details on ?mvOU.}
\item{param }{List of model fit parameters (optimization, method, model, number of parameters...).}
\item{llik }{The log-likelihood function evaluated in the model fit "$llik(par, root.mle=TRUE)".}

%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Bartoszek K., Pienaar J., Mostad P., Andersson S., Hansen T.F. 2012. A phylogenetic comparative method for studying multivariate adaptation. J. Theor. Biol. 314:204-215.

Clavel J., Escarguel G., Merceron G. 2015. mvMORPH: an R package for fitting multivariate evolutionary models to morphometric data. Methods Ecol. Evol. 6(11):1311-1319.

Hunt G., Bell M.A., Travis M.P. 2008. Evolution toward a new adaptive optimum: phenotypic evolution in a fossil stickleback lineage. Evolution 62(3):700-710.

Pinheiro J.C., Bates D.M. 1996. Unconstrained parameterizations for variance-covariance matrices. Stat. Comput. 6:289-296.

Sy J.P., Taylor J.M.G., Cumberland W.G. 1997. A stochastic model for the analysis of bivariate longitudinal AIDS data. Biometrics. 53:542-555.
%% ~put references to the literature/web site here ~
}
\author{
Julien Clavel
%%  ~~who you are~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{mvMORPH}}
\code{\link{halflife}}
\code{\link{stationary}}
\code{\link{mvOU}}
\code{\link{mvRWTS}}
\code{\link{mvBM}}
\code{\link{mvEB}}
\code{\link{mvSHIFT}}
\code{\link{mvSIM}}
\code{\link{LRT}}
\code{\link{optim}}

%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
\donttest{
# Simulate the time series
set.seed(14)
timeseries <- 0:49
# Parameters with general alpha matrix on two competitive species (or two traits)
# asymetric (drift) matrix with intervention from the lowest layer
alpha <- matrix(c(0.15,0,0.1,0.1),2,2) 
# scatter matrix
sigma <- matrix(c(0.01,0.005,0.005,0.01),2)
# ancestral states and long term optimum expectation
theta <- matrix(c(0,1,0,.5),2) # columns=traits

# Simulate the data
traits <- mvSIM(timeseries, model="OUTS", param=list(theta=theta, alpha=alpha, sigma=sigma))

# Plot the time series
matplot(traits,type="o",pch=1,  xlab="Time (relative)")

fit1 <- mvOUTS(timeseries, traits, param=list(decomp="qr"))

fit2 <- mvOUTS(timeseries, traits, param=list(decomp="eigen"))

fit3 <- mvOUTS(timeseries, traits, param=list(decomp="diagonal"))

results <- list(fit1,fit2,fit3)
aicw(results)

# Simulate under the MLE
traits2 <- simulate(fit1,tree=timeseries)
matplot(traits2, type="o", pch=1, xlab="Time (relative)")

mvOUTS(timeseries, traits2, param=list(decomp="eigen"))
mvOUTS(timeseries, traits2, param=list(decomp="diagonal"))
mvOUTS(timeseries, traits2, param=list(decomp="upper"))
mvOUTS(timeseries, traits2, param=list(decomp="lower"))


# try user defined constraints
set.seed(100)
ts <- 49
timeseries <- 1:ts

sigma <- matrix(c(0.01,0.005,0.003,0.005,0.01,0.003,0.003,0.003,0.01),3)
# upper triangular matrix with effect of trait 2 on trait 1.
alpha <- matrix(c(0.4,0,0,-0.5,0.3,0,0,0,0.2),3,3) 
theta <- matrix(c(0,0,0,1,0.5,0.5),byrow=TRUE, ncol=3); root=TRUE

data <- mvSIM(timeseries, model="OUTS", param=list(alpha=alpha, 
              sigma=sigma, theta=theta, root=root, 
              names_traits=c("sp 1", "sp 2", "sp 3")))

# plot
matplot(data, type="o", pch=1, xlab="Time (relative)")
legend("bottomright", inset=.05, legend=colnames(data), pch=19, col=c(1,2,3), horiz=TRUE)

# define an user constrained drift matrix
indice <- matrix(NA,3,3)
diag(indice) <- c(1,2,3)
indice[1,2] <- 4

# fit the model
fit_1 <- mvOUTS(timeseries, data, param=list(vcv="fixedRoot", decomp=indice))
fit_2 <- mvOUTS(timeseries, data, param=list(vcv="fixedRoot", decomp="diagonal"))

LRT(fit_1, fit_2)
}



}% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Ornstein Uhlenbeck }
\keyword{ Time series }
\keyword{ Hessian }
\keyword{ OU }% __ONLY ONE__ keyword per line
