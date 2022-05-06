\name{mvols}
\alias{mvols}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit linear model using Ordinary Least Squares to multivariate (high-dimensional) data sets
%%  ~~function to do ... ~~
}
\description{
This function uses maximum likelihood (or restricted likelihood) and penalized likelihood approaches to fit linear models with independent observations (this is the multivariate (and penalized) counterpart to the base \code{lm} function). \code{mvols} uses a penalized-likelihood (PL) approach (see descriptions in Clavel et al. 2019) to fit linear models to high-dimensional data sets (where the number of variables \emph{p} is approaching or is larger than the number of observations \emph{n}). The PL approach generally provides improved estimates compared to ML. OLS is a special case  of GLS linear models (and a wrapper of \code{mvgls}) and can be used with all the package functions working on \code{mvgls} class objects.

%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mvols(formula, data, method=c("PL-LOOCV","LL"),
		REML=TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{ 

\item{formula}{
An object of class "\code{formula}" (a two-sided linear formula describing the model to be fitted. See for instance ?\code{lm})
}

\item{data}{
An optional list, data.frame or environment containing the variables in the model. If not found in \emph{data} the variables are taken from the current environment. Prefer \code{list} for blocks of multivariate responses unless you're specifying the response variables by their names using \code{cbind} with data.frame.
}

\item{method}{
The method used to fit the model. "PL-LOOCV" (or equivalently just "LOOCV") is the nominal leave one out cross-validation of the penalized log-likelihood, "LL" is the log-likelihood (used in the conventional ML and REML estimation). Two approximated LOOCV methods are also available: "H&L" and "Mahalanobis". The method "H&L" is a fast LOOCV approach based on Hoffbeck and Landgrebe (1996) tricks, and "Mahalanobis" is an approximation of the LOOCV score proposed by Theiler (2012). Both "H&L" and "Mahalanobis" work only with the "RidgeArch" penalty and for intercept only models (i.e. of the form Y~1, see also details).
}

\item{REML}{
Use REML (default) or ML for estimating the parameters.
}

\item{...}{
Options to be passed through. For instance the type of penalization: 
\code{penalty="RidgeArch"} (default), \code{penalty="RidgeAlt"}, or \code{penalty="LASSO"}. The target matrices used by "RidgeArch" and "RidgeAlt" penalizations: \code{target="unitVariance"}, \code{target="Variance"} or \code{target="null"}... etc. (see details). One can also define contrasts options as for the \code{lm} function.
}

}
\details{
\code{mvols} allows fitting multivariate linear models to multivariate (possibly high-dimensional, i.e. where the number of variables \emph{p} is larger than \emph{n}) datasets. Models estimated using penalized likelihood (e.g., method="PL-LOOCV") are generally more accurate than those estimated by maximum likelihood methods (method="LL") when the number of traits approach the number of observations. PL is the only solution when \emph{p}>\emph{n}.  Models fit can be compared using the GIC or EIC criterion (see ?\code{GIC} and ?\code{EIC}) and hypothesis testing can be performed using the \code{manova.gls} function.


The various \emph{arguments} that can be passed through \bold{"..."}:

\bold{"penalty"} - The "penalty" argument allows specifying the type of penalization used for regularization (described in Clavel et al. 2019). The various penalizations are: \code{penalty="RidgeArch"} (the default), \code{penalty="RidgeAlt"} and \code{penalty="LASSO"}. The "RidgeArch" penalization shrink linearly the "sample"" covariance matrix toward a given target matrix with a specific structure (see below for \code{target}). This penalization is generally fast and the tuning parameter is bounded between 0 and 1 (see van Wieringen & Peeters 2016, Clavel et al. 2019). The "RidgeAlt" penalization scheme uses a quadratic ridge penalty to shrink the covariance matrix toward a specified target matrix (see \code{target} below and also see van Wieringen & Peeters 2016). Finally, the "LASSO" regularize the covariance matrix by estimating a sparse estimate of its inverse - the precision matrix (Friedman et al. 2008). Solving the LASSO penalization is computationally intensive. Moreover, this penalization scheme is not invariant to arbitrary rotations of the data.

\bold{"target"} - This argument allows specifying the target matrix toward which the covariance matrix is shrunk for "Ridge" penalties. \code{target="unitVariance"} (for a diagonal target matrix proportional to the identity) and \code{target="Variance"} (for a diagonal matrix with unequal variance) can be used with both "RidgeArch" and "RidgeAlt" penalties. \code{target="null"} (a null target matrix) is only available for "RidgeAlt". Penalization with the "Variance" target shrinks the eigenvectors of the covariance matrix and is therefore not rotation invariant. See details on the various target properties in Clavel et al. (2019).

\bold{"weights"} - A (named) vector of weights (variances) for all the observations. If provided, a weighted least squares (WLS) rather than OLS fit is performed.

\bold{"echo"} - Whether the results must be returned or not.

\bold{"grid_search"} - A logical indicating whether or not a preliminary grid search must be performed to find the best starting values for optimizing the log-likelihood (or penalized log-likelihood). User-specified starting values can be provided through the \bold{start} argument. Default is \code{TRUE}.

\bold{"tol"} - Minimum value for the regularization parameter. Singularities can occur with a zero value in high-dimensional cases. (default is \code{NULL})


%%  ~~ If necessary, more details than the description above ~~
}

\note{
	This function is a wrapper to the \code{mvgls} function (it uses gls with a diagonal covariance). For these reasons, the function can be used with all the methods working with \code{mvgls} class objects.
	}
	
	
\value{

An object of class '\code{mvols}'. It contains a list including the same components as the \code{mvgls} function (see ?mvgls).

}
\references{

Clavel, J., Aristide, L., Morlon, H., 2019. A Penalized Likelihood framework for high-dimensional phylogenetic comparative methods and an application to new-world monkeys brain evolution. Systematic Biology 68(1): 93-116.

Clavel, J., Morlon, H. 2020. Reliable phylogenetic regressions for multivariate comparative data: illustration with the MANOVA and application to the effect of diet on mandible morphology in phyllostomid bats. Systematic Biology 69(5): 927-943.

Friedman J., Hastie T., Tibshirani R. 2008. Sparse inverse covariance estimation with the graphical lasso. Biostatistics. 9:432-441.

Hoffbeck J.P., Landgrebe D.A. 1996. Covariance matrix estimation and classification with limited training data. IEEE Trans. Pattern Anal. Mach. Intell. 18:763-767.

Theiler J. 2012. The incredible shrinking covariance estimator. In: Automatic Target Recognition XXII. Proc. SPIE 8391, Baltimore, p. 83910P.

van Wieringen W.N., Peeters C.F.W. 2016. Ridge estimation of inverse covariance matrices from high-dimensional data. Comput. Stat. Data Anal. 103:284-303.

%% ~put references to the literature/web site here ~
}
\author{
Julien Clavel
%%  ~~who you are~~
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{manova.gls}}
\code{\link{mvgls}}
\code{\link{EIC}}
\code{\link{GIC}}
\code{\link{mvgls.pca}}
\code{\link{fitted.mvgls}}
\code{\link{residuals.mvgls}}
\code{\link{coef.mvgls}}
\code{\link{vcov.mvgls}}
\code{\link{predict.mvgls}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}

\examples{
\donttest{
set.seed(1)
n <- 32 # number of species
p <- 50 # number of traits (p>n)

tree <- pbtree(n=n, scale=1) # phylogenetic tree
R <- crossprod(matrix(runif(p*p), ncol=p)) # a random covariance matrix
# simulate a BM dataset
Y <- mvSIM(tree, model="BM1", nsim=1, param=list(sigma=R, theta=rep(0,p))) 
data=list(Y=Y)

fit1 <- mvgls(Y~1, data=data, tree, model="BM", penalty="RidgeArch")

# compare to OLS?
fit2 <- mvols(Y~1, data=data, penalty="RidgeArch")

GIC(fit1); GIC(fit2); 

## Fit a model by Maximum Likelihood (rather than Penalized likelihood) when p<<n
fit_ml <- mvols(Y[,1:2]~1, data=data, method="LL")
summary(fit_ml)


}

}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ GIC }
\keyword{ OLS }
\keyword{ GLS }
\keyword{ Regularization }
\keyword{ Penalized likelihood }
\keyword{ Multivariate Linear Models }
\keyword{ High dimensions }% __ONLY ONE__ keyword per line
