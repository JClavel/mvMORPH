\name{mvgls}
\alias{mvgls}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit linear model using Generalized Least Squares to multivariate (high-dimensional) data sets.
%%  ~~function to do ... ~~
}
\description{
This function use maximum likelihood (or restricted likelihood) and penalized likelihood approaches to fit linear models where the errors are allowed to be correlated (i.e. a GLS model for serially correlated phylogenetic and time-series data). \code{mvgls} uses a penalized-likelihood (PL) approach (see descriptions in Clavel et al. 2019) to fit linear models to high-dimensional data sets (where the number of variables \emph{p} is approaching or is larger than the number of observations \emph{n}). The PL approach generally provides improved estimates compared to ML.

%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mvgls(formula, data, tree, model, method=c("PL-LOOCV","LL"),
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

\item{tree}{
Phylogenetic tree (an object of class "\code{phylo}") or a time-series object (not yet available).
}

\item{model}{
The evolutionary model: "BM" is Brownian Motion, "OU" is Ornstein-Uhlenbeck, "EB" is Early Burst, and "lambda" is Pagel's lambda transformation.

}

\item{method}{
The method used to fit the model. "PL-LOOCV" (or equivalently just "LOOCV") is the nominal leave one out cross-validation of the penalized log-likelihood, "LL" is the log-likelihood (used in the conventional ML and REML estimation). Two approximated LOOCV methods are also available: "H&L" and "Mahalanobis". The method "H&L" is a fast LOOCV approach based on Hoffbeck and Landgrebe (1996) tricks, and "Mahalanobis" is an approximation of the LOOCV score proposed by Theiler (2012). Both "H&L" and "Mahalanobis" work only with the "RidgeArch" penalty and for intercept only models (see details).
}

\item{REML}{
Use REML (default) or ML for estimating the parameters.
}

\item{...}{
Options to be passed through. For instance the type of penalization: 
\code{penalty="RidgeArch"} (default), \code{penalty="RidgeAlt"}, or \code{penalty="LASSO"}. The target matrices used by "RidgeArch" and "RidgeAlt" penalizations: \code{target="unitVariance"}, \code{target="Variance"} or \code{target="null"}... etc. (see details)
}

}
\details{
\code{mvgls} allows fitting various multivariate linear models to high-dimensional datasets (i.e. where the number of variables \emph{p} is larger than \emph{n}) for which the residuals have a correlated structure (e.g. evolutionary models such as BM and OU). Models estimates are generally more accurate than maximum likelihood methods. Models fit can be compared using the GIC or EIC criterion (see ?\code{GIC} and ?\code{EIC}) and hypothesis testing can be performed using the \code{manova.gls} function.

The tree is assumed to be fully dichotomic and in "postorder", otherwise the functions \code{multi2di} and \code{reorder.phylo} are used internally.

The various \emph{arguments} that can be passed through \bold{"..."}:

\bold{"penalty"} - The "penalty" argument allow specifying the type of penalization used for regularization (described in Clavel et al. 2019). The various penalizations are: \code{penalty="RidgeArch"} (the default), \code{penalty="RidgeAlt"} and \code{penalty="LASSO"}. The "RidgeArch" penalization shrink linearly the covariance matrix toward a target structure (see below for \code{target}). This penalization is generally fast and the tuning parameter is bounded between 0 and 1 (see van Wieringen & Peeters 2016). The "RidgeAlt" penalization scheme uses a quadratic ridge penalty to shrink the covariance matrix toward a specified target matrix (see \code{target} below and also see van Wieringen & Peeters 2016). Finally, the "LASSO" regularize the covariance matrix by estimating a sparse estimate of its inverse (Friedman et al. 2008). The computation of the solution for this penalization is computationally intensive. Moreover, this penalization scheme is not invariant to arbitrary rotations of the data.

\bold{"target"} - This argument allows specifying the target matrix toward which the covariance matrix is shrunk for "Ridge" penalties. \code{target="unitVariance"} (for a diagonal target proportional to the identity) and \code{target="Variance"} (for a diagonal unequal variance target) can be used with both "RidgeArch" and "RidgeAlt" penalties. \code{target="null"} (a null target matrix) is only available for "RidgeAlt". Penalization with the "Variance" target shrink the eigenvectors of the covariance matrix and is therefore not rotation invariant. See details on the various target properties in Clavel et al. (2019).

\bold{"error"} - If \code{TRUE} the measurement error (or intra-specific variance) is estimated from the data as a nuisance parameter (like in mixed models). It should probably be systematically used with empirical data. See also Housworth et al. 2004 and Clavel et al. 2019 for details on the proposed implementation.

\bold{"scale.height"} - Whether the tree should be scaled to unit length or not.

\bold{"echo"} - Whether the results must be returned or not.

\bold{"grid_search"} - A logical indicating whether or not a preliminary grid search must be performed to find the best starting values for optimizing the log-likelihood (or penalized log-likelihood). User-specified starting values can be provided through the \bold{start} argument. Default is \code{TRUE}.

\bold{"upper"} - The upper bound for the parameter search with the "\code{L-BFGS-B}" method. See \code{optim} for details.

\bold{"lower"} - The lower bound for the parameter search with the "\code{L-BFGS-B}" method. See \code{optim} for details.

\bold{"tol"} - Minimum value for the regularization parameter. Singularities can occur with a zero value in high-dimensional cases. (default is \code{NULL})


%%  ~~ If necessary, more details than the description above ~~
}
\value{

An object of class '\code{mvgls}'. It contains a list including the following components:

\item{coefficients  }{a named vector of coefficients}
\item{residuals     }{the residuals ("raw") of the model. That is response minus fitted values. Uses the \code{residuals(x, type="normalized")} function to obtain the normalized residuals.}
\item{fitted.values }{the fitted mean values}
\item{variables     }{the variables used for model fit}
\item{sigma         }{the estimated covariance (Pi) and precision (P) matrix, as well as the sample estimate (S)}
\item{model         }{the evolutionary model. But more generally, the model used to specify the structure within the residuals}
\item{logLik	    	}{either the (negative) log-likelihood when \code{method="LL"} or the cross-validated penalized likelihood}
\item{param			    }{the (evolutionary) model parameter estimates}
\item{tuning		    }{the regularization/tuning parameter of the penalized likelihood}
\item{mserr		    	}{the estimated standard error when \code{error=TRUE}}
\item{start_values  }{the starting parameters used for the optimization of the LL or PL}
\item{corrSt        }{a list including the transformed tree, the determinant obtained from its covariance matrix and the normalized variables (by the inverse square root of the covariance matrix of the phylogenetic tree or the time-series)}
\item{penalty       }{the penalty used for the penalized likelihood approach}
\item{target        }{the target used with the "RidgeArch" or "RidgeAlt" penalized likelihood approaches}
\item{REML	    		}{logical indicating if the REML (\code{TRUE}) or ML (\code{FALSE}) method has been used}
\item{opt	    	  	}{optimizing function output. See \code{optim}}


%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{

Clavel, J., Aristide, L., Morlon, H., 2019. A Penalized Likelihood framework for high-dimensional phylogenetic comparative methods and an application to new-world monkeys brain evolution. Systematic Biology 68(1): 93-116.

Friedman J., Hastie T., Tibshirani R. 2008. Sparse inverse covariance estimation with the graphical lasso. Biostatistics. 9:432-441.

Hoffbeck J.P., Landgrebe D.A. 1996. Covariance matrix estimation and classification with limited training data. IEEE Trans. Pattern Anal. Mach. Intell. 18:763-767.

Housworth E.A., Martins E.P., Lynch M. 2004. The phylogenetic mixed model. Am. Nat. 163:84-96.

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
\code{\link{EIC}}
\code{\link{GIC}}
\code{\link{mvgls.pca}}
\code{\link{fitted.mvgls}}
\code{\link{residuals.mvgls}}
\code{\link{coef.mvgls}}
\code{\link{vcov.mvgls}}
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
fit2 <- mvgls(Y~1, data=data, tree, model="OU", penalty="RidgeArch")
fit3 <- mvgls(Y~1, data=data, tree, model="EB", penalty="RidgeArch")

GIC(fit1); GIC(fit2); GIC(fit3) # BM have the lowest GIC value

# Testing for phylogenetic signal with model fit
signal <- mvgls(Y~1, data=data, tree, model="lambda", penalty="RidgeArch")
summary(signal)


# A High-dimensional dataset
p <- 200 # number of traits (p>n)

R <- crossprod(matrix(runif(p*p), ncol=p)) # a random symmetric matrix (covariance)
# simulate a BM dataset
Y <- mvSIM(tree, model="BM1", nsim=1, param=list(sigma=R, theta=rep(0,p))) 
data=list(Y=Y)

# Fast LOOCV using "H&L" with RidgeArch penalization
summary(mvgls(Y~1, data=data, tree, model="BM", penalty="RidgeArch", method="H&L"))

}

}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ GIC }
\keyword{ GLS }
\keyword{ Regularization }
\keyword{ Penalized likelihood }
\keyword{ Multivariate Linear Models }
\keyword{ High dimensions }% __ONLY ONE__ keyword per line
