\name{manova.gls}
\alias{manova.gls}

\title{
Multivariate Analysis of Variance
}
\description{
Performs a Multivariate Analysis of Variance (MANOVA) on an object fitted by the \code{mvgls} or the \code{mvols} function. With the regularized methods by penalized likelihood implemented in \code{mvgls} and \code{mvols} (ridgeArch penalty), this function can be used to compare model fit on high-dimensional datasets (where the number of variables is larger than the number of observations). When model fit is performed by maximum likelihood (\code{method="LL"}), both parametric and permutation tests are possible.}
\usage{

manova.gls(object, test=c("Pillai", "Wilks", "Hotelling-Lawley", "Roy"),
            type=c("I","II","III"), nperm=1000L, L=NULL, ...)
  
}
\arguments{
  \item{object}{
A model fit obtained by the \code{mvgls} or \code{mvols} function.}

  \item{test}{
The multivariate test statistic to compute - "Wilks", "Pillai", "Hotelling-Lawley", or "Roy"}

  \item{type}{
The type of test (sums of squares and cross-products) - "I", "II", or "III"}

  \item{nperm}{
The number of permutations used for building the null distribution of the chosen statistic. Permutation is the only available approach for high-dimensional PL models, but either permutations or parametric tests can be used with maximum likelihood (method "LL" in \code{mvgls} and \code{mvols})}

  \item{L}{
A (contrasts) matrix or a vector giving linear combinations of the coefficients rows.}

 \item{...}{
Further arguments to be passed through. (e.g., \code{nbcores=2L} to provide the number of cores used for parallel calculus; \code{parametric=FALSE} to obtain permutation instead of parametric tests for maximum likelihood fit; \code{verbose=TRUE} to display a progress bar during permutations; \code{rhs=0} the "right-hand-side" vector for general linear hypothesis testing; \code{P} can be used to specify a matrix of contrasts giving linear combinations of the coefficient columns. See details)}
}


\value{
	An object of class 'manova.mvgls' which is usually printed. It contains a list including the following components:
  
  \item{test}{the multivariate test statistic used}
  \item{type}{the type of tests used to compute the SSCP matrices}
  \item{stat}{the statistic calculated for each terms in the model}
  \item{pvalue}{the pvalues calculated for each terms in the model}
 
  }

\details{
\code{manova.gls} allows performing multivariate tests (e.g. Pillai's, Wilks, Hotelling-Lawley and Roy largest root) on generalized least squares (GLS) linear model (objects of class "mvgls", or OLS with objects of class "mvols") fit by either maximum likelihood (\code{method="LL"}) or penalized likelihood (\code{method="PL-LOO"}) using the \code{mvgls} or \code{mvols} function.

General Linear Hypothesis of the form:

\deqn{\bold{LB=O}} or \deqn{\bold{LBP=O}}

Where \bold{L} is a matrix specifying linear combinations of the model coefficients (\bold{B}) can be provided through the argument \code{L}. This type of "contrasts" matrix allows testing specific hypotheses (for instance pairwise differences - see \code{?pairwise.contrasts} and \code{?pairwise.glh}). The right-hand-side matrix \bold{O} is a constant matrix (of zeros by default) that can be provided through the argument \code{rhs}. \bold{P} is a matrix specifying linear combinations of the model coefficients (\bold{B}) estimated for each responses (usually used in repeated measures designs or for testing linear, quadratic, etc. relationships between successive responses).

Permutations on high-dimensional datasets is time consuming. You can use the option \code{nbcores} to parallelize the calculus over several cores using forking in UNIX platforms (default is \code{nbcores=1L}. Estimated time to completion is displayed when \code{verbose=TRUE}.
}

\note{
	For PL methods, only the "RidgeArch" penalty is allowed for now.
	A tutorial is available from Dryad: https://doi.org/10.5061/dryad.jsxksn052
	}
	
\references{
Clavel, J., Aristide, L., Morlon, H., 2019. A Penalized Likelihood framework for high-dimensional phylogenetic comparative methods and an application to new-world monkeys brain evolution. Systematic Biology 68(1): 93-116.

Clavel, J., Morlon, H. 2020. Reliable phylogenetic regressions for multivariate comparative data: illustration with the MANOVA and application to the effect of diet on mandible morphology in phyllostomid bats. Systematic Biology 69(5): 927-943.
}

\author{J. Clavel}

\seealso{
\code{\link{mvgls}},
\code{\link{mvols}},
\code{\link{pairwise.glh}},
\code{\link{GIC}},
\code{\link{EIC}}
}

\examples{
\donttest{

# ---------------------------------------------------- #
# Multivariate regression tests (continuous predictor) #
# ---------------------------------------------------- #

set.seed(1)
n <- 32 # number of species
p <- 30 # number of traits

tree <- pbtree(n=n) # phylogenetic tree
R <- crossprod(matrix(runif(p*p),p))  # a random symmetric matrix (covariance)

# simulate a dataset
Y <- mvSIM(tree, model="BM1", nsim=1, param=list(sigma=R))
X <- rnorm(n) # continuous
grp <- rep(1:2, each=n/2)
dataset <- list(y=Y, x=X, grp=as.factor(grp))

# Model fit
model1 <- mvgls(y~x, data=dataset, tree=tree, model="BM", method="LOO")

# Multivariate test
(multivariate_test <- manova.gls(model1, nperm=999, test="Pillai"))

# ---------------------------------------------------- #
# Multivariate regression tests (discrete predictor)   #
# ---------------------------------------------------- #

# MANOVA on a binary predictor
model2 <- mvgls(y~grp, data=dataset, tree=tree, model="lambda", method="LOO")

# Multivariate test
(multivariate_test <- manova.gls(model2, nperm=999, test="Pillai", verbose=TRUE))

# ---------------------------------------------------- #
# Parametric MANOVA tests                              #
# ---------------------------------------------------- #

# When p<n we can use non-penalized approaches and parametric tests
model2b <- mvgls(y[,1:2]~grp, data=dataset, tree=tree, model="lambda", method="LL")
(multivariate_test2b <- manova.gls(model2b, test="Pillai"))

# ---------------------------------------------------- #
# Multivariate contrasts tests                         #
# ---------------------------------------------------- #

# Multivariate contrasts allow testing specific hypotheses 
# (see also ?pairwise.glh and ?pairwise.contrasts)

# We can replicate the above result by testing if the
# group means are different using the following contrast:
L = matrix(c(0,1), ncol=2)
(manova.gls(model2b, test="Pillai", L=L))

# ---------------------------------------------------- #
# Repeated measures design tests                       #
# ---------------------------------------------------- #

# Contrasts can be used also to test if there's differences
# between repeated measures (responses variables)
# For instance, for comparing y[,1] and y[,2], define the contrast:
P = matrix(c(1,-1), nrow=2)
(manova.gls(model2b, test="Pillai", P=P, L=L))

}
}

\keyword{ GIC }
\keyword{ GLS }
\keyword{ OLS }
\keyword{ MANOVA }
\keyword{ Multivariate tests }
\keyword{ General Linear Hypothesis }
\keyword{ Models comparison }
\keyword{ Regularization }
\keyword{ Penalized likelihood }
\keyword{ High dimensions }% __ONLY ONE__ keyword per line
