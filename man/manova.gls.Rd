\name{manova.gls}
\alias{manova.gls}

\title{
Multivariate Analysis of Variance
}
\description{
Performs a Multivariate Analysis of Variance (MANOVA) on an object fitted by the \code{mvgls} function. With the regularized methods by penalized likelihood implemented in \code{mvgls} (ridgeArch penalty), this function can be used to compare model fit on high-dimensional datasets (where the number of variable is larger than the number of observations).}
\usage{

manova.gls(object, test=c("Pillai", "Wilks", "Hotelling-Lawley", "Roy"),
            type=c("I","II","III"), nperm=100L, L=NULL, ...)
  
}
\arguments{
  \item{object}{
A model fit obtained by the \code{mvgls} function.
}
  \item{test}{
The multivariate test statistic to compute - "Wilks", "Pillai", "Hotelling-Lawley", or "Roy"}

  \item{type}{
The type of test (sums of squares and cross-products) - "I", "II", or "III"}

  \item{nperm}{
The number of permutations used for building the null distribution of the chosen statistic. Permutation is the only available approach for high-dimensional PL models, but either permutations or parametric tests can be used with maximum likelihood (method "LL")}

  \item{L}{
A (contrasts) matrix or a vector giving linear combinations of the coefficients rows.}

 \item{...}{
Further arguments to be passed through. (e.g., \code{nbcores=2L} to provide the number of cores used for parallel calculus; \code{parametric=FALSE} to obtain permutation instead of parametric tests for maximum likelihood fit; \code{verbose=TRUE} to display a progress bar during permutations; \code{rhs=0} the "right-hand-side" vector for general linear hypothesis testing. See details)}
}


\value{
	An object of class 'manova.mvgls' which is usually printed. It contains a list including the following components:
  
  \item{test}{the multivariate test statistic used}
  \item{type}{the type of tests used to compute the SSCP matrices}
  \item{stat}{the statistic calculated for each terms in the model}
  \item{pvalue}{the pvalues calculated for each terms in the model}
 
  }

\details{
\code{manova.gls} allows performing multivariate tests (e.g. Pillai's, Wilks, Hotelling-Lawley and Roy largest root) on generalized least squares (GLS) linear model (objects of class "mvgls") fit by either maximum likelihood (\code{method="LL"}) or penalized likelihood (\code{method="LOO"}) using the \code{mvgls} function. 

General Linear Hypothesis of the form:

\bold{LB=O}

an \bold{L} matrix specifying linear combinations of the model coefficients (\bold{B}) can be provided through the argument \code{L}. The right-hand-side matrix \bold{O} is a constant matrix (of zeros by default) that can be provided through the argument \code{rhs}.

Permutations on high-dimensional datasets is time consuming. You can use the option \code{nbcores} to parallelize the calculus over several cores using forking in UNIX platforms (default is \code{nbcores=1L}. Estimated time to completion is displayed when \code{verbose=TRUE}.
}

\note{
	For PL methods, only the "RidgeArch" penalty is allowed for now.}
	
\references{
Clavel, J., Aristide, L., Morlon, H., 2019. A Penalized Likelihood framework for high-dimensional phylogenetic comparative methods and an application to new-world monkeys brain evolution. Systematic Biology 68(1): 93-116.

Clavel, J., Morlon, H. (submitted). Reliable phylogenetic regressions for multivariate comparative data: illustration with the MANOVA and application to the effect of diet on mandible morphology in phyllostomid bats.
}

\author{J. Clavel}

\seealso{
\code{\link{mvgls}},
\code{\link{GIC}}
\code{\link{EIC}}
}

\examples{
\donttest{
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

# MANOVA on a binary predictor
model2 <- mvgls(y~grp, data=dataset, tree=tree, model="lambda", method="LOO")

# Multivariate test
(multivariate_test <- manova.gls(model2, nperm=999, test="Pillai", verbose=TRUE))
}
}

\keyword{ GIC }
\keyword{ GLS }
\keyword{ MANOVA }
\keyword{ Multivariate tests }
\keyword{ General Linear Hypothesis }
\keyword{ Models comparison }
\keyword{ Regularization }
\keyword{ Penalized likelihood }
\keyword{ High dimensions }% __ONLY ONE__ keyword per line
