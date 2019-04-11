\name{EIC}
\alias{EIC}

\title{
Extended Information Criterion (EIC) to compare models fit with \code{mvgls} by Maximum Likelihood (ML) or Penalized Likelihood (PL).
}

\description{
The EIC (Kitagawa & Konishi 2010), uses bootstrap to estimate the bias term of the Extended Information Criterion. This criterion allows comparing models fit by Maximum Likelihood (ML) or Penalized Likelihood (PL).
}

\usage{

EIC(object, nboot=100L, nbcores=1L, ...)
  
}

\arguments{
  \item{object}{
An object of class 'mvgls'. See ?mvgls}
  \item{nboot}{
The number of boostrap replicates used for estimating the EIC.}
  \item{nbcores}{
The number of cores used to speed-up the computations (uses the 'parallel' package)}
  \item{...}{
Options to be passed through.}
}


\value{
		a list with the following components
  
  \item{LogLikelihood}{the log-likelihood estimated for the model with estimated parameters}
  \item{EIC}{the EIC criterion}
  \item{se}{the standard error of the bias term estimated by bootstrap}
  \item{bias}{the values of the bias term estimated from the boostrapped replicates to compute the EIC}

}

\details{
The Extended Information Criterion (\code{EIC}) allows comparing the fit of various models estimated by Penalized Likelihood or Maximum Likelihood (see ?\code{mvgls}). The estimation of the criterion is based on bootstrap simulations rather than analytical solutions or approximations (see for instance ?\code{GIC}). With small sample sizes, the variability around the bootstrap estimates is expected to be high, and one must increase the number of bootstrap replicates. Parallel computation allows to speed-up the computations.
}

	
\references{
Clavel, J., Aristide, L., Morlon, H., 2019. A Penalized Likelihood framework for high-dimensional phylogenetic comparative methods and an application to new-world monkeys brain evolution. Syst. Biol. 68:93-116.

Kitagawa G., Konishi S. 2010. Bias and variance reduction techniques for bootstrap information criterion. Ann. Inst. Stat. Math. 62:209-234.

Konishi S., Kitagawa G. 1996. Generalised information criteria in model selection. Biometrika. 83:875-890.

}

\author{J. Clavel}

\seealso{
\code{\link{GIC}}
\code{\link{mvgls}}
\code{\link{mvgls.pca}}
}

\examples{
\donttest{

set.seed(1)
n <- 32 # number of species
p <- 50 # number of traits

tree <- pbtree(n=n) # phylogenetic tree
R <- crossprod(matrix(runif(p*p), ncol=p)) # a random symmetric matrix (covariance)
# simulate a dataset
Y <- mvSIM(tree, model="BM1", nsim=1, param=list(sigma=R))

fit1 <- mvgls(Y~1, tree=tree, model="BM", method="H&L")
fit2 <- mvgls(Y~1, tree=tree, model="OU", method="H&L")


EIC(fit1); EIC(fit2)
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Model comparison }
\keyword{ GLS }
\keyword{ High dimensions }% __ONLY ONE__ keyword per line