\name{GIC}
\alias{GIC}

\title{
Generalized Information Criterion (GIC) to compare models fit with \code{mvgls} by Maximum Likelihood (ML) or Penalized Likelihood (PL).
}
\description{
The GIC (Konishi & Kitagawa 1996) allows comparing models fit by Maximum Likelihood (ML) or Penalized Likelihood (PL).
}
\usage{


GIC(object, ...)
  
  
}
\arguments{
  \item{object}{
An object of class 'mvgls'. See ?mvgls}
  \item{...}{
Options to be passed through.}
}


\value{
		a list with the following components
  
  \item{LogLikelihood}{the log-likelihood estimated for the model with estimated parameters}
  \item{GIC}{the GIC criterion}
  \item{bias}{the value of the bias term estimated to compute the GIC}

}

\details{
The Generalized Information Criterion (\code{GIC}) allows comparing the fit of various models estimated by Penalized Likelihood (see ?\code{mvgls}). See also the \code{gic_criterion} function in the RPANDA package.
Note that the current implementation of the criterion has not been tested for multiple predictors comparison. Prefer simulation based comparisons instead.
}

	
\references{
Clavel, J., Aristide, L., Morlon, H., 2018. A Penalized Likelihood framework for high-dimensional phylogenetic comparative methods and an application to new-world monkeys brain evolution. Systematic Biology DOI:10.1093/sysbio/syy045.

Konishi S., Kitagawa G. 1996. Generalised information criteria in model selection. Biometrika. 83:875-890.


}

\author{J. Clavel}

\seealso{
\code{\link{mvgls}}
\code{\link{mvgls.pca}}
}

\examples{
\donttest{

set.seed(1)
n <- 32 # number of species
p <- 50 # number of traits

tree <- pbtree(n=n) # phylogenetic tree
R <- crossprod(matrix(rnorm(p*p), ncol=p)) # a random symmetric matrix (covariance)
# simulate a dataset
Y <- mvSIM(tree, model="BM1", nsim=1, param=list(sigma=R))

fit1 <- mvgls(Y~1, tree=tree, model="BM", method="H&L")
fit2 <- mvgls(Y~1, tree=tree, model="OU", method="H&L")


GIC(fit1); GIC(fit2)
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Model comparison }
\keyword{ GLS }
\keyword{ High dimensions }% __ONLY ONE__ keyword per line
