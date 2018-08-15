\name{mvgls.pca}
\alias{mvgls.pca}

\title{
Principal Component Analysis (PCA) based on GLS estimate of the traits variance-covariance matrix (possibly regularized).
}
\description{
Performs a principal component analysis (PCA) on a regularized variance-covariance matrix obtained using the \code{mvgls} function. With "evolutionary" models, this performs the so-called phylogenetic PCA.}
\usage{

mvgls.pca(object, plot=TRUE, ...)
  
}
\arguments{
  \item{object}{
A model fit obtained by the \code{mvgls} function.
}
  \item{plot}{
Plot of the PC's axes. Default is TRUE (see details).'}

 \item{...}{
Options to be passed through. (e.g., \code{axes=c(1,2)}, \code{col}, \code{pch}, \code{cex}, \code{mode="cov"} or \code{"corr"}, etc.)}
}


\value{
	a list with the following components
  
  \item{scores}{the PC scores}
  \item{values}{the eigenvalues of the variance-covariance matrix estimated by mvgls}
  \item{vectors}{the eigenvectors of the variance-covariance matrix estimated by mvgls}
  \item{rank}{the rank of the estimated variance-covariance matrix}
  }

\details{
\code{mvgls.pca} allows computing a principal component analysis based on a GLS estimate of the covariance matrix (see \code{mvgls}). The phylogenetic PCA (following Revell 2009) is a special case obtained from the (possibly regularized) evolutionary variance-covariance matrix (see also the \code{phyl.pca_pl} function in RPANDA). In the high-dimensional case the contribution of the firsts PC axes tend to be overestimated with traditional maximum likelihood approaches. Penalized/regularized model fits reduce this bias and allow the estimation of various error structure models (see Clavel et al. 2018).
Ploting options, the number of axes to display (\code{axes=c(1,2)} is the default), and whether the covariance (\code{mode="cov"}) or correlation (\code{mode="corr"}) should be used can be specified through the ellipsis "\code{...}" argument.
}

\note{
	Contrary to conventional PCA, the principal axes of the gls PCA are not orthogonal, they represent the main axes of independent (according to a given phylogenetic or time-series model) changes.}
	
\references{
Clavel, J., Aristide, L., Morlon, H., 2018. A Penalized Likelihood framework for high-dimensional phylogenetic comparative methods and an application to new-world monkeys brain evolution. Systematic Biology DOI:10.1093/sysbio/syy045.

Revell, L.J., 2009. Size-correction and principal components for intraspecific comparative studies. Evolution, 63:3258-3268.
}

\author{J. Clavel}

\seealso{
\code{\link{mvgls}},
\code{\link{GIC}}
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

# The conventional phylogenetic PCA
phylo_pca <- mvgls(Y~1, tree=tree, model="BM", method="LL")
mvgls.pca(phylo_pca, plot=T) 




# fit a multivariate Pagel lambda model with Penalized likelihood
fit <- mvgls(Y~1, tree=tree, model="lambda", method="LOO", penalty="RidgeAlt")

# Perform a regularized phylogenetic PCA using the model fit (Pagel lambda model)
pca_results <- mvgls.pca(fit, plot=T) 

# retrieve the scores
head(pca_results$scores)
}
}

\keyword{ GIC }
\keyword{ GLS }
\keyword{ PCA }
\keyword{ Regularization }
\keyword{ Penalized likelihood }
\keyword{ High dimensions }% __ONLY ONE__ keyword per line
