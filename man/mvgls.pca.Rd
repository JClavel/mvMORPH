\name{mvgls.pca}
\alias{mvgls.pca}

\title{
Principal Component Analysis (PCA) based on GLS (or OLS) estimate of the traits variance-covariance matrix (possibly regularized)
}
\description{
Performs a principal component analysis (PCA) on a regularized variance-covariance matrix obtained using the \code{mvgls} or the \code{mvols} function. With "evolutionary" models in \code{mvgls}, this performs the so-called phylogenetic PCA.}
\usage{

mvgls.pca(object, plot=TRUE, ...)
  
}
\arguments{
  \item{object}{
A model fit obtained by the \code{mvgls} or \code{mvols} function.
}
  \item{plot}{
Plot of the PC's axes. Default is TRUE (see details).'}

 \item{...}{
Options to be passed through. (e.g., \code{axes=c(1,2)}, \code{col}, \code{pch}, \code{cex}, \code{mode="cov"} or \code{"corr"}, etc.)}
}


\value{
	a list with the following components
  
  \item{scores}{the PC scores}
  \item{values}{the eigenvalues of the variance-covariance matrix estimated by mvgls or mvols}
  \item{vectors}{the eigenvectors of the variance-covariance matrix estimated by mvgls or mvols}
  \item{rank}{the rank of the estimated variance-covariance matrix}
  }

\details{
\code{mvgls.pca} allows computing a principal component analysis based on a GLS (or OLS) estimate of the covariance matrix (see \code{mvgls} and \code{mvols}). The phylogenetic PCA (following Revell 2009) is a special case obtained from the (possibly regularized) evolutionary variance-covariance matrix (see also the \code{phyl.pca_pl} function in RPANDA). In the high-dimensional case the contribution of the firsts PC axes tend to be overestimated with traditional maximum likelihood approaches. Penalized/regularized model fit reduce this bias and allow incorporating various residuals structures (see Clavel et al. 2019).
Plotting options, the number of axes to display (\code{axes=c(1,2)} is the default), and whether the covariance (\code{mode="cov"}) or correlation (\code{mode="corr"}) should be used can be specified through the ellipsis "\code{...}" argument.
}

\note{
	Contrary to conventional PCA (for instance using \code{mvols} with "LL" method), the principal axes of the gls PCA are not orthogonal, they represent the main axes of independent (according to a given phylogenetic or time-series model) evolutionary changes.}
	
\references{
Clavel, J., Aristide, L., Morlon, H., 2019. A Penalized Likelihood framework for high-dimensional phylogenetic comparative methods and an application to new-world monkeys brain evolution. Systematic Biology 68(1): 93-116.

Revell, L.J., 2009. Size-correction and principal components for intraspecific comparative studies. Evolution, 63:3258-3268.
}

\author{J. Clavel}

\seealso{
\code{\link{mvgls}},
\code{\link{mvols}},
\code{\link{GIC}},
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

# The conventional phylogenetic PCA
phylo_pca <- mvgls(Y~1, tree=tree, model="BM", method="LL")
mvgls.pca(phylo_pca, plot=TRUE) 




# fit a multivariate Pagel lambda model with Penalized likelihood
fit <- mvgls(Y~1, tree=tree, model="lambda", method="LOO", penalty="RidgeAlt")

# Perform a regularized phylogenetic PCA using the model fit (Pagel lambda model)
pca_results <- mvgls.pca(fit, plot=TRUE) 

# retrieve the scores
head(pca_results$scores)
}
}

\keyword{ GIC }
\keyword{ GLS }
\keyword{ OLS }
\keyword{ PCA }
\keyword{ Regularization }
\keyword{ Penalized likelihood }
\keyword{ High dimensions }% __ONLY ONE__ keyword per line
