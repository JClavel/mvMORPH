\name{mvgls.dfa}
\alias{mvgls.dfa}

\title{
Discriminant Function Analysis (DFA) - also called Linear Discriminant Analysis (LDA) or Canonical Variate Analysis (CVA) - based on multivariate GLS (or OLS) model fit 
}
\description{
Performs a discriminant analysis (DFA) on a regularized variance-covariance matrix obtained using either the \code{mvgls} or \code{mvols} function.
}
\usage{

mvgls.dfa(object, ...)

}
\arguments{
  \item{object}{
A model fit obtained by the \code{mvgls} or the \code{mvols} function.
}

 \item{...}{
Options to be passed through. (e.g., \code{term="the term corresponding to the factor of interest"}, \code{type="I"} for the type of decomposition of the hypothesis matrix (see also manova.gls) , etc.)
}
}


\value{
	a list with the following components
  
  \item{coeffs}{a matrix containing the raw discriminants}
  \item{coeffs.std}{a matrix containing the standardized discriminants}
  \item{scores}{a matrix containing the discriminant scores [residuals X coeffs]}
  \item{residuals}{the centered [with GLS or OLS] response variables}
  \item{H}{the hypothesis (or between group model matrix)}
  \item{E}{the error (or residual model matrix)}
  \item{rank}{the rank of HE^{-1}}
  \item{pct}{the percentage of the discriminant functions}

  }

\details{
\code{mvgls.dfa} allows computing a discriminant analysis based on GLS (or OLS) estimates from a regression model (see \code{mvgls} and \code{mvols}). Discriminant functions can be used for dimensionality reduction, to follow up a MANOVA analysis to describe group separation, or for group prediction.
}

\note{
	Still in development, may not handle special designs. }

\references{
Clavel, J., Aristide, L., Morlon, H., 2019. A Penalized Likelihood framework for high-dimensional phylogenetic comparative methods and an application to new-world monkeys brain evolution. Systematic Biology 68(1): 93-116.

Clavel, J., Morlon, H., 2020. Reliable phylogenetic regressions for multivariate comparative data: illustration with the MANOVA and application to the effect of diet on mandible morphology in Phyllostomid bats. Systematic Biology 69(5): 927-943.

}

\author{J. Clavel}

\seealso{
\code{\link{mvgls}},
\code{\link{mvols}},
\code{\link{manova.gls}},
\code{\link{mvgls.pca}},
\code{\link{predict.mvgls.dfa}},
}

\examples{
\donttest{
library(mvMORPH)
n=64
p=4

tree <- pbtree(n=n)
sigma <- crossprod(matrix(runif(p*p),p,p))
resid <- mvSIM(tree, model="BM1", param=list(sigma=sigma))
Y <- rep(c(0,1.5), each=n/2) + resid
grp <- as.factor(rep(c("gp1","gp2"),each=n/2))
names(grp) = rownames(Y)
data <- list(Y=Y, grp=grp)
mod <- mvgls(Y~grp, data=data, tree=tree, model="BM")

# fda
da1 <- mvgls.dfa(mod)

plot(da1)
}
}

\keyword{ LDA }
\keyword{ CVA }
\keyword{ DFA }
\keyword{ Discriminant }
\keyword{ Regularization }
\keyword{ Penalized likelihood }
\keyword{ High dimensions }% __ONLY ONE__ keyword per line
