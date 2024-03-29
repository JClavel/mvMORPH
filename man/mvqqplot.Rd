\name{mvqqplot}
\alias{mvqqplot}

\title{
Quantile-Quantile plots for multivariate models fit with \code{mvgls} or \code{mvols}
}
\description{
The quantile-quantile plots of the Chi square distribution is used to assess multivariate normality and detect outliers using the squared Mahalanobis distances from the models residuals.
}
\usage{

mvqqplot(object, conf=0.95, ...)
  
}
\arguments{
  \item{object}{
A model fit obtained by the \code{mvgls} or \code{mvols} function.
}
  \item{conf}{
Confidence interval for the approximate envelope. Default is 0.95.}

 \item{...}{
Graphical options.}
}

\details{
The empirical quantiles of standardized Mahalanobis distances (Caroni 1987) estimated from models fit by \code{mvgls} (or \code{mvols}) are compared to the quantiles of a Chi square distribution with 'p' degrees of freedom (where 'p' is the number of dimensions) when models are fit by maximum likelihood (\code{method='LL'}). For penalized likelihood model fit (regularized covariance), a matching moments method is used to map the standardized Mahalanobis distances to the Chi square distribution (Clavel, in prep.). This last option is experimental and still under development.
}

\value{
	a list with components
  
  \item{squared_dist}{the squared Mahalanobis distances (standardized)}
  \item{chi2q}{the chi squared quantiles}
  }


\note{
	Chi square Q-Q plots may be outperformed by F based Q-Q plots for identifying outliers (Hardin & Rocke 2005). The function is still under development.}
	
\references{
Caroni, C. 1987. Residuals and Influence in the multivariate linear model. Journal of the Royal Statistical Society 36(4): 365-370.

Clavel, J., Aristide, L., Morlon, H., 2019. A Penalized Likelihood framework for high-dimensional phylogenetic comparative methods and an application to new-world monkeys brain evolution. Systematic Biology 68(1): 93-116.

Clavel, J., Morlon, H. 2020. Reliable phylogenetic regressions for multivariate comparative data: illustration with the MANOVA and application to the effect of diet on mandible morphology in phyllostomid bats. Systematic Biology 69(5): 927-943.
}

\author{J. Clavel}

\seealso{
\code{\link{mvgls}},
\code{\link{mvols}},
\code{\link{manova.gls}}
}

\examples{
\donttest{
data(phyllostomid)

# Fit a linear model by PL
fit <- mvgls(mandible~grp1, data=phyllostomid, phyllostomid$tree, model="lambda", method="PL") 

# QQ plots
mvqqplot(fit, lty=2, conf=0.99)
}
}

\keyword{ QQ plots }
\keyword{ GLS }
\keyword{ OLS }
\keyword{ Mahalanobis }
\keyword{ Regularization }
\keyword{ Penalized likelihood }
\keyword{ High dimensions }% __ONLY ONE__ keyword per line
