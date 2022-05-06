\name{pairwise.glh}
\alias{pairwise.glh}

\title{
Pairwise multivariate tests between levels of a factor
}
\description{
Performs pairwise multivariate tests (e.g. "Pillai") on levels of a factor in a model fitted by the \code{mvgls} or \code{mvols} function. This is achieved by evaluating all the pairwise contrasts using generalized linear hypothesis tests (see also ?manova.gls).}

\usage{

pairwise.glh(object, term=1, test=c("Pillai", "Wilks", "Hotelling-Lawley", "Roy"),
            adjust="holm", nperm=1000L, ...)
  
}
\arguments{
  \item{object}{
A model fit obtained by the \code{mvgls} or \code{mvols} function.
}
  \item{term}{
The factor term in the "object" model fit on which the pairwise tests should be evaluated.
}
  \item{test}{
The multivariate test statistic to compute - "Wilks", "Pillai", "Hotelling-Lawley", or "Roy"}

  \item{adjust}{
The multiple comparison adjustment. See \code{?p.adjust}.}

  \item{nperm}{
The number of permutations used for building the null distribution of the chosen statistic. Permutation is the only available approach for high-dimensional PL models, but either permutations or parametric tests can be used with maximum likelihood (method "LL" in \code{mvgls} and \code{mvols})}

 \item{...}{
Further arguments to be passed through. (e.g., \code{nbcores=2L} to provide the number of cores used for parallel calculus; \code{parametric=FALSE} to obtain permutation instead of parametric tests for maximum likelihood fit; \code{verbose=TRUE} to display a progress bar during permutations; \code{rhs=0} the "right-hand-side" vector for general linear hypothesis testing. See details)}
}


\value{
	An object of class 'pairs.mvgls' which is usually printed. It contains a list including the following components:
  
  \item{test}{the multivariate test statistic used}
  \item{L}{the contrasts used for all the pairwise tests}
  \item{stat}{the statistic calculated for each pairwise comparisons}
  \item{pvalue}{the p-values calculated for each pairwise comparisons}
  \item{adjust}{the adjusted (for multiple comparisons) p-values calculated for each pairwise comparisons}
 
  }

\details{
\code{pairwise.glh} allows performing multivariate tests (e.g. Pillai's, Wilks, Hotelling-Lawley and Roy largest root) on generalized least squares (GLS) linear model (objects of class "mvgls") fit by either maximum likelihood (\code{method="LL"}) or penalized likelihood (\code{method="PL-LOO"}) using the \code{mvgls} or \code{mvols} function. 

General Linear Hypothesis of the form:

\deqn{\bold{LB=O}}

is used internally with an \bold{L} matrix specifying linear combinations ("contrasts") of the model coefficients (\bold{B}) for each pairwise comparisons. The right-hand-side matrix \bold{O} is a constant matrix (of zeros by default) that can be provided through the argument \code{rhs} (to test specific values for instance).

Permutations on high-dimensional datasets is time consuming. You can use the option \code{nbcores} to parallelize the computations over several cores using forking in UNIX platforms (default is \code{nbcores=1L}). Estimated time to completion is displayed when \code{verbose=TRUE}.
}

\note{
	For PL methods, only the "RidgeArch" penalty is allowed for now. Due to corrections for multiple comparisons, one should ensure that the number of permutations is large enough.}
	
\references{
Clavel, J., Aristide, L., Morlon, H., 2019. A Penalized Likelihood framework for high-dimensional phylogenetic comparative methods and an application to new-world monkeys brain evolution. Systematic Biology 68(1): 93-116.

Clavel, J., Morlon, H. 2020. Reliable phylogenetic regressions for multivariate comparative data: illustration with the MANOVA and application to the effect of diet on mandible morphology in phyllostomid bats. Systematic Biology 69(5): 927-943.
}

\author{J. Clavel}

\seealso{
\code{\link{mvgls}},
\code{\link{mvols}},
\code{\link{pairwise.contrasts}},
\code{\link{manova.gls}}
}

\examples{
\donttest{
data("phyllostomid")

# model fit with mandible~"grp2"
fit <- mvgls(mandible~grp2, data=phyllostomid, phyllostomid$tree, model="lambda", method="PL") 

# pairwise tests 
pairwise.glh(fit, term="grp2", test="Pillai", adjust="holm", nperm=1000, verbose=TRUE)

# fit the model by ML (p<n) and use parametric tests
fitb <- mvgls(mandible[,1:5]~grp2, data=phyllostomid, 
              phyllostomid$tree, model="lambda", method="LL") 
pairwise.glh(fitb, term="grp2", test="Pillai", adjust="holm", verbose=TRUE)

# use permutations on ML fit
pairwise.glh(fitb, term="grp2", test="Pillai", adjust="holm", nperm=1000,  parametric=FALSE)
}
}

\keyword{ Contrasts }
\keyword{ GLS }
\keyword{ OLS }
\keyword{ Pairwise }
\keyword{ MANOVA }
\keyword{ Multivariate tests }
\keyword{ General Linear Hypothesis }
\keyword{ Models comparison }
\keyword{ Regularization }
\keyword{ Penalized likelihood }
\keyword{ High dimensions }% __ONLY ONE__ keyword per line
