\name{pairwise.contrasts}
\alias{pairwise.contrasts}

\title{
Pairwise contrasts
}
\description{
Generates pairwise contrasts (for factor levels) from an object fitted by the \code{mvgls} function. This function is used internally in \code{pairwise.glh} for the generalized linear hypothesis tests (see also ?manova.gls).}

\usage{

pairwise.contrasts(object, term=1, ...)
  
}
\arguments{
  \item{object}{
A model fit obtained by the \code{mvgls} function.
}
  \item{term}{
The factor term in the "object" model fit for which all the pairwise contrasts are built.
}
 
 \item{...}{
Further arguments to be passed through. Not used.}
}


\value{
	Returns a matrix of contrasts for all the pairwise comparisons between levels of "term".
  }


\note{
	The function assumes "effect" dummy coding for the factor levels (see ?contr.treatment)}
	
\references{
Clavel, J., Morlon, H. 2020. Reliable phylogenetic regressions for multivariate comparative data: illustration with the MANOVA and application to the effect of diet on mandible morphology in phyllostomid bats. Systematic Biology 69(5): 927-943.
}

\author{J. Clavel}

\seealso{
\code{\link{mvgls}},
\code{\link{manova.gls}},
\code{\link{pairwise.glh}}
}

\examples{
\donttest{
data("phyllostomid")

# model fit with mandible~"grp2"
fit <- mvgls(mandible~grp2, data=phyllostomid, phyllostomid$tree, model="lambda", method="PL") 

# pairwise tests 
pairwise.contrasts(fit, term="grp2")
}
}

\keyword{ Contrasts }
\keyword{ GLS }
\keyword{ Pairwise }
\keyword{ MANOVA }
\keyword{ Multivariate tests }
\keyword{ General Linear Hypothesis }
\keyword{ Models comparison }
\keyword{ Regularization }
\keyword{ Penalized likelihood }
\keyword{ High dimensions }% __ONLY ONE__ keyword per line
