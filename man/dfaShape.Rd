\name{dfaShape}
\alias{dfaShape}

\title{
Projection of 2D and 3D shapes (from geometric morphometric datasets) on Discriminant axes
}
\description{
The function extracts the shape changes along discriminant axes computed by a DFA (\code{mvgls.dfa}).
}

\usage{

dfaShape(object, reference, axis=1, ndim=3, spp=NULL, scaling=1, plot=FALSE, ...)  

}
\arguments{
  \item{object}{
A discriminant analysis obtained by the \code{mvgls.dfa} function.
}
  \item{reference}{
The reference shape used to compare the deformations. Usually the mean shape.}

 \item{axis}{
The discriminant axis on which morphological changes are projected.}

 \item{ndim}{
The number of dimensions of the GMM data set (2 for 2D and 3 for 3D).}

 \item{spp}{
Names of the species (should match names in the dataset) shape to project onto the PC axis. If null, the two extreme shapes along \code{axis} are reported.  
}

 \item{scaling}{
An arbitrary factor used to multiply the effects (for better visualization)  
}

 \item{plot}{
Should the projected landmarks be plotted?  
}

 \item{...}{
Further options.}
}

\details{
The function will project the shape changes along discriminant axes obtained from a DFA by \code{mvgls.dfa}.  This can be used to display morphological changes (for 2D and 3D geometric morphometric data) that best separate individuals from distinct groups.
}

\value{
	a list with 2D or 3D coordinates for the shape projected on the selected PC axis.
  }


	
\references{

Claude, J., 2008. Morphometrics with R. Springer Science.

Clavel, J., Aristide, L., Morlon, H., 2019. A Penalized Likelihood framework for high-dimensional phylogenetic comparative methods and an application to new-world monkeys brain evolution. Systematic Biology 68(1): 93-116.

}

\author{J. Clavel}

\seealso{
\code{\link{mvgls}},
\code{\link{mvols}},
\code{\link{mvgls.dfa}}
\code{\link{pcaShape}}
}

\examples{
\donttest{
data(phyllostomid)

# Fit a linear model by PL
fit <- mvgls(mandible[,-1]~grp1, data=phyllostomid, phyllostomid$tree, model="lambda", method="PL") 

# Discriminant analysis
da <- mvgls.dfa(fit)

# Project the mandible shape extremes on the first discriminant axe
proj_shape <- dfaShape(da, reference=coef(fit), axis=1, ndim=2, plot=TRUE)

polygon(proj_shape$min)
polygon(proj_shape$max, border="red")

}
}

\keyword{ DFA projections }
\keyword{ GLS }
\keyword{ OLS }
\keyword{ Geometric Morphometrics }
\keyword{ High dimensions }% __ONLY ONE__ keyword per line
