\name{pcaShape}
\alias{pcaShape}

\title{
Projection of 2D and 3D shapes (from geometric morphometric datasets) on Principal Component Axes (PCA)
}
\description{
The function extracts the shape changes along PCA (or phylogenetic PCA) axes computed from a linear model obtained by \code{mvgls} or \code{mvols}.
}

\usage{

pcaShape(object, axis=1, ndim=3, spp=NULL, plot=FALSE, ...)
  
}
\arguments{
  \item{object}{
A model fit obtained by the \code{mvgls} or \code{mvols} function.
}
  \item{axis}{
PC axis on which the shape data will be projected.}

 \item{ndim}{
The number of dimensions of the GMM data set (2 for 2D and 3 for 3D).
}

 \item{spp}{
Names of the species (should match names in the dataset) shape to project onto the PC axis. If null, the two extreme shapes along \code{axis} are reported. 
}

 \item{plot}{
Should the projected landmarks be plotted?
}

 \item{...}{
Further options.}
}

\details{
The function will project the shape changes along PC axes obtained from models fit by \code{mvgls} (or \code{mvols}).  This can be used to display main (evolutionary) morphological changes for 2D and 3D geometric morphometric data. See for instance Clavel et al. 2019.
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
\code{\link{mvgls.pca}}
\code{\link{dfaShape}}

}

\examples{
\donttest{
data(phyllostomid)

# Fit a linear model by PL
fit <- mvgls(mandible[,-1]~grp1, data=phyllostomid, phyllostomid$tree, model="lambda", method="PL") 

# Project the mandible shape of Ametrida bat genus 
proj_shape <- pcaShape(fit, axis=1, ndim=2, spp="Ametrida", plot=TRUE)

polygon(proj_shape$Ametrida)

}
}

\keyword{ PCA projections }
\keyword{ GLS }
\keyword{ OLS }
\keyword{ Geometric Morphometrics }
\keyword{ High dimensions }% __ONLY ONE__ keyword per line
