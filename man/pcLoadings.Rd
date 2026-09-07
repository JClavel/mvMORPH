\name{pcaLoadings}
\alias{pcaLoadings}

\title{
Plot the loadings from Principal Component Analysis (PCA) or Probabilistic Phylogenetic Principal Component Analysis (P3CA)
}
\description{
The function extracts and plot the PCs loadings obtained by \code{p3ca} or \code{mvgls.pca} to help interpreting P3CA or phylogenetic PCA.
}

\usage{

pcaLoadings(object, ...)
  
}
\arguments{
  \item{object}{
A PCA fit obtained by the \code{p3ca} or \code{mvgls.pca} function.
}

 \item{...}{
Further options.  For instance, \code{color} is used to define the color ramp palette used for the plot (default is \code{c('red','white','blue4')}).
The argument \code{scale} is used to scale the bubble plot (default is 3). The argument \code{horizontal=TRUE} (default) plot the PCs on the 'x' axis.}
}

\details{
The function will plot (using a bubble plot) the correlations (loadings) of each variable used in a P3CA to each of the "q" P3Cs axes obtained by \code{p3ca} function. 
This can be used to display main (evolutionary) changes associated with each P3Cs axes (see Montoya et al. 2026; Joseph et al. 2026).
}

\value{
	a plot showing the correlation (loadings) of each variable with each PC axis.
  }


	
\references{
Montoya P., Joseph J., Goswami A., Morlon H., Clavel J. 2026. A probabilistic and phylogenetic principal component analysis for modelling high-dimensional trait evolution. doi.org/10.64898/2026.05.27.728209.

Joseph J., Montoya P., Baudat F., Ruiz-Herrera A., Clavel J. 2026. Explaining the rapid evolution of mammalian meiotic recombination proteins. doi.org/10.64898/2026.06.01.729254. 
}

\author{J. Clavel}

\seealso{
\code{\link{p3ca}},
\code{\link{pcaShape}},
\code{\link{mvgls.pca}}
}

\examples{
\donttest{

set.seed(2508)

# Loading the data
data(phyllostomid)
phyllos_data = phyllostomid$mandible[,-1]
phyllos_tree = phyllostomid$tree

# Perfoming the P3CA on the complete dataset - Analytical solution
p3ca_phyllos = p3ca(phyllos_data, phyllos_tree, q=5, model='BM')

# plot the loadings
pcaLoadings(p3ca_phyllos)
}
}

\keyword{ PCA loadings }
\keyword{ P3CA }% __ONLY ONE__ keyword per line
