\name{phyllostomid}
\alias{phyllostomid}
\docType{data}
\title{
Phylogeny and trait data for a sample of Phyllostomid bats
%%   ~~ data name/kind ... ~~
}
\description{
Phylogeny, diet, and morphological variables for 49 species of Phyllostomid bats.
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("phyllostomid")}

\details{
Illustrative phylogeny and morphological data (73 variables composed of the superimposed procrustes 2D-coordinates for the mandible and the condylobasal length) of 49 species of Phyllostomid bats from Monteiro & Nogueira (2011). The firsts 22 coordinates represent anatomical landmarks and the last 50 coordinates are semilandmarks. The four grouping factor variables are the adaptive regime models for association between mandible morphology and diet considered in Monteiro & Nogueira (2011).
%%  ~~ If necessary, more details than the __description__ above ~~
}

\references{
Monteiro L.R., Nogueira M.R. 2011. Evolutionary patterns and processes in the radiation of phyllostomid bats. BMC Evolutionary Biology. 11:1-23.

Clavel, J., Morlon, H. (submitted). Reliable phylogenetic regressions for multivariate comparative data: illustration with the MANOVA and application to the effect of diet on mandible morphology in phyllostomid bats.

%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(phyllostomid)
plot(tree)
head(data$mandible)

\donttest{
# Fit a linear model by PL
fit1 <- mvgls(mandible~grp1, data=data, tree, model="lambda", method="LOO") 

# regularized MANOVA test
(manova.gls(fit1, test="Wilks", verbose=TRUE))
}
}

\keyword{mvgls}
\keyword{manova.gls}
\keyword{datasets}
\keyword{bats}