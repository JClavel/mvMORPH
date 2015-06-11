# mvMORPH
mvMORPH: an R package for fitting multivariate evolutionary models to morphometric data    

This package allows the fitting of multivariate evolutionary models (Ornstein-Uhlenbeck, Brownian motion, Early burst, Shift models).
It also provides functions to compute log-likelihood of users specified models with fast methods (*e.g.*, for Bayesian approaches or customized comparative methods), simulates correlated traits under various models, constrain various parts of multivariate models...

The package is designed to handle ultrametric and non-ultrametric trees (*i.e.* with fossil species) and missing data in multivariate datasets (NA values), SIMMAP mapping of discrete traits, measurement error, etc...

**mvMORPH 1.0.6**

1. This is the developmental version 1.0.6:
  + Allows missing values in multivariate datasets (NA)
  + Package vignette

2. _TODO_:
  + Incorporation of a tests-suite

The current stable version of the mvMORPH package (1.0.5) is on the CRAN repository.
[http://cran.r-project.org/web/packages/mvMORPH/index.html](http://cran.r-project.org/web/packages/mvMORPH/index.html)

##**Package Installation**

You can download the current beta version (1.0.6) binaries for Windows and Mac OS X from the [release page](https://github.com/JClavel/mvMORPH/releases)

##**Report an issue**
Any bugs encountered when using the package can be reported [here](https://github.com/JClavel/mvMORPH/issues)

##**Package citation**

**Clavel, J., Escarguel, G., Merceron, G. 2015.** mvMORPH: an R package for fitting multivariate evolutionary models to morphometric data. Methods in Ecology and Evolution, DOI: 10.1111/2041-210X.12420
