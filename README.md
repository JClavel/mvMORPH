# mvMORPH
mvMORPH: an R package for fitting multivariate evolutionary models to morphometric data    

This package allows the fitting of multivariate evolutionary models (Ornstein-Uhlenbeck, Brownian motion, Early burst, Shift models) on species trees and time series.
It also provides functions to compute log-likelihood of users specified models with fast methods (*e.g.*, for Bayesian approaches or customized comparative methods), simulates correlated traits under various models, constrain various parts of multivariate models...

The package is designed to handle ultrametric and non-ultrametric trees (*i.e.* with fossil species) and missing data in multivariate datasets (NA values), SIMMAP mapping of discrete traits, measurement error, etc...

See the packages vignettes for details and examples: browseVignettes("mvMORPH").

**mvMORPH 1.0.7**

1. This is the version 1.0.7:
  + Allows estimating the missing cases (NA)
  + Allows estimating trends
  + User defined constrained models and parameterizations  
  + Return the log-likelihood function
  + Package vignette 2
  + Multivariate models for time-series (TS)
  + Partial implementation of a tests-suite

2. _TODO_:
  + Incorporation of a tests-suite
  + Implement the sampler (upcomming mvMORPH 1.0.8) 
  + Code improvements
  + Extend the shift model to TS
  + Formula option for independent variables

The current stable version of the mvMORPH package (1.0.6) is on the CRAN repository.
[https://cran.r-project.org/package=mvMORPH](https://cran.r-project.org/package=mvMORPH)

##**Package Installation**

You can download the current beta version (1.0.7) binaries for Windows and Mac OS X from the [release page](https://github.com/JClavel/mvMORPH/releases)

You can also install it directly from gitHub through devtools:

library(devtools)

install_github("JClavel/mvMORPH", build_vignettes = TRUE)

(The installation may crash if your dependencies are not up to date. Note that you may also need to install Rtools to compile the C codes included in the package. For [Windows] (https://cran.r-project.org/bin/windows/Rtools/) and for [Mac] (http://r.research.att.com) (and [Tools] (https://r.research.att.com/tools/) )

##**Report an issue**
Any bugs encountered when using the package can be reported [here](https://github.com/JClavel/mvMORPH/issues)

##**Package citation**

**Clavel, J., Escarguel, G., Merceron, G. 2015.** mvMORPH: an R package for fitting multivariate evolutionary models to morphometric data. Methods in Ecology and Evolution, 6(11):1311-1319.    DOI: 10.1111/2041-210X.12420

[Download version with appended supplementary material.](http://www.researchgate.net/publication/277711429_mvMORPH_an_R_package_for_fitting_multivariate_evolutionary_models_to_morphometric_data)
