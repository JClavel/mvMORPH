# mvMORPH
mvMORPH: an R package for fitting multivariate evolutionary models to morphometric data    

This package allows the fitting of multivariate evolutionary models (Ornstein-Uhlenbeck, Brownian motion, Early burst, Shift models) on species trees and time series.
It also provides functions to compute log-likelihood of users specified models with fast methods (*e.g.*, for Bayesian approaches or customized comparative methods), simulates correlated traits under various models, constrain various parts of multivariate models...

The package implement now efficient methods for high-dimensional multivariate comparative methods (mvgls) based on Penalized likelihood as well as associated tests (Wilks, Pillai...)

The package is designed to handle ultrametric and non-ultrametric trees (*i.e.* with fossil species) and missing data in multivariate datasets (NA values), SIMMAP mapping of discrete traits, measurement error, etc...

See the packages vignettes for details and examples: browseVignettes("mvMORPH").

**mvMORPH 1.1.5**

1. This is the version 1.1.5:
  + mvqqplot for model diagnostics
  + multivariate association
  + pairwise comparison
  + repeated measures design
  + DFA
  

2. _TODO_:
  + Incorporation of a tests-suite
  + Implement the sampler (upcomming mvMORPH) 
  + Code improvements
  + Extend the shift model to TS
  + Improved mvOU model
  + Threshold model for categorical data

The current stable version of the mvMORPH package (1.1.4) is on the CRAN repository.
[https://cran.r-project.org/package=mvMORPH](https://cran.r-project.org/package=mvMORPH)

## **Package Installation**

You can install the package directly from gitHub through devtools:

```
library(devtools)

install_github("JClavel/mvMORPH", build_vignettes = TRUE)

```


(The installation may crash if your dependencies are not up to date. Note that you may also need to install Rtools to compile the C codes included in the package. For [Windows] (https://cran.r-project.org/bin/windows/Rtools/) and for [Mac] (https://mac.r-project.org/) (and [Tools] (https://mac.r-project.org/tools/) )

## **Report an issue**
Any bugs encountered when using the package can be reported [here](https://github.com/JClavel/mvMORPH/issues)

## **Package citation**

**Clavel, J., Escarguel, G., Merceron, G. 2015.** mvMORPH: an R package for fitting multivariate evolutionary models to morphometric data. Methods in Ecology and Evolution, 6(11):1311-1319.

