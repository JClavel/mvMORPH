## mvMORPH 1.1.4
    + model BMM in mvgls
    + "predict" function for mvgls
    + plot function for mvgls objects
    + dfa on gls fit (beta)
    + fix typo in mvSHIFT output print 
    + fix error in estimation of OU with n-ultrametric trees in mvgls + new adjusted bounds for parameters search
    + fix typo in summary.mvgls print option in the calculation of the AIC.
    + handling of design matrices with deficient ranks in regressions
## mvMORPH 1.1.3
    + CRAN request to remove the export statements for S3 classes
    + replace is.binary.tree to is.binary.phylo, the former being deprecated from "ape".
    + fix error in mvSHIFT. The wrong values were returned (but not printed) for "beta" in BMEB models.
## mvMORPH 1.1.2
    + replaced "F" by "FALSE" in example files to follow CRAN policies
    + optimization of some diagonal matrices computations
    + bug fix in mvSIM (SHIFT model without simmap tree provided)
## mvMORPH 1.1.1
    +  update help pages of various functions
    +  Bugs fixes for bounds in the parameter search in "mvgls", and missing values estimation in "estim" with OU1 model.
    +  manova.gls: MANOVA methods and multivariate statistics based on ML (n<p) and PL (any p) model fit.
    +  Phyllostomid dataset from Monteiro & Nogueira (2011)
    +  EIC: function to compute the EIC score for both ML and PL techniques (class mvgls)
## mvMORPH 1.1.0
    + mvgls: a multivariate GLS function based on ML (n<p) and PL (any p).
    + mvgls.pca: a function to perform PCA on serially correlated data.
    + GIC: function to compute the GIC score for both ML and PL techniques (class mvgls)
    + suite of functions to use with mvgls
    + Bug fix for "mvSIM" with sparse methods and speed-up for "BM1" model.
    + Typo fix in testLRT
    + Typo fix in weight-matrix-mvmorph.c (indexing of variable)
    + add options for the pruning algorithm.
    + fix error with NA values for the mvOU function.
## mvMORPH 1.0.9
    + Registration of compiled codes
    + New pruning algorithm to compute the determinant and square root matrix
    + Bugs fixes in "estim" for ancestral state reconstructions with OU models and "smean=FALSE" option.
    + Bug fix for "estim" with univariate data.
    + Corrected bug in mvSIM.
    + Fix bug with matrix parametrization "equaldiagonal".
    + add the "svd" option for simulating traits in mvSIM.
    + Bug fix in "mvEB".
## mvMORPH 1.0.8
    + Allows estimating the missing cases (NA)
    + Allows estimating trends
    + User defined constrained models and parameterizations  
    + Return the log-likelihood function
    + Simulating traits on package vignette
    + Multivariate models for time-series (TS)
    + Partial implementation of a tests-suite
