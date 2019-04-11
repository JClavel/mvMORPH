## mvMORPH 1.1.1
   + manova.gls: MANOVA methods and multivariate statistics based on ML (n<p) and PL (p>n) model fit.
   + EIC: function to compute the EIC score for both ML and PL techniques (class mvgls)
   + update help pages of various functions
## mvMORPH 1.1.0
    + mvgls: a multivariate GLS function based on ML (n<p) and PL (p>n).
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
