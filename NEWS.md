## mvMORPH 1.1.9
    + reverse change for zero branch lengths in "penalized" file. A typo was introduced. Now a check is performed first in the "mvgls" function
## mvMORPH 1.1.8
    + add the pcaShape and dfaShape functions
    + add wrapper "simulate()" to simulate from 'mvgls' or 'mvols' fit
    + fix error in "estim" for BM with trend
    + fix issues with zero branch lengths in the pruning algorithm used in 'mvgls'
    + fix typo on factor labels in DFA predictions
    + fix CRAN requests (typo on Rd file and arguments in error function in C code)
## mvMORPH 1.1.7
    + fix error in handling pairwise glh tests with the "effectSize" function
    + fix error in mvols with EIC - no use of the pruning algorithm now
    + add ploting option for pairwise tests
    + add AIC extractor
## mvMORPH 1.1.6
    + fix error in the estimation of GIC with "BMM" model in mvgls. Now the 'mvgls' function uses a more robust parameterization of BMM which ease the computation of GIC.
    + mvols function - wrapper to mvgls to fit OLS (or WLS) multivariate models (possibly regularized)
## mvMORPH 1.1.5
    + mvqqplot (multivariate normality and outliers assessment - beta)
    + effectsize (multivariate measures of association - beta)
    + pairs.contrasts (build a matrix of pairwise contrasts)
    + pairwise.glh (performs multivariate pairwises tests)
    + predict.mvgls.dfa (predict option for DFA - beta)
    + manova.gls (implements now contrasts for repeated measures designs)
    + fix error in estim
    + add "ancestral" function to estimate ancestral states for mvgls objects
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
