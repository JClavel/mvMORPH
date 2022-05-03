\name{mvOU}
\alias{mvOU}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multivariate Ornstein-Uhlenbeck model of continuous traits evolution
%%  ~~function to do ... ~~
}
\description{
This function allows the fitting of a multivariate Ornstein-Uhlenbeck (OU1) model with possibly multiple optima (OUM) for different "selective regimes". A "phylo" object with SIMMAP-like mapping of ancestral states is required to subdivise the tree (or branches) into "selective regimes".
Species measurement errors or dispersions can also be included in the model.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mvOU(tree, data, error = NULL, model = c("OUM", "OU1"), param = list(sigma = NULL,
    alpha = NULL, vcv = "fixedRoot", decomp = c("cholesky","spherical","eigen","qr",
    "diagonal","upper","lower")), method = c("rpf", "sparse", "inverse",
    "pseudoinverse", "univarpf"), scale.height = FALSE, optimization = c("L-BFGS-B",
    "Nelder-Mead", "subplex"), control = list(maxit = 20000), precalc = NULL, 
    diagnostic = TRUE, echo = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree}{
Phylogenetic tree with mapped ancestral states in phytools' SIMMAP format for "OUM" model.
(See make.simmap, paintBranches, paintSubTree, and make.era.map functions from the phytools package). A "phylo" object can be used with model "OU1".
%%     ~~Describe \code{tree} here~~
}
  \item{data}{
Matrix or data frame with species in rows and continuous traits in columns. NA values are allowed with the "rpf", "inverse", and "pseudoinverse" methods.
%%     ~~Describe \code{data} here~~
}
  \item{error}{
Matrix or data frame with species in rows and continuous trait sampling variance (squared standard errors) in columns.
%%     ~~Describe \code{error} here~~
}

  \item{model}{
Choose between "OUM" for a multiple selective regime model, or "OU1" for a unique selective regime for the whole tree.
%%     ~~Describe \code{model} here~~
}
  \item{param}{
List of arguments to be passed to the function. See details below.
%%     ~~Describe \code{param} here~~
}
  \item{method}{
Choose between "rpf", "sparse", "inverse", "pseudoinverse", or "univarpf" for computing the log-likelihood during the fitting process. See details below.
%%     ~~Describe \code{method} here~~
}
  \item{scale.height}{
Whether the tree should be scaled to unit length or not.
%%     ~~Describe \code{scale.height} here~~
}
  \item{optimization}{
Methods used by the optimization routines (see ?optim and ?subplex for details). The "fixed" method returns the log-likelihood function only.
%%     ~~Describe \code{optimization} here~~
}
  \item{control}{
Max. bound for the number of iteration of the optimizer; other options can be fixed in the list. (See ?optim or ?subplex for details).
%%     ~~Describe \code{control} here~~
}
  \item{precalc}{
Optional. Precalculation of fixed parameters. See ?mvmorph.Precalc for details.
%%     ~~Describe \code{precalc} here~~
}
  \item{diagnostic}{
 Whether the convergence diagnostics should be returned or not.
%%     ~~Describe \code{diagnostic} here~~
}
  \item{echo}{
Whether the results must be returned or not.
%%     ~~Describe \code{echo} here~~
}
}
\details{
The mvOU function fits a multivariate model of evolution according to an Ornstein-Uhlenbeck process:

\deqn{dX(t) = A(\Theta - X(t))dt + \Sigma^{1/2}dW(t)}

The user can incorporate measurement errors and uses SIMMAP-like mapping of ancestral states (phytools objects of class "simmap"). SIMMAP mapping allows one to assign parts of branches to different selective regimes for estimating regime-specific evolutionary optima. See the package vignette: browseVignettes("mvMORPH").

Mapping of ancestral states can be done using the "make.simmap", "make.era.map" or "paintSubTree" functions from the "phytools" package.

The "method" argument allows the user to try different algorithms for computing the log-likelihood. The \code{"rpf"}, \code{"univarpf"} (for univariate analysis) and \code{"sparse"} methods use fast GLS algorithms based on factorization to avoid explicit computation of the inverse of the variance-covariance matrix and its determinant during log-likelihood estimation. The \code{"inverse"} approach uses a "stable" (based on QR decomposition) explicit computation of the inverse and determinant of the matrix and is therefore slower. The \code{"pseudoinverse"} method uses a generalized inverse that is safer for matrix near singularity but highly time consuming. See ?mvLL for details.

Arguments in the \bold{"param"} \code{list} are:

\bold{"sigma"} or \bold{"alpha"} - Starting values for the likelihood search can be specified through the "alpha" and "sigma" arguments in the param list. It is also possible to test for the significance of the off-diagonal sigma (scatter) and alpha (drift) matrix in the full model by making comparison with a constrained model (using sigma="constraint", or alpha="constraint") in the "param" argument list. You can also provide starting values for the constrained model. For instance, for two traits use sigma=list("constraint", c(0.5,0.5))  (or alpha=list("constraint", c(0.5,0.5))).

\bold{"decomp"} - You can further constrain the alpha matrix by specifying the decomposition of the matrix through the "decomp" argument in the "param" list. The multivariate Ornstein-Uhlenbeck model is described by the spectral decomposition of the alpha matrix. It is possible to parameterize the alpha matrix to be decomposable using various parameterizations (e.g., on its eigenvalues with different biological interpretations; Sy et al. 1997, Bartoszek et al. 2012; Clavel et al. 2015). For a symmetric matrix parameterization, the user can choose the \code{"cholesky"}, \code{"eigen"}, or \code{"spherical"} option. 

For general square (non-symmetric) matrices the \code{"svd"}, \code{"qr"} and \code{"schur"} parameterizations can be used. The \code{"schur"} parameterization constrains the eigenvalues of the alpha matrix to be real numbers. The \code{"svd+"}, \code{"qr+"} or \code{"eigen+"} options force the eigenvalues to be positives by taking their logarithm. It is also possible to specify \code{"diagonal"} which is similar to the use of the "constraint" argument for "alpha" argument, or to use \code{"equal"} and \code{"equaldiagonal"}. Finally, one can specify that the alpha matrix is \code{"upper"} or \code{"lower"} triangular (i.e., one process affect the other unilateraly). Details can be found in the package vignette: browseVignettes("mvMORPH") and on Clavel et al. 2015.

\bold{"decompSigma"} - The sigma matrix is parameterized by various methods to ensure its positive definiteness (Pinheiro and Bates, 1996). These methods can be accessed through the "decompSigma" argument and are the \code{"cholesky"}, \code{"eigen+"}, and \code{"spherical"} parameterization. The sigma matrix can also be forced to be diagonal using \code{"diagonal"} or \code{"equaldiagonal"} and forced to have the same variances using \code{"equal"}. Details can be found in the package vignette: browseVignettes("mvMORPH").

\bold{"vcv"} - It is possible to specify in the "param" list what kind of variance-covariance matrix to use with the "vcv" argument, depending on how the root is treated. 
The \code{vcv="randomRoot"} option assumes that the value at the root is a random variable with the stationary distribution of the process. It cannot be used with the "sparse" method to speed up the computations. The \code{vcv="fixedRoot"} option assumes that the root is a fixed parameter. On ultrametric trees both approaches should converge on the same results when the OU process is stationary.

\bold{"root"} - This argument allows the user to specify if the ancestral state at the root (theta 0) should be estimated (\code{root=TRUE}), or assumed to be at the oldest regime state stationary distribution (\code{root=FALSE}). An alternative is to follow Beaulieu et al. (2012) and explicitly drop the root state influence (\code{root="stationary"}). The first option should be used with non-ultrametric trees (i.e., with fossil species; e.g., Hansen 1997) where information on the ancestral state is directly available from the data. Note that estimating shifts from the ancestral state to the optimum(s) from extant species can be problematic and it can be preferable to assume each regime optimum(s) to be at the stationary distribution.

For the \bold{"decomp"} and \bold{"decompSigma} arguments, a user-defined matrix with integer values taken as indices of the parameters to be estimated can be provided. See ?mvBM and ?mvRWTS.


Note on the returned Hessian matrix in the result list (param$opt$hessian):

The hessian is the matrix of second order partial derivatives of the likelihood function with respect to the maximum likelihood parameter values. This matrix provides a measure of the steepness of the likelihood surface in the vicinity of the optimum. The eigen-decomposition of the hessian matrix allows assessing the reliability of the model fit (even if the optimizer has converged).
When the optimization function does not converge on a stable result, the user may consider increasing the "maxit" argument in the "control" option, or try a simpler model with fewer parameters to estimate. Changing the starting values ("alpha" and "sigma" options in the param list) as well as the optimizing method ("optimization" option) may helps sometimes (e.g., alpha=runif(3) for a two-trait analysis with random starting values - i.e., the lower triangular alpha matrix). Note that the number of starting values to provide depends on the matrix decomposition chosen for the alpha parameter (p*(p+1)/2 values for symmetric alpha matrix, but p*p values for non-symmetric ones - with p the number of traits).

%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{LogLik }{The log-likelihood of the optimal model.}
\item{AIC }{Akaike Information Criterion for the optimal model.}
\item{AICc }{Sample size-corrected AIC.}
\item{theta }{Estimated ancestral states.}
\item{alpha }{Matrix of estimated alpha values (strength of selection).}
\item{sigma }{Evolutionary rate matrix (drift).}
\item{convergence }{Convergence status of the optimizing function; "0" indicates convergence. (see ?optim for details).}
\item{hess.values }{Reliability of the likelihood estimates calculated through the eigen-decomposition of the hessian matrix. "0" means that a reliable estimate has been reached. See details above.}
\item{param }{List of model fit parameters (optimization, method, model, number of parameters...).}
\item{llik }{The log-likelihood function evaluated in the model fit "$llik(par, root.mle=TRUE)".}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Bartoszek K., Pienaar J., Mostad P., Andersson S., Hansen T.F. 2012. A phylogenetic comparative method for studying multivariate adaptation. J. Theor. Biol. 314:204-215.

Beaulieu J.M., Jhwueng D.-C., Boettiger C., O'Meara B.C. 2012. Modeling stabilizing selection: Expanding the Ornstein-Uhlenbeck model of adaptive evolution. Evolution. 66:2369-2389.

Butler M.A., King A.A. 2004. Phylogenetic comparative analysis: a modeling approach for adaptive evolution. Am. Nat. 164:683-695.

Clavel J., Escarguel G., Merceron G. 2015. mvMORPH: an R package for fitting multivariate evolutionary models to morphometric data. Methods Ecol. Evol. 6(11):1311-1319.

Hansen T.F. 1997. Stabilizing selection and the comparative analysis of adaptation. Evolution. 51:1341-1351.

Pinheiro J.C., Bates D.M. 1996. Unconstrained parameterizations for variance-covariance matrices. Stat. Comput. 6:289-296.

Sy J.P., Taylor J.M.G., Cumberland W.G. 1997. A stochastic model for the analysis of bivariate longitudinal AIDS data. Biometrics. 53:542-555.
%% ~put references to the literature/web site here ~
}
\author{
Julien Clavel
%%  ~~who you are~~
}
\note{
This function partly uses a modified version of the C code from the "OUCH" package built by Aaron King, as well as a C code which is part of the "ape" package by Emmanuel Paradis.  I kindly thank those authors for sharing their sources.
Note that Bartoszek et al. (2012) proposed the mvSLOUCH package also dedicated to multivariate Ornstein-Uhlenbeck processes, which allows fitting regression models with randomly evolving predictor variables.

The "symmetric", "nsymmetric", "symmetricPositive", and "nsymPositive" options for the "decomp" argument are deprecated.
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{mvMORPH}}
\code{\link{mvgls}}
\code{\link{halflife}}
\code{\link{stationary}}
\code{\link{mvBM}}
\code{\link{mvEB}}
\code{\link{mvSHIFT}}
\code{\link{mvOUTS}}
\code{\link{mvRWTS}}
\code{\link{mvSIM}}
\code{\link{LRT}}
\code{\link{optim}}
\code{\link{make.simmap}}
\code{\link{make.era.map}}
\code{\link{paintSubTree}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# Simulated dataset
set.seed(14)
# Generating a random tree
tree<-pbtree(n=50)

# Setting the regime states of tip species
sta<-as.vector(c(rep("Forest",20),rep("Savannah",30))); names(sta)<-tree$tip.label

# Making the simmap tree with mapped states
tree<-make.simmap(tree,sta , model="ER", nsim=1)
col<-c("blue","orange"); names(col)<-c("Forest","Savannah")

# Plot of the phylogeny for illustration
plotSimmap(tree,col,fsize=0.6,node.numbers=FALSE,lwd=3, pts=FALSE)

# Simulate the traits
alpha<-matrix(c(2,0.5,0.5,1),2)
sigma<-matrix(c(0.1,0.05,0.05,0.1),2)
theta<-c(2,3,1,1.3)
data<-mvSIM(tree, param=list(sigma=sigma, alpha=alpha, ntraits=2, theta=theta,
            names_traits=c("head.size","mouth.size")), model="OUM", nsim=1)

## Fitting the models
\donttest{
# OUM - Analysis with multiple optima
 mvOU(tree, data)

# OU1 - Analysis with a unique optimum
 mvOU(tree, data, model="OU1", method="sparse")

# various options
mvOU(tree, data, model="OUM", method="sparse", scale.height=FALSE,
    param=list(decomp="svd", root="stationary"))# non-symmetric alpha
mvOU(tree, data, model="OUM", method="sparse", scale.height=FALSE,
    param=list(decomp="qr", root=TRUE)) # non-symmetric alpha
mvOU(tree, data, model="OUM", method="sparse", scale.height=FALSE,
    param=list(decomp="cholesky", root=TRUE)) # symmetric-positive
# OUCH setting
mvOU(tree, data, model="OUM", method="rpf", scale.height=FALSE,
    param=list(decomp="cholesky", root=FALSE, vcv="ouch"))

## Univariate case - FAST with RPF
 set.seed(14)
 tree<-pbtree(n=500)

# Setting the regime states of tip species
 sta<-as.vector(c(rep("Forest",200),rep("Savannah",300))); names(sta)<-tree$tip.label

# Making the simmap tree with mapped states
 tree<-make.simmap(tree,sta , model="ER", nsim=1)
 col<-c("blue","orange"); names(col)<-c("Forest","Savannah")

# Plot of the phylogeny for illustration
 plotSimmap(tree,col,fsize=0.6,node.numbers=FALSE,lwd=3, pts=FALSE)

# Parameters
 alpha<-2.5
 sigma<-0.1
 theta<-c(0,2)
 data<-mvSIM(tree, param=list(sigma=sigma, alpha=alpha, ntraits=1, theta=theta,
             names_traits=c("body_size")), model="OUM", nsim=1)

# Fit the model
 system.time(mvOU(tree, data, model="OUM", method="univarpf",
                param=list(root="stationary")))
 system.time(mvOU(tree, data, model="OU1", method="univarpf",
                param=list(root="stationary")))

# Add measurement error
 error=rnorm(500,sd=0.1)
 mvOU(tree, data, error=error^2, model="OUM", method="univarpf",
    param=list(root="stationary"))

	}
}% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Ornstein Uhlenbeck }
\keyword{ Hessian }
\keyword{ OU }% __ONLY ONE__ keyword per line
