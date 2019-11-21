\name{mvLL}
\alias{mvLL}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multivariate (and univariate) algorithms for log-likelihood estimation of arbitrary covariance matrix/trees
%%  ~~function to do ... ~~
}
\description{
This function allows computing the log-likelihood and estimating ancestral states of an arbitrary tree or variance-covariance matrix with differents algorithms based on GLS (Generalized Least Squares) or Independent Contrasts. Works for univariate or multivariate models. Can be wrapped for maximizing the log-likelihood of user-defined models.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mvLL(tree, data, error = NULL, method = c("pic", "rpf", "sparse", "inverse",
    "pseudoinverse"), param = list(estim = TRUE, mu = 0, sigma = 0, D = NULL,
    check = TRUE), precalc = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree}{
A phylogenetic tree of class "phylo" or a variance-covariance matrix (vcv) of that tree (or time-series).
%%     ~~Describe \code{tree} here~~
}
  \item{data}{
Matrix or data frame with species in rows and continuous traits in columns. NA values are allowed with the "rpf", "inverse" and "pseudoinverse" methods.
%%     ~~Describe \code{data} here~~
}
  \item{error}{
Matrix or data frame with species in rows and continuous trait sampling variance (squared standard errors) in columns.
%%     ~~Describe \code{error} here~~
}
  \item{method}{
Method used for computing the log-likelihood. Could be "pic", "sparse", "rpf", "inverse", or "pseudoinverse". See details below.
%%     ~~Describe \code{method} here~~
}
  \item{param}{
List of additional arguments to be passed through the function. The "estim", "mu" and "sigma" arguments are only used with the "pic" method. The "D" argument is used with the others to specify the design matrix. See details below.
%%     ~~Describe \code{param} here~~
}
  \item{precalc}{
Optional. Object of class "precalc.mvmorph". See ?mv.Precalc for details.
%%     ~~Describe \code{precalc} here~~
}
}
\details{
The mvLL function computes the log-likelihood and the ancestral states (mean at the root-theta) for an arbitrary variance-covariance matrix (or trees for the prunning algorithm based on independent contrasts "pic") provided by the user. This function can be wrapped for optimizing various multivariate models of trait evolution (by transforming the branch lengths of a tree for the "pic" method, or feeding it with variance-covariance and design matrices for the other methods).

Five methods are proposed to compute the log-likelihood:

-"pic" is a very fast prunning algorithm based on independent contrasts which should be used with strictly dichotomic trees (i.e., no polytomies). This method can neither be used with measurement errors nor for multiple ancestral states estimation (theta values).

-"rpf" is a GLS algorithm using the rectangular packed format Cholesky factorization for solving the linear system without computing the inverse of the variance-covariance matrix and its determinant (normally used in the loglikelihood estimation). This algorithm uses fast BLAS 3 routines with half storage in packed format for computing the Cholesky upper factor. This method is more efficient than the "inverse" method and can be used with dense matrices (no zero entries).

-"sparse" is a GLS algorithm using Cholesky factorization for sparse matrices (including zero entries). The matrices are stored in the "old Yale sparse format" internally. Depending on the sparsity structure of the variance-covariance matrix this algorithm can be more efficient than the "rpf" method.

-"inverse" is a GLS algorithm that uses explicit inversion of the variance-covariance matrix (through QR decomposition) as well as computation of its determinant in the log-likelihood estimation. This is the "textbook" method, that is computationally more intensive than the previous approaches.

-"pseudoinverse" is a GLS method that uses a generalized inverse (through SVD) for computing the log-likelihood. This method is safer when the matrix is near singularity, but it is the most time-consuming.

The user must provide a variance-covariance matrix (e.g., vcv.phylo(tree)) or a multivariate variance-covariance matrix (e.g., kronecker(matrix(c(2,1,1,1.5),2),vcv.phylo(tree)) as well as a design matrix (or multivariate design matrix) with the "rpf", "sparse", "inverse", and "pseudoinverse" methods.

Use the "param" list of arguments to define whether or not the brownian rate should be estimated and returned (estim=TRUE) with the "pic" method. Otherwise, the rate parameter (also called sigma) is fixed to 1. The arguments "mu" and "sigma" can be used to specify (e.g., in a MCMC setting) the mean at the root and the (squared) brownian rate, respectively.

You can choose to provide differently scaled trees for multivariate data with the "pic" method. In such a case, the trees (one per trait) should be embedded within a list() object. See example below.

%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{logl }{Estimated log-likelihood for the data with the given matrix/tree.}
\item{theta }{Estimated ancestral states at the root. They are defined by the design matrix (D) for all the methods but "pic".}
\item{sigma }{Estimated (squared) rate parameters. Only when param$estim=TRUE with the "pic" method.}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Andersen B. S., Wasniewski J., Gustavson F. G. 2001. A recursive formulation of Cholesky
factorization of a matrix in packed storage. ACM Trans. Math. Soft. 27:214-244.

Clavel J., Escarguel G., Merceron G. 2015. mvMORPH: an R package for fitting multivariate evolutionary models to morphometric data. Methods Ecol. Evol. 6(11):1311-1319.

Freckleton R.P. 2012. Fast likelihood calculations for comparative analyses. Methods Ecol. Evol. 3:940-947.

Golub G.H., Van Loan C.F. 2013. Matrix computations. Baltimore: The John Hopkins University Press.

Gustavson, F.G., Wasniewski, J., Dongarra, J.J., Langou, J. 2010. Rectangular full packed format for Cholesky's algorithm: factorization, solution and inversion. ACM Trans. Math. Soft., 37:1-33.

%% ~put references to the literature/web site here ~
}
\author{
Julien Clavel
%%  ~~who you are~~
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{mvMORPH}}
\code{\link{mvgls}}
\code{\link{mvOU}}
\code{\link{mvEB}}
\code{\link{mvBM}}
\code{\link{mvSHIFT}}
\code{\link{mvSIM}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

## Simulated dataset
set.seed(14)
# Generating a random tree with 50 tips
n=50
tree<-pbtree(n=n)

# Simulated trait
data=rTraitCont(tree)

# Design matrix
D=matrix(rep(1,n),ncol=1)

## Compute the log-likelihood
# Inverse
mvLL(vcv.phylo(tree),data,method="inverse",param=list(D=D))

# Pseudoinverse
mvLL(vcv.phylo(tree),data,method="pseudoinverse",param=list(D=D))

# Sparse
mvLL(vcv.phylo(tree),data,method="sparse",param=list(D=D))

# RPF
mvLL(vcv.phylo(tree),data,method="rpf",param=list(D=D))

# Pic
mvLL(tree,data,method="pic",param=list(estim=TRUE))

# Pic with arbitrary values
mvLL(tree,data,method="pic",param=list(estim=FALSE, mu=0, sigma=1))
mvLL(tree,data,method="pic",param=list(estim=FALSE))
mvLL(tree,data,method="pic",param=list(estim=FALSE, sigma=1)) # similar to mu=NULL

# Arbitrary value for mu with other methods (similar to mu=0 and sigma=1 with "pic")
mvLL(vcv.phylo(tree),data,method="rpf",param=list(D=D, estim=FALSE, mu=0)) 

## Multivariate cases
# Simulate traits
data2<-mvSIM(tree,nsim=1,model="BM1",param=list(sigma=diag(2),theta=c(0,0),ntraits=2))
# Design matrix
D<-cbind(rep(c(1,0),each=50),rep(c(0,1),each=50))

# RPF
mvLL(kronecker(diag(2),vcv.phylo(tree)),data2,method="rpf", param=list(D=D))

# Inverse (with default design matrix if not provided)
mvLL(kronecker(diag(2),vcv.phylo(tree)),data2,method="inverse")

# Pic
mvLL(tree,data2,method="pic")
# NB: The trees in the list could be differently scaled for each traits...
mvLL(list(tree,tree),data2,method="pic")

## VERY FAST COMPUTATION FOR LARGE TREES (take few seconds)


\donttest{
# Big tree (1,000,000 species) - It's the time consuming part...
 tree2<-rtree(1000000)
# Simulate trait with a Brownian motion process
 trait<-rTraitCont(tree2)
 system.time(mvLL(tree2,trait,method="pic",param=list(estim=FALSE, sigma=1)))

 precal<-mv.Precalc(tree2,nb.traits=1, param=list(method="pic"))
 system.time(mvLL(tree2,trait,method="pic",param=list(estim=FALSE, sigma=1),
   precalc=precal))

# Check=FALSE !! Your tree should be in post-order !!
 tr2<-reorder(tree2,"postorder")
 system.time(mvLL(tr2,trait,method="pic",param=list(estim=FALSE, sigma=1, check=FALSE)))
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Loglikelihood }
\keyword{ Methods }
\keyword{ GLS }
\keyword{ Independent contrasts }% __ONLY ONE__ keyword per line
