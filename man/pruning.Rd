\name{pruning}
\alias{pruning}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Pruning algorithm to compute the square root of the phylogenetic covariance matrix and it's determinant.
%%  ~~function to do ... ~~
}
\description{
This function use the pruning algorithm (Felsenstein 1973) to efficiently compute the determinant of the phylogenetic covariance matrix as well as the square root of this matrix (or it's inverse; Stone 2011, Khabbazian et al. 2016). This algorithm is faster than using "eigen" or "cholesky" function to compute the determinant or the square root (see e.g., Clavel et al. 2015) and can be used to compute the log-likelihood of a model in linear time.

%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
pruning(tree, inv=TRUE, check=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{ 

\item{tree}{
Phylogenetic tree (an object of class "phylo" or "simmap").
}

\item{inv}{
Return the transpose of the matrix square root of either the covariance matrix (inv=FALSE) or it's inverse (inv=TRUE, the default).
}

\item{check}{
Check if the input tree is dichotomous and in "postorder" (see ?is.binary.tree and ?reorder.phylo).
}

}
\details{
The tree is assumed to be fully dichotomic and in "postorder", otherwise the functions \emph{multi2di} and \emph{reorder.phylo} are used internally when \emph{check=TRUE}. Note that the transpose of the square root is returned.

%%  ~~ If necessary, more details than the description above ~~
}
\value{

\item{sqrtMat  }{The matrix square root }
\item{varNode  }{Variance associated to each node values (similar to "contrasts" variance)}
\item{varRoot  }{Variance associated to the root value (similar to the ancestral state variance)}
\item{det      }{Log-determinant of the phylogenetic covariance of the tree}

%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Clavel J., Escarguel G., Merceron G. 2015. mvMORPH: an r package for fitting multivariate evolutionary models to morphometric data. Methods Ecol. Evol. 6:1311-1319.

Felsenstein J. 1973. Maximum-likelihood estimation of evolutionary trees from continuous characters. Am. J. Hum. Genet. 25:471-492.

Khabbazian M., Kriebel R., Rohe K., Ane C. 2016. Fast and accurate detection of evolutionary shifts in Ornstein-Uhlenbeck models. Methods Ecol. Evol. 7:811-824.

Stone E.A. 2011. Why the phylogenetic regression appears robust to tree misspecification. Syst. Biol. 60:245-260

%% ~put references to the literature/web site here ~
}
\author{
Julien Clavel
%%  ~~who you are~~
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{mvLL}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}

\examples{

## Simulated dataset
set.seed(14)
# Generating a random tree
tree<-pbtree(n=50)
Y <- mvSIM(tree, model="BM1", param=list(sigma=1, theta=0)) # trait
X <- matrix(1, nrow=Ntip(tree), ncol=1) # design matrix

## Use the GLS trick
# Compute the matrix square root
C <- vcv.phylo(tree)
D <- chol(C)
Cinv <- solve(C)
Di <- chol(Cinv)

# transform the traits
Xi <- Di\%*\%X
Yi <- Di\%*\%Y

# Compute the GLS estimate and determinant (see Clavel et al. 2015)
# GLS estimate for the root
print(pseudoinverse(Xi)\%*\%Yi)

# Determinant of the phylogenetic covariance matrix
print(sum(log(diag(D)^2)))    


## Use the pruning algorithm (much faster)

M <- pruning(tree, inv=TRUE)

Xi <- t(M$sqrtMat)\%*\%X
Yi <- t(M$sqrtMat)\%*\%Y

# GLS estimate
print(pseudoinverse(Xi)\%*\%Yi)

# determinant
print(M$det)

## REML determinant (without variance of the root state; see Felsenstein 1973)
# full REML
log(det(C)) + log(det(t(X)\%*\%Cinv\%*\%X))

# pruning REML
sum(log(M$varNode))

}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Matrix square root }
\keyword{ Determinant }
\keyword{ GLS }
\keyword{ Independent contrasts }% __ONLY ONE__ keyword per line
