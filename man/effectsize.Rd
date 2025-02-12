\name{effectsize}
\alias{effectsize}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multivariate measure of association/effect size for objects of class "manova.gls"
%%  ~~function to do ... ~~
}
\description{
This function estimate the multivariate effectsize for all the outcomes variables of a multivariate analysis of variance
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
effectsize(x,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{ 

 \item{x}{
An object of class "manova.gls"
%%     ~~Describe \code{model2} here~~
}
  \item{...}{
One can specify \code{adjusted=TRUE} to obtain Serlin' adjustment to Pillai trace effect size, or Tatsuoka' adjustment for Wilks' lambda. These adjustments are correcting positive bias with increasing number of variables.
%%     ~~Describe \code{model1} here~~
}
}
\details{
This function allows estimating multivariate effect size / multivariate measure of association for the four multivariate statistics implemented in \code{manova.gls} (Pillai, Wilks, Roy, Hotelling-Lawley). These multivariate measures are common generalizations of univariate measures such as the squared-multiple correlation, with values ranging between 0 and 1 (see for instance details in Rencher 2002). Note that these measures are known to be upwardly biased with increased dimensionality (Kim and Olejnik 2005), and several adjustments were proposed (e.g., Serlin 1982, Tatsuoka 1973; see  \code{adjusted=TRUE} argument above). For models fit by PL, the adjustments is done by estimating the bias from the permuted statistics under the null. Note that adjusted measures might be sligthly negative, in such a case the measure should be interpreted as virtually 0 (accordingly, statistical tests from the \code{manova.gls} function should be non-significant in these situations).
%%  ~~ If necessary, more details than the description above ~~
}

\value{
Return the effect size for all the terms of the MANOVA or pairwise tests.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}

\note{
	This function is still under development.}
	
\references{

Clavel et al. in prep.

Kim, S., Olejnik S., 2005. Bias and precision of measures of association for a fixed-effect multivariate analysis of variance model. Multivariate Behavioral Research, 40(4):401-421. 

Rencher, 2002. Methods of Multivariate Analysis (Second Edition). Wiley and Sons. pp. 705.

Serlin, R. C., 1982. A multivariate measure of association based on the Pillai-Bartlett procedure. Psycological Bulletin 91(2):413-417.

Tatsuoka, M. M. 1973. An examination of the statistical properties of a multivariate measure of strength of relationship (final report). Project No. 2-E-020, Grant No. OEG-5-72-0027(509).

}
	
\author{
Julien Clavel
%%  ~~who you are~~
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
\code{\link{manova.gls}}
\code{\link{mvgls}}
\code{\link{mvols}}
\code{\link{pairwise.glh}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
set.seed(123)
n <- 32 # number of species
p <- 3  # number of traits
tree <- pbtree(n=n) # phylogenetic tree
R <- crossprod(matrix(runif(p*p),p))  # a random symmetric matrix (covariance)

# simulate a dataset
Y <- mvSIM(tree, model="BM1", nsim=1, param=list(sigma=R))
X <- rnorm(n) # continuous
grp <- rep(1:2, each=n/2)
dataset <- list(y=Y, x=X, grp=as.factor(grp))

# Model fit
model1 <- mvgls(y~x+grp, data=dataset, tree=tree, model="BM", method="LL")

# Multivariate test
(multivariate_test <- manova.gls(model1, test="Pillai"))
effectsize(multivariate_test)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Effect size }
\keyword{ R-squared }
\keyword{ Multivariate measure of association }% __ONLY ONE__ keyword per line