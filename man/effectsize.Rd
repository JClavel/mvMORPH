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
This function allows estimating multivariate effect size for the four multivariate statistics implemented in \code{manova.gls} (Pillai, Wilks, Roy, Hotelling-Lawley). For models fit by PL, a multivariate measure of effect size is estimated from the permuted data. Interpret only relatively.
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
	
\author{
Julien Clavel
%%  ~~who you are~~
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
\code{\link{manova.gls}}
\code{\link{mvgls}}
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