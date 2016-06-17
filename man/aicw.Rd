\name{aicw}
\alias{aicw}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Akaike weights
%%  ~~function to do ... ~~
}
\description{
This function return the Akaike weights for a set of fitted models.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
aicw(x,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A list with the fitted objects or a list/vector of AIC
%%     ~~Describe \code{tree} here~~
	}
  \item{...}{
Options to be passed through; e.g. aicc=TRUE when a list of fitted objects is provided.
%%     ~~Describe \code{tree} here~~
	}
}
 
\details{
This function compute the Akaike weights for a set of model AIC or AICc. Akaike weights can be used for model comparison and model averaging.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{models }{List of models}
\item{AIC }{Akaike Information Criterion}
\item{diff }{AIC difference with the best fit model}
\item{wi }{Absolute weight}
\item{aicweights }{Akaike weights (relative weights)}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Burnham K.P., Anderson D.R. 2002. Model selection and multi-model inference: a practical information-theoric approach. New York: Springer-Verlag.

%% ~put references to the literature/web site here ~
}
\author{
Julien Clavel
%%  ~~who you are~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{AIC}}
\code{\link{mvMORPH}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
set.seed(1)
# Generating a random tree
tree<-pbtree(n=50)

#simulate the traits
sigma <- matrix(c(0.01,0.005,0.003,0.005,0.01,0.003,0.003,0.003,0.01),3)
theta<-c(0,0,0)
data<-mvSIM(tree, model="BM1", nsim=1, param=list(sigma=sigma, theta=theta))

## Fitting the models
# BM1 - General structure
fit1 <- mvBM(tree, data, model="BM1", method="pic")

# BM1 - No covariations
fit2 <- mvBM(tree, data, model="BM1", method="pic", param=list(constraint="diagonal"))

# BM1 - Equal variances/rates
fit3 <- mvBM(tree, data, model="BM1", method="pic", param=list(constraint="equal"))

results <- list(fit1,fit2,fit3)

# or
# results <- c(AIC(fit1), AIC(fit2), AIC(fit3))

# Akaike weights
aicw(results)

# AICc weights
aicw(results, aicc=TRUE)

# we can compare the MSE...
# mean((fit1$sigma-sigma)^2)
# mean((fit3$sigma-sigma)^2)

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ AIC }
\keyword{ Akaike weights }
% __ONLY ONE__ keyword per line
