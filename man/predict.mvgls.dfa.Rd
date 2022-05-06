\name{predict.mvgls.dfa}
\alias{predict.mvgls.dfa}

\title{
Predictions from Discriminant analysis conducted with a mvgls model fit
}
\description{
Returns the prediction(s) of DFA of class 'mvgls.dfa'.
}
\usage{

\method{predict}{mvgls.dfa}(object, newdata, prior, ...)

}
\arguments{
  \item{object}{
an object of class 'mvgls' obtained from a \code{mvgls} or \code{mvols} fit.
}
  \item{newdata}{
a matrix with new observation(s) for the response variables. 
Note: the predictions are performed on fitted values if "newdata" is not provided.
}
 \item{prior}{
the group priors. If not provided, assumes equal prior.
}
  \item{...}{
further arguments for this generic function. If \code{tree} is provided (with tip name(s) matching rownames in newdata and in the training sample (model fit)), then the best unbiased linear prediction (BLUP) for the model is returned. Otherwise the GLS coefficients are used to predict "newdata", in this condition classification might be less optimal than performing a regular DFA (see \code{lda} from MASS or mvgls.dfa on a \code{mvols} fit).
}

}


\value{

\item{class }{The class assigned to each new observations}
\item{posterior }{The posterior probabilities used to classify each new observations}
\item{prior }{The prior used to classify each new observations to each categories}
}


\author{J. Clavel}

\references{
Duhamel A. et al. in prep.
}

\seealso{

\code{\link{mvgls.dfa}} 
\code{\link{predict.mvgls}} 
\code{\link{fitted.mvgls}} 
\code{\link{vcov.mvgls}} 
\code{\link{residuals.mvgls}}
\code{\link{coef.mvgls}}
\code{\link{mvgls}}
\code{\link{mvols}}
}

\examples{
\donttest{
library(mvMORPH)
n=64
p=4

tree <- pbtree(n=n)
sigma <- crossprod(matrix(runif(p*p),p,p))
resid <- mvSIM(tree, model="BM1", param=list(sigma=sigma))
Y <- rep(c(0,1.5), each=n/2) + resid
grp <- as.factor(rep(c("gp1","gp2"),each=n/2))
names(grp) = rownames(Y)
data <- list(Y=Y, grp=grp)
mod <- mvgls(Y~grp, data=data, tree=tree, model="BM")

# fda
da1 <- mvgls.dfa(mod)

}
}

% R documentation directory.
\keyword{ DFA }
\keyword{ LDA }
\keyword{ Prediction }
\keyword{ CVA }% __ONLY ONE__ keyword per line