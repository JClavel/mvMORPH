\name{LRT}
\alias{LRT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Likelihood Ratio Test
%%  ~~function to do ... ~~
}
\description{
This function compares the fit of two nested models of trait evolution with a loglikelihood-ratio statistic.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
LRT(model1, model2, echo = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model1}{
The most parameterized model. A fitted object from an mvMORPH model.
%%     ~~Describe \code{model1} here~~
}
  \item{model2}{
The second model under comparison (fitted object).
%%     ~~Describe \code{model2} here~~
}
  \item{echo}{
Whether to return the result or not.
%%     ~~Describe \code{echo} here~~
}
}
\details{
The LRT function extracts the log-likelihood of two nested models to compute the loglikelihood-ratio statistic which is compared to a Chi-square distribution. Note that if the models are not nested, the LRT is not an appropriate test and you should rely instead on Information criteria, evidence ratios, or simulated distributions (e.g., Lewis et al. 2011).
%%  ~~ If necessary, more details than the description above ~~
}
\value{

\item{pval }{The p-value of the LRT test (comparison with Chi-square distribution).}
\item{ratio }{The LRT (Loglikelihood-ratio test) statistic.}
\item{ddf   }{The number of degrees of freedom between the two models.}
\item{model1   }{Name of the first model.}
\item{model2   }{Name of the second model.}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Neyman J., Pearson E.S. 1933. On the problem of the most efficient tests of statistical hypotheses. Philos. Trans. R. Soc. A. 231:289-337.

Lewis F., Butler A., Gilbert L. 2011. A unified approach to model selection using the likelihood ratio test. Meth. Ecol. Evol. 2:155-162.

%% ~put references to the literature/web site here ~
}
\author{
Julien Clavel
%%  ~~who you are~~
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{mvMORPH}}
\code{\link{mvOU}}
\code{\link{mvEB}}
\code{\link{mvBM}}
\code{\link{mvSHIFT}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## Simulated dataset
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

# Simulate two correlated traits evolving along the phylogeny
traits<-mvSIM(tree,nsim=1, model="BMM", param=list(sigma=list(matrix(c(2,1,1,1.5),2,2),
         matrix(c(4,1,1,4),2,2)), ntraits=2, names_traits=c("head.size","mouth.size")))

# Fit of model 1
mod1<-mvBM(tree,traits,model="BMM")

# Fit of model 2
mod2<-mvBM(tree,traits,model="BM1")

# comparing the fit using LRT...
LRT(mod1,mod2)

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ LRT }
\keyword{ Loglikelihood ratio test }
\keyword{ Models comparison }% __ONLY ONE__ keyword per line
