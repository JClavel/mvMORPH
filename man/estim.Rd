\name{estim}
\alias{estim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Ancestral states reconstructions and missing value imputation with phylogenetic/time-series models
%%  ~~function to do ... ~~
}
\description{
This function imputes the missing cases (NA values) according to a given phylogenetic model (object of class "mvmorph"); it can also do ancestral state reconstruction.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
estim(tree, data, object, error=NULL, asr=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{ 

 \item{tree}{
Phylogenetic tree (an object of class "phylo" or "simmap") or a time-series.
%%     ~~Describe \code{model2} here~~
}
  \item{data}{
Matrix or data frame with species in rows and continuous traits with missing cases (NA values) in columns (preferentially with names and in the same order than in the tree).
%%     ~~Describe \code{model1} here~~
}
  \item{object}{
A fitted object from an mvMORPH model (class "mvmorph").
%%     ~~Describe \code{echo} here~~
}
  \item{error}{
Matrix or data frame with species in rows and continuous traits sampling variance (squared standard errors) in columns.
%%     ~~Describe \code{echo} here~~
}

\item{asr}{
If asr=TRUE, the ancestral states are estimated instead of the missing cases.
%%     ~~Describe \code{echo} here~~
}
}
\details{
Missing observations for species in a phylogenetic tree are estimated according to a given evolutionary model (and parameters). Multivariate models are useful to recover the variance and covariance structure of the dataset to be imputed.

When \emph{asr=TRUE}, the estimates, their variances and standard errors are those of the ancestral states at each node of the tree (this option is not available for the time-series). Note that if there are missing cases, they are first imputed before estimating the ancestral states.


%%  ~~ If necessary, more details than the description above ~~
}
\value{

\item{estimates}{The imputed dataset }
\item{var      }{Variance of the estimates}
\item{se       }{Standard error of the estimates}
\item{NA_index }{Position of the missing cases in the dataset}

%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Clavel J., Merceron G., Escarguel G. 2014. Missing Data Estimation in Morphometrics: How Much is Too Much? Syst. Biol. 63:203-218.

Cunningham C.W., Omland K.E., Oakley T.H. 1998. Reconstructing ancestral character states: a critical reappraisal. Trends Ecol. Evol. 13:361-366.
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
         matrix(c(4,1,1,4),2,2)), names_traits=c("head.size","mouth.size")))

# Introduce some missing cases (NA values)
data<-traits
data[8,2]<-NA
data[25,1]<-NA

# Fit of model 1
fit<-mvBM(tree,data,model="BMM")

# Estimate the missing cases
imp<-estim(tree, data, fit)

# Check the imputed data
imp$estim[1:10,]

## We want the ancestral states values at each nodes:
nodelabels() # To see where the nodes are situated

imp2<-estim(tree, data, fit, asr=TRUE)

# Check the 10 firsts ancestral states
imp2$estim[1:10,]
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Estim }
\keyword{ Imputation }
\keyword{ Missing values }% __ONLY ONE__ keyword per line
