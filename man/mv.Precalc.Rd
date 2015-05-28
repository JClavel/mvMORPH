\name{mv.Precalc}
\alias{mv.Precalc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Model parameterization for the various mvMORPH functions
%%  ~~function to do ... ~~
}
\description{
This function allows computing the fixed parameters or objects needed by the mvMORPH functions. This could be useful in bootstrap-like computations (see exemple)
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mv.Precalc(tree, nb.traits = 1, scale.height = FALSE, param = list(pivot = "MMD",
            method = c("sparse"), smean = TRUE, model = "OUM"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree}{
A "phylo" (or SIMMAP like) object representing the tree for which we want to precalculate parameters.
%%     ~~Describe \code{tree} here~~
}
  \item{nb.traits}{
The number of traits involved in the subsequent analysis.
%%     ~~Describe \code{nb.traits} here~~
}
  \item{scale.height}{
Whether the tree should be scaled to unit length.
%%     ~~Describe \code{scale.height} here~~
}
  \item{param}{
A list of parameters used in the computations (see details)
%%     ~~Describe \code{param} here~~
}
}
\details{
The mv.Precalc function allows the pre-computation of the fixed parameters required by the different mvMORPH models (e.g., the design matrix, the vcv matrix, the sparsity structure...).
In the "param" list you should provide the details about the model fit:

-model name (e.g., "OUM", "OU1")

-method (which kind of algorithm is used for computing the log-likelihood).

-smean (whether there is one ancestral state per traits or per selective regimes - for mvBM only).

Additional parameters can be fixed:
-root (estimation of the ancestral state for the Ornstein-Uhlenbeck model, see ?mvOU).
-pivot (pivot method used by the "sparse" matrix method for computing the log-likelihood, see ?spam).
%%  ~~ If necessary, more details than the description above ~~
}
\value{
An object of class "mvmorph.precalc" which can be used in the "precalc" argument of the various mvMORPH functions.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}

\author{
Julien Clavel
%%  ~~who you are~~
}
\note{
This function is mainly used internally; it is still in development. A misuse of this functions can result in a crash of the R session.
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{mvMORPH}}
\code{\link{mvOU}}
\code{\link{mvEB}}
\code{\link{mvBM}}
\code{\link{mvSHIFT}}
\code{\link{mvLL}}

%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
set.seed(14)
# Generating a random tree
tree<-pbtree(n=50)

# Simulate two correlated traits evolving along the phylogeny according to a
# Ornstein-Uhlenbeck process
alpha<-matrix(c(2,1,1,1.3),2,2)
sigma<-matrix(c(1,0.5,0.5,0.8),2,2)
theta<-c(3,1)
nsim<-50
simul<-mvSIM(tree,param=list(sigma=sigma, alpha=alpha, ntraits=2, mu=theta,
             names_traits=c("head.size","mouth.size")), model="OU1", nsim=nsim)

# Do the pre-calculations
precal<-mv.Precalc(tree,nb.traits=2, param=list(method="sparse",model="OU1", root=FALSE))

mvOU(tree, simul[[1]], method="sparse", model="OU1", precalc=precal,
    param=list(decomp="symmetric"))

### Bootstrap
# Fit the model to the "nsim" simulated datasets
# results<-lapply(1:nsim,function(x){
# mvOU(tree, simul[[x]], method="sparse", model="OU1", precalc=precal,
#    param=list(decomp="symmetric"))
# })

### Use parallel package (not run)
# library(parallel)
# number_of_cores<-2L
# results<-mclapply(1:nsim, function(x){
#    mvOU(tree, simul[[x]], method="sparse", model="OU1", precalc=precal,
#       param=list(decomp="symmetric"))
# }, mc.cores = getOption("mc.cores", number_of_cores))


# Summarize the likelihood
# loglik<-sapply(1:nsim,function(x){logLik(results[[x]])})
# hist(loglik)

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ precalculation }
\keyword{ parameters }% __ONLY ONE__ keyword per line
