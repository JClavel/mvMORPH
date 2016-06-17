\name{halflife}
\alias{halflife}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The phylogenetic half-life for an Ornstein-Uhlenbeck process
%%  ~~function to do ... ~~
}
\description{
This function returns the phylogenetic half-life for an Ornstein-Uhlenbeck process (object of class "ou").
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
halflife(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
Object fitted with the "mvOU" function.
%%     ~~Describe \code{tree} here~~
	}
}
 
\details{
The phylogenetic half-life describes the time to move halfway from the ancestral state to the primary optimum (Hansen, 1997).
The multivariate counterpart is computed on the eigenvalues of the "selection" matrix (Bartoszeck et al. 2012).
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The phylogenetic half-life computed from each eigenvalues (or alpha for the univariate case)
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Bartoszek K., Pienaar J., Mostad P., Andersson S., Hansen T.F. 2012. A phylogenetic comparative method for studying multivariate adaptation. J. Theor. Biol. 314:204-215.

Hansen T.F. 1997. Stabilizing selection and the comparative analysis of adaptation. Evolution. 51:1341-1351.

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
\code{\link{stationary}}
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
# OUM - Analysis with multiple optima
result<-mvOU(tree, data)

halflife(result)

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Ornstein Uhlenbeck }
\keyword{ half-life }
\keyword{ OU }% __ONLY ONE__ keyword per line
