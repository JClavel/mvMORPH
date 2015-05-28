\name{stationary}
\alias{stationary}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The stationary variance of an Ornstein-Uhlenbeck process
%%  ~~function to do ... ~~
}
\description{
This function return the stationary variance for an Ornstein-Uhlenbeck process (object of class "ou").
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
stationary(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
Object fitted with the "mvOU" function.
%%     ~~Describe \code{tree} here~~
}
}
 
\details{
This function compute the dispersion parameter of the Ornstein-Uhlenbeck process (i.e., the expected variance when the process is stationary).
The multivariate normal stationary distribution of the Ornstein-Uhlenbeck process is computed following Bartoszeck et al. (2012).
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The stationary variance-covariance matrix of the OU process
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Bartoszek K., Pienaar J., Mostad P., Andersson S., Hansen T.F. 2012. A phylogenetic comparative method for studying multivariate adaptation. J. Theor. Biol. 314:204-215.

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
\code{\link{halflife}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## Toy example

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
data<-mvSIM(tree, param=list(sigma=sigma, alpha=alpha, ntraits=2, mu=theta,
            names_traits=c("head.size","mouth.size")), model="OUM", nsim=1)

## Fitting the models
# OUM - Analysis with multiple optima
result<-mvOU(tree, data)

stationary(result)

# Expected values when the process is stationary
expected<-list(alpha=alpha,sigma=sigma)
class(expected)<-c("mvmorph","ou")
stationary(expected)

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Ornstein Uhlenbeck }
\keyword{ stationary }
\keyword{ OU }% __ONLY ONE__ keyword per line
