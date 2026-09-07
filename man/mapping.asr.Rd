\name{mapping.asr}
\alias{mapping.asr}

\title{
Mapping ancestral character estimates onto a phylogenetic tree
}
\description{
The function maps ancestral state estimates for discrete characters onto a phylogenetic tree. Ancestral character estimates should be obtained prior to use the function.
}
\usage{

mapping.asr(tree, ancestral, tips)
  
}
\arguments{
  \item{tree}{
Phylogenetic tree (class \code{phylo}).
}
  \item{ancestral}{
Ancestral state estimates. Either an object of class \code{describe.simmap} or \code{ace}. 
}
 \item{tips}{
Character states at the tips. Either an object of class \code{character} or \code{factor}, which elements are named and ordered according to the phylogenetic tree.
}
}

\details{
The function follows the structure of the function \code{make.simmap()} from the package \code{phytools}, adding two elements describing the time spent in each state per edge. A phylogenetic tree with the ancestral state estimate mapped onto it is used to model trait evolution under multiple selective regimes in functions such as \code{mvgls()}.  }

\value{
	A phylogenetic tree of class \code{simmap}. Two elements are added:
  
  \item{mapped.edge}{\code{matrix} describing the duration time of each state per edge}
  \item{maps}{\code{list}. Its elements indicate the time spent in each state per branch, in chronological order }
  }
  

\author{Julien Clavel}

\seealso{
\code{\link{mvgls}}
\code{\link[ape]{ace}}
\code{\link[phytools]{make.simmap}}
}

\examples{
\donttest{
#### Example starts ####

library(mvMORPH)
library(ape)

set.seed(2508)

# Loading the data
data(phyllostomid)
phyllos_tree = phyllostomid$tree
phyllos_trait = phyllostomid$grp1

# The states at the tips should follow the order of the phylogenetic tree
phyllos_trait = phyllos_trait[phyllos_tree$tip.label]

# Estimating the ancestral character states
anc_state = ace(phyllos_trait, phyllos_tree, type="discrete",
    model="SYM")

# Mapping the reconstruction on the phylogenetic tree
new_tree = mapping.asr(tree=phyllos_tree, ancestral=anc_state, tips=phyllos_trait)
plot(new_tree)

#### Example ends ####
}
}

\keyword{ simmap }
\keyword{ ancestral state reconstruction }
\keyword{ BMM }
\keyword{ OUM }% __ONLY ONE__ keyword per line
