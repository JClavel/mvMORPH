\name{coef.mvgls}
\alias{coef.mvgls}

\title{
Extract multivariate gls model coefficients
}
\description{
Returns the coefficients of a linear model fit of class 'mvgls'.
}
\usage{

\method{coef}{mvgls}(object, ...)

}
\arguments{
  \item{object}{
an object of class 'mvgls' obtained from a \code{mvgls} fit.
}
  
  \item{...}{
other arguments (not used).
}

}

\note{
For an intercept only model with phylogenetic structure this correspond to the ancestral states.
}


\value{
	The coefficients extracted from the model.
}


\author{J. Clavel}

\seealso{
%% add later  \code{\link{simulate_t_env}}  
\code{\link{vcov.mvgls}} 
\code{\link{residuals.mvgls}}
\code{\link{fitted.mvgls}}
\code{\link{mvgls}}
}

