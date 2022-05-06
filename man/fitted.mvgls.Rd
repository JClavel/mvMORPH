\name{fitted}
\alias{fitted.mvgls}

\title{
Extract multivariate gls (or ols) model fitted values
}
\description{
Returns the fitted values of a linear model of class 'mvgls'.
}
\usage{

\method{fitted}{mvgls}(object, ...)

}
\arguments{
  \item{object}{
an object of class 'mvgls' obtained from a \code{mvgls} or \code{mvols} fit.
}
  
  \item{...}{
other arguments (not used).
}

}

\value{
	The fitted values extracted from the model.
}


\author{J. Clavel}

\seealso{
%% add later  \code{\link{simulate_t_env}}  
\code{\link{vcov.mvgls}} 
\code{\link{residuals.mvgls}}
\code{\link{coef.mvgls}}
\code{\link{mvgls}}
\code{\link{mvols}}
}

