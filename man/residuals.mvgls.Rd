\name{residuals.mvgls}
\alias{residuals.mvgls}

\title{
Extract gls model residuals
}
\description{
Returns the residuals of a linear model of class 'mvgls'.
}
\usage{

\method{residuals}{mvgls}(object, type, ...)

}
\arguments{
  \item{object}{
an object of class 'mvgls' obtained from a \code{mvgls} fit.
}
  \item{type}{
an optional character string specifying the type of residuals to be used. To match conventions used in the \emph{nlme} package: if "\code{response}", the "raw" residuals (observed-fitted) are used; else, if "normalized", the normalized residuals (the residuals pre-multiplied by the inverse square-root factor of the estimated (between observations) covariance matrix) are used. Note however that there are still between variables correlations with both types.
}

  \item{...}{
other arguments for this generic function (not used).
}

}


\value{
	A matrix with the residuals for the linear model fitted by \code{mvgls}.
}


\author{J. Clavel}

\seealso{
%% add later  \code{\link{simulate_t_env}}  
\code{\link{vcov.mvgls}} 
\code{\link{residuals.mvgls}}
\code{\link{coef.mvgls}}
\code{\link{mvgls}}
}

