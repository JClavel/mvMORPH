\name{vcov}
\alias{vcov.mvgls}

\title{
Calculate variance-covariance matrix for a fitted object of class 'mvgls'
}
\description{
Returns the variance-covariance matrix of the coefficients or the traits.
}
\usage{

\method{vcov}{mvgls}(object, ...)

}
\arguments{
  \item{object}{
an object of class 'mvgls' obtained from a \code{mvgls} or \code{mvols} fit.
}
  
  \item{...}{
additional arguments for methods function. See \emph{details} below.
}

}

\details{
The \code{vcov} function returns by default the variance-covariance matrix of the main parameters of a fitted model object. The main parameters are the coefficients (this correspond to the argument \code{type="coef"}; see also \code{coef.mvgls}). With \code{type="covariance"}, the \code{vcov.mvgls} function returns the estimated traits covariance matrix (possibly regularized for PL approaches) while \code{type="precision"} return the precision matrix (i.e. the inverse of the covariance).
}


\value{
	A matrix of the estimated covariances between the parameter estimates (of type "coef", "covariance", or "precision").
  
}


\author{J. Clavel}

\seealso{
%% add later 
\code{\link{coef.mvgls}} 
\code{\link{residuals.mvgls}}
\code{\link{fitted.mvgls}}
\code{\link{mvgls}}
\code{\link{mvols}}
}

