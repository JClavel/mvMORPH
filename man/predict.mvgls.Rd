\name{predict.mvgls}
\alias{predict.mvgls}

\title{
Predictions from (multivariate) gls model fit
}
\description{
Returns the prediction(s) of a linear model of class 'mvgls'.
}
\usage{

\method{predict}{mvgls}(object, newdata, ...)

}
\arguments{
  \item{object}{
an object of class 'mvgls' obtained from a \code{mvgls} fit.
}
  \item{newdata}{
a dataframe with new observation(s). The column names must match the names of the predictors in the model fit object. The type (e.g. factors, numeric) must also match the type of the predictors in the model fit object.
Note: the fitted values are simply returned if "newdata" is not provided.
}

  \item{...}{
other arguments for this generic function. If \code{tree} is provided (with tip name(s) matching rowname(s) in newdata), then the best unbiased linear prediction for the model is returned. Otherwise the GLS coefficients are used to predict "newdata".
}

}


\value{
	A matrix with the predictions for the linear model fitted by \code{mvgls}.
}


\author{J. Clavel}

\seealso{

\code{\link{fitted.mvgls}} 
\code{\link{vcov.mvgls}} 
\code{\link{residuals.mvgls}}
\code{\link{coef.mvgls}}
\code{\link{mvgls}}
}

