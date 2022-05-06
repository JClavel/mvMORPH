\name{predict}
\alias{predict.mvgls}

\title{
Predictions from (multivariate) gls or ols model fit
}
\description{
Returns the prediction(s) of a linear model of class 'mvgls'.
}
\usage{

\method{predict}{mvgls}(object, newdata, ...)

}
\arguments{
  \item{object}{
an object of class 'mvgls' obtained from a \code{mvgls} or \code{mvols} fit.
}
  \item{newdata}{
a dataframe with new observation(s). The column names must match the names of the predictors in the model fit object. The type (e.g. factors, numeric) must also match the type of the predictors in the model fit object.
Note: the fitted values are simply returned if "newdata" is not provided.
}

  \item{...}{
further arguments for this generic function. For models fit by \code{mvgls}, if \code{tree} is provided (with tip name(s) matching rowname(s) in newdata and in the training (model fit) dataset), then the best unbiased linear prediction (BLUB) for the model is returned. Otherwise the GLS coefficients are used to predict "newdata".
}

}


\value{
	A matrix with the predictions for the linear model fitted by \code{mvgls} or \code{mvols}.
}


\author{J. Clavel}

\seealso{

\code{\link{fitted.mvgls}} 
\code{\link{vcov.mvgls}} 
\code{\link{residuals.mvgls}}
\code{\link{coef.mvgls}}
\code{\link{mvgls}}
\code{\link{mvols}}
}

