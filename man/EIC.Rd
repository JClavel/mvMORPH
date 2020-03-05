\name{EIC}
\alias{EIC}

\title{
Extended Information Criterion (EIC) to compare models fit with \code{mvgls} by Maximum Likelihood (ML) or Penalized Likelihood (PL).
}

\description{
The EIC (Ishiguro et al. 1997, Kitagawa & Konishi 2010), uses bootstrap to estimate the bias term of the Extended Information Criterion. This criterion allows comparing models fit by Maximum Likelihood (ML) or Penalized Likelihood (PL).
}

\usage{

EIC(object, nboot=100L, nbcores=1L, ...)
  
}

\arguments{
  \item{object}{
An object of class 'mvgls'. See ?mvgls}
  \item{nboot}{
The number of boostrap replicates used for estimating the EIC.}
  \item{nbcores}{
The number of cores used to speed-up the computations (uses the 'parallel' package)}
  \item{...}{
Options to be passed through (see details).}
}


\value{
		a list with the following components
  
  \item{LogLikelihood}{the log-likelihood estimated for the model with estimated parameters}
  \item{EIC}{the EIC criterion}
  \item{se}{the standard error of the bias term estimated by bootstrap}
  \item{bias}{the values of the bias term estimated from the boostrapped replicates to compute the EIC}

}

\details{
The Extended Information Criterion (\code{EIC}) allows comparing the fit of various models estimated by Penalized Likelihood or Maximum Likelihood (see ?\code{mvgls}). Similar to the GIC or the more common AIC, the EIC has the form:

\emph{EIC = -2*(Likelihood) + 2*bias}

Where \emph{Likelihood} corresponds to either the full or the restricted likelihood (see the note below), and the bias term is estimated by (semi-parametric) bootstrap simulations rather than by using analytical or approximations solutions (see for instance ?\code{GIC}). The smaller the EIC, the better is the model. With small sample sizes, the variability around the bootstrap estimates is expected to be high, and one must increase the number of bootstrap replicates. Parallel computation (argument \code{nbcores}) allows to speed-up the computations.

Note: for models estimated by REML, it is generally not possible to compare the restricted likelihood when they have different fixed effects. However, it is possible to compare models with different fixed effects by using the full likelihood evaluated with the REML estimates (see e.g. Yafune et al. 2006, Verbyla 2019). Both options - evaluating the restricted likelihood or the full likelihood with parameters estimated by REML - are available through the \code{REML} argument in the \code{EIC} function. The default has been set to \code{REML=FALSE} to allow the comparison of models with different fixed effects using the full likelihood evaluated with the REML estimates (see Verbyla 2019).
}

	
\references{
Clavel J., Aristide L., Morlon H., 2019. A Penalized Likelihood framework for high-dimensional phylogenetic comparative methods and an application to new-world monkeys brain evolution. Syst. Biol. 68:93-116.

Ishiguro M., Sakamoto Y., Kitagawa G., 1997. Bootstrapping log likelihood and EIC, an extension of AIC. Ann. Inst. Stat. Math. 49:411-434.

Kitagawa G., Konishi S., 2010. Bias and variance reduction techniques for bootstrap information criterion. Ann. Inst. Stat. Math. 62:209-234.

Konishi S., Kitagawa G., 1996. Generalised information criteria in model selection. Biometrika. 83:875-890.

Verbyla A. P., 2019. A note on model selection using information criteria for general linear models estimated using REML. Aust. N. Z. J. Stat. 61:39-50.

Yafune A., Funatogawa T., Ishiguro M., 2005. Extended information criterion (EIC) approach for linear mixed effects models under restricted maximum likelihood (REML) estimation. Statist. Med. 24:3417-3429.
}

\author{J. Clavel}

\seealso{
\code{\link{GIC}}
\code{\link{mvgls}}
\code{\link{mvgls.pca}}
}

\examples{
\donttest{

set.seed(1)
n <- 32 # number of species
p <- 50 # number of traits

tree <- pbtree(n=n) # phylogenetic tree
R <- crossprod(matrix(runif(p*p), ncol=p)) # a random symmetric matrix (covariance)
# simulate a dataset
Y <- mvSIM(tree, model="BM1", nsim=1, param=list(sigma=R))

fit1 <- mvgls(Y~1, tree=tree, model="BM", method="H&L")
fit2 <- mvgls(Y~1, tree=tree, model="OU", method="H&L")


EIC(fit1); EIC(fit2)
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Model comparison }
\keyword{ GLS }
\keyword{ High dimensions }% __ONLY ONE__ keyword per line