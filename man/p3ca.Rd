\name{p3ca}
\alias{p3ca}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{

Modelling multivariate trait evolution using a Probabilistic and Phylogenetic Principal Component Analysis

%%  ~~function to do ... ~~
}
\description{
The function performs a probabilistic and phylogenetic PCA, a maximum likelihood-based approach which models the evolution of multiple traits (describing one or more phenotypes) using a 
fewer set of variables (called latent variables) than in the original dataset. This approach incorporates the dependency arising from the evolutionary shared history of comparative data in the probabilistic PCA (Tipping and Bishop, 1999; Roweis, 1997; Li et al., 2009). It includes an Expectation-Maximisation (EM) algorithm that allows to make parameter inference while accounting for missing values. A complete description of the approach can be found in Montoya et al., 2026.

%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
p3ca(Y, tree, q=4, model="lambda", REML=TRUE, tol=1e-4, 
     EM=FALSE, plot=TRUE, ...)
      
      
}
%- maybe also 'usage' for other objects documented here.
\arguments{ 

\item{Y}{
Traits in the form of \code{data frame} or \code{matrix} , where every row should contain the traits for one species (number of columns = number of traits; number of rows = number of species). The row names should match the tip names in the phylogenetic tree. Only quantitative measurements can be included as data. Missing values are allowed and should be noted as NA.
}

\item{tree}{
Phylogenetic tree (object of class \code{phylo}) including only the species in Y. The tips names should match with the row names in Y.
}

\item{q}{
Number of latent variables describing the reduced space (\code{integer}). \emph{q} should be lower than the minimum between the number of traits (\emph{p}) and the number of species (\emph{n}). The accuracy in the inference reduces as \emph{q} approaches \emph{p} or \emph{n}. 
}

\item{model}{
Trait evolution model. Only "\code{lambda}" for Pagels' lambda and "\code{BM}" for Brownian motion, are currently available.
}

\item{REML}{
When \code{TRUE} (default) the inference is performed using the Restricted Maximum Likelihood. Otherwise, the maximum likelihood (ML) is used.
}

\item{tol}{
Tolerance value used to assess convergence in the EM algorithm. 
}

\item{EM}{
When \code{TRUE}, the Expected-Maximisation (EM) algorithm is used for parameter inference. See Details. 
}

\item{plot}{
When \code{TRUE}, plot the two firsts axes of the P3CA. The "axes" argument can be provided through the ellipsis (e.g., \code{axes=c(1,3)}). 
}

 \item{...}{
Further arguments to be passed through. For instance, When \code{verbose = TRUE} (default), the function prints the number of iterations of the EM algorithm at convergence.
\code{maxit=5000} (default) is the maximum number of iterations of the EM algorithm (\code{integer}).
}


}
\details{
The probabilistic and phylogenetic PCA (P3CA) performs parameter inference by maximum likelihood, using either analytical solutions or an Expectation-Maximisation algorithm (see Tipping and Bishop, 1999; Montoya et al., 2026). 
Unlike the analytical solution, the EM algorithm allows fitting the P3CA with missing values. At convergence, the EM is expected to reach the ML solution.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%% TODO the latent variable and the PC axes are different. The loadings are the correlations between the PCs (i.e. the rotated latent variables) and the data.
A \code{list} containing the following elements :

\item{par       }{Parameter describing the trait evolution model. Lambda for Pagels lambda, and NA when the BM model is used.}
\item{loglik    }{Log-likelihood. It corresponds to the restricted maximum likelihood when \code{REML=TRUE}.}
\item{W         }{\code{matrix} mapping the observed space (number of traits) to the reduced space (number of latent variables).}
\item{sigma2    }{Estimated parameter describing the noise of the latent variable model underlying the P3CA (see Tipping and Bishop, 1999; Montoya et al., 2026).}
\item{L         }{\code{matrix} containing the loadings, the correlations between each trait and the latent variables. The stronger the loadings, the greater the contribution of the trait to the latent variable.}
\item{scores    }{PCs scores (\code{matrix}). It describes the position of each species in the reduced space. }
\item{vectors   }{Eigenvectors of the trait covariance matrix described from the latent variable model underlying the P3CA (see Montoya et al., 2026). This is the rotation matrix that was used to calculate PC scores.}
\item{eigenval  }{Eigenvalues of the trait covariance matrix described from the latent variable model underlying the P3CA (see Montoya et al., 2026).}
\item{varExp    }{Variance explained by each PCs (rotated latent variables) describing the reduced space.}
\item{coef      }{Ancestral values estimated for each trait (i.e., mean used to center the data), using the generalized least squares solution. When the dataset includes missing values, the reconstruction is made considering only the observed values.}
\item{model     }{Trait evolutionary model fitted by the function.}
\item{imputed   }{Trait \code{matrix} including the estimated values for the missing cases.}
\item{count     }{Number of iterations required for reaching convergence in the EM algorithm.}
\item{tol       }{Tolerance used for assessing convergence in the EM algorithm.}
\item{maxit     }{Maximum number of iterations used in the EM algorithm.} %% is it useful ? This is an argument of the original function
\item{prop_missing}{Relative proportion of missing values if any.}


%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{

Li W.-J., Yeung D.-Y., Zhang Z. 2009. Probabilistic Relational PCA. Advances in neural information processing systems. 22.

Montoya P., Joseph J., Goswami A., Morlon H., Clavel J. 2026. A probabilistic and phylogenetic principal component analysis for modelling high-dimensional trait evolution. doi.org/10.64898/2026.05.27.728209.

Roweis S. 1997. EM Algorithms for PCA and SPCA. Advances in neural information processing systems. 10.

Tipping M.E., Bishop C.M. 1999. Probabilistic Principal Component Analysis. Journal of the Royal Statistical Society Series B: Statistical Methodology. 61:611-622.

%% ~put references to the literature/web site here ~
}
\author{
Paola Montoya and Julien Clavel
%%  ~~who you are~~
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{mvgls}}
\code{\link{mvgls.pca}}
\code{\link{pcaShape}}
\code{\link{pcaLoadings}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}

\examples{
\donttest{

### Example starts ###
library(mvMORPH)

set.seed(2508)

# Loading the data
data(phyllostomid)
phyllos_data = phyllostomid$mandible[,-1]
phyllos_tree = phyllostomid$tree

# Perfoming the P3CA on the complete dataset - Analytical solution
p3ca_phyllos = p3ca(phyllos_data, phyllos_tree, q=5, model='lambda')
p3ca_phyllos$logLik
p3ca_phyllos$par

# Dataset including missing values
id_for_nas = sample(1:length(phyllos_data),100)
phyllos_data_nas = phyllos_data
phyllos_data_nas[id_for_nas]=NA

# Perfoming the P3CA on the incomplete dataset - EM algorithm
p3ca_phyllos_nas = p3ca(phyllos_data_nas, phyllos_tree, q=5, model='lambda')
p3ca_phyllos_nas$logLik
p3ca_phyllos_nas$par

# Visualising and comparing the reduced spaces (by hand) obtained on the
# original and imputed datasets
par(mfrow=c(1,2))
axes=1:2
plot(p3ca_phyllos$scores[,axes], pch=19, col=c('dodgerblue3','forestgreen')[phyllostomid$grp1],
     xlab = paste('P3C',axes[1],' ',round(p3ca_phyllos$varExp[axes[1]],1),'\%', sep=''),
     ylab = paste('P3C',axes[2],' ',round(p3ca_phyllos$varExp[axes[2]],1),'\%', sep=''),
     main = paste('P3CA (AS) - ', p3ca_phyllos$model, ' ', 
     round(p3ca_phyllos$par,1), sep=''))

plot(p3ca_phyllos_nas$scores[,axes], pch=19, col=c('dodgerblue3','forestgreen')[phyllostomid$grp1],
     xlab = paste('P3C',axes[1],' ',round(p3ca_phyllos_nas$varExp[axes[1]],1),'\%', sep=''),
     ylab = paste('P3C',axes[2],' ',round(p3ca_phyllos_nas$varExp[axes[2]],1),'\%', sep=''),
     main = paste('P3CA (EM) - ', p3ca_phyllos_nas$model, ' ',
     round(p3ca_phyllos_nas$par,1), sep=''))

### Example ends ###

}

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Dimensionality reduction technique }
\keyword{ High dimension }
\keyword{ Phylogenetic Factor Analysis }
\keyword{ Rotation invariance }% __ONLY ONE__ keyword per line
