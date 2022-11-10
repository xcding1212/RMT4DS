\name{GeneralMPLaw}
\alias{qgmp}
\alias{rgmp}
\alias{pgmp}
\alias{dgmp}
\title{General Marchenko-Pastur Distribution}
\description{
    Density, distribution function, quantile function and random generation for
    the general Marchenko-Pastur distribution, the limiting distribution
    of empirical spectral measure for large Wishart matrices.
}
\usage{
qgmp(p, ndf=NULL, pdim=NULL, svr=ndf/pdim, eigens=NULL, lower.tail=TRUE,
    log.p=FALSE, m=500)
rgmp(n, ndf=NULL, pdim=NULL, svr=ndf/pdim, eigens=NULL, m=500)
pgmp(q, ndf=NULL, pdim=NULL, svr=ndf/pdim, eigens=NULL, lower.tail=TRUE,
    log.p=FALSE, m=500)
dgmp(x, ndf=NULL, pdim=NULL, svr=ndf/pdim, eigens=NULL, log.p=FALSE, m=500)
}
\arguments{
    \item{x,q}{vector of quantiles.}
    \item{p}{vector of probabilities.}
    \item{n}{number of observation.}
    \item{m}{number of points used in estimating density.}
    \item{ndf}{the number of degrees of freedom for the Wishart matrix.}
    \item{pdim}{the number of dimensions (variables) for the Wishart matrix.}
    \item{svr}{samples to variables ratio; the number of degrees of freedom per dimension.}
    \item{log, log.p}{logical; if TRUE, probabilities p are given as log(p).}
    \item{lower.tail}{logical; if TRUE (default), probabilities are
    \eqn{P[X \le x]}, otherwise, \eqn{P[X > x]}.}
    \item{eigens}{input eigenvalues of population covariance matrix.}
}
\value{
    \code{dgmp} gives the density,

    \code{pgmp} gives the distribution function,

    \code{qgmp} gives the quantile function,

    \code{rgmp} generates random deviates,
}
\details{
    Those functions work only for non-spiked part.

    To achieve high accuracy of estimation, \code{eigens} should be large, like
    larger than 500.

    In general Marchenko Pastur distributions, the support of density is the union
    of one or more intervals.
}
\source{
    If \code{eigens} is missing, functions from package
    \code{RMTstat} will be used to compute classical Marchenko-Pastur
    distribution.
}
\examples{
N = 1000
M = 300
d = c(rep(3.8,M/3),rep(1.25,M/3),rep(0.25,M/3))
qgmp(0.5, ndf=N, pdim=M, eigens=d)
pgmp(3, ndf=N, pdim=M, eigens=d)
dgmp(2, ndf=N, pdim=M, eigens=d)
rgmp(2, ndf=N, pdim=M, eigens=d)
}
\references{
  [1] Knowles, A., & Yin, J. (2017). Anisotropic local laws for random matrices. Probability Theory and Related Fields, 169(1), 257-352.

  [2] Bai, Z., & Yao, J. (2012). On sample eigenvalues in a generalized spiked population model. Journal of Multivariate Analysis, 106, 167-177.

  [3] Ding, X. (2021). Spiked sample covariance matrices with possibly multiple bulk components. Random Matrices: Theory and Applications, 10(01), 2150014.

  [4] Ding, X., & Trogdon, T. (2021). A Riemann--Hilbert approach to the perturbation theory for orthogonal polynomials: Applications to numerical linear algebra and random matrix theory. arXiv preprint arXiv:2112.12354.
}
\author{Xiucai Ding, Yichen Hu}
\keyword{distribution}
