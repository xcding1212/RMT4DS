\name{CovEst}
\alias{MPEst}
\alias{MomentEst}
\title{Estimation of the Spectrum of Population Covariance Matrix}
\description{
    Estimation of the eigenvalues of population covariance matrix given samples.
}
\usage{
MPEst(X, n=nrow(X), k=1, num=NULL, penalty=FALSE, n_spike=0)
MomentEst(X, n=nrow(X), k=1, n_spike=0)
}
\arguments{
    \item{X}{n by p data matrix.}
    \item{n}{sample size.}
    \item{k}{repeated times in estimation. If \code{k>1}, estimation will be
        the average.}
    \item{num}{numbers of mass points chosen in estimation.}
    \item{penalty}{whether to implement L-1 penalty in inverting Marchenko-Pastur
        law}
    \item{n_spike}{number of spikes in population spectral.}
}
\value{
    \code{MPEst} and \code{MomentEst} give estimation of the spectrum of population covariance matrix and corresponding spectral density.
}
\details{
    Given \eqn{E(X)=0} and \eqn{Cov(X)=\Sigma} with \eqn{\Sigma} unknown and fourth moment of X exists, we want to estimate spectrum of \eqn{\Sigma} from sample covariance matrix \eqn{X'X/n}.

    \code{MPEst} estimates spectrum by inverting Marchenko-Pastur law while \code{MomentEst} estimates spectrum
    by estimating the moment of population spectral density.

    Those two functions give estimates of the eigenvalues by \code{d} and estimates of spectral density by \code{xs} and \code{cdf}.
}
\examples{
require(MASS)
n = 500
p = 250
X = mvrnorm(n, rep(0,p), diag(c(rep(2,p/2),rep(1,p/2))))
MPEst(X, n)$d
MomentEst(X, n)$d
}

\references{
    [1] El Karoui, N. (2008). Spectrum estimation for large dimensional covariance matrices using random matrix theory. The Annals of Statistics, 36(6), 2757-2790.

    [2] Kong, W., & Valiant, G. (2017). Spectrum estimation from samples. The Annals of Statistics, 45(5), 2218-2247.
}
\author{Xiucai Ding, Yichen Hu}
\keyword{estimation}
