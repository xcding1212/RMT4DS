\name{GeneralWM}
\alias{dWishartMax}
\alias{pWishartMax}
\alias{qWishartMax}
\alias{rWishartMax}
\title{The Wishart Maximum Eigenvalue Distribution}
\description{
    Density, distribution function, quantile function and random generation for the maximum eigenvalue from a general non-spiked Wishart matrix
    (sample covariance matrix) with \code{ndf} degrees of freedom,
    \code{pdim} dimensions, and order
    parameter \code{beta}.
}
\usage{
dWishartMax(x, eigens, ndf, pdim, beta, log = FALSE)
pWishartMax(q, eigens, ndf, pdim, beta, lower.tail = TRUE, log.p = FALSE)
qWishartMax(p, eigens, ndf, pdim, beta, lower.tail = TRUE, log.p = FALSE)
rWishartMax(n, eigens, ndf, pdim, beta)
}
\arguments{
    \item{x,q}{vector of quantiles.}
    \item{p}{vector of probabilities.}
    \item{n}{number of observations.}
    \item{eigens}{eigenvalues of population covariance matrix.}
    \item{ndf}{the number of degrees of freedom for the Wishart matrix}
    \item{pdim}{the number of dimensions (variables) for the Wishart matrix}
    \item{beta}{the order parameter. 1 for real Wishart and 2 for complex Wishart.}
    \item{log, log.p}{logical; if TRUE, probabilities p are given as log(p).}
    \item{lower.tail}{logical; if TRUE (default), probabilities are
        \eqn{P[X \le x]}, otherwise, \eqn{P[X > x]}.}
}
\value{
    \code{dWishartMax} gives the density,
    
    \code{pWishartMax} gives the distribution function,
    
    \code{qWishartMax} gives the quantile function,
    
    \code{rWishartMax} generates random deviates.  
}
\details{
    A real Wishart matrix is equal in distribution to 
    \eqn{X^T X/n}, where \eqn{X} are 
    \eqn{n\times p} real matrix with elements of mean zero and 
    covariance matrix \eqn{\Sigma}.
    A complex Wishart matrix is equal in distribution to 
    \eqn{X^* X/n}, where both real and imagety part of \eqn{X} are 
    \eqn{n\times p} complex matrice with elements of mean zero and 
    covariance matrix \eqn{\Sigma/2}. \code{eigens} are eigenvalues of 
    \eqn{\Sigma}. These functions give the limiting distribution of the largest 
    eigenvalue from the such a matrix when \code{ndf} and \code{pdim} both tend to 
    infinity.
}
\examples{
n = 500
p = 100
eigens = c(rep(2,p/2), rep(1, p/2))
beta = 2
rWishartMax(5, eigens, n, p, beta=beta)
qWishartMax(0.5, eigens, n, p, beta=beta)
pWishartMax(3.5, eigens, n, p, beta=beta)
dWishartMax(3.5, eigens, n, p, beta=beta)
}
\references{
    [1] El Karoui, N. (2007). Tracy–Widom limit for the largest eigenvalue of a large class of complex sample covariance matrices. The Annals of Probability, 35(2), 663-714.
    
    [2] Lee, J. O., & Schnelli, K. (2016). Tracy–Widom distribution for the largest eigenvalue of real sample covariance matrices with general population. The Annals of Applied Probability, 26(6), 3786-3839.
}
\author{Xiucai Ding, Yichen Hu}
\keyword{distribution}