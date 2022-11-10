\name{CovTest}
\alias{CovarianceTest}
\alias{OneSampleCovTest}
\alias{TwoSampleCovTest}
\alias{MultiSampleCovTest}
\title{High-dimensional Covariance Test}
\description{
    Test of given population covariance matrix, test of equal covariance of two or more
    samples.
}
\usage{
OneSampleCovTest(X, mean=NULL, S=NULL)
TwoSampleCovTest(X1, X2, mean=NULL)
MultiSampleCovTest(..., input=NULL)
}
\arguments{
    \item{X, X1, X2}{input samples in the form n by p where p is the dimension.}
    \item{mean}{population mean of samples. If it is missing,
        sample mean will be used.}
    \item{S}{covariance matrix to be tested. If it is missing, test of
        identity covariance will be performed.}
    \item{...}{any samples to be tested.}
    \item{input}{list of samples to be tested. Please choose either \code{...}
        or \code{input} as input form.}
}
\value{
    \code{OneSampleCovTest} tests given covariance matrix of one sample,

    \code{TwoSampleCovTest} tests equal covariance matrices of two samples,

    \code{MultiSampleCovTest} tests equal covariance matrices of multiple samples.
}
\source{
    Maximal likelihood tests fail in high-dimensional settings, so corrections are
    made. Note all tests are one-sided. Large statistics indicate violation of
    null hypothesis.
}
\examples{
require(MASS)
n = 500
p = 100
S1 = diag(rep(1,p))
S2 = diag(sample(c(1,4),p,replace=TRUE))
OneSampleCovTest(mvrnorm(n,rep(0,p),S2), S=S1)
TwoSampleCovTest(mvrnorm(n,rep(0,p),S1), mvrnorm(n,rep(0,p),S2))
MultiSampleCovTest(mvrnorm(n,rep(0,p),S1), mvrnorm(n,rep(0,p),S2))
}

\references{
    [1] Zheng, S., Bai, Z., & Yao, J. (2015). Substitution principle for CLT of linear spectral statistics of high-dimensional sample covariance matrices with applications to hypothesis testing. The Annals of Statistics, 43(2), 546-591.
}
\author{Xiucai Ding, Yichen Hu}
\keyword{test}
