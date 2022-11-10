\name{Limits}
\alias{MP_vector_dist}
\alias{cov_spike}
\alias{quadratic}
\title{Limits in High-dimensional Sample Covariance}
\description{
    Some limits of eigenvalues and eigenvectors in high-dimensional sample covariance.
}
\usage{
MP_vector_dist(k, v, ndf=NULL, pdim, svr=ndf/pdim, cov=NULL)
cov_spike(spikes, eigens, ndf, svr)
quadratic(k, cov, svr, spikes, type=1)
}
\arguments{
    \item{k}{k-th eigenvector. In \code{MP_vector_dist}, k can be a serie.}
    \item{v}{vector to be projected on.}
    \item{ndf}{the number of degrees of freedom for the Wishart matrix.}
    \item{pdim}{the number of dimensions (variables) for the Wishart matrix.}
    \item{svr}{samples to variables ratio; the number of degrees of freedom per dimension.}
    \item{cov}{population covariace matrix. If it is null, it will be regarded as identity.}
    \item{eigens}{input eigenvalues of population covariance matrix without spikes.}
    \item{spikes}{spikes in population covariance matrix.}
    \item{type}{transformation of eigenvalues. n for n-th power. 0 for logarithm.}
}
\value{
    \code{MP_vector_dist} gives asymptotic variance of projection of eigenvectors of non-spiked Wishart matrix,
    
    \code{cov_spike} gives spikes in sample covariance matrix and their asymptotic variance.
  
    \code{quadratic} gives mean of certain quadratic forms of k-th sample eigenvector of spiked models. Note k should be within the spikes.
}
\details{
    In \code{MP_vector_dist}, the variance computed is for \eqn{\sqrt{\code{pdim}}u_k^T v}, where \eqn{u_k} is the k-th eigenvector.
    
    Note in \code{quadratic}, k should be within the spikes.
}
\references{
    [1] Knowles, A., & Yin, J. (2017). Anisotropic local laws for random matrices. Probability Theory and Related Fields, 169(1), 257-352.
    
    [2] Jolliffe, I. (2005). Principal component analysis. Encyclopedia of statistics in behavioral science.
}
\examples{
k = 1
n = 200
p = 100
v = runif(p)
v = v/sqrt(sum(v^2))
MP_vector_dist(k,v,n,p,cov=diag(p))
cov_spike(c(10),rep(1,p),n,n/p)
quadratic(k,diag(p),n/p,c(30))
}
\author{Xiucai Ding, Yichen Hu}
\keyword{limit}
