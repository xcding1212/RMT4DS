\name{SignalNoise}
\alias{StepWiseSVD}
\alias{ScreeNot}
\alias{GetRank}
\alias{signal_value}
\alias{signal_vector}
\title{Signal-Plus-Noise Models}
\description{
    Estimation of signals, rank of signals.
}
\usage{
StepWiseSVD(Y, threshold=NULL, B=1000, level=0.02, methods='kmeans',
    u_threshold=NULL, v_threshold=NULL, sparse=TRUE)
ScreeNot(Y, r1)
GetRank(Y, r1, type=c("1","2"), level=0.1, B=500)
signal_value(d, svr)
signal_vector(k1, k2, d1, d2, svr, left=TRUE)
}
\arguments{
    \item{Y}{matrix to be denoised.}
    \item{B}{repeat time of simulations.}
    \item{threshold}{threshold used in determining rank of signal.}
    \item{level}{significance level in determing ranks.}
    \item{methods}{methods used in determining sparse structure.}
    \item{u_threshold, v_threshold}{thresholds used in determining sparse
        structure if kmeans is not used.}
    \item{sparse}{whether signals have sparse structure.}
    \item{r1}{upper bound of rank.}
    \item{type}{type of test.}
    \item{k1, k2}{k-th eigenvector.}
    \item{d, d1, d2}{eigenvalues of corresponding signal matrix}
    \item{left}{whether to use left singular vectors.}
    \item{svr}{ndf/ndim of Y.}
}
\value{
    \code{StepWiseSVD} performs step-wise SVD to denoise and returns decomposed strcuture,

    \code{ScreeNot} performs ScreeNot to denoise and returns decomposed strcuture,

    \code{GetRank} gives rank of signals.

    \code{signal_value} gives corrected signal eigenvalue from SVD result,

    \code{signal_vector} gives limiting inner product between signal vector and corresponding signal-plus-noise vector.
}
\details{
    \code{StepWiseSVD} works well in sparse setting and requires i.i.d normal noise and a lot simulation time.\code{SreeNot} is to pick the best TSVD result so works well in general setting.

    When using signal-plus-noise related limits, make sure they are limits of signal-related values or vectors.
}
\references{
    [1] Ding, X. (2020). High dimensional deformed rectangular matrices with applications in matrix denoising. Bernoulli, 26(1), 387-417.

    [2] Donoho, D. L., Gavish, M., & Romanov, E. (2020). Screenot: Exact mse-optimal singular value thresholding in correlated noise. arXiv preprint arXiv:2009.12297.

    [3] Ding, X., & Yang, F. (2022). Tracy-Widom distribution for heterogeneous Gram matrices with applications in signal detection. IEEE Transactions on Information Theory, vol. 68, no. 10, pp. 6682-6715.
}
\author{Xiucai Ding, Yichen Hu}
\keyword{estimation}
