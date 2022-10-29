#' @include GeneralMP.R


scaling = function(eigens, beta, d){
    ls = eigens
    # Use binary search to solve the equation
    f1 = function(x){
        mean((x*ls/(1-x*ls))^2)-1/d
    }
    epsilon = 1e-7
    high = 1/max(ls)-epsilon
    low = 0
    while(abs(high-low)>1e-10){
        if(f1(low)*f1((low+high)/2)<=0){
            high = (low+high)/2
        } else{
            low = (low+high)/2
        }
    }
    c = (low+high)/2
    # Compute the mean and variance of largest eigenvalue.
    mu = 1/c*(1+d*mean(c*ls/(1-c*ls)))
    sigma = (1/c^3*(1+d*mean((c*ls/(1-c*ls))^3)))^(1/3)
    list("mu"=mu, "sigma"=sigma)
}


#' @export
rWishartMax = function(n, eigens, ndf, pdim, beta){
    # generate largest eigenvalues of sample covariance based on population covariance
    out = scaling(eigens, beta, pdim/ndf)
    RMTstat::rtw(n, beta=beta)/ndf^(2/3)*out$sigma + out$mu
}


#' @export
qWishartMax = function(p, eigens, ndf, pdim, beta, lower.tail=TRUE, log.p=FALSE){
    # quantile of largest eigenvalues of sample covariance based on population covariance
    out = scaling(eigens, beta, pdim/ndf)
    mu = out$mu
    sigma = out$sigma
    if(log.p==TRUE){
        p = exp(p)
    }
    RMTstat::qtw(p, beta=beta, lower.tail=lower.tail)/ndf^(2/3)*sigma + mu
}


#' @export
pWishartMax = function(q, eigens, ndf, pdim, beta, lower.tail=TRUE, log.p=FALSE){
    # cdf of largest eigenvalues of sample covariance based on population covariance
    out = scaling(eigens, beta, pdim/ndf)
    mu = out$mu
    sigma = out$sigma
    q = (q-mu)/sigma*ndf^(2/3)
    RMTstat::ptw(q, beta=beta, lower.tail=lower.tail, log.p=log.p)
}


#' @export
dWishartMax = function(x, eigens, ndf, pdim, beta, log=FALSE){
    # pdf of largest eigenvalues of sample covariance based on population covariance
    out = scaling(eigens, beta, pdim/ndf)
    mu = out$mu
    sigma = out$sigma
    x = (x-mu)/sigma*ndf^(2/3)
    RMTstat::dtw(x, beta=beta, log=log)
}
