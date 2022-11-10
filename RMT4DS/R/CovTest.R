#' @include CovEst.R


d2 = function(y1, y2) {
    (y1+y2-y1*y2)/(y1*y2)*log((y1+y2)/(y1+y2-y1*y2))+
        y1*(1-y2)/(y2*(y1+y2))*log(1-y2)+y2*(1-y1)/(y1*(y1+y2))*log(1-y1)
}
mu2 = function(y1, y2) 1/2*log((y1+y2-y1*y2)/(y1+y2))-(y1*log(1-y2)+y2*log(1-y1))/(y1+y2)
sigma2_2 = function(y1, y2) {
    -(2*y1^2*log(1-y2)+2*y2^2*log(1-y1))/(y1+y2)^2-
        2*log((y1+y2)/(y1+y2-y1*y2))
}
d2_hat = function(y1, y2){
    d2(y1,y2)-y1/(y1+y2)*log(y1/(y1+y2))-y2/(y1+y2)*log(y2/(y1+y2))
}
d3 = function(y1, y2, alpha, beta, h){
    a = (1-y2)/y2*log(alpha/(1-y2))-(y1+y2)/(y1*y2)*log((h*alpha-y2*beta)/(h*(1-y2)))
    if(y1>=1){
        a+(y1-1)/y1*log((alpha-beta/h)/(1-y2))
    }else{
        a+(1-y1)/y1*log((alpha-h*beta)/(1-y2))
    }
}

#' @export
OneSampleCovTest = function(X, mean=NULL, S=NULL){
    # X: input data matrix n by p
    # mu
    # S: covariance matrix for test, default identity
    n = dim(X)[1]
    p = dim(X)[2]
    y = p/n
    N = n-1
    yN = p/N
    if(!is.null(S)){
        S_half = chol(S)
        X = X%*%solve(S_half)
    }
    if(is.null(mean)){
        X = X - tcrossprod(rep(1,n), colMeans(X))
        S = t(X)%*%X/N
    } else{
        X = X - tcrossprod(rep(1,n), mean)
        S = t(X)%*%X/n
    }
    lrt = sum(diag(S))-log(det(S))-p
    mu1 = -1/2*log(1-y)
    sigma1 = -2*log(1-y)-2*y
    z_value = (lrt-p*(1+(1-yN)/yN*log(1-yN))-mu1)/sqrt(sigma1)
    p_value = stats::pnorm(z_value, lower.tail=FALSE)
    list("p_value"=p_value, "z_value"=z_value, "lrt"=lrt)
}


#' @export
TwoSampleCovTest = function(X1, X2, mean=NULL){
    # X1, X2: input data matrix n1 by p and n2 by p
    if(dim(X1)[2]!=dim(X2)[2]) stop('Input data should have the same dimension')
    n1 = dim(X1)[1]
    n2 = dim(X2)[1]
    if(!is.null(mean)){
        N1 = n1
        N2 = n2
        X1 = X1 - tcrossprod(rep(1,n1), mean)
        X2 = X2 - tcrossprod(rep(1,n2), mean)
    } else{
        N1 = n1-1
        N2 = n2-1
        X1 = X1 - tcrossprod(rep(1,n1), colMeans(X1))
        X2 = X2 - tcrossprod(rep(1,n2), colMeans(X2))
    }
    N = N1+N2
    c1 = N1/N
    c2 = N2/N
    p = dim(X1)[2]
    yN1 = p/N1
    yN2 = p/N2
    S1 = t(X1)%*%X1/N1
    S2 = t(X2)%*%X2/N2

    log_V1_ = log(det(S1%*%solve(S2)))*(N1/2)-
        log(det(c1*S1%*%solve(S2)+c2*diag(rep(1,p))))*(N/2)
    z_value = (-2/N*log_V1_-p*d2(yN1,yN2)-mu2(yN1,yN2))/sqrt(sigma2_2(yN1,yN2))
    p_value = stats::pnorm(z_value, lower.tail=FALSE)

    list("p_value"=p_value, "z_value"=z_value, "lrt"=-2*log_V1_)
}


#' @export
MultiSampleCovTest = function(..., input=NULL){
    # Xi: input data matrix n1 by p and n2 by p
    if(is.null(input)){
        matrices = list(...)
    } else{
        matrices = input
    }
    ps = c()
    Ns = c()
    q = length(matrices)
    for(i in 1:q){
        ps = c(ps, dim(matrices[[i]])[2])
        Ns = c(Ns, dim(matrices[[i]])[1]-1)
    }
    n0 = max(Ns)
    if(mean(ps)!=ps[1]) stop('Input data should have the same dimension')
    p = ps[1]
    As = array(0, dim=c(q,p,p))
    for(i in 1:q){
        As[i,,] = stats::cov(matrices[[i]])*Ns[i]/n0
    }
    y_n1 = c()
    y_n2 = c()
    logV1 = c()
    for(i in 2:q){
        y_n1 = c(y_n1, p/sum(Ns[1:(i-1)]))
        y_n2 = c(y_n2, p/Ns[i])
        if(i==2){
            logV1 = c(logV1, Ns[1]/2*log(det(As[1,,]))+
                          Ns[2]/2*log(det(As[2,,]))-
                          sum(Ns[1:2])/2*log(det(colSums(As[1:2,,],dims=1))))
        } else {
            logV1 = c(logV1, sum(Ns[1:(i-1)])/2*log(det(colSums(As[1:(i-1),,],dims=1)))+
                          Ns[i]/2*log(det(As[i,,]))-
                          sum(Ns[1:i])/2*log(det(colSums(As[1:i,,],dims=1))))
        }
    }
    test = 0
    sigma2 = 0
    mu = 0
    lrt = 0

    for(i in 2:q){
        test = test-2/sum(Ns[1:i])*logV1[i-1]-p*d2_hat(y_n1[i-1],y_n2[i-1])
        lrt = lrt-2/sum(Ns[1:i])*logV1[i-1]
        mu = mu + mu2(y_n1[i-1],y_n2[i-1])
        sigma2 = sigma2 + sigma2_2(y_n1[i-1],y_n2[i-1])
    }
    z_value = (test-mu)/sqrt(sigma2)
    p_value = stats::pnorm(z_value, lower.tail=FALSE)
    list("p_value"=p_value, "z_value"=z_value, "lrt"=lrt)
}
