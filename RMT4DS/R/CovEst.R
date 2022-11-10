#' @include Limits.R


# invert m-p equation
MPEstOne = function(X, sigma=1, penalty=FALSE, num=NULL){
    n = nrow(X)
    p = ncol(X)
    c = p/n
    S = t(X)%*%X/sigma**2/n
    sample_eigen = svd(S)$d
    l_max = max(sample_eigen)
    l_min = min(sample_eigen)
    # l_min<=lambda_-1<=lambda_1<=l_max, due to convexity
    start_point = l_min / l_max
    if(is.null(num)){
        ts = seq(start_point, 1, by=0.005)*l_max
    } else {
        ts = seq(start_point, 1, by=(1-start_point)/(num-1))*l_max
    }

    # pick v(z_j) and then solve for corresponding z_j
    vt = c()
    xs = seq(0.02, 1, 0.02)
    for(i in 1:length(xs)){
        vt = c(vt, xs[i]+0.01*1i, xs[i]-0.001*1i)
    }
    # convert complex computation into real and image parts
    v_F = function(z, t){
        x = z[1]
        y = z[2]
        real = mean((sample_eigen-x)/((sample_eigen-x)**2+y**2))*c - (1-c)*x/(x**2+y**2) - t[1]
        image = mean(y/((sample_eigen-x)**2+y**2))*c - (1-c)*(-y)/(x**2+y**2) - t[2]
        c(real, image)
    }
    v_f = function(z){
        mean(1/(sample_eigen-z))*c-(1-c)/z
    }

    zs = c()
    for(i in 1:length(vt)){
        sol = nleqslv::nleqslv(c(1,1), v_F, t=c(Re(vt[i]),Im(vt[i])))$x
        zs = c(zs, sol[1]+sol[2]*1i)
    }

    z_t = matrix(0, nrow=length(zs), ncol=length(ts))
    for(j in 1:length(vt)){
        for(k in 1:length(ts)){
            z_t[j,k] = c*ts[k]/(1+v_f(zs[j])*ts[k])
        }
    }

    z_t_re = Re(z_t)
    z_t_im = Im(z_t)

    const = c()
    for(i in 1:length(vt)){
        const = c(const, 1/v_f(zs[i])+zs[i])
    }
    # const = 1/vt+zs
    const_re = Re(const)
    const_im = Im(const)

    # perform first order penalty
    f.obj = c(1, rep(0, length(ts)))
    if(penalty==TRUE){
        f.rhs = c(const_re, -const_re, const_im, -const_im, 1,-1,
                  mean(sample_eigen), -mean(sample_eigen))
    } else{
        f.rhs = c(const_re, -const_re, const_im, -const_im, 1,-1)
    }
    f.con = cbind(-1, z_t_re)
    f.con = rbind(f.con, cbind(-1, -z_t_re))
    f.con = rbind(f.con, cbind(-1, z_t_im))
    f.con = rbind(f.con, cbind(-1, -z_t_im))
    f.con = rbind(f.con, c(0, rep(1, length(ts))))
    f.con = rbind(f.con, c(0, rep(-1, length(ts))))
    if(penalty==TRUE){
        f.con = rbind(f.con, c(0,ts))
        f.con = rbind(f.con, -c(0,ts))
    }
    f.dir = rep("<=", dim(f.con)[1])

    lp_obj = lpSolve::lp("min", f.obj, f.con, f.dir, f.rhs)
    values = lp_obj$solution[2:length(lp_obj$solution)]
    list("xs"=ts, "pdf"=values)
}


#' @export
MPEst = function(X, n=nrow(X), k=1, num=NULL, penalty=FALSE, n_spike=0){
    X1 = X[sample(1:nrow(X), n, replace=FALSE),]
    decompose = svd(t(X1)%*%X1/n)
    u = decompose$u
    d1 = decompose$d

    if(n_spike > 0){
        X_spike = X[,1:n_spike]
        X = X[,(n_spike+1):ncol(X)]
    }
    xs = c()
    pdf = c()
    for(i in 1:k){
        X1 = X[sample(1:nrow(X), n, replace=FALSE),]
        obj = MPEstOne(X, num=num, sigma=1, penalty=penalty)
        xs = c(xs, obj$xs)
        pdf = c(pdf, obj$pdf/k)
    }
    pdf = pdf[order(xs, decreasing=FALSE)]
    xs = xs[order(xs, decreasing=FALSE)]
    cdf = cumsum(pdf)
    p = ncol(X)
    d = c()
    for(i in 1:p){
        q = (p-i+1/2)/p
        d = c(d, xs[which(cdf>=q)[1]])
    }
    if(n_spike > 0){
        lamdas = d1[(n_spike+1):length(d1)]
        d_spike = c()
        for(i in 1:n_spike){
            sigma_i = -1/(p/n*mean(1/(lamdas-d1[i]))+(1-p/n)*1/(0-d1[i]))
            d_spike = c(d_spike, sigma_i)
        }
        d = c(d_spike, d)
    }
    list("d"=d, "xs"=xs, "cdf"=cdf)
}


# method of moment
MomentEstOne = function(Y, k=7){
    n = nrow(Y)
    d = ncol(Y)
    moments = c()
    A = Y%*%t(Y)
    G = A
    G[lower.tri(G,diag=TRUE)] = 0
    Gk = diag(n)
    xs = seq(0, rARPACK::svds(t(Y)%*%Y/n, k=1, nu=0, nv=0)$d[1],1/max(d,n))
    V = matrix(0, nrow=k, ncol=length(xs))

    k_moment = function(n, d, k, A, Gk){
        sum(diag(Gk%*%A))/(d*choose(n,k))
    }

    for(i in 1:k){
        moments = c(moments, k_moment(n, d, i, A, Gk))
        V[i,] = xs^i
        if(i<k){
            Gk = Gk%*%G
        }
    }
    f.obj = c(rep(0, length(xs)), rep(1,k))
    f.con = c(rep(1, length(xs)), rep(0,k))
    f.con = rbind(f.con, c(rep(-1, length(xs)), rep(0,k)))
    f.con = rbind(f.con, cbind(V, -diag(k)))
    f.con = rbind(f.con, cbind(-V, -diag(k)))
    f.rhs = c(1,-1,moments,-moments)
    f.dir = rep("<=", dim(f.con)[1])
    lp_obj = lpSolve::lp("min", f.obj, f.con, f.dir, f.rhs)
    pdf = lp_obj$solution[1:length(xs)]
    list("xs"=xs, "cdf"=cumsum(pdf), "pdf"=pdf)
}

#' @export
MomentEst = function(X, n=nrow(X), k=1, n_spike=0){
    X1 = X[sample(1:nrow(X), n, replace=FALSE),]
    decompose = svd(t(X1)%*%X1/n)
    u = decompose$u
    d1 = decompose$d

    if(n_spike > 0){
        X_spike = X[,1:n_spike]
        X = X[,(n_spike+1):ncol(X)]
    }
    xs = c()
    pdf = c()
    for(i in 1:k){
        obj = MomentEstOne(X[sample(1:nrow(X), n, replace=FALSE),])
        xs = c(xs, obj$xs)
        pdf = c(pdf, obj$pdf/k)
    }
    pdf = pdf[order(xs, decreasing=FALSE)]
    xs = xs[order(xs, decreasing=FALSE)]
    cdf = cumsum(pdf)

    p = ncol(X)
    X1 = X[sample(1:nrow(X), n, replace=FALSE),]
    u = svd(t(X1)%*%X1/n)$u
    d = c()
    for(i in 1:p){
        q = (p-i+1/2)/p
        d = c(d, xs[which(cdf>=q)[1]])
    }
    if(n_spike > 0){
        lamdas = d1[(n_spike+1):length(d1)]
        d_spike = c()
        for(i in 1:n_spike){
            sigma_i = -1/(p/n*mean(1/(lamdas-d1[i]))+(1-p/n)*1/(0-d1[i]))
            d_spike = c(d_spike, sigma_i)
        }
        d = c(d_spike, d)
    }
    list("d"=d, "xs"=xs, "cdf"=cdf)
}
