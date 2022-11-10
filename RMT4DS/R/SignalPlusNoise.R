#' @export
StepWiseSVD = function(Y, threshold=NULL, B=1000, level=0.02,
                       methods='kmeans', u_threshold=NULL, v_threshold=NULL,
                       sparse=TRUE){
    n = nrow(Y)
    p = ncol(Y)
    C = n/p

    simulation = function(x){
        Matrix_sim = MASS::mvrnorm(x[1], rep(0,x[2]), diag(x[2])/x[2])
        eigs = rARPACK::svds(Matrix_sim, k=2, nu=0, nv=0)$d
        (eigs[1]/eigs[2])^2
    }

    if(is.null(threshold)){
        # generate B random matrices and record the ratios between their
        # first and second eigenvalues
        ratios = c()
        for(i in 1:B){
            ratios = c(ratios, simulation(c(n,p)))
        }

        # using 1-level percentile to get the threshold to pick q
        threshold = stats::quantile(ratios, probs=1-level)
        threshold = unname(threshold)
    }

    # perform SVD to determine q
    mu = (svd(Y)$d)**2
    mu_ratio = mu[1:(min(n,p)-1)]/mu[2:min(n,p)]
    q = utils::tail(which(mu_ratio[1:floor(min(n,p)/2)]>threshold), 1)

    p1 = function (d, a, c) (d**2+1)*(d**2+1/c)/d**2-a**2
    a1 = function(d, c) (d**4-1/c)/(d**2*(d**2+1/c))
    a2 = function(d, c) (d**4-1/c)/(d**2*(d**2+1))

    us = matrix(nrow=n, ncol=0)
    vs = matrix(nrow=p, ncol=0)
    ds = c()

    # main loop
    for(i in 1:q){
        # just get the first component
        decomp = rARPACK::svds(Y, k=1, nu=1, nv=1)
        t = decomp$d
        u = decomp$u
        v = decomp$v
        # compute p^(-1)(t_j^2)
        d = nleqslv::nleqslv(t, p1, a=t, c=C)$x

        # incorporate both sparse and non-sparse case
        if(sparse==TRUE){
            if(methods=='kmeans'){
                kmeans_u = stats::kmeans(abs(u), 2)
                kmeans_v = stats::kmeans(abs(v), 2)
                if(kmeans_u$centers[1]>1/sqrt(n)){
                    u_cluster = 1
                } else {
                    u_cluster = 2
                }
                if(kmeans_v$centers[1]>1/sqrt(p)){
                    v_cluster = 1
                } else {
                    v_cluster = 2
                }
                u_idx = which(kmeans_u$cluster==u_cluster)
                v_idx = which(kmeans_v$cluster==v_cluster)
            } else{
                u_idx = which(abs(u)>=u_threshold)
                v_idx = which(abs(v)>=v_threshold)
            }

            # obtain submatrix to correct eigenvectors
            submatrix = as.matrix(Y[u_idx, v_idx])
            decomp_sub = svd(submatrix, nu=1, nv=1)
            u_sub = decomp_sub$u
            v_sub = decomp_sub$v
            u_new = rep(0, n)
            v_new = rep(0, p)
            u_new[u_idx] = u_sub
            v_new[v_idx] = v_sub

            ds = c(ds, d)
            us = cbind(us, u_new)
            vs = cbind(vs, v_new)
            Y = Y - d * u_new %*% t(v_new)
        } else{
            eta = d*a1(d,C)*a2(d,C)
            ds = c(ds, eta)
            us = cbind(us, u)
            vs = cbind(vs, v)
            Y = Y - t * u %*% t(v)
        }
    }
    dimnames(us) = NULL
    dimnames(vs) = NULL
    list("r"=q, "d"=ds, "u"=us, "v"=vs, "threshold"=threshold)
}


#' @export
ScreeNot = function(Y, r1){
    n = nrow(Y)
    p = ncol(Y)
    gamma = p/n
    decompose = svd(Y, nu=r1, nv=r1)
    fY = decompose$d
    fZ = fY
    for(i in 1:r1){
        fZ[i] = fY[r1+1]+(fY[r1+1]-fY[2*r1+1])*(1-((i-1)/r1)**(2/3))/(2**(2/3)-1)
    }
    psi = function(y){
        phi = mean(y/(y**2-fZ**2))
        phid = -mean((y**2+fZ**2)/(y**2-fZ**2)**2)
        D = gamma*phi+(1-gamma)/y
        Dd = gamma*phid-(1-gamma)/y**2
        y*(phid/phi+Dd/D) + 4
    }
    thre = abs(nleqslv::nleqslv(fY[2], psi)$x)
    d = fY[fY>thre]
    r = length(d)
    u = decompose$u[,1:r]
    v = decompose$v[,1:r]
    list("r"=r, "threshold"=thre, "u"=u, "v"=v, "d"=d)
}


TestRank = function(Y, r1, r0, type=c("1","2"), level=0.1, threshold=NULL, B=5000){
    n = nrow(Y)
    p = ncol(Y)

    simulation = function(x, r1){
        Matrix_sim = matrix(stats::rnorm(x[1]*x[2], mean=0, sd=sqrt(1/x[2])), nrow=x[1], byrow=TRUE)
        rate = (rARPACK::svds(Matrix_sim, k=(r1+2))$d)^2
        rate1 = (rate[1:r1]-rate[2:(1+r1)])/(rate[2:(1+r1)]-rate[3:(2+r1)])
        rate2 = (rate[1]-rate[2])/(rate[1:r1]-rate[2:(1+r1)])
        list("rate1"=rate1, "rate2"=rate2)
    }

    if(is.null(threshold)){
        rates1 = matrix(0, nrow=B, ncol=r1)
        rates2 = matrix(0, nrow=B, ncol=r1)
        for(i in 1:B){
            obj = simulation(c(p,n), r1)
            rates1[i,] = obj$rate1
            rates2[i,] = obj$rate2
        }
    }
    # test given rank of signal
    ds = (rARPACK::svds(Y, k=(r1+2))$d)^2
    if(type=="1"){
        a = max((ds[(r0+1):r1]-ds[(r0+2):(r1+1)])/(ds[(r0+2):(r1+1)]-ds[(r0+3):(r1+2)]))
        if(is.null(threshold)){
            if(r0==(r1-1)){
                b = stats::quantile(rates1[,1], 1-level)
            } else{
                b = stats::quantile(apply(rates1[,1:(r1-r0)],1,max), 1-level)
            }
        } else{
            b = threshold
        }
        if(a<b){
            result=1
        }else{
            result=0
        }
    }
    if(type=="2"){
        a = max((ds[r0+1]-ds[r0+2])/(ds[r1+1]-ds[r1+2]))
        if(is.null(threshold)){
            if(r0==(r1-1)){
                b = stats::quantile(rates2[,1], 1-level)
            } else{
                b = stats::quantile(apply(rates2[,1:(r1-r0)],1,max), 1-level)
            }
        } else{
            b = threshold
        }
        if(a<b){
            result=1
        }else{
            result=0
        }
    }
    attr(r0, "names") = "rank of signal"
    out = list(method="Signal Plus Noise Rank Test",
               null.value=r0,
               data.name="ratio",
               p.value=result,
               parameter=NULL,
               statistics=a,
               alternative="greater")
    class(out) = "htest"
    out
}


#' @export
GetRank = function(Y, r1, type=c("1","2"), level=0.1, B=500){
    # using test of rank to get rank
    for(r0 in 0:(r1-1)){
        if(TestRank(Y=Y, r1=r1, r0=r0, type=type, level=level, B=B)$p.value==1){
            break
        }
    }
    list("rank"=r0)
}


#' @export
signal_value = function(d, svr){
    p = function (d, a, c) (d^2+1)*(d^2+1/c)/d^2-a^2
    nleqslv::nleqslv(d, p, a=d, c=svr)$x
}


#' @export
signal_vector = function(k1, k2, d1, d2, svr, left=TRUE){
    if(k1!=k2){
        0
    }else{
        if(left==TRUE){
            (d1^4-1/svr)/(d1^2*(d1^2+1/svr))
        }else{
            (d1^4-1/svr)/(d1^2*(d1^2+1))
        }
    }
}
