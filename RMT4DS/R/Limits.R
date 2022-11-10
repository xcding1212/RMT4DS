#' @export
MP_vector_dist = function(k, v, ndf=NULL, pdim, svr=ndf/pdim, cov=NULL){
    if(is.null(cov)){
        sigma2 = rep(sum(v^2),length(k))
    }else{
        decompose = svd(cov)
        eigens = decompose$d
        U = decompose$u
        v = t(U)%*%v

        if(length(eigens)<1000){
            d_ = sample(eigens, 1000, replace=TRUE)
        } else{
            d_ = eigens
        }
        gmp = GeneralMarchenkoPasturPar(d=d_, phi=1/svr, m=10*pdim)
        dcdf = diff(gmp$cdf)
        xs = (gmp$Xs[1:(length(gmp$Xs)-1)]+gmp$Xs[2:length(gmp$Xs)])/2

        get_quantile = function(x){
            # Our computation is from left to right so x==1 is the special case.
            if(x==1){
                return(utils::tail(gmp$Xs, 1))
            }

            # The quantile we want will fall into interval [Xs[idx],Xs[idx+1])
            idx = utils::tail(which(gmp$cdf<=x),1)

            # Compute the approximate location between two quantiles
            s0 = x - gmp$cdf[idx]
            # Simple quadratic equation using the root formula.
            a = gmp$Ys[idx]
            b = gmp$Ys[idx+1]
            c = gmp$Xs[idx+1] - gmp$Xs[idx]
            if(b==a){
                k = s0/(a*c)
            }else{
                k = (-a*c+sqrt(a^2*c^2+2*s0*c*(b-a)))/(b-a)
            }
            gmp$Xs[idx] + c*k
        }
        sigma2 = c()
        for(i in 1:length(k)){
            eta = min(1/10/sum(eigens^(1.5))*k[i],0.005)
            gamma = get_quantile(1-k[i]/pdim) + eta*1i
            m_gamma = sum(1/(xs-gamma)*dcdf)/svr-(1-1/svr)/Re(gamma)
            sigma2 = c(sigma2, sum(1/svr*v^2*eigens/(Re(gamma)*Mod(1+m_gamma*eigens)^2)))
        }
    }
    list("variance"=sigma2)
}


#' @export
cov_spike = function(spikes, eigens, ndf, svr){
    M = length(eigens)
    N = ndf
    phi = M/N
    X = matrix(stats::rnorm(M*N,sd=sqrt(1/N)),nrow=M)
    T = diag(sqrt(eigens))
    if(M<N){
        sample_cov = T%*%X%*%t(X)%*%t(T)
    }else{
        sample_cov = t(X)%*%diag(eigens)%*%X
    }
    lamdas = eigen(sample_cov)$value

    # First we try to get the boundaries of sample spectral density
    d = eigens
    counter = as.data.frame(table(d))
    counter$d = as.numeric(as.character(counter$d))
    counter$Freq = counter$Freq/M
    counter = counter[order(counter$d, decreasing=TRUE),]

    # The y values of critical points of f are the boundaries of sample spectral density
    f = function(x){
        out = -1/x
        for(i in 1:nrow(counter)){
            out = out + phi*counter$Freq[i]/(x+1/counter$d[i])
        }
        out
    }
    # Derivative of f
    fd = function(x){
        out = 1/x^2
        for(i in 1:nrow(counter)){
            out = out - phi*counter$Freq[i]/(x+1/counter$d[i])^2
        }
        out
    }
    divides = sort(-1/counter$d)
    edges = c()
    # Introduce a epsilon term to avoid 1/0 in computation.
    epsilon = 1e-7
    solution = rootSolve::uniroot.all(fd, lower=-100+epsilon, upper=divides[1]-epsilon)
    if(length(solution)>0){
        edges = c(edges, solution)
    }
    if(length(divides)>1){
        for(i in 2:(nrow(counter))){
            solution = rootSolve::uniroot.all(fd, lower=divides[i-1]+epsilon, upper=divides[i]-epsilon)
            if(length(solution)>0){
                edges = c(edges, solution)
            }
        }
    } else{
        i = 1
    }
    solution = rootSolve::uniroot.all(fd, lower=divides[i]+epsilon, upper=-epsilon)
    if(length(solution)>0){
        edges = c(edges, solution)
    }
    solution = rootSolve::uniroot.all(fd, lower=0+epsilon, upper=10)
    if(length(solution)>0){
        edges = c(edges, solution)
    }

    # compute spikes
    SPK = c()
    SPK_cov = c()
    for(i in 1:length(spikes)){
        spk = spikes[i]
        idx = which(edges>(-1/spk))[1]
        if(is.na(idx)){
            if(-1/spk>=(utils::tail(edges,1)+0*N^(-1/3))){
                SPK = c(SPK, rep(f(-1/spk)))
                SPK_cov = c(SPK_cov, 2*f(-1/spk)^2/(N-1))
            }
        }else{
            if(idx%%2==1|(-1/spk>=(edges[idx-1]+0*N^(-1/3)))){
                SPK = c(SPK, rep(f(-1/spk)))
                SPK_cov = c(SPK_cov,2*f(-1/spk)^2/(N-1))
            }
        }
    }
    list("spikes"=SPK, "variance"=SPK_cov)
}


#' @export
quadratic = function(k, cov, svr, spikes, type=1){
    n_spike = length(spikes)
    c = 1/svr
    p = nrow(cov)
    N = p/c
    decompose = svd(cov)
    eigens = decompose$d
    d = c()
    if(n_spike > 0){
        for(i in 1:n_spike){
            d = c(d, spikes[i]/eigens[i]-1)
        }
    }
    f = function(x){
        -1/x + c*mean(1/(x+1/eigens))
    }
    fd = function(x){
        1/x^2 - c*mean(1/(x+1/eigens)^2)
    }
    m = function(x){
        mean(1/(eigens-x))
    }
    md = function(x){
        mean(1/(eigens-x)^2)
    }
    g = function(x){
        if(type==0){
            log(x)
        }else{
            x^type
        }
    }
    if(k<=n_spike){
        ai = f(-1/spikes[k])
        bi = fd(-1/spikes[k])/f(-1/spikes[k])
        di = d[k]
        ci = di^2/spikes[k]^2*bi
        psai = g(spikes[k])/spikes[k]*bi+ci/(1/eigens[k]+m(ai))^2*
            1/N*sum(g(eigens[(1+n_spike):p])/eigens[(1+n_spike):p]*md(ai)/
                        (1/eigens[(1+n_spike):p]+m(ai))^2)
        psai
    } else{
        stop("Only spiked part is supported.")
    }
}
