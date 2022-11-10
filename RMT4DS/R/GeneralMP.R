GeneralMarchenkoPasturPar = function(d, phi, m=500){
    M = length(d)
    N = M/phi
    if(M==N) stop("N and M cannot be the same!")
    X = matrix(stats::rnorm(M*N,sd=sqrt(1/N)),nrow=M)
    T = diag(sqrt(d))
    if(M<N){
        sample_cov = T%*%X%*%t(X)%*%t(T)
    }else{
        sample_cov = t(X)%*%diag(d)%*%X
    }
    lamdas = eigen(sample_cov)$value

    # First we try to get the boundaries of sample spectral density
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
    edge=sort(f(edges))
    as = edge[c(TRUE, FALSE)]
    bs = edge[c(FALSE, TRUE)]

    if(length(as)==1){
        l = 12
    }else{
        l = 4
    }

    # Then we try to use orthogonal polynomials to estimate density.
    # lim_{eta->0}r(x+eta*i)
    # We approximate r(z), and this is equivalent to a quadratic programming.
    k = 20
    g = length(as)
    # In construction of z's, we set imagery part to be 0.1 and uniformly choose m
    # numbers in every interval of the support as real part.

    # xk's are choose from the roots of k-th second type of chebyshev polynomial.
    # Ek is for restriction.
    Xks = suppressMessages(rootSolve::uniroot.all(as.function(mpoly::chebyshev(k, "u")),
                                                  lower=-1, upper=1))
    funcs = list()
    for(i in 1:l){
        funcs = suppressMessages(append(funcs, list(as.function(mpoly::chebyshev(i-1, "u")))))
    }
    Ek = matrix(0, nrow=k, ncol=l)
    for(i in 1:k){
        for(j in 1:l){
            Ek[i,j] = funcs[[j]](Xks[i])
        }
    }


    # CM is the basis matrix for our fit.
    C_M = function(z,a,b,j){
        z1 = (z-(a+b)/2)/((b-a)/2)
        -2*(z1-sqrt(z1-1)*sqrt(z1+1))^(j+1)*(b-a)*pi/4
    }

    # Computing r(z)
    # r(z)'s are approximation of stieljes transformations.
    stj = function(z){
        mean(1/(lamdas-z))
    }

    zs = c()
    for(i in 1:g){
        zs = c(zs, seq(-1,1,2/(m-1))*(bs[i]-as[i])/2+(bs[i]+as[i])/2+0.1i)
    }
    rz = unlist(lapply(zs, stj))


    CM = matrix(0, nrow=m*g, ncol=l*g)
    for(j in 1:g){
        for(i in 1:(m*g)){
            for(t in 1:l){
                CM[i,(t+(j-1)*l)] = C_M(zs[i], as[j], bs[j], t-1)
            }
        }
    }

    # We want to minimize the norm of the difference. It is more convenient to
    # convert the loss function into quadratic loss and we can drop the constant term.
    A1 = Re(CM)
    A2 = Im(CM)
    y1 = Re(rz)
    y2 = Im(rz)

    # Constraints Ek*dj>=0 <=> -Ek*dj<=0
    A = matrix(0, nrow=g*k, ncol=g*l)
    for(i in 1:g){
        A[((i-1)*k+1):(i*k),((i-1)*l+1):(i*l)] = -Ek
    }

    # Solve the quadratic programming.
    ds = pracma::quadprog(t(A1)%*%A1+t(A2)%*%A2, -as.vector(t(A1)%*%y1+t(A2)%*%y2),
                          A, rep(0,k*g))$xmin

    # Then we set the imagery part to be 0 and compute the fit.
    # Notice due to numerical instability, small negative results may be produced,
    # so it is better to take the absolute value.
    # But that is not a big issue according to our experiments since it only influence
    # a very small part on the boundaries whose density is almost zero.
    Xs = Re(zs)
    CM1 = matrix(0, nrow=m*g, ncol=l*g)
    for(j in 1:g){
        for(i in 1:(m*g)){
            for(t in 1:l){
                CM1[i,(t+(j-1)*l)] = C_M(Xs[i]+0i, as[j], bs[j], t-1)
            }
        }
    }
    Ys = abs(Im(CM1%*%ds)/pi)

    # Finally we try to compute CDF of the density.
    # We compute the trapezoid area as the density of every small interval.
    # Then we normalize the pdf to make sure its area is 1.
    delta = diff(Xs)
    areas = (Ys[1:(g*m-1)]+Ys[2:(g*m)])*delta/2
    Ys = Ys/sum(areas)
    areas = areas/sum(areas)
    cdf = c(0, cumsum(areas))

    list("edge"=edge, "Xs"=Xs, "Ys"=Ys, "cdf"=cdf)
}


#' @export
qgmp = function(p, ndf=NULL, pdim=NULL, svr=ndf/pdim, eigens=NULL,
                lower.tail=TRUE, log.p=FALSE, m=500){
    if(log.p){
        p = log(p)
    }
    if(!lower.tail){
        p = 1-p
    }
    if(is.null(eigens)){
        q = RMTstat::qmp(p, svr=svr)
    }else{
        general_mp = GeneralMarchenkoPasturPar(d=eigens, phi=1/svr, m=m)
        get_quantile = function(x){
            # Our computation is from left to right so x==1 is the special case.
            if(x==1){
                return(utils::tail(general_mp$Xs, 1))
            }

            # The quantile we want will fall into interval [Xs[idx],Xs[idx+1])
            idx = utils::tail(which(general_mp$cdf<=x),1)

            # Compute the approximate location between two quantiles
            s0 = x - general_mp$cdf[idx]
            # Simple quadratic equation using the root formula.
            a = general_mp$Ys[idx]
            b = general_mp$Ys[idx+1]
            c = general_mp$Xs[idx+1] - general_mp$Xs[idx]
            if(b==a){
                k = s0/(a*c)
            }else{
                k = (-a*c+sqrt(a^2*c^2+2*s0*c*(b-a)))/(b-a)
            }
            general_mp$Xs[idx] + c*k
        }
        q = unlist(lapply(p, get_quantile))
    }
    q
}


#' @export
rgmp = function(n, ndf=NULL, pdim=NULL, svr=ndf/pdim, eigens=NULL, m=500){
    if(is.null(eigens)){
        RMTstat::rmp(n, var=1, svr=svr)
    }else{
        u = stats::runif(n)
        qgmp(u, svr=1/svr, eigens=eigens, m=m)
    }
}


#' @export
pgmp = function(q, ndf=NULL, pdim=NULL, svr=ndf/pdim, eigens=NULL,
                lower.tail=TRUE, log.p=FALSE, m=500){
    if(is.null(eigens)){
        p = RMTstat::pmp(q, svr=svr)
    } else{
        general_mp = GeneralMarchenkoPasturPar(d=eigens, phi=1/svr, m=500)
        get_cdf = function(x){
            # Our computation is from left to right so x==1 is the special case.
            if(x>=utils::tail(general_mp$Xs,1)){
                return(1)
            }
            if(x<general_mp$Xs[1]){
                return(0)
            }
            # The cdf we want will fall into interval [cdf[idx],cdf[idx+1])
            idx = utils::tail(which(general_mp$Xs<=x),1)
            a = general_mp$Ys[idx]
            b = general_mp$Ys[idx+1]
            c = general_mp$Xs[idx+1] - general_mp$Xs[idx]
            k = (x-general_mp$Xs[idx])/c
            s0 = (2*a+(b-a)*k)*k*c/2
            general_mp$cdf[idx] + s0
        }
        p = unlist(lapply(q, get_cdf))
    }
    if(log.p){
        p = log(p)
    }
    if(!lower.tail){
        p = 1-p
    }
    p
}


#' @export
dgmp = function(x, ndf=NULL, pdim=NULL, svr=ndf/pdim, eigens=NULL,
                log.p=FALSE, m=500){
    if(is.null(eigens)){
        d = RMTstat::dmp(x, svr=svr, log=log.p)
    } else{
        general_mp = GeneralMarchenkoPasturPar(d=eigens, phi=1/svr, m=m)
        get_pdf = function(x){
            # Our computation should be restricted to the support
            if(x>=utils::tail(general_mp$Xs, 1)|x<=general_mp$Xs[1]){
                return(0)
            }
            # The pdf we want will fall into interval [Ys[idx],Ys[idx+1])
            idx = utils::tail(which(general_mp$Xs<=x),1)
            a = general_mp$Ys[idx]
            b = general_mp$Ys[idx+1]
            c = general_mp$Xs[idx+1] - general_mp$Xs[idx]
            k = (x-general_mp$Xs[idx])/c
            general_mp$Ys[idx] + k*(b-a)
        }
        d = unlist(lapply(x, get_pdf))
    }
    if(log.p){
        d = log(d)
    }
    d
}
