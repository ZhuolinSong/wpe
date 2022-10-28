#' Covariance estimation by basis expansion
#' @param Lt a list (for irregular design) or a vector (for regular design)
#' @param Ly a list (for irregular design) or a matrix (for regular design). If \code{Ly} is a matrix, then \code{ncol(Ly)} must be equal to \code{length(Lt)}
#' @param delta the snippet parameter, only used for irregular design
#' @keywords internal
cov.be <- function(Lt,Ly,p=NULL,lam=NULL,ext=NULL,newt=NULL,mu=NULL,tuning='lle',weig=NULL,maxIt=3,domain=c(0,1), delta=NULL){
    if(is.null(p)) p <- seq(5,7)
    if(is.null(lam)) lam <- 0
    if(is.null(ext)) ext <- c(0,0.05,0.1)


    if(is.list(Ly)){
        n <- length(Ly)
        datatype <- 'irregular'
        if(is.null(weig)) weig <- get.weig.cov(Lt,Ly,'OBS')
    }else{
        n <- nrow(Ly)
        datatype <- 'regular'
        weig <- rep(1/n,n)
    }


    if(is.null(mu)) mu <- meanfunc(Lt,Ly,method='FOURIER', domain=domain, newt=newt)

    if(is.null(domain)){
        tmp <- unique(unlist(Lt))
        domain <- c(min(tmp),max(tmp))
    }

    if(domain[1] != 0 || domain[2] != 1){# standardize tobs
        if(datatype=='irregular')
            Lt <- lapply(Lt,function(x) (x-domain[1])/(domain[2]-domain[1]))
        else
            Lt <- (Lt-domain[1])/(domain[2]-domain[1])
    }

    R <- list(C=NULL,p=p,lam=lam,ext=ext,mu=mu,datatype=datatype,domain=domain,
              maxIt=maxIt,Lt=Lt,Ly=Ly,weig=weig,method='FOURIER')
    class(R) <- 'covfunc'


    
    if(is.list(Ly) && is.list(Lt)){ # irregular design
        delta <- estimate.delta(Lt)

        tmp_grid <- sort(unique(unlist(Lt)))
        xcov <- cov.pace(Lt, Ly, bw=NULL, newt=NULL, mu=mu, kernel='epanechnikov',delta=delta)

        Lr <- lapply(1:n,function(i){
            predict.covfunc(xcov,Lt[[i]])
        })

        R$Lr <- Lr
        R$delta <- delta

    }else if(is.vector(Lt) && is.matrix(Ly)){ # regular design
        mui <- predict(mu,Lt)
        Lr <- lapply(1:n,function(i){
            (Ly[i,]-mui) %*% t(Ly[i,]-mui)
        })

        R$Lr <- Lr
        R$delta <- 1
    }else stop('unsupported data type')

    if(length(p) > 1 || length(lam) > 1 || length(ext) > 1){
        if(tuning=='CV') stop('unsupported tuning method')
        else if(tuning=='GCV') stop('unsupported tuning method')
        else{
            R <- cov.be.lle(R,p,lam,ext)$R # select tuning parameters
        }
    }

    lam <- R$lam
    p <- R$p
    ext <- R$ext

    if(datatype=='irregular'){
        if(ext > 0) Ls <- shrink(Lt,ext)
        else Ls <- Lt
        auxmat <- comp.aux.mat(p,Ls)
        C <- estimate.cov.coef(Ls,Lr,auxmat,lam,weig,maxIt)
    }else{
        Lt <- list(Lt)
        if(ext > 0) Ls <- shrink(Lt,ext)
        else Ls <- Lt
        auxmat <- comp.aux.mat(p,Ls)

        Lr <- list(Reduce('+',Lr)/length(Lr))
        C <- estimate.cov.coef(Ls,Lr,auxmat,lam,weig,maxIt)
    }


    R$C <- C

    if(!is.null(newt)) R$fitted <- predict.covfunc(R,newt)

    return(R)

}

cov.be.lle <- function(R,p,lam,ext){

    if(R$datatype == 'irregular'){
        xcov <- cov.pace(R$Lt,R$Ly,bw=NULL,newt=NULL,mu=R$mu,
                         tuning='GCV',weig='SUBJ',kernel='epanechnikov',
                         delta=R$delta)
        newt <- xcov$newt
        tmp <- pracma::meshgrid(newt)
        x1 <- as.vector(tmp$X)
        x2 <- as.vector(tmp$Y)
        idx <- abs(x1-x2) < R$delta
        idx <- idx & (x1!=x2)
        y <- as.vector(xcov$fitted)[idx]

        Lr <- R$Lr
        Lt <- R$Lt
        Ly <- R$Ly
    }else{
        newt <- R$Lt
        n <- length(R$Lr)
        tr <- sample(n,ceiling(0.75*n))
        te <- setdiff(1:n,tr)

        tmp <- R$Lr[te]
        y <- Reduce('+',tmp) / length(te)
        y <- as.vector(y)
        idx <- rep(TRUE,length(y))

        Lr <- list(Reduce('+',R$Lr[tr])/length(tr))
        Lt <- list(R$Lt)

    }

    # find p and lambda first, then ext
    R$ext <- 0
    #system.time({ 
    if(length(p) > 1 || length(lam) > 1){
        err <- matrix(Inf,length(p),length(lam))
        for(u in 1:length(p)){
            auxmat <- comp.aux.mat(p[u],Lt)
            R$p <- p[u]
            for(v in 1:length(lam)){
                C <- estimate.cov.coef(Lt,Lr,auxmat,lam[v],R$weig,maxIt=R$maxIt, method='naive')
                R$C <- C
                R$lam <- lam[v]
                covhat <- predict.covfunc(R,newt)
                yhat <- as.vector(covhat)
                yhat <- yhat[idx]
                err[u,v] <- mean((y-yhat)^2) 
            }
        }
        min.idx <- which(err==min(err),arr.ind=T)
        p <- p[min.idx[1]]
        lam <- lam[min.idx[2]]
    }
    #})
    R$p <- p
    R$lam <- lam

    if(length(ext) > 1){
        err <- rep(0,length(ext))
        for(u in 1:length(ext)){
            Ls <- shrink(Lt,ext[u])
            auxmax <- comp.aux.mat(p,Ls)
            C <- estimate.cov.coef(Lt,Lr,auxmat,lam,R$weig,maxIt=R$maxIt)
            R$ext <- ext[u]
            covhat <- predict.covfunc(R,newt)
            yhat <- as.vector(covhat)
            yhat <- yhat[idx]
            err[u] <- mean((y-yhat)^2)
        }
        ext <- ext[which.min(err)]
    }
    else{
        R$ext <- ext
    }

    return(list(R=R,ext=ext,p=p,lam=lam))
}
