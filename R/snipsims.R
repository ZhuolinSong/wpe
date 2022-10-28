# Required libraries
#library("synfd")
#library("fields")
#library("fdapace")
#library("fda")
#library(splines2)

# Helper functions


# 'fourier'
C.fourier.coef2 = function(l=5) {
  temp = matrix(NA, nrow=l, ncol=l)
  for(i in seq(l)){
    for(j in seq(l)){
      if (i==j) temp[i,j] = 1.5^(1-j)
      if (i!=j) temp[i,j] = 2^(-abs(i-j)-5/2)
    }
  }

  temp
}
C.fourier.offd = function(j,k)  exp(-abs(j-k))/5
C.fourier.diag = function(j) 1.5^(1-j)
C.fourier.coef = function(l=5) {
  temp = matrix(NA, nrow=l, ncol=l)
  for(i in seq(l)){
    for(j in seq(l)){
      temp[i,j] = C.fourier.offd(i, j)
    }
  }
  #diag(temp) = C.fourier.diag(seq(l))

  temp
}
C.fourier = function(ts, l=5) {
  if(l%%2 == 0) l = l+1
  f.coef = C.fourier.coef(l)
  #f.basis = create.fourier.basis(rangeval = c(0,1),nbasis = l)
  #mat = eval.basis(ts,basisobj=f.basis)
  mat = evaluate.basis(l,grid=ts,type = 'FOURIER') 
  if(length(ts)<2){
    return(crossprod(mat, f.coef)%*%mat)
  }
  mat %*% tcrossprod(f.coef, mat)

}


# 'matern'
C.matern.sig = function(t) sqrt(sqrt(t) * exp(-(t-0.1)^2/10) + 1)
C.matern.stf = function(s,t) C.matern.sig(s) * matern(s,t,nu=0.5,rho=1)* C.matern.sig(t)
C.matern = function(ts){
  sig = C.matern.sig(ts)
  mat = matern(ts,nu=1,rho=1)
  mat = tcrossprod(sig) * mat
  return(mat)
}



# nonperiodic and nonsmooth covariance in Lin
f_rho = function(s,t) exp(-(s-t)^2)

f_v = function(t){
  sqrt((1 + t + max(t-2/9,0) + max(t-4/9,0) + max(t-6/9, 0) + max(t-8/9, 0))/sqrt(2))
}

C.2.coef = function(ts) {
  n=length(ts)
  diff = ts[2]-ts[1]
  temp = matrix(1, nrow=n, ncol=n)
  for (i in seq(n-1)) temp[abs(row(temp) - col(temp)) == i] <- exp(-(i*diff)^2)
  temp
}

C.2= function(ts) {
  tcrossprod(sapply(ts, f_v)) * C.2.coef(ts)
}


# 'general'
C3stf = function(j,k,l) exp(-abs(j-k))
C3.diag = function(j,l) exp((1-j))
C3.coef = function(l=10) {
  temp = matrix(NA, nrow=l, ncol=l)
  for(i in seq(l)){
    for(j in seq(l)){
      temp[i,j] = C3stf(i, j, l)
    }
  }
  diag(temp) = C3.diag(seq(l), 1)
  r = eigen(temp)
  V = r$vectors
  lam = r$values
  lam[which(lam < 0)] = 0
  V %*% diag(lam) %*% solve(V)
}
C3 = function(ts, l=10) {
  if(l%%2 != 0) l = l+1
  f.coef = C3.coef(l)
  f.basis = create.fourier.basis(rangeval = c(0,1),nbasis = (l+1))
  mat = eval.basis(ts,basisobj=f.basis)[,-1]
  if(length(ts)<2){
    return(crossprod(mat, f.coef)%*%mat)
  }
  return(mat %*% tcrossprod(f.coef, mat))
}




# mean function mu0 
mu0f = function(t) 0*t



sim.truth = function(params){
  grid = params[['grid']]
  mu = params[['mu']]
  sigx = params[['sigx']]

  truecovgrid = sigx(grid)
  truemugrid = mu(grid)

  return(list(truemugrid=truemugrid, truecovgrid=truecovgrid))

}


# Sim 1 and Sim 2 (regular design)
sim.regular <- function(params){
  n = params[['n']]
  sige2 = params[['sige2']]
  delta = params[['delta']]
  m_avg = params[['m_avg']]
  grid = params[['grid']]
  mu = params[['mu']]
  sigx = params[['sigx']]

  # fix the grid
  t_num = length(grid)
  span = ceiling(delta*t_num)
  width = ceiling(delta/2*t_num)

  m = rpois(n,m_avg)
  m[m>span] = span
  m[m<2] = 2
  if (length(width:(t_num-width))==1) Oi = rep(width, n)
  if (length(width:(t_num-width))>1) Oi = sample(width:(t_num-width), n, replace = T)
  t = lapply(1:n, function(i) sort(grid[sample((Oi[i] - width + 1) : (Oi[i] + width), m[i])]))

  truecovgrid = as.vector(sigx(grid))
  extract_sigx = function(t){
    tmp <- pracma::meshgrid(which(grid %in% t))
    x1 <- as.vector(tmp$X)
    x2 <- as.vector(tmp$Y)
    idx = (x1-1)*t_num + x2
    matrix(truecovgrid[idx], length(t), length(t))
  }

  Xs = lapply(1:n,function(i) extract_sigx(t[[i]]))
  mu = lapply(1:n,function(i) mu(t[[i]]))

  X = lapply(1:n, function(i) c(mvtnorm::rmvnorm(n=1, sigma=Xs[[i]])))
  eps = lapply(1:n, function(i) rnorm(m[i], 0, sqrt(sige2)))
  Y = lapply(1:n, function(i) mu[[i]] + X[[i]] + eps[[i]])
  ts = unlist(t)
  ys = unlist(Y)
  id = c()
  for(i in seq_along(t)){id = c(id, rep(i, length(t[[i]])))}
  data=data.frame(argvals=ts, subj=id, y=ys)
  cache = list(m=m, grid=grid, O=Oi)
  return(list(data=data,t=t, y=Y, cache=cache))
}

# Sim 3 (random design)
sim.random <- function(params){
  n = params[['n']]
  sige2 = params[['sige2']]
  delta = params[['delta']]
  m_avg = params[['m_avg']]
  grid = params[['grid']]
  mu = params[['mu']]
  sigx = params[['sigx']]

  # Generate # obs and time points
  m = rpois(n,m_avg)
  m[m<2] = 2
  Oi = runif(n,delta/2, 1-delta/2)
  t = lapply(1:n, function(i) sort(runif(m[i], Oi[i] - delta/2, Oi[i] + delta/2)))

  # Generate Mean and Cov function
  mu = lapply(1:n,function(i) mu(t[[i]]))
  Xs = lapply(1:n,function(i) sigx(t[[i]]))

  X = lapply(1:n, function(i) c(mvtnorm::rmvnorm(n=1, sigma=Xs[[i]])))
  eps = lapply(1:n, function(i) rnorm(m[i], 0, sqrt(sige2)))
  Y = lapply(1:n, function(i) mu[[i]] + X[[i]] + eps[[i]])
  ts = unlist(t)
  ys = unlist(Y)
  id = c()
  for(i in seq_along(t)){id = c(id, rep(i, length(t[[i]])))}
  data=data.frame(argvals=ts, subj=id, y=ys)
  cache = list(m=m, grid=grid, O=Oi)
  return(list(data=data,t=t, y=Y, cache=cache))
}




# Sim full (regular design with full curves)
sim.regular.full <- function(params){
  n = params[['n']]
  sige2 = params[['sige2']]
  delta = params[['delta']]
  m_avg = params[['m_avg']]
  grid = params[['grid']]
  mu = params[['mu']]
  sigx = params[['sigx']]

  # fix the grid
  t_num = length(grid)
  span = ceiling(delta*t_num)
  width = ceiling(delta/2*t_num)

  # generate the full functional curves
  # full.sims = sim.regular(list(delta=1, m_avg=2*t_num, n = n, sige2 = sige2, grid=grid, mu=mu, sigx=sigx))
  Xs = sigx(grid)
  mu = mu(grid)
  X = lapply(1:n, function(i) c(mvtnorm::rmvnorm(n=1, sigma=Xs)))
  eps = lapply(1:n, function(i) rnorm(t_num, 0, sqrt(sige2)))
  full.Y = lapply(1:n, function(i) mu + X[[i]] + eps[[i]])

  # generate the random curves
  m = rpois(n,m_avg)
  m[m>span] = span
  m[m<2] = 2
  if (length(width:(t_num-width))==1) Oi = rep(width, n)
  if (length(width:(t_num-width))>1) Oi = sample(width:(t_num-width), n, replace = T)
  t = lapply(1:n, function(i) sort(grid[sample((Oi[i] - width + 1) : (Oi[i] + width), m[i])]))
  
  #selected curves
  Y = lapply(1:n, function(i) full.Y[[i]][which(grid%in%t[[i]])] )

  id = c()
  for(i in seq_along(t)){id = c(id, rep(i, length(t[[i]])))}
  data=data.frame(argvals=unlist(t), subj=id, y=unlist(Y))
  data.full = data.frame(argvals=rep(grid, times=n) , subj=rep(seq(n), each=t_num), y=unlist(full.Y))
  cache = list(m=m, grid=grid, O=Oi)
  return(list(data=data, data.full = data.full, t=t, y=Y, cache = cache))
}

# Sim full (random design with full curves)
sim.random.full <- function(params){
  n = params[['n']]
  sige2 = params[['sige2']]
  delta = params[['delta']]
  m_avg = params[['m_avg']]
  grid = params[['grid']]
  mu = params[['mu']]
  sigx = params[['sigx']]

  # Generate # obs and time points
  m = rpois(n,m_avg)
  m[m<2] = 2
  Oi = runif(n,delta/2, 1-delta/2)
  t = lapply(1:n, function(i) sort(runif(m[i], Oi[i] - delta/2, Oi[i] + delta/2)))
  grid = sort(unique(unlist(t)))

  # fix the grid
  t_num = length(grid)

  # generate the full functional curves
  # full.sims = sim.random(list(delta=1, m_avg=2*t_num, n = n, sige2 = sige2, grid=grid, mu=mu, sigx=sigx))
  Xs = sigx(grid)
  mu = mu(grid)
  X = lapply(1:n, function(i) c(mvtnorm::rmvnorm(n=1, sigma=Xs)))
  eps = lapply(1:n, function(i) rnorm(t_num, 0, sqrt(sige2)))
  full.Y = lapply(1:n, function(i) mu + X[[i]] + eps[[i]])

  # Generate random curves
  # sims = sim.random(params)
  # data.full = data.frame(argvals=rep(grid, times=n) , subj=rep(seq(n), each=t_num), y=unlist(full.Y))
  
  #selected curves
  Y = lapply(1:n, function(i) full.Y[[i]][which(grid %in% t[[i]])] )
  id = c()
  for(i in seq_along(t)){id = c(id, rep(i, length(t[[i]])))}
  data=data.frame(argvals=unlist(t), subj=id, y=unlist(Y))
  data.full = data.frame(argvals=rep(grid, times=n) , subj=rep(seq(n), each=t_num), y=unlist(full.Y))
  cache = list(m=m, grid=grid, O=Oi)
  sims = list(data=data, y=Y, t=t, cache=cache)
  
  return(list(data=sims$data, data.full = data.full, t=sims$t, y=sims$Y, cache = sims$cache))
}



# Sim 2 (sparse regular grid)
sim2 <- function(n,sige2,delta=.2, m_avg = 4, t_num = 100,
mu=c(1,2), sigx=c('fourier', "matern", "bspline", "general")){
  grid <- regular.grid()
  # fix the grid number as t_num
  t_grid <- seq(0.005, 1-0.005,length.out=t_num)
  diff <- t_grid[2] - t_grid[1]
  span <- delta / diff

  m = rpois(n,m_avg) + 1
  m[m>span] = span
  m[m<2] = 2
  Oi <- sample(ceiling(span/2):(t_num-ceiling(span/2)), n, replace = T)
  t = lapply(1:n, function(i) sort(t_grid[sample((Oi[i] - ceiling(span/2) + 1) : (Oi[i] + ceiling(span/2)), m[i])]))

  if (sigx=='fourier'){
    sigx = lapply(1:n,function(i) C.fourier(t[[i]]))
    truecovgrid = C.fourier(grid)
  }else if (sigx=='matern'){
    sigx = lapply(1:n,function(i) C.matern(t[[i]]))
    truecovgrid = C.matern(grid)
  }else if (sigx=='bspline'){
    sigx = lapply(1:n,function(i) C.bspline(t[[i]]))
    truecovgrid = C.bspline(grid)
  }else if (sigx=='general'){
    sigx = lapply(1:n,function(i) C3(t[[i]]))
    truecovgrid = C3(grid)
  }

  if (mu==1){
    truemugrid = mu1f(grid)
    mu = lapply(1:n,function(i) mu1f(t[[i]]))
  } else if (mu==2){
    truemugrid = mu2f(grid)
    mu = lapply(1:n,function(i) mu2f(t[[i]]))
  }


  X = lapply(1:n, function(i) c(mvtnorm::rmvnorm(n=1, sigma=sigx[[i]])))
  eps = lapply(1:n, function(i) rnorm(m[i], 0, sqrt(sige2)))
  Y = lapply(1:n, function(i) mu[[i]] + X[[i]] + eps[[i]])
  ts = unlist(t)
  ys = unlist(Y)
  Zs = Cs = id = c()
  for(i in 1:length(t)){id = c(id, rep(i, length(t[[i]])))}
  data=data.frame(argvals=ts, subj=id, y=ys)
  return(list(data=data, m = m, grid=t_grid, t=t, y=Y, truemugrid=truemugrid, truecovgrid=truecovgrid))
}





# extra helpers

# Identity
C.identity = function(ts) diag(length(ts))

# FULL
C.ones = function(ts) matrix(1, nrow = length(ts), ncol=length(ts))

# full 1 coef
C1.comp = function(l){
  
  # test inputs
  temp = matrix(pi/2, l, l)
  diag(temp) = pi/2
  for (i in seq(l-1)) temp[abs(row(temp) - col(temp)) == i] <- temp[abs(row(temp) - col(temp)) == i]*1.1^(-i)
  
  temp
}

# 'b-spline'

C.bspline.offd = function(j,k,l) exp(-abs(j-k)/l)
C.bspline.diag = function(j, l) pi^(-j^2/l+1)
C.bspline.coef = function(l=10) C1.comp(l)

C.bspline = function(ts, l=10) {
  f.coef = C.bspline.coef(l)
  knots <- seq(0,1,length=l-2)
  knots <- knots[-c(1, length(knots))]
  mat = bSpline(ts, knots = knots, degree = 3, intercept = T, Boundary.knots = c(0,1))
  return(mat %*% tcrossprod(f.coef, mat))
}

# nonperiodic 1 in Delaigle
# phi1 = function(t) rep(1,length(t))
# phi2 = function(t) (2*t - 1)*sqrt(3)
# phi3 = function(t) (6*t^2 - 6*t + 1)*sqrt(5)
# phi4 = function(t) (20*t^3 - 30*t^2 + 12*t - 1)*sqrt(7)
# 
# C.2= function(ts) {
#   tcrossprod(phi1(ts)) + 0.5*tcrossprod(phi2(ts)) + 0.5^2*tcrossprod(phi3(ts)) + 0.5^3*tcrossprod(phi4(ts))
# }

# mean function mu1 
mu1f = function(t) 2 * t^2 * cos(2*pi*t)

# mean function mu2
mu2f = function(t) exp(t)/2
