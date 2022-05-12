# Required libraries
#library("synfd")
#library("fields")
#library("fdapace")
#library("fda")
#library(splines2)

# Helper functions
C.fourier.offd = function(j,k) {
  2^(-abs(j-k)-5/2)
}
C.fourier.diag = function(j) {
  1.5^(1-j)
}
C.fourier.coef = function(l=5) {
  temp = matrix(NA, nrow=l, ncol=l)
  for(i in seq(l)){
    for(j in seq(l)){
      temp[i,j] = C.fourier.offd(i, j)
    }
  }
  diag(temp) = C.fourier.diag(seq(l))
  (temp+t(temp))/2
}
C.fourier = function(ts, l=50) {
  if(l%%2 != 0) l = l+1
  f.coef = C.fourier.coef(l)
  f.basis = create.fourier.basis(rangeval = c(0,1),nbasis = (l+1))
  mat = eval.basis(ts,basisobj=f.basis)[,-1]
  if(length(ts)<2){
    return(crossprod(mat, f.coef)%*%mat)
  }
  out = mat %*% tcrossprod(f.coef, mat)
  return((out+t(out))/2)
}

C.matern.var = function(t){
  sqrt(sqrt(t) * exp((-(t-.1)^2)/10) + 1)
}
C1stf = function(s,t){
  C.matern.var(s) * matern(s,t,nu=0.5,rho=1)* C.matern.var(t)
}
C.matern = function(ts){
  mat = matrix(NA,length(ts),length(ts))
  for(i in 1:length(ts)){
    for(j in 1:length(ts)){
      s = ts[i]
      t = ts[j]
      mat[i,j] = C1stf(s, t)
    }
  }
  return(mat)
}

C.bspline.offd = function(j,k){
  2^(-abs(j-k)-5/2)
}
C.bspline.diag = function(j) {
  1.5^(1-j)
}
C.bspline.coef = function(l=50) {
  temp = matrix(NA, nrow=l, ncol=l)
  for(i in seq(l)){
    for(j in seq(l)){
      temp[i,j] = C.bspline.offd(i, j)
    }
  }
  diag(temp) = C.bspline.diag(seq(l))
  temp
}
C.bspline = function(ts, l=50) {
  f.coef = C.bspline.coef(l)
  knots <- seq(0,1,length=l-2)
  knots <- knots[-c(1, length(knots))]
  mat = bSpline(ts, knots = knots, degree = 3, intercept = T, Boundary.knots = c(0,1))
  return(mat %*% tcrossprod(f.coef, mat))
}


C3stf = function(j,k,l) exp(-abs(j-k)/l)
C3.coef = function(l=10) {
  temp = matrix(NA, nrow=l, ncol=l)
  for(i in seq(l)){
    for(j in seq(l)){
      temp[i,j] = C3stf(i, j, l)
    }
  }
  diag(temp) = C.fourier.diag(seq(l))
  temp
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

# mean function mu1 
mu1f = function(t){
  2 * t^2 * cos(2*pi*t)
}
# mean function mu2
mu2f = function(t){
  exp(t)/2
}




# Sim 1 (dense regular grid)
sim1 <- function(n,sige2,delta=.2, m_avg = 50, t_num=100,
mu=c(1,2), sigx=c('fourier', "matern", "bspline", "general")){
  grid <- regular.grid()
  # fix the grid number as t_num
  t_grid <- seq(0.005, 1-0.005,length.out=t_num)
  diff <- t_grid[2] - t_grid[1]
  span <- delta / diff

  m = rep(m_avg,n)
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

# Sim 3 (random design)
sim3 <- function(n,sige2,delta=.2, m_avg = 4,
mu=c(1,2), sigx=c('fourier', "matern", "bspline", "general")){
  grid <- regular.grid()
  
  m = rpois(n,m_avg) + 1
  m[m<2] = 2
  Oi = runif(n,delta/2, 1-delta/2)
  t = lapply(1:n, function(i) sort(runif(m[i], Oi[i] - delta/2, Oi[i] + delta/2)))

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
  #ord = order(ts)
  #mus = unlist(mu)
  ys = unlist(Y)
  Zs = Cs = id = c()
  for(i in 1:length(t)){id = c(id, rep(i, length(t[[i]])))}
  data=data.frame(argvals=ts, subj=id, y=ys)
  return(list(data=data, m = m, grid=grid, t=t, y=Y, truemugrid=truemugrid, truecovgrid=truecovgrid))
}

# Sim full (fixed grid with full curves)
sim.full <- function(n,sige2,delta=.2, m_avg = 4, t_num = 100,
mu=c(1,2), sigx=c("matern", 'fourier', "bspline", "general")){
  grid <- regular.grid()

  m = rpois(n,m_avg) + 1
  # fix the grid number as t_num
  t_grid <- seq(0.005, 1-0.005,length.out=t_num)
  diff <- t_grid[2] - t_grid[1]
  span <- delta / diff
  m[m>span] = span
  Oi <- sample(ceiling(span/2):(t_num-ceiling(span/2)), n, replace = T)
  t = lapply(1:n, function(i) sort(t_grid[sample((Oi[i] - ceiling(span/2) + 1) : (Oi[i] + ceiling(span/2)), m[i])]))

  #full curves
  if (sigx=='fourier'){
    sigx = C.fourier(t_grid)
    truecovgrid = C.fourier(grid)
  }else if (sigx=='matern'){
    sigx = C.matern(t_grid)
    truecovgrid = C.matern(grid)
  }else if (sigx=='bspline'){
    sigx = C.bspline(t_grid)
    truecovgrid = C.bspline(grid)
  }else if (sigx=='general'){
    sigx = C3(t_grid)
    truecovgrid = C3(grid)
  }

  if (mu==1){
    truemugrid = mu1f(grid)
    mu = mu1f(t_grid)
  } else if (mu==2){
    truemugrid = mu2f(grid)
    mu = mu2f(t_grid)
  }

  X = lapply(1:n, function(i) c(mvtnorm::rmvnorm(n=1, sigma=sigx)))
  eps = lapply(1:n, function(i) rnorm(t_num, 0, sqrt(sige2)))
  full.Y = lapply(1:n, function(i) mu + X[[i]] + eps[[i]])

  #selected curves
  Y = lapply(1:n, function(i) full.Y[[i]][which(t_grid%in%t[[i]])] )
  ts = unlist(t)

  ys = unlist(Y)
  id = c()
  for(i in 1:length(t)){id = c(id, rep(i, length(t[[i]])))}
  data=data.frame(argvals=ts, subj=id, y=ys)
  data.full = data.frame(argvals=rep(t_grid, n), subj=rep(seq(n),each=t_num), y=unlist(full.Y))
  return(list(data=data, m = m, grid=t_grid, t=t, y=Y, data.full = data.full, truemugrid=truemugrid, truecovgrid=truecovgrid))
}





