devtools::install_github("linulysses/mcfda")

library(mcfda)
library(synfd)
library(fields)
library(fdapace)

# mean function
mu <- function(s) 2*s^2 *cos(2*pi*s)

# cov function
sig2 = function(t){
  sqrt(sqrt(t) * exp((-(t-.1)^2)/10) + 1)
}
C1stf = function(s,t){
  sig2(s) * matern(s,t,nu=0.5,rho=1)* sig2(t)
}
C1f = function(ts){
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

# generate data with delta = 0.25
set.seed(3)
grid <- regular.grid()
cov0 <- C1f(grid)
D <- synfd::irreg.fd(mu=mu, X=synfd::gaussian.process(C1f), n=50, m=4, delta = 0.25, sig = 0.25, snr = NULL)

# fit with 'PACE' without delta = 0.25 (error)
error <- try(cov.obj <- covfunc(D$t,D$y,newt=NULL,mu=NULL,method='PACE', tuning='GCV',weig=NULL,kernel='epanechnikov',delta=0), silent = T)
print(error)

# fit with 'PACE' with delta = 0.25
cov.obj <- covfunc(D$t,D$y,newt=NULL,mu=NULL,method='PACE', tuning='GCV',weig=NULL,kernel='epanechnikov',delta=0.25)
cov.hat <- predict(cov.obj,grid)
mean((cov.hat-cov0)^2) / mean(cov0^2)

# fit with 'FOURIER basis' (no error)
cov.obj <- covfunc(D$t,D$y,newt=NULL,mu=NULL,method='FOURIER', tuning='lle',weig=NULL)
cov.hat <- predict(cov.obj,grid)
mean((cov.hat-cov0)^2) / mean(cov0^2)
