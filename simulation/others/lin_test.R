devtools::load_all()
library(synfd)
library(fields)
library(fdapace)
library(ggplot2)

# mu <- function(s) 2*s^2 *cos(2*pi*s)
# cov <- synfd::matern
# D <- synfd::irreg.fd(mu=mu, X=synfd::gaussian.process(C1f), n=50, m=4, delta = 0.25, sig = 0.25)
# cov0 <- cov(grid)

C1stf = function(s,t){
  sig2(s) * matern(s,t,nu=0.5,rho=1)* sig2(t)
}
sig2 = function(t){
  sqrt(sqrt(t) * exp((-(t-.1)^2)/10) + 1)
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

mu1f = function(t){
  2 * t^2 * cos(2*pi*t)
}

# cov
grid <- regular.grid()
sim <- function(n_obj, sige, del=.25, m_avg = 4, snr = NULL){
  D <- synfd::irreg.fd(mu=mu1f, X=synfd::gaussian.process(C1f), n=n_obj, m=m_avg, sig = sqrt(sige), delta = del, snr = snr)
  ts = unlist(D$t)
  ys = unlist(D$y)
  truemugrid = mu1f(grid)
  truecovgrid = C1f(grid)
  Zs = Cs = id = c()
  for(i in seq_len(n_obj)){id = c(id, rep(i, length(D$t[[i]])))}
  data=data.frame(argvals=ts, subj=id, y=ys)
  return(list(data=data, t = D$t, y = D$y, D=D, truemugrid=truemugrid, truecovgrid=truecovgrid, truesig2 = attr(D,'sig')^2 ))
}
cov_aggf = function(u){
  mumseface = mean((u$muface - u$truemugrid)^2)
  mumsesnip = mean((u$musnip4 - u$truemugrid)^2)
  covmseface = mean((u$Cface - u$truecovgrid)^2)
  #covmsesnip2 = mean((u$Csnip2 - u$truecovgrid)^2)
  covmsesnip3 = mean((u$Csnip3 - u$truecovgrid)^2)
  covmsesnip4 = mean((u$Csnip4 - u$truecovgrid)^2)
  return(c(mumseface,mumsesnip,covmseface,covmsesnip3,covmsesnip4))
}
cov_fit = function(n_obj, sige, del=.25, m_avg = 4){
  ss <- sim(n_obj, sige, del, m_avg)
  fit1 = face::face.sparse(ss$data,argvals.new = grid)
  #fit2 = covfunc(ss$t,ss$y,newt=NULL,mu=NULL,method='PACE', tuning='GCV',weig=NULL,kernel='epanechnikov',delta=del)
  fit3 = covfunc(ss$t,ss$y,newt=NULL,mu=NULL,method='FOURIER', tuning='lle',weig=NULL, domain=c(0,1))
  fit4 = covfunc(ss$t,ss$y,,newt=NULL,method='SP',weig=NULL, domain=c(0,1))
  #cov.hat2 <- predict(fit2, grid)
  cov.hat3 <- predict(fit3, grid)
  mu.hat4 = predict(fit4$mu, grid)
  cov.hat4 <- predict(fit4, grid)
  return(list(sim=ss,
            Cface = fit1$Chat.new, muface = fit1$mu.new, #Csnip2 = cov.hat2,
            Csnip3 = cov.hat3,musnip4 = mu.hat4, Csnip4 = cov.hat4,
            truemugrid=ss$truemugrid, truecovgrid=ss$truecovgrid))
}

set.seed(20220323)
n_times = 200
cov_sims <- lapply(1:n_times, function(i){cov_fit(50, 0.25, del=.25, m_avg = 4)})
agg1 = lapply(cov_sims, cov_aggf)
A = do.call(rbind, agg1)
colnames(A) = c("mumseface","mumsesnip","covmseface","covmsebasis","covmsesp")
round(apply(A, 2, median), 3)
round(apply(A, 2, IQR), 3)
round(apply(A, 2, mean), 3)
round(apply(A, 2, sd), 3)

# infinity covariance investigate: reason is optim('bfgs') sometimes output Inf (both cov.sp and cov.basis)
which(A[,5] == Inf) #169
A[169,]

ss_inf <- cov_sims[[169]]$sim
dat_inf <- ss_inf$data
dat_inf <- within(dat_inf, {  subj <- factor(subj)})

fit3 = covfunc(ss_inf$t,ss_inf$y,newt=NULL,mu=NULL,method='FOURIER', tuning='lle',weig=NULL, domain=c(0,1))
predict(fit3, grid)[1,]
cov_sims[[169]]$Csnip3[1,]

inf.fit = covfunc(ss_inf$t,ss_inf$y,,newt=NULL,method='SP',weig=NULL, domain=c(0,1))
cov.inf <- predict(inf.fit, grid)
cov_sims[[169]]$Csnip4[1,]



# wide simulation
library(R.matlab)
set.seed(20220323)
long.data <- sim(50, 0.25, del=.25, m_avg = 4)$data
long.data <- long.data[order(long.data$argvals),]
wide.data <- as.matrix(reshape(long.data, v.names="y",idvar="subj",
                     timevar="argvals",direction="wide"))[,-1]

filename <- paste(tempfile(), ".mat", sep = "")
writeMat(con=filename, x=wide.data)

filename <- paste(tempfile(), ".mat", sep = "")
writeMat(con=filename, t=sort(long.data$argvals))

filename <- paste(tempfile(), ".mat", sep = "")
writeMat(con=filename, grid=grid)
# mean
mu_aggf = function(u){
  mumseface = mean((u$muface - u$truemugrid)^2)
  mumsebasis = mean((u$mu_basis - u$truemugrid)^2)
  mumsepace = mean((u$mu_pace - u$truemugrid)^2)
  sigmseface = mean((u$sig_face - u$truesig2)^2)
  sigmsemcfda = mean((u$sig_mcfda - u$truesig2)^2)
  return(c(mumseface,mumsebasis,mumsepace,sigmseface,sigmsemcfda))
}
mu_fit = function(n_obj, sige, del=.25, m_avg = 4){
  ss <- sim(n_obj, sige, del, m_avg)
  fit1 = face::face.sparse(ss$data,argvals.new = grid)
  mu.basis = meanfunc(ss$t,ss$y,newt=grid,method='FOURIER', tuning='cv',weig=NULL, domain=c(0,1))
  mu.pace = meanfunc(ss$t,ss$y,newt=grid,method='PACE', tuning='cv',weig=NULL,kernel='epanechnikov',deg=1)

  sig.face <- fit1$sigma2
  sig.mcfda <- mcfda::sigma2(ss$t,ss$y)

  return(list(sim=ss,muface = fit1$mu.new, mu_basis = mu.basis$fitted,
   mu_pace = mu.pace$fitted, sig_face = sig.face, sig_mcfda = sig.mcfda,
    truemugrid=ss$truemugrid,  truesig2 = ss$truesig2))
}


set.seed(2022317)
n_times = 200
mu_sims <- lapply(1:n_times, function(i){mu_fit(n_obj=50, sige=.25, del=.25, m_avg = 4)})

agg2 = lapply(mu_sims, mu_aggf)
B = do.call(rbind, agg2)
colnames(B) = c("mumseface","mumsebasis","mumsepace","sigmseface","sigmsemcfda")
round(apply(B, 2, median), 3)
round(apply(B, 2, IQR), 3)
round(apply(B, 2, mean), 3)
round(apply(B, 2, sd), 3)







# mean
## mu.obj <- meanfunc(D$t,D$y,newt=NULL,method='FOURIER', tuning='cv',weig=NULL, domain=c(0,1))
## plot(mu.obj)
## lines(regular.grid(),mu(regular.grid()))
## 
## mu.obj <- meanfunc(D$t,D$y,newt=NULL,method='PACE', tuning='cv',weig=NULL,kernel='epanechnikov',deg=1)
## plot(mu.obj)
## lines(regular.grid(),mu(regular.grid()))

# cov

## pace
sim <- sim(50, 0.25, del=.25, m_avg = 4)
cov.obj <- covfunc(sim$t,sim$y,newt=NULL,mu=NULL,method='PACE', tuning='GCV',weig=NULL,kernel='epanechnikov',delta=0.25)
cov.hat <- predict(cov.obj,grid)
mean((cov.hat-sim$truecovgrid)^2) / mean(sim$truecovgrid^2)

## basis
cov.obj <- covfunc(sim$t,sim$y,newt=NULL,mu=NULL,method='FOURIER', tuning='lle',weig=NULL)
cov.hat <- predict(cov.obj,grid)
mean((cov.hat-sim$truecovgrid)^2) / mean(sim$truecovgrid^2)

## SP
cov.obj <- covfunc(sim$t,sim$y,,newt=NULL,method='SP',weig=NULL)
cov.hat <- predict(cov.obj,grid)
mean((cov.hat-sim$truecovgrid)^2) / mean(sim$truecovgrid^2)

# Error var
sig2 <- mcfda::sigma2(sim$t,sim$y)
sig2
abs(sig2 - attr(sim$D,'sig')^2)
