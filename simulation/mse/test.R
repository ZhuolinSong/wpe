#install.packages("synfd", repos = "https://cloud.r-project.org")
#install.packages("fdapace", repos = "https://cloud.r-project.org")
#install.packages("fda", repos = "https://cloud.r-project.org")
devtools::load_all()


aggf = function(u){
  mumseface = mean((u$muface - u$truemugrid)^2)
  covmseface = mean((u$Cface - u$truecovgrid)^2)

  mumsepace = mean((u$mupace - u$truemugrid)^2)
  covmsepace = mean((u$Cpace - u$truecovgrid)^2)

  mumsefour = mean((u$mufour - u$truemugrid)^2)
  covmsefour = mean((u$Cfour - u$truecovgrid)^2)

  mumsesnpt = mean((u$musnpt - u$truemugrid)^2)
  covmsesnpt = mean((u$Csnpt - u$truecovgrid)^2)

  covmsesamc = mean((u$Csamc - u$truecovgrid)^2)
  return(c(mumseface,covmseface,
           mumsepace,covmsepace,
           mumsefour,covmsefour,
           mumsesnpt,covmsesnpt,
           covmsesamc))
}

fittings = function(sim){
  grid <- regular.grid()
  fit1 = face::face.sparse(sim$data,argvals.new = grid)
  fit2 = covfunc(sim$t, sim$y, method='PACE', newt=grid, bw=sqrt(2)/2*(1-0.2), delta=1)
  fit3 = covfunc(sim$t, sim$y, method='FOURIER', newt=grid, ext=0.1, p=10, domain=c(0,1))
  fit4 = covfunc(sim$t, sim$y, method='SP', newt=grid, domain=c(0,1))

  X = t(to.wide(sim$data,sim$grid)); p=length(sim$grid)
  A=getA1_new_cv(X, p, dl=0.9*100*0.5, incre=0.1*100*0.5, bw = NA, kernel = 'epan', sigma2hat= NA)$A
  fit5 = tcrossprod(A)

  return(list(data=sim$data, truemugrid=sim$truemugrid, truecovgrid=sim$truecovgrid,
  muface = fit1$mu.new, Cface = fit1$Chat.new,
  mupace = fit2$mu$fitted, Cpace = fit2$fitted,
  mufour = fit3$mu$fitted, Cfour = fit3$fitted,
  musnpt = fit4$mu$fitted, Csnpt = fit4$fitted,
  Csamc = fit5))
}



set.seed(20200507)
n_sims <- 1
sim_dat <- list()
sims1 <- list()
i = 1
for (i in 1:n_sims){
  sim_dat[[i]] <- sim1(50, 0.2, delta=.2, m_avg = 4, sigx='fourier', mu=1)
  sims1[[i]] <- fittings(sim_dat[[i]])
}
sim = sim_dat[[1]]



sim_dat[[i]]$m
sim_dat[[i]]$grid


agg1 = lapply(sims1, aggf)
A = do.call(rbind, agg1)
colnames(A) = c("mumseface","covmseface",
                "mumsepace","covmsepace",
                "mumsefour","covmsefour",
                "mumsesnpt","covmsesnpt",
                "covmsezhang")
apply(A, 2, median)
apply(A, 2, IQR)

for(tmu in c(1,2)){
  for(tsigx in c('fourier', "matern", "bspline")){
    temp=sim1(10,0,0.5,mu=tmu,sigx=tsigx)
    print(mean(diag(temp$truecovgrid^2)/4))
  }
}

v_sige2 <- c(6.292991, 0.6890763, 0.001649281)
