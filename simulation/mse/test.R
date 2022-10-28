#install.packages("synfd", repos = "https://cloud.r-project.org")
#install.packages("fdapace", repos = "https://cloud.r-project.org")
#install.packages("fda", repos = "https://cloud.r-project.org")
#install.packages("face", repos = "https://cloud.r-project.org")
#devtools::install_github("linulysses/synfd") 

devtools::load_all()
library(plotly)

set.seed(13579)


regular.illustrate = function(params){
grid = params[['grid']]
example = sim.regular(params)
truth = sim.truth(params)

start_time <- proc.time()
fit1 = covfunc(example$t, example$y, method='BE', lam=0, domain=c(0,1), mu=mu0f, newt=grid)
print(proc.time() - start_time)

start_time <- proc.time()
fit2 = covfunc(example$t, example$y, method='FOURIER', domain=c(0,1), mu=mu0f, newt=grid)
print(proc.time() - start_time)

start_time <- proc.time()
fit3 = covfunc(example$t, example$y, method='SP', domain=c(0,1), mu=mu0f, newt=grid)
print(proc.time() - start_time)

start_time <- proc.time()
fit4 = covfunc(example$t, example$y, method='PACE', kernel='gauss', mu=mu0f, newt=grid)
print(proc.time() - start_time)

start_time <- proc.time()
fit5 = face::face.sparse(example$data, argvals.new = grid, center = FALSE)
print(proc.time() - start_time)

start_time <- proc.time()
fit6 = tcrossprod(getA1_new_eig(Lt=example$t, Ly=example$y, r=10, mu=mu0f, newt=grid))
print(proc.time() - start_time)

print(c(mean((fit1$fitted - truth$truecovgrid)^2),
mean((fit2$fitted - truth$truecovgrid)^2),
mean((fit3$fitted - truth$truecovgrid)^2),
mean((fit4$fitted - truth$truecovgrid)^2),
mean((fit5$Chat.new - truth$truecovgrid)^2),
mean(((fit6) - truth$truecovgrid)^2)))

# print(c(
# norm(fit1$fitted - truth$truecovgrid, "F")^2,
# norm(fit2$fitted - truth$truecovgrid, "F")^2,
# norm(fit3$fitted - truth$truecovgrid, "F")^2,
# norm(fit4$fitted - truth$truecovgrid, "F")^2,
# norm(fit5$Chat.new - truth$truecovgrid, "F")^2,
# norm(fit6 - truth$truecovgrid, "F")^2))

print(mean(truth$truecovgrid^2))

print(plot_ly(z = ~truth$truecovgrid, x = ~grid, y = ~grid, type='surface'))
print(plot_ly(z = ~fit1$fitted, x = ~grid, y = ~grid, type='surface'))
print(plot_ly(z = ~fit2$fitted, x = ~grid, y = ~grid, type='surface'))
print(plot_ly(z = ~fit3$fitted, x = ~grid, y = ~grid, type='surface'))
print(plot_ly(z = ~fit4$fitted, x = ~grid, y = ~grid, type='surface'))
print(plot_ly(z = ~fit5$Chat.new, x = ~grid, y = ~grid, type='surface'))
print(plot_ly(z = ~fit6, x = ~grid, y = ~grid, type='surface'))
}
random.illustrate = function(params){
example = sim.random(params)
truth = sim.truth(params)
grid = params[['grid']]


start_time <- proc.time()
fit1 = covfunc(example$t, example$y, method='BE', lam=0, domain=c(0,1), mu=mu0f, newt=grid)
print(proc.time() - start_time)

start_time <- proc.time()
fit2 = covfunc(example$t, example$y, method='FOURIER', domain=c(0,1), mu=mu0f, newt=grid)
print(proc.time() - start_time)

start_time <- proc.time()
fit3 = covfunc(example$t, example$y, method='SP', domain=c(0,1), mu=mu0f, newt=grid)
print(proc.time() - start_time)

start_time <- proc.time()
fit4 = covfunc(example$t, example$y, method='PACE', kernel='gauss', mu=mu0f, newt=grid)
print(proc.time() - start_time)

start_time <- proc.time()
fit5 = face::face.sparse(example$data, argvals.new = grid, center = FALSE)
print(proc.time() - start_time)

start_time <- proc.time()
fit6 = tcrossprod(getA1_new_eig(Lt=example$t, Ly=example$y, r=NULL, mu=mu0f, newt=grid))
print(proc.time() - start_time)


print(c(mean((fit1$fitted - truth$truecovgrid)^2),
mean((fit2$fitted - truth$truecovgrid)^2),
mean((fit3$fitted - truth$truecovgrid)^2),
mean((fit4$fitted - truth$truecovgrid)^2),
mean((fit5$Chat.new - truth$truecovgrid)^2),
mean(((fit6) - truth$truecovgrid)^2)))

print(mean(truth$truecovgrid^2))

print(plot_ly(z = ~truth$truecovgrid, x = ~grid, y = ~grid, type='surface'))
print(plot_ly(z = ~fit1$fitted, x = ~grid, y = ~grid, type='surface'))
print(plot_ly(z = ~fit2$fitted, x = ~grid, y = ~grid, type='surface'))
print(plot_ly(z = ~fit3$fitted, x = ~grid, y = ~grid, type='surface'))
print(plot_ly(z = ~fit4$fitted, x = ~grid, y = ~grid, type='surface'))
print(plot_ly(z = ~fit5$Chat.new, x = ~grid, y = ~grid, type='surface'))
print(plot_ly(z = ~fit6, x = ~grid, y = ~grid, type='surface'))
}
system.time(regular.illustrate(params))
system.time(random.illustrate(params))

params = list(n=100, delta=0.2, m_avg=4, sige2=0.2722285, sigx=C.fourier, mu=mu0f, grid=regular.grid())
params = list(n=100, delta=0.2, m_avg=4, sige2=0.5916025, sigx=C.2, mu=mu0f, grid=regular.grid())
params = list(n=50, delta=0.2, m_avg=4, sige2=0.6890763, sigx=C.matern, mu=mu0f, grid=regular.grid())

grid = regular.grid()
delta = params[['delta']]
truth = sim.truth(params)
example = sim.regular(params)
fit1 = covfunc(example$t, example$y, method='FOURIER', ext=0.1, p=10, domain=c(0,1), mu=mu0f, newt=grid)
fit2 = covfunc(example$t, example$y, method='SP', domain=c(0,1), mu=mu0f, newt=regular.grid())
fit3 = face::face.sparse(example$data, argvals.new = grid, center = TRUE)
fit4 = covfunc(example$t, example$y, method='PACE', kernel='gauss', mu=mu0f, newt=grid)
fit5 = tcrossprod(getA1_new_eig(Lt=example$t, Ly=example$y, newt=grid, r=10, mu=mu0f))
fit6 = covfunc(example$t, example$y, method='BE', lam=0, ext=0.1, p=10, domain=c(0,1), mu=mu0f, newt=grid)

print(plot_ly(z = ~truth$truecovgrid, x = ~grid, y = ~grid, type='surface'))
print(plot_ly(z = ~fit1$fitted, x = ~grid, y = ~grid, type='surface'))
print(plot_ly(z = ~fit2$fitted, x = ~grid, y = ~grid, type='surface'))
print(plot_ly(z = ~fit3$Chat.new, x = ~grid, y = ~grid, type='surface'))
print(plot_ly(z = ~fit4$fitted, x = ~grid, y = ~grid, type='surface'))
print(plot_ly(z = ~fit5, x = ~grid, y = ~grid, type='surface'))
print(plot_ly(z = ~fit6$fitted, x = ~grid, y = ~grid, type='surface'))

r = NULL

devtools::load_all()
for(mu in c(mu0f)){
for(sigx in c(C.fourier, C.2, C.matern)){
    temp=sim.truth(list(sigx=sigx, mu=mu, grid=regular.grid()))
    print(c(mean(temp$truecovgrid^2), mean(diag(temp$truecovgrid^2)/4)))
  }
}

v_sige2 <- c(0.2722285, 0.5916025, 0.6890763)


# random
fit2 = covfunc(summary[[1]][[2]][[1]][[1]][[1]][[1]]$cache$ts, summary[[1]][[2]][[1]][[1]][[1]][[1]]$cache$ys, method='SP', newt=grid, domain=c(0,1))
fit2$fitted

tmp = pracma::meshgrid(grid)
x1 = tmp$X
x2 = tmp$Y
idx_Sdelta = (which(abs(x1 - x2) <= delta))
idx_Sc = (which(abs(x1 - x2) > delta))
Sdelta_truecovgrid = truth$truecovgrid
Sc_truecovgrid = truth$truecovgrid
Sdelta_truecovgrid[idx_Sc] = NaN
Sc_truecovgrid[idx_Sdelta] = NaN

fig = plot_ly(x = ~grid, y = ~grid, type='surface')
fig = fig %>% add_surface(z = ~Sc_truecovgrid)
fig %>% add_surface(z= ~Sdelta_truecovgrid)