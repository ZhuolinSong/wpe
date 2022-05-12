devtools::load_all()

library(parallel)
RNGkind("L'Ecuyer-CMRG")
s_seed <- 999983
set.seed(s_seed)


v_n <- c(50, 200)
v_del <- c(0.2, 0.5)
v_C <- c("fourier", "matern", "bspline");v_sige2 <- c(6.292991, 0.6890763, 0.001649281)
v_mu <- c(1)
k <- 100

sparse_four <- mclapply(v_n, n_loop <- function(n){
    mclapply(v_del, del_loop <- function(del){
    mclapply(v_C, C_loop <- function(sigx){
    mclapply(v_mu, mu_loop <- function(mu){
    mclapply(1:k, sim_loop <- function(i){
        sige = v_sige2[match(sigx, v_C)]
        sim <- sim2(n, sige, del, m_avg = 4, sigx=sigx, mu=mu)
        grid <- regular.grid()
        fit3 = covfunc(sim$t, sim$y, method='FOURIER', newt=grid, ext=0.1, p=10, domain=c(0,1))

        return(list(data=sim$data, truemugrid=sim$truemugrid, truecovgrid=sim$truecovgrid,
                    mufour = fit3$mu$fitted, Cfour = fit3$fitted))
    }, mc.cores = 16)
    }, mc.cores = 1)
    }, mc.cores = 1)
    }, mc.cores = 1)
}, mc.cores = 1)
save(sparse_four, file = "sparse_four.RData")


random_four <- mclapply(v_n, n_loop <- function(n){
    mclapply(v_del, del_loop <- function(del){
    mclapply(v_C, C_loop <- function(sigx){
    mclapply(v_mu, mu_loop <- function(mu){
    mclapply(1:k, sim_loop <- function(i){
        sige = v_sige2[match(sigx, v_C)]
        sim <- sim3(n, sige, del, m_avg = 4, sigx=sigx, mu=mu)
        grid <- regular.grid()
        fit3 = covfunc(sim$t, sim$y, method='FOURIER', newt=grid, ext=0.1, p=10, domain=c(0,1))

        return(list(data=sim$data, truemugrid=sim$truemugrid, truecovgrid=sim$truecovgrid,
                    mufour = fit3$mu$fitted, Cfour = fit3$fitted))
    }, mc.cores = 16)
    }, mc.cores = 1)
    }, mc.cores = 1)
    }, mc.cores = 1)
}, mc.cores = 1)
save(random_four, file = "random_four.RData")