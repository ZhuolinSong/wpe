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

dense_zhang <- mclapply(v_n, n_loop <- function(n){
    mclapply(v_del, del_loop <- function(del){
    mclapply(v_C, C_loop <- function(sigx){
    mclapply(v_mu, mu_loop <- function(mu){
    mclapply(1:k, sim_loop <- function(i){
        sige = v_sige2[match(sigx, v_C)]
        sim <- sim1(n, sige, del, m_avg = 20, sigx=sigx, mu=mu)
        grid <- regular.grid()
        X = t(to.wide(sim$data,sim$grid)); p=length(sim$grid)
        A=getA1_new_cv(X, p, dl=floor(0.9*100*del), incre=floor(0.1*100*del), bw = NA, kernel = 'epan', sigma2hat= NA)$A
        fit5 = tcrossprod(A)

        return(list(data=sim$data, truemugrid=sim$truemugrid, truecovgrid=sim$truecovgrid,
                    Czhang = fit5))
    }, mc.cores = 8)
    }, mc.cores = 1)
    }, mc.cores = 1)
    }, mc.cores = 1)
}, mc.cores = 1)
save(dense_zhang, file = "dense_zhang.RData")