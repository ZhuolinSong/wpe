# wide simulation
#library(R.matlab)
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

dense_mat <- mclapply(v_n, n_loop <- function(n){
    mclapply(v_del, del_loop <- function(del){
    mclapply(v_C, C_loop <- function(sigx){
    mclapply(v_mu, mu_loop <- function(mu){
    mclapply(1:k, sim_loop <- function(i){
        sige = v_sige2[match(sigx, v_C)]
        sim <- sim1(n, sige, del, m_avg = 50, sigx=sigx, mu=mu)
        return(to.wide(sim$data, sim$grid))

    }, mc.cores = 8)
    }, mc.cores = 1)
    }, mc.cores = 1)
    }, mc.cores = 1)
}, mc.cores = 1)
save(dense_mat, file = "dense_mat.RData")


sparse_mat <- mclapply(v_n, n_loop <- function(n){
    mclapply(v_del, del_loop <- function(del){
    mclapply(v_C, C_loop <- function(sigx){
    mclapply(v_mu, mu_loop <- function(mu){
    mclapply(1:k, sim_loop <- function(i){
        sige = v_sige2[match(sigx, v_C)]
        sim <- sim2(n, sige, del, m_avg = 4, sigx=sigx, mu=mu)
        return(to.wide(sim$data, sim$grid))

    }, mc.cores = 8)
    }, mc.cores = 1)
    }, mc.cores = 1)
    }, mc.cores = 1)
}, mc.cores = 1)
save(sparse_mat, file = "sparse_mat.RData")


random_mat <- mclapply(v_n, n_loop <- function(n){
    mclapply(v_del, del_loop <- function(del){
    mclapply(v_C, C_loop <- function(sigx){
    mclapply(v_mu, mu_loop <- function(mu){
    mclapply(1:k, sim_loop <- function(i){
        sige = v_sige2[match(sigx, v_C)]
        sim <- sim3(n, sige, del, m_avg = 4, sigx=sigx, mu=mu)
        return(to.wide(sim$data, sim$grid))

    }, mc.cores = 8)
    }, mc.cores = 1)
    }, mc.cores = 1)
    }, mc.cores = 1)
}, mc.cores = 1)
save(random_mat, file = "random_mat.RData")



# filename <- paste(tempfile(), ".mat", sep = "")
# writeMat(con=filename, x=wide.data)
# 
# filename <- paste(tempfile(), ".mat", sep = "")
# writeMat(con=filename, t=sort(thess$grid))
# 
# filename <- paste(tempfile(), ".mat", sep = "")
# writeMat(con=filename, grid=grid)