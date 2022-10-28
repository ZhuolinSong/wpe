devtools::load_all()

library(parallel)
RNGkind("L'Ecuyer-CMRG")
s_seed <- 13579+40
set.seed(s_seed)


v_n <- c(100, 200)
v_delta <- c(0.2, 0.4, 0.6)
v_m = c(20)
v_C <- c(C.fourier, C.2, C.matern);v_sige2 <- c(0.2722285, 0.5916025, 0.6890763)
v_i = seq_along(v_C)
v_mu <- c(mu0f)
k <- 50
grid = regular.grid()

truth_set <- mclapply(v_n, n_loop <- function(n){
    lapply(v_delta, FUN = function(delta){
    lapply(v_m, FUN = function(m){
    lapply(v_i, FUN = function(idx){
    lapply(v_mu, FUN = function(mu){
        params = list(n=n, delta=delta, m_avg=m, sige2=v_sige2[idx], sigx=v_C[[idx]], mu=mu, grid=grid)
        truth = sim.truth(params)
        return(list(truemugrid = truth$truemugrid, truecovgrid=truth$truecovgrid, params = params))
})
})
})
})}, mc.cores = 16, mc.set.seed = TRUE)

dense_regular_design <-  mclapply(41:60, FUN = function(i){
    lapply(seq_along(v_n), FUN = function(n){
    lapply(seq_along(v_delta), FUN = function(delta){
    lapply(seq_along(v_m), FUN = function(m){
    lapply(seq_along(v_i), FUN = function(idx){
    lapply(seq_along(v_mu), FUN = function(mu){
        truemugrid = truth_set[[n]][[delta]][[m]][[idx]][[mu]]$truemugrid
        truecovgrid = truth_set[[n]][[delta]][[m]][[idx]][[mu]]$truecovgrid
        params = truth_set[[n]][[delta]][[m]][[idx]][[mu]]$params

        sim <- sim.regular(params)

        fit1 = covfunc(sim$t, sim$y, method='BE', lam=0, ext=0, domain=c(0,1), mu=mu0f, newt=grid)
        fit2 = covfunc(sim$t, sim$y, method='FOURIER', domain=c(0,1), mu=mu0f, newt=grid)
        fit3 = covfunc(sim$t, sim$y, method='SP', domain=c(0,1), mu=mu0f, newt=grid)
        fit4 = covfunc(sim$t, sim$y, method='PACE', kernel='gauss', mu=mu0f, newt=grid)
        fit5 = face::face.sparse(sim$data, argvals.new = grid, center = FALSE)

        if (idx==3){
            fit6 = getA1_new_eig(sim$t, sim$y, r=NULL, newt=grid, mu=mu0f)
        }else if(idx==2){
            fit6 = getA1_new_eig(sim$t, sim$y, r=NULL, newt=grid, mu=mu0f)
        }else if(idx==1){
            fit6 = getA1_new_eig(sim$t, sim$y, r=5, newt=grid, mu=mu0f)
        }

        cov1 = fit1$fitted
        cov1.ise = mean((truecovgrid-cov1)^2)

        cov2 = fit2$fitted
        cov2.ise = mean((truecovgrid-cov2)^2)

        cov3 = fit3$fitted
        cov3.ise = mean((truecovgrid-cov3)^2)

        cov4 = fit4$fitted
        cov4.ise = mean((truecovgrid-cov4)^2)

        cov5 = fit5$Chat.new
        cov5.ise = mean((truecovgrid-cov5)^2)

        cov6 = tcrossprod(fit6)
        cov6.ise = mean((truecovgrid-cov6)^2)

        return(list(cache = list(truecovgrid=truecovgrid, ts = sim$t, ys = sim$y),
                    cov1 = cov1, cov1.ise = cov1.ise,
                    cov2 = cov2, cov2.ise = cov2.ise,
                    cov3 = cov3, cov3.ise = cov3.ise,
                    cov4 = cov4, cov4.ise = cov4.ise,
                    cov5 = cov5, cov5.ise = cov5.ise,
                    cov6 = cov6, cov6.ise = cov6.ise
                    ))
})
})
})
})
})
}, mc.cores = 16, mc.set.seed = TRUE)
save(dense_regular_design, file = "dense_regular_design20c.RData")