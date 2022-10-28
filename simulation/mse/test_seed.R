devtools::load_all()

library(parallel)
RNGkind("L'Ecuyer-CMRG")
s_seed <- 999983
set.seed(s_seed)

RNGkind("L'Ecuyer-CMRG")

rng1 <- unlist(
    mclapply(X = 1:3, FUN = function(x){
        mclapply(X = 1:3, FUN = function(x) runif(1)
        , mc.cores = 1, mc.set.seed = TRUE)
    }, mc.cores = 2, mc.set.seed = TRUE)
)

rng1



rng2 <- unlist(
    mclapply(X = 1:3, FUN = function(x){
        mclapply(X = 1:3, FUN = function(x) runif(1)
        , mc.cores = 2, mc.set.seed = F)
    }, mc.cores = 2, mc.set.seed = F)
)

rng2

# > [1] 0.67681994 0.54730337 0.05398847 0.19480448 0.94954659 0.35727778
#   [7] 0.17057359 0.83029494 0.37063552 0.24445617