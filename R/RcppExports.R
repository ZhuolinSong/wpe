# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

RgradQ <- function(Lt, Lr, L, B, U, V, W, lam, weight) {
    .Call('_wpe_RgradQ', PACKAGE = 'wpe', Lt, Lr, L, B, U, V, W, lam, weight)
}

RhessianQ <- function(Lt, Lr, S, B, U, V, W, lam, weight, X) {
    .Call('_wpe_RhessianQ', PACKAGE = 'wpe', Lt, Lr, S, B, U, V, W, lam, weight, X)
}

crawcov <- function(Lt, Ly, weig) {
    .Call('_wpe_crawcov', PACKAGE = 'wpe', Lt, Ly, weig)
}

csmoothcov <- function(h, kernel, xy, z, w, xgrid, ygrid, delta) {
    .Call('_wpe_csmoothcov', PACKAGE = 'wpe', h, kernel, xy, z, w, xgrid, ygrid, delta)
}

csmoothmean <- function(x, z, w, h, kernel, d, newx) {
    .Call('_wpe_csmoothmean', PACKAGE = 'wpe', x, z, w, h, kernel, d, newx)
}

