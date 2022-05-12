library(xtable)

v_n <- c(50, 200)
v_del <- c(0.2, 0.5)
v_C <- c("fourier", "matern", "bspline")
v_mu <- c(1)
k <- 100
summarize = function(input){
    size=length(v_n)*length(v_del)*length(v_C)
    result = matrix(NA, nrow=size, ncol=2)
    for(n in seq_along(v_n)){
        for(del in seq_along(v_del)){
        for(C in seq_along(v_C)){
        for(mu in seq_along(v_mu)){
            temp = lapply(input[[n]][[del]][[C]][[mu]], agg)
            A = do.call(rbind, temp)
            A = A[which(!is.na(A)),]
            A = as.matrix(A)
            temp1 = apply(A, 2, mean)
            temp2 = apply(A, 2, sd)
            print(length(A))
            idx = (C-1)*4+(n-1)*2+del
            result[idx,1] = temp1[1]; result[idx,2] = temp2[1]
    }}}}
    #return(result)
    print(xtable(result, digits=2))
}

# face
agg = function(u){
  covmseface = mean((u$Cface - u$truecovgrid)^2)/mean((u$truecovgrid)^2)
  return(c(covmseface))
}
load("simulation/result/dense_face.RData")
summarize(dense_face)
load("simulation/result/sparse_face.RData")
summarize(sparse_face)
load("simulation/result/random_face.RData")
summarize(random_face)


# pace
agg = function(u){
  if(is.atomic(u)) return(NA)
  covmseface = mean((u$Cpace - u$truecovgrid)^2)/mean((u$truecovgrid)^2)
  return(c(covmseface))
}
load("simulation/result/dense_pace.RData")
summarize(dense_pace)
load("simulation/result/sparse_pace.RData")
summarize(sparse_pace)
load("simulation/result/random_pace.RData")
summarize(random_pace)


# four
agg = function(u){
  if(is.atomic(u)) return(NA)
  covmseface = mean((u$Cfour - u$truecovgrid)^2)/mean((u$truecovgrid)^2)
  return(c(covmseface))
}
load("simulation/result/dense_four.RData")
summarize(dense_four)
load("simulation/result/sparse_four.RData")
summarize(sparse_four)
load("simulation/result/random_four.RData")
summarize(random_four)


# sp
agg = function(u){
  if(is.atomic(u)) return(NA)
  covmseface = mean((u$Csnpt - u$truecovgrid)^2)/mean((u$truecovgrid)^2)
  return(c(covmseface))
}
load("simulation/result/dense_snpt.RData")
summarize(dense_snpt)
load("simulation/result/sparse_snpt.RData")
summarize(sparse_snpt)
load("simulation/result/random_snpt.RData")
summarize(random_snpt)


# zhang
agg = function(u){
  covmseface = mean((u$Czhang - u$truecovgrid)^2)/mean((u$truecovgrid)^2)
  return(c(covmseface))
}
load("simulation/result/dense_zhang.RData")
summarize(dense_zhang)
load("simulation/result/sparse_zhang.RData")
summarize(sparse_zhang)
