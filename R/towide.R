to.wide = function(long.data, grid){
  grididx <- !grid %in% long.data[long.data$subj == 1,]$argvals
  long.data <- data.frame(argvals=c(long.data$argvals, grid[grididx]), subj=c(long.data$subj, rep(1, sum(grididx))), y=c(long.data$y, rep(NA, sum(grididx))))
  long.data <- long.data[order(long.data$argvals),]
  wide.data <- as.matrix(reshape(long.data, v.names="y",idvar="subj",
                      timevar="argvals",direction="wide"))[,-1]
  return(wide.data)
}