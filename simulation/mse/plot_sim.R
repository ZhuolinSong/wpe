devtools::load_all()
library(synfd)
library(fields)
library(fdapace)
library(ggplot2)

grid <- regular.grid()

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

C1stf = function(s,t){
  sig2(s) * matern(s,t,nu=0.5,rho=1)* sig2(t)
}
sig2 = function(t){
  sqrt(sqrt(t) * exp((-(t-.1)^2)/10) + 1)
}
C1f = function(ts){
  mat = matrix(NA,length(ts),length(ts))
  for(i in 1:length(ts)){
    for(j in 1:length(ts)){
      s = ts[i]
      t = ts[j]
      mat[i,j] = C1stf(s, t)
    }
  }
  return(mat)
}
mu1f = function(t){
  2 * t^2 * cos(2*pi*t)
}


# Caleb (sphagetti plot)
sim1 <- function(n,sige2,delta=.25, m_avg = 4){
  m = rpois(n,m_avg) + 1
  m[m<2] = 2
  Oi = runif(n,delta/2, 1-delta/2)
  t = lapply(1:n, function(i) sort(runif(m[i], Oi[i] - delta/2, Oi[i] + delta/2)))
  sigx = lapply(1:n,function(i) C1f(t[[i]]))
  truemugrid = mu1f(grid)
  truecovgrid = C1f(grid)
  mu = lapply(1:n,function(i) mu1f(t[[i]]))
  X = lapply(1:n, function(i) c(mvtnorm::rmvnorm(n=1, sigma=sigx[[i]])))
  eps = lapply(1:n, function(i) rnorm(m[i], 0, sqrt(sige2)))
  Y = lapply(1:n, function(i) mu[[i]] + X[[i]] + eps[[i]])
  ts = unlist(t)
  #ord = order(ts)
  #mus = unlist(mu)
  ys = unlist(Y)
  Zs = Cs = id = c()
  for(i in 1:length(t)){id = c(id, rep(i, length(t[[i]])))}
  data=data.frame(argvals=ts, subj=id, y=ys)
  return(list(data=data, m = m, grid=grid, t=t, y=Y, truemugrid=truemugrid, truecovgrid=truecovgrid))
}
# Zhenhua Lin
sim2 <- function(n_obj, sige2, del=.25, m_avg = 4, snr = NULL){
  D <- synfd::irreg.fd(mu=mu1f, X=synfd::gaussian.process(C1f), n=n_obj, m=m_avg, sig = sqrt(sige2), delta = del, snr = snr)
  ts = unlist(D$t)
  ys = unlist(D$y)
  truemugrid = mu1f(grid)
  truecovgrid = C1f(grid)
  Zs = Cs = id = c()
  for(i in seq_len(n_obj)){id = c(id, rep(i, length(D$t[[i]])))}
  data=data.frame(argvals=ts, subj=id, y=ys)
  return(list(data=data, t = D$t, y = D$y, D=D, truemugrid=truemugrid, truecovgrid=truecovgrid, truesig2 = attr(D,'sig')^2 ))
}


# generate data
set.seed(20220317)
ss1 <- sim1(10,.25,delta=.25, m_avg = 4)
df1 <- ss1$data
set.seed(20220317)
ss2 <- sim2(10,.25,del=.25, m_avg = 4, snr = NULL)
df2 <- ss2$data
dat <- data.frame(argvals=c(df1$argvals, df2$argvals), y=c(df1$y, df2$y), subj=c(df1$subj, df2$subj), sim=c(rep(0, length(df1$y)), rep(1, length(df2$y))))
dat <- within(dat, {  subj <- factor(subj)
    sim <- factor(sim,levels=0:1,labels=c("Caleb","Zhenhua"))
})


pdf("plot_sim.pdf",width=8)
pp <- ggplot(data=dat,aes(x=argvals,y=y,group=subj)) +  geom_line(color='gray') +
geom_smooth(aes(group = sim), size = 1.5, method="loess", se=FALSE) + 
geom_point(shape=18) + facet_grid(. ~ sim) + xlim(0, 1) + xlab("t") + ylab("y") +
geom_line(data = data.frame(argvals=grid, y=ss$truemugrid, subj = 1), size = 1.5)
print(pp)
dev.off()


ss1$t
ss2$t

ss1$y
ss2$y

fit2 = covfunc(ss1$t, ss1$y, method='FOURIER', tuning='lle',weig=NULL, domain=c(0,1))
fit2 = covfunc(ss2$t, ss2$y, method='FOURIER', tuning='lle',weig=NULL, domain=c(0,1))

