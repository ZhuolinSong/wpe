devtools::load_all()
library(ggplot2)



illustration <- function(thess){
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
  long.data <- thess$data
  long.full <- thess$data.full

  # spaghetti plot
  long10 = long.data[long.data$subj %in% seq(10),]
  full10 = long.full[long.full$subj %in% seq(10),]
  pp <- ggplot(full10,aes(x=argvals,y=y, group=subj))
  qq1 <- pp + geom_point(color='grey', size=0.1) +
    geom_line(data=long10, color='black') +
    xlim(0, 1) + xlab("t") + ylab("X(t)") + theme_classic()
  print(qq1)


  grididx <- !thess$grid %in% long.data[long.data$subj == 1,]$argvals
  long.data <- data.frame(argvals=c(long.data$argvals, thess$grid[grididx]), subj=c(long.data$subj, rep(1, sum(grididx))), y=c(long.data$y, rep(NA, sum(grididx))))
  long.data <- long.data[order(long.data$argvals),]
  wide.data <- as.matrix(reshape(long.data, v.names="y",idvar="subj",
                      timevar="argvals",direction="wide"))[,-1]

  temp = wide.data
  temp[is.na(temp)] = 0
  temp = crossprod(temp)

  idx = which(temp != 0)
  x = c(matrix(rep(thess$grid, 100), 100, 100))[idx]
  y = c(matrix(rep(thess$grid, 100), 100, 100, byrow = T))[idx]
  pp <- ggplot(data.frame(x=x, y=y), aes(x=x,y=y))
  qq2 <- pp + geom_point(shape=18)+
    xlim(0, 1) + xlab("t") + ylab("t") + theme_classic()
  print(qq2)
  plot(x,y, xlim=c(1,0), ylim=c(1,0), pch=43)

  pdf("illustration.pdf",width=8, height = 4)
  multiplot(qq1,qq2,cols=2)
  dev.off()
}
set.seed(20220327)
thess <- sim.full(n=100,sige2=0,delta=.1, m_avg = 4, sigx="matern", mu=1)
illustration(thess)

illustration <- function(thess){
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
  long.data <- thess$data
  long.full <- thess$data.full

  # spaghetti plot
  long10 = long.data[long.data$subj %in% seq(10),]
  full10 = long.full[long.full$subj %in% seq(10),]
  pp <- ggplot(full10,aes(x=argvals,y=y, group=subj))
  qq1 <- pp + geom_point(color='grey', size=0.1) +
    geom_point(data=long10, color='black', size=0.5) +
    xlim(0, 1) + xlab("t") + ylab("X(t)") + theme_classic()
  print(qq1)


  grididx <- !thess$grid %in% long.data[long.data$subj == 1,]$argvals
  long.data <- data.frame(argvals=c(long.data$argvals, thess$grid[grididx]), subj=c(long.data$subj, rep(1, sum(grididx))), y=c(long.data$y, rep(NA, sum(grididx))))
  long.data <- long.data[order(long.data$argvals),]
  wide.data <- as.matrix(reshape(long.data, v.names="y",idvar="subj",
                      timevar="argvals",direction="wide"))[,-1]

  temp = wide.data
  temp[is.na(temp)] = 0
  temp = crossprod(temp)

  idx = which(temp != 0)
  x = c(matrix(rep(thess$grid, 100), 100, 100))[idx]
  y = c(matrix(rep(thess$grid, 100), 100, 100, byrow = T))[idx]
  pp <- ggplot(data.frame(x=x, y=y), aes(x=x,y=y))
  qq2 <- pp + geom_point(shape=18)+
    xlim(0, 1) + xlab("t") + ylab("t") + theme_classic()
  print(qq2)
  plot(x,y, xlim=c(1,0), ylim=c(1,0), pch=43)

  pdf("illustration2.pdf",width=8, height = 4)
  multiplot(qq1,qq2,cols=2)
  dev.off()
}
set.seed(20220327)
thess <- sim.full(n=100,sige2=0,delta=0.9, m_avg = 3, sigx="matern", mu=1)
illustration(thess)


require(plot3D)
persp3D(z = temp, theta = 0)

# 2 D local polynomial kernel smoother
out <- fdapace::Lwls2D(1,kern="epan",cbind(x,y),temp[temp!=0],win = NULL,
xout1 = NULL,xout2 = NULL,xout = cbind(x,y),subset = NULL,crosscov = FALSE)
length(out)
smooth.temp = temp
smooth.temp[temp != 0] = out

persp3D(z = smooth.temp, theta = 0)

Lwls2D(0.5, kernel, xin=rcov$tPairs, yin=rcov$cxxn,
                        xout1=cutRegGrid, xout2=cutRegGrid, win=covobj$weig,
                        delta=1)


fdapace::Lwls2D(0.4,kern="epan",xin=rcov$tPairs,yin=rcov$cxxn,win=covobj$weig,
xout1=cutRegGrid, xout2=cutRegGrid)




bmd = read.table('simulation/bone.data', header=T)
bmd$spnbmd
bmd <- within(bmd, 
{  idnum <- factor(idnum)
    gender <- factor(gender)
})
range(bmd$age)
bmd.female = bmd[which(bmd$gender=="female"),]
# spaghetti plot
pp <- ggplot(bmd.female,aes(x=age,y=spnbmd, group=idnum))
qq <- pp +  geom_line(color='gray') + geom_point(shape=18)+
  xlim(9, 26) + xlab("age") + ylab("spnbmd") + theme_classic()
print(qq)
