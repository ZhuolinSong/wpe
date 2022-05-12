# wpc
Written Preliminary Exam or Literature Review: Methods for Functional Snippets. (modified from the [mcfda](https://github.com/linulysses/mcfda) package)

## Install
devtools::install_github("zhuolinsong/wpe")


## Demonstration
library(wpe)

set.seed(1)

### set grid and simulate functional snippet dataset
grid <- regular.grid()

sim <- sim1(50, 0.1, delta=.2, m_avg = 4, sigx='matern', mu=1)

### spaghetti plot
library(ggplot2)

pp <- ggplot(data=sim$data,aes(x=argvals,y=y,group=subj)) +  
geom_line(color='gray') + 
geom_point(shape=18) + xlim(0, 1) + xlab("t") + ylab("Y(t)") +
geom_line(data = data.frame(argvals=grid, y=sim$truemugrid, subj = 1), size = 1.5)

print(pp)

### estimate by 'FACE' method
fit1 = face::face.sparse(sim$data,argvals.new = grid)

### estimate by 'PACE' method
fit2 = covfunc(sim$t, sim$y, method='PACE', newt=grid, bw=sqrt(2)/2*(1-0.2), delta=1)

### estimate by 'FOURIER' method
fit3 = covfunc(sim$t, sim$y, method='FOURIER', newt=grid, ext=0.1, p=10, domain=c(0,1))

### estimate by 'SP' method
fit4 = covfunc(sim$t, sim$y, method='SP', newt=grid, domain=c(0,1))

### estimate by 'SAMC' method
X = t(to.wide(sim$data,sim$grid)); p=length(sim$grid)
A=getA1_new_cv(X, p, dl=0.9*100*0.5, incre=0.1*100*0.5, bw = NA, kernel = 'epan',sigma2hat= NA)$A
fit5 = tcrossprod(A)

### estract mean and cov estimation over grid
sim1 = list(data=sim$data, truemugrid=sim$truemugrid, truecovgrid=sim$truecovgrid,
  muface = fit1$mu.new, Cface = fit1$Chat.new,
  mupace = fit2$mu$fitted, Cpace = fit2$fitted,
  mufour = fit3$mu$fitted, Cfour = fit3$fitted,
  musnpt = fit4$mu$fitted, Csnpt = fit4$fitted,
  Csamc = fit5)

### calculate relative error in mise
aggf = function(u){
  mumseface = mean((u$muface - u$truemugrid)^2)/mean(u$truemugrid^2)
  covmseface = mean((u$Cface - u$truecovgrid)^2)/mean(u$truemugrid^2)

  mumsepace = mean((u$mupace - u$truemugrid)^2)/mean(u$truemugrid^2)
  covmsepace = mean((u$Cpace - u$truecovgrid)^2)/mean(u$truemugrid^2)

  mumsefour = mean((u$mufour - u$truemugrid)^2)/mean(u$truemugrid^2)
  covmsefour = mean((u$Cfour - u$truecovgrid)^2)/mean(u$truemugrid^2)

  mumsesnpt = mean((u$musnpt - u$truemugrid)^2)/mean(u$truemugrid^2)
  covmsesnpt = mean((u$Csnpt - u$truecovgrid)^2)/mean(u$truemugrid^2)

  covmsesamc = mean((u$Csamc - u$truecovgrid)^2)/mean(u$truemugrid^2)
  return(c(mumseface,covmseface,
           mumsepace,covmsepace,
           mumsefour,covmsefour,
           mumsesnpt,covmsesnpt,
           covmsesamc))
}

agg1 = aggf(sim1)
A = matrix(agg1, 1)
colnames(A) = c("mumseface","covmseface",
                "mumsepace","covmsepace",
                "mumsefour","covmsefour",
                "mumsesnpt","covmsesnpt",
                "covmsesamc")
print(A)


## References
Lin, Z. and Wang, J.-L. (2020+). [Mean and covariance estimation for functional snippets](https://www.tandfonline.com/doi/full/10.1080/01621459.2020.1777138). Journal of the American Statistical Association. Volume 117, 2022 - Issue 537

Lin, Z., Wang, J.-L. and Zhong, Q. (2021). [Basis expansions for functional snippets](https://academic.oup.com/biomet/article-abstract/108/3/709/5937818?redirectedFrom=fulltext#no-access-message). Biometrika, Volume 108, Issue 3, September 2021, Pages 709–726

Zhang, A., and Chen, K. (2022). [Nonparametric covariance estimation for mixed longitudinal studies, with applications in midlife women's health](http://www3.stat.sinica.edu.tw/LatestART/SS-2019-0219_fp.pdf). Statistica Sinica 32 (2022), 1-21.


Yao, F., Müller, H.-G. and Wang, J.-L. (2005). [Functional Data Analysis for Sparse Longitudinal Data](https://www.tandfonline.com/doi/abs/10.1198/016214504000001745). Journal of the American Statistical Association. 100(470): 577-590.
