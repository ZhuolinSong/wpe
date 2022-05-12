# The algorithm with cross-validation, possible smoothing, and incomplete samples
getA1_new_cv <- function(X, p, dl, incre, bw = NA, kernel = 'epan', sigma2hat= NA){
  # Estimation of A with cross-validation
  # Cross validation to select the rank r
  # Range for candidate r's: 2 to (dl-incre)
  # 5-fold cross validation, i.e. 4/5 samples for estimation, 1/5 for validation
  # Return: hat_A and hat_r
  n = dim(X)[2]; n_train = round(n * 4/5); n_test = n - n_train
  lmax = 10  # repetition time
  
  X0 = X; X0[is.na(X0)] = 0; # observations with NA's replaced by zeros
  M = 1 - is.na(X) # indicator of observable entries 
  cross_error = rep(0, dl-incre)
  
  for (l in 1:lmax){
    samp_train = sample(n, n_train)
    samp_test = setdiff((1:n), samp_train)
    X_train = X[, samp_train]
    X_test = X0[, samp_test]
    # Idea for Constructing testing covariance: use all entries of Sigma_test which appears for at least 4 times.
    # then compare the entries of Sigma_train and Sigma_test
    Sigma_test = X_test %*% t(X_test) / pmax(1, M[, samp_test] %*% t(M[, samp_test]))
    M_test = ((M[, samp_test] %*% t(M[, samp_test]))>=4)*(1-diag(p))
    for (r in 2:(dl-incre)){
      A_train <- getA1_new_eig(X, p, dl, incre, r, bw, kernel, sigma2hat) 
      Sigma_train = A_train %*% t(A_train) # Training covariance
      cross_error[r] = cross_error[r] + norm((Sigma_train-Sigma_test)*M_test, "F")^2
    }
  }
  hat_r = which.min(cross_error[2:(dl-incre)]) + 1
  A <- getA1_new_eig(X, p, dl, incre, hat_r, bw, kernel, sigma2hat)
  return( list(A = A, r = hat_r))
}



getA1_new_eig <- function(X, p, dl, incre, r, bw = NA, kernel= 'epan', sigma2hat = NA){ 
  # getA1_new_eig use all, including incomplete, samples to calculate each piece of the A matrix, with or without smoothing.
  # Input: 
  # X: data matrix, with NaNs as unobservable entries
  # p: dimension of the samples
  # dl: size of block
  # incre: increment of concatenation in each step
  # r: estimated rank
  # bw: a vector with length r if apply additional smoothing, no smoothing by default
  # Output: hat_A 
  
  X0 = X; X0[is.na(X0)] = 0;
  M = 1 - is.na(X)
  
  piece = ceiling((p-dl)/incre) # number of incremental intervals we need to construct A
  
  #n.piece = rep(0, piece) # number of samples within each inverval
  A = matrix(0, ncol=r, nrow=p)
  
  #A_pre_aggre = matrix(NA, nrow=p, ncol=r*(piece+1))
  A_pre_aggre = array(NA, dim=c(p, r, piece+1))
  # A pre-aggregation 3-d array
  # A_pre_aggre contains all A.hat information: A_1.hat, A_2.hat ......
  for (l in 1:piece){
    # Construction of sample covariance matrix for [(l-1)*incre+1: ((l-1)*incre+dl)]
    n_count = (M[((l-1)*incre+1): ((l-1)*incre+dl),] %*% t(M[((l-1)*incre+1): ((l-1)*incre+dl),]))
    XX_trans = X0[((l-1)*incre+1): ((l-1)*incre+dl),] %*% t(X0[((l-1)*incre+1): ((l-1)*incre+dl),])
    sample.cov = XX_trans/pmax(1, n_count)
    # svd.cov is the sample covariance for Sigma[((l-1)*dl/2+1): ((l+1)*dl/2), ((l-1)*dl/2+1): ((l+1)*dl/2)]
    eig.cov = eigen(sample.cov) 
    if(!is.na(sigma2hat)){
      hat.d = sigma2hat
    }else{ 
      hat.d = max(mean(eig.cov$values[-(1:r)]), 0)
    }
    A.hat = eig.cov$vectors[,1:r] %*% diag(sqrt(pmax(0, eig.cov$values[1:r]-hat.d))) # A.tilde is the factorization of sample.cov
    if(l == 1){
      A_pre_aggre[1:dl, 1:r, l] = A.hat
    }else{
      # the next two lines calculate the best angle to embed A.tilde to A
      svd.rota = svd(t(A.hat[1:(dl-incre), ]) %*% A_pre_aggre[((l-1)*incre+1):((l-2)*incre+dl), 1:r, l-1])
      rota = svd.rota$u %*% t(svd.rota$v)
      A_pre_aggre[((l-1)*incre+1):((l-1)*incre+dl), 1:r, l] = A.hat %*% rota
    }
  }
  # Construction Last piece
  l = piece+1
  #  sample.cov = X0[((l-1)*incre+1): p,] %*% t(X0[((l-1)*incre+1): p,])/(M[((l-1)*incre+1): p,] %*% t(M[((l-1)*incre+1): p,]))
  ind.sample = (colSums(is.na(X[((l-1)*incre+1): p, ])) == 0) 
  # search for all samples without NaN's in [((l-1)*dl/2+1): (l*dl/2)]
  # ind.sample is the set of samples which can be used in calculating the corresponding part of A: A[((l-1)*dl/2+1): (l*dl/2),]
  n_count = (M[((l-1)*incre+1): p,] %*% t(M[((l-1)*incre+1): p,]))
  XX_trans = X0[((l-1)*incre+1): p,] %*% t(X0[((l-1)*incre+1): p,])
  sample.cov = XX_trans/pmax(1, n_count)
  
  # svd.cov is the sample covariance for Sigma[((l-1)*dl/2+1): ((l+1)*dl/2), ((l-1)*dl/2+1): ((l+1)*dl/2)]
  eig.cov = eigen(sample.cov) 
  if(!is.na(sigma2hat)){
    hat.d = sigma2hat
  }else{ 
    hat.d = max(mean(eig.cov$values[-(1:r)]), 0)
  }
  A.hat = eig.cov$vectors[,1:r] %*% diag(sqrt(pmax(0, eig.cov$values[1:r]-hat.d))) # A.tilde is the factorization of sample.cov
  svd.rota = svd(t(A.hat[1:(dl-incre), ]) %*% A_pre_aggre[((l-1)*incre+1):((l-2)*incre+dl), 1:r, l-1])
  rota = svd.rota$u %*% t(svd.rota$v)
  #A_pre_aggre[((l-1)*incre+1): p, 1:r, l] = A.hat[(dl-incre+1): dim(A.hat)[1], ] %*% rota 
  A_pre_aggre[((l-1)*incre+1): p, 1:r, l] = A.hat %*% rota 
  
  #A = f.aggre(A_pre_aggre) #Just average without smoothing
  A = f.aggre.smooth(A_pre_aggre, p, r, bw, kernel) # Average with/without smoothing
  return(A)
}


f.aggre <- function(A_pre_aggre, p, dl, incre, r) {#Aggregation function
  A = apply(A_pre_aggre, c(1,2), mean.rmna) # no smoothing, just the average of all hat.A's
  return(A)
}

f.aggre.smooth <- function(A_pre_aggre, p, r, bw, kernel = 'epan'){
	if (is.na(bw)){
	   A = apply(A_pre_aggre, c(1,2), mean.rmna)# no smoothing, just the average of all hat.A's
	   return(A)
	} else {
		if (length(bw) ==1){
			bw = rep(bw, r)
		}
		q = dim(A_pre_aggre)[3]
		A.smooth = matrix(0, p, r)
		for (i in 1:r){
			yin = matrix(A_pre_aggre[,i,], 1, p*q)
			xin = matrix(1:p, 1, p*q)
			win = (!is.na(yin))
			xou = 1:p
			A.smooth[,i]  = lwls(bw[i], xin = xin, yin = yin, win = win, xou = xou, kernel = kernel)$mu
		}
		return(A.smooth)
	}
}

# f.aggre.smooth2 <- function(A_pre_aggre, p, r, bw, kernel = 'epan'){
# 	if (is.na(bw)){
# 	   A = apply(A_pre_aggre, c(1,2), mean.rmna)# no smoothing, just the average of all hat.A's
# 	   return(A)
# 	} else {
# 		if (length(bw) ==1){
# 			bw = rep(bw, r)
# 		}
# 		A = apply(A_pre_aggre, c(1,2), mean.rmna)# no smoothing, just the average of all hat.A's
# 		A.smooth = matrix(0, p, r)
# 		for (i in 1:r){
# 			yin = A[,i]
# 			xin = 1:p
# 			win = (!is.na(yin))
# 			xou = 1:p
# 			A.smooth[,i]  = lwls(bw[i], xin = xin, yin = yin, win = win, xou = xou, kernel = kernel)$mu
# 		}
# 		return(A.smooth)
# 	}
# }

mean.rmna <- function(v){
  return(mean(v, na.rm=TRUE))
}

lwls <- function(bw, kernel = c("epan","gauss", "rect", "quar"), npoly = nder+1, nder = 0, xin, yin, win, xou){

   if(npoly < nder)
      stop("Degree of polynomial should be no less than the order of derivative!")
   
   kernel = kernel[1]
   require(MASS)
   actobs = which(win != 0)
   xin = xin[actobs]
   yin = yin[actobs]
   win = win[actobs]
   invalid = 0

   aa = 1
   mu = numeric(length(xou))
   gap = mu
   bw = rep(bw,length(xou)) 
  
   
   #LWLS with different weight functions
   for(i in 1:length(xou)){

       #(3-1) Locating local window
       if(kernel != "gauss"){
          idx = xin <= xou[i] + aa*bw[i] & xin >= xou[i]-aa*bw[i]       
       }else{
          idx = 1:length(xin)
       }
       lx = xin[idx]
       ly = yin[idx]
       lw = win[idx]

       if(length(unique(lx)) >= (npoly+1)){

          #Sepcify weight matrix
          llx = (lx-xou[i])/bw[i]
          
          if(kernel == "epan"){
              w = lw*(1-llx^2)*0.75
          }else if(kernel == "rect"){
              w = lw
          }else if(kernel == "gauss"){
              w = lw*dnorm(llx)
          }else if(kernel == "quar"){
              w = lw*(15/16)*(1-llx^2)^2
          }else{
              cat("Invalid kernel, Epanechnikov kernel is used!\n")
              w = lw*(1-llx^2)*0.75
          }
          W = diag(w, length(w), length(w))
          # Define design matrix
          dx = matrix(1,length(lx),npoly+1)
          for(j in 1:npoly){
            dx[,j+1] = (xou[i]-lx)^j
          }       
          
          p = ginv(t(dx)%*%W%*%dx)%*%t(dx)%*%W%*%ly  

          #Find estimate
          mu[i] = p[(nder+1)*gamma(nder+1)*((-1)^nder)]          
          
       }else{
          gap[i] = 1
          invalid = 1
       }

   }
   indx = which(gap == 0)
   if((length(indx)>= 0.9*length(xou))&& (length(indx)<length(xou))){
        mu1 = mu[indx]
        rr = myunique(xou[indx])
        mu=interp1(rr$out1,mu1[rr$id],xou);
   }else if(length(indx) < 0.9*length(xou)){
        mu = NULL
        cat("Too many gaps, please increase bandwidth!\n")
        invalid = 1
   }
   return(list(invalid = invalid, mu = mu))
}



getA2_new_eig <- function(X, p, dl, incre, r, bw = NA, kernel= 'epan', sigma2hat = NA){ 
  # use the possibly incomplete samples that has sufficient coverage to calculate each piece of the A matrix, with/without smoothing.
  # Input: 
  # X: data matrix, with NaNs as unobservable entries
  # p: dimension of the samples
  # dl: size of block
  # incre: increment of concatenation in each step
  # r: estimated rank
  # bw: a vector with length r if apply additional smoothing, no smoothing by default
  # Output: hat_A 
  X0 = X; X0[is.na(X0)] = 0;
  M = 1 - is.na(X)
  
  n = dim(X)[2]
  piece = ceiling((p-dl)/incre) # number of incremental intervals we need to construct A
  
  #n.piece = rep(0, piece) # number of samples within each inverval
  A = matrix(0, ncol=r, nrow=p)
  
  #A_pre_aggre = matrix(NA, nrow=p, ncol=r*(piece+1))
  A_pre_aggre = array(NA, dim=c(p, r, piece+1))
  # A pre-aggregation 3-d array
  # A_pre_aggre contains all A.hat information: A_1.hat, A_2.hat ......
  indlist = matrix(0, piece+1, n)
  for (v in 1:n){
    temp = which(M[,v] == 1)
    tempmin = ceiling((temp[1]-1)/incre+1); tempmax = floor((temp[length(temp)]-dl)/incre+1)
    if(tempmax >= tempmin)
      indlist[tempmin:min(tempmax, p),v] = 1
  }
  for (l in 1:piece){
    # Construction of sample covariance matrix for [(l-1)*incre+1: ((l-1)*incre+dl)]
    ind = which(indlist[l,]==1)
    n_count = M[((l-1)*incre+1): ((l-1)*incre+dl),ind] %*% t(M[((l-1)*incre+1): ((l-1)*incre+dl),ind])
    XX_trans = X0[((l-1)*incre+1): ((l-1)*incre+dl),ind] %*% t(X0[((l-1)*incre+1): ((l-1)*incre+dl),ind])
    sample.cov = XX_trans/pmax(1, n_count)
    # svd.cov is the sample covariance for Sigma[((l-1)*dl/2+1): ((l+1)*dl/2), ((l-1)*dl/2+1): ((l+1)*dl/2)]
    eig.cov = eigen(sample.cov) 
    if(!is.na(sigma2hat)){
      hat.d = sigma2hat
    }else{ 
      hat.d = max(mean(eig.cov$values[-(1:r)]), 0)
    }
    A.hat = eig.cov$vectors[,1:r] %*% diag(sqrt(pmax(0, eig.cov$values[1:r]-hat.d))) # A.tilde is the factorization of sample.cov
    if(l == 1){
      A_pre_aggre[1:dl, 1:r, l] = A.hat
    }else{
      # the next two lines calculate the best angle to embed A.tilde to A
      svd.rota = svd(t(A.hat[1:(dl-incre), ]) %*% A_pre_aggre[((l-1)*incre+1):((l-2)*incre+dl), 1:r, l-1])
      rota = svd.rota$u %*% t(svd.rota$v)
      A_pre_aggre[((l-1)*incre+1):((l-1)*incre+dl), 1:r, l] = A.hat %*% rota
    }
  }
  # Construct the Last piece
  l = piece+1
  #  sample.cov = X0[((l-1)*incre+1): p,] %*% t(X0[((l-1)*incre+1): p,])/(M[((l-1)*incre+1): p,] %*% t(M[((l-1)*incre+1): p,]))
  ind = which(indlist[l,]==1)
  
  n_count = M[((l-1)*incre+1): p,ind] %*% t(M[((l-1)*incre+1): p,ind])
  XX_trans = X0[((l-1)*incre+1): p,ind] %*% t(X0[((l-1)*incre+1): p,ind])
  sample.cov = XX_trans/pmax(1, n_count)
  
  # svd.cov is the sample covariance for Sigma[((l-1)*dl/2+1): ((l+1)*dl/2), ((l-1)*dl/2+1): ((l+1)*dl/2)]
  eig.cov = eigen(sample.cov) 
  if(!is.na(sigma2hat)){
    hat.d = sigma2hat
  }else{ 
    hat.d = max(mean(eig.cov$values[-(1:r)]), 0)
  }
  A.hat = eig.cov$vectors[,1:r] %*% diag(sqrt(pmax(0, eig.cov$values[1:r]-hat.d))) # A.tilde is the factorization of sample.cov
  svd.rota = svd(t(A.hat[1:(dl-incre), ]) %*% A_pre_aggre[((l-1)*incre+1):((l-2)*incre+dl), 1:r, l-1])
  rota = svd.rota$u %*% t(svd.rota$v)
  #A_pre_aggre[((l-1)*incre+1): p, 1:r, l] = A.hat[(dl-incre+1): dim(A.hat)[1], ] %*% rota 
  A_pre_aggre[((l-1)*incre+1): p, 1:r, l] = A.hat %*% rota 
  
  #A = f.aggre(A_pre_aggre) #Just average without smoothing
  A = f.aggre.smooth(A_pre_aggre, p, r, bw, kernel) # Average with/without smoothing
  return(A)
}



getA2_new_cv <- function(X, p, dl, incre, bw = NA, kernel = 'epan', sigma2hat = NA){
  # Estimation of A with cross-validation
  # Cross validation to select the rank r
  # Range for candidate r's: 2 to (dl-incre)
  # 5-fold cross validation, i.e. 4/5 samples for estimation, 1/5 for validation
  # Return: hat_A and hat_r
  n = dim(X)[2]; n_train = round(n * 4/5); n_test = n - n_train
  lmax = 10  # repetition time
  
  X0 = X; X0[is.na(X0)] = 0; # observations with NA's replaced by zeros
  M = 1 - is.na(X) # indicator of observable entries
  
  cross_error = rep(0, dl-incre)
  
  for (l in 1:lmax){
    samp_train = sample(n, n_train)
    samp_test = setdiff((1:n), samp_train)
    X_train = X[, samp_train]
    X_test = X0[, samp_test]
    # Idea for Constructing testing covariance: use all entries of Sigma_test which appears for at least 4 times.
    # then compare the entries of Sigma_train and Sigma_test
    Sigma_test = X_test %*% t(X_test) / pmax(1, M[, samp_test] %*% t(M[, samp_test]))
    M_test = ((M[, samp_test] %*% t(M[, samp_test]))>=4)*(1-diag(p))
    for (r in 2:(dl-incre)){
      A_train <- getA2_new_eig(X, p, dl, incre, r, bw, kernel,sigma2hat) 
      Sigma_train = A_train %*% t(A_train) # Training covariance
      cross_error[r] = cross_error[r] + norm((Sigma_train-Sigma_test)*M_test, "F")^2
    }
  }
  hat_r = which.min(cross_error[2:(dl-incre)]) + 1
  A <- getA2_new_eig(X, p, dl, incre, hat_r, bw, kernel,sigma2hat)
  return( list(A =A, r= hat_r))
}
