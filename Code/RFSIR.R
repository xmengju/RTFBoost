



j.0.5 <- function(x){
  z <- svd(x)
  S <- z$v %*% diag(sqrt(z$d)) %*% t(z$v)
  S
}

# spatial median 
SpMed=function(x){
  
  n=nrow(x)
  p=ncol(x)
  t=seq(0,1,length=p)
  
  # similar as pairwise distance 
  A=0.5*x[,(1:(p-1))]%*%diag(t[2:p]-t[1:(p-1)])%*%t(x[,(1:(p-1))])+0.5*x[,(2:p)]%*%diag(t[2:p]-t[1:(p-1)])%*%t(x[,(2:p)])
  w=matrix(c(rep(1,n)),n,1) # equal weight 
  norms=sqrt(diag(A)+ as.numeric(t(w)%*%A%*%(w))-2*A%*%w) #last equation on page 233
  f=sum(norms)
  err=1
  
  iter=0
  while(err>1e-0& iter<50){
    iter=iter+1
    f0=f
    if(sum(norms<2.2204e-016)>0){
      i0=which(norms<2.2204e-016)
      w=matrix(0,n,1)
      w(i0)=1/length(i0)
    }else{
      w=1/norms # downweight large norm 
      w=w/sum(w)
    }
    norms=sqrt(diag(A)+as.numeric(t(w)%*%A%*%(w))-2*A%*%w)
    f=sum(norms)
    err=abs(f/f0-1)
  }
  
  med=t(w)%*%x
  res=list(med=med,norms=norms,w=w)
  res
}

fte<-function(x,q,alpha,beta){
  
  n=nrow(x)
  p=ncol(x)
  t=seq(0,1,length=p)
  na=ceiling(alpha*n)  # percentage to keep 
  nb=n-floor(n*beta)
  r=matrix(0,n,1)
  for (i in 1:n){
    d=sqrt(apply((x-t(matrix(x[i,],p,n)))^2,1,sum))  # distance to the ith observation 
    d=sort(d)
    r[i]=d[na]  # radius
  }
  
  w=matrix(0,n,1)
  d=sort(r)  # sort all radius 
  rnkn=apply(matrix(rep(1,n),n,1,byrow=T)%*%t(r)<=r%*%matrix(rep(1,n),1,n),2,sum)/n  # rank of the radius/n    
  a=0.5   # tau_1
  b=1-beta # tau_2
  w[rnkn<=a]=1   # hard weight 
  ind=which(rnkn>a&rnkn<=b)
  w[ind]=(rnkn[ind]-b)*((1/(a-b))+(rnkn[ind]-a)*(2*rnkn[ind]-(a+b))/(b-a)^3)
  mu=t(w)%*%x/sum(w)  # weighted average to calculate mean
  xw=diag(sqrt(as.vector(w)))%*%(x-t(matrix(mu,p,n)))
  G=matrix(0,n,n)
  for(i in 1:n){
    for(j in 1:n){
      G[i,j]=crossprod(xw[i,],xw[j,])  # robust covariance function 
    }
  }
  eig=eigen(G,symmetric=T)
  u=eig$vectors
  dd=eig$values
  l=dd[1:q]
  U=u[,1:q]
  pc=matrix(0,q,p)
  pc=t(U%*%diag(1/sqrt(l)))%*%xw
  l=l/sum(w)
  res<-list(mu=mu,l=l,pc=pc)
  res
}

invgen<-function(A){
  ###### From Herve Cardot's thesis.######
  spa <- eigen(A, symmetric = TRUE)
  p <- ncol(A)
  valp <- spa$values
  ind <- (1:p)[valp > 1e-08]
  U <- spa$vectors
  U <- U[, ind]
  if(length(ind) > 1) {
    B <- U %*% diag(1/valp[ind]) %*% t(U)
  }else{
    B <- (1/valp[ind]) * as.matrix(U) %*% t(as.matrix(U))
  }
  B
}

robust_mspe <- function(x){
  
  return(median(x)^2 + mscale(x)^2)
}
# q: number of indices 
# robust mean, covariance, and PCs
# q: number of PCs
# alpha: percentage of points (to get the radius)
# beta: tau_2
RFSIR <- function(x_train, y_train, x_val, y_val, x_test, y_test, k, qs = 1:3, make_prediction){
  
  cr=fte(x_train,50,0.1,0.3)
  cr_mu=as.vector(cr$mu)
  pc=cr$pc
  l=cr$l
  GG=t(pc)%*%diag(l)%*%pc
  
  n=nrow(x_train)
  p=ncol(x_train)
  mx <- t(matrix(cr_mu, p, n))
  x_train_tmp <- (x_train - mx)  # center x 
  sir = matrix(0, p, p)
  h = rep(0, k + 1)
  k1 = floor(n/k)
  a = sort(y_train)
  h[1] = a[1]
  
  for(i in 2:k) {
    h[i] = a[(i - 1) * k1 + 1]  # threshold of y to make slides 
  }
  
  h[k + 1] = a[n] + 0.001
  ai = rep(0, k)
  mu = array(0, dim = c(1, p, k))
  
  for(j in 1:k) { 
    indictor=which(y_train>=h[j]&y_train<h[j+1])
    cr=SpMed(x_train_tmp[indictor,])
    ai[j]=sum(cr$w) 
    mu[,  , j]=cr$med
  }
  
  for(j in 1:k){
    sir = sir + (ai[j]/sum(ai)) * mu[,  , j] %*% t(mu[,  , j])
  }
  
  G=j.0.5(sir) # gamma_e^{1/2}
  gama=GG  # robust covariance 
  spa <- eigen(gama, symmetric = TRUE)
  valp <- spa$values
  U <- spa$vectors
  
  x_val_tmp <- x_val  - t(matrix(cr_mu, p, nrow(x_val)))
  
  err_vals <- rep(NA, length(qs))
  
  for(kk in 1:length(qs)){
    
    print(kk)
    q <- qs[kk]
    U_tmp <- U[, 1:q]
    pk <- U_tmp %*% t(U_tmp) 
    GE <- G %*% invgen(pk %*% gama %*% pk) %*% G
    
    et <- eigen(invgen(GE), symmetric = TRUE)
    val <- et$values[1:q]
    ett <- et$vectors[, 1:q]
    
    DDR<- t(val * t(invgen(pk %*% gama %*% pk) %*% G %*% ett))
    proj <- x_train_tmp%*%DDR
    dat2smooth <- data.frame(proj = x_train_tmp%*%DDR, y = y_train)
    
    
    if(ncol(proj) == 1){
      ff <- as.formula("y~proj")
    }else{
      ff <- as.formula(paste("y~", paste("proj.", 1:q, sep = "", collapse= "+")))
    }
    
    bw.subset <- npregbw(formula =ff,  data =dat2smooth)
    model.np <- npreg(bws = bw.subset)
    
    dat2smooth.val <- data.frame(proj = x_val_tmp%*%DDR)
    pred_val <- predict(model.np, newdata = dat2smooth.val) 
    
    err_vals[kk] <- robust_mspe(pred_val - y_val)
 
  }
 
  J <- which.min(err_vals)
  Q <- qs[J]
  
  # fit the real model 
  U_tmp <- U[, 1:Q]
  pk <- U_tmp %*% t(U_tmp) 
  GE <- G %*% invgen(pk %*% gama %*% pk) %*% G
  
  et <- eigen(invgen(GE), symmetric = TRUE)
  val <- et$values[1:Q]
  ett <- et$vectors[, 1:Q]
  
  DDR<- t(val * t(invgen(pk %*% gama %*% pk) %*% G %*% ett))
  proj <- x_train_tmp%*%DDR
  dat2smooth <- data.frame(proj = x_train_tmp%*%DDR, y = y_train)
  
  
  if(ncol(proj) == 1){
    ff <- as.formula("y~proj")
  }else{
    ff <- as.formula(paste("y~", paste("proj.", 1:Q, sep = "", collapse= "+")))
  }
  
  bw.subset <- npregbw(formula =ff,  data =dat2smooth)
  model.np <- npreg(bws = bw.subset)
  pred_train <- predict(model.np)
  err_train <- mean( (pred_train - y_train)^2)
  
  dat2return <- list(J = J,  err_vals = err_vals, bw=bw.subset$bw, coefs = DDR,  pred_train = pred_train ,err_train = err_train )                                               
  
  if(make_prediction){
    x_test_tmp <- x_test  - t(matrix(cr_mu, p, nrow(x_test)))
    dat2smooth.new <- data.frame(proj = x_test_tmp%*%DDR)
    pred_test <- predict(model.np, newdata = dat2smooth.new) 
    dat2return <- c(dat2return, list(pred_test = pred_test))
    
    if(!missing(y_test)){
      err_test <- mean( (pred_test- y_test)^2)
      dat2return <- c(dat2return, list(err_test = err_test)) 
    }
  }
  
  return(dat2return)
}

