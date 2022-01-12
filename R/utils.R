#' @import rpart stats lme4 splines graphics grid
mse <- function(x){
  
  return(mean(x^2))
}

aad <- function(x){
  return(mean(abs(x)))
}

tmse <- function(trim_prop = NULL, trim_c = NULL, x){
  if(length(trim_c) != 0){
    idx <- (x < (median(x) + trim_c*mad(x))) & x > (median(x) - trim_c*mad(x))
  }else{
    if(length(trim_prop) != 0){
      idx <-  (abs(x) < quantile(abs(x), 1-trim_prop))
    }
  }
  return(list(tmse = mse(x[idx]), idx = idx))
}


func.l2 <- function(x, cc = NULL) {
  return((x)^2/2)
}

func.l2.grad<- function(x, cc = NULL) {
  return(x)
}

func.lad <- function(x, cc = NULL) {
  return(abs(x))
}

func.lad.grad<- function(x, cc = NULL) {
  return(sign(x))
}


func.tukey <- function(r, cc= 4.685) {
  w <- as.numeric(abs(r) <= cc)
  v <- w*(1 - (1 - (r/cc)^2)^3)  +(1-w)*1
  return(v)
}


func.tukey.grad <- function(r, cc = 4.685) {
  w <- as.numeric(abs(r) <= cc )
  gg <- w*6*r*(1 - (r/cc)^2)^2/(cc^2)  +(1-w)*0
  return(gg)
}

func.huber <- function(r, cc= 0.98) {
  res <- (r^2)/2
  res[abs(r) >cc] <- cc*abs(r)[abs(r) >cc] - (cc^2)/2
  return (res)
}



func.huber.grad <- function(r, cc = 0.98) {
  res <- r
  res[abs(r) > cc] = sign(r)[abs(r) > cc]*cc
  return(res)
}

func.huber.grad.prime <- function(r, cc = 0.98) {
  return(0 + (abs(r) <=cc)*1)
}

cal.efficiency <- function(e, psi,psi.prime) {
   tmp_1 <- integrate(function(a, cc) (psi(a, cc)^2)*dnorm(a), cc=e, lower=-Inf, upper=+Inf)$value
   tmp_2 <- integrate(function(a, cc) psi.prime(a, cc)*dnorm(a), cc=e, lower=-Inf, upper=+Inf)$value
   return( 1/(tmp_1/tmp_2^2) )
 }

cal.efficiency <- function(e, psi,psi.prime) {
  tmp_1 <- integrate(function(a, cc) (psi(a, cc)^2)*dnorm(a), cc=e, lower=-Inf, upper=+Inf)$value
  tmp_2 <- integrate(function(a, cc) psi.prime(a, cc)*dnorm(a), cc=e, lower=-Inf, upper=+Inf)$value
  
  return( 1/(tmp_1/tmp_2^2) )
}


init.boosting <- function(type, init.type = "mean"){
  switch (type,
          "L2" = {
            func <- func.l2
            func.grad <- func.l2.grad
          },
          "LAD" = {
            func <- func.lad
            func.grad <- func.lad.grad
          },
          "LAD-M" = {
            func <- func.lad
            func.grad <- func.lad.grad
          },
          "RR" = {
            func <- func.tukey
            func.grad <- func.tukey.grad
          }
  )
  
  switch(init.type,
          mean = {
            init.func <- function(x){
              mean(x)
            }
          },
          median = {
            init.func <- function(x){
              median(x)
            }
          },
         LADTree = {
           init.func <- function(x){
           }
         }
  )
  return (list(func = func, func.grad = func.grad, init.func = init.func))
}



cal.neggrad <- function(type, y.train, f.train, func.grad, init.status, ss, cc){
  
  if(type == "RR" & init.status == 0){
    tmp_numerator <- -func.grad((f.train - y.train)/ss, cc = cc)/ss
    tmp_denominator <- -sum(func.grad((f.train - y.train)/ss, cc = cc)*(y.train - f.train))/(ss^2)
    neg_grad <- tmp_numerator/(tmp_denominator)
  }
  if(type == "RR" & init.status == 1) {
    neg_grad <- - func.grad((f.train - y.train)/ss, cc = cc)/ss
  }
  if(type == "LAD-M" & init.status == 1) {
    neg_grad <- - func.grad((f.train - y.train)/ss, cc = cc)/ss
  }
  if(type == "LAD-M" & init.status == 0) {
    neg_grad <- -func.grad(f.train - y.train)
  }
  if(type %in% c("L2","LAD")){
    neg_grad <- -func.grad(f.train - y.train)
  }
  return(neg_grad)
}


cal.ss.rr <- function(f.train, y.train,  cc, bb) {
  
  ss <- RobStatTM::mscale(f.train - y.train, delta = bb)
  return(ss)
}


cal.alpha <- function(f.train, h.train, y.train, func, type, init.status, ss, bb, cc){
  
  if(type == "L2"){
    return(1)
  }
  
  if( (type == "LAD") ||  (type == "LAD-M" & init.status == 0) ){
    ff = function(a,r,h){
      return(mean(func(r - a*h)))
    }
    obj_val <- Inf
    min_val <- 0
    for(upper in c(1,5,10)){
      tmp <- optimize(ff, lower = -1, upper = upper, r = y.train - f.train, h = h.train)
      if(tmp$objective < obj_val){
        obj_val <- tmp$objective
        min_val <- tmp$minimum
      }
    }
    return(min_val)
  }
  
  if(type == "RR" & init.status == 0) {
    ff3 <- function(a, r, h) return(RobStatTM::mscale(r - a*h, delta = bb))
    upper_region = c(0.5,10,100,300,500,800,1000,3000,5000)
    tmp <-  tmp_val <- rep(NA, length(upper_region))
    for(i in 1:length(upper_region)){
      val = optimize(ff3, lower = -1, upper = upper_region[i], r = y.train - f.train, h = h.train)
      tmp[i] <- val$minimum
      tmp_val[i] <- val$objective
    }
    idx <- min(which(tmp_val == min(tmp_val)))
    order_val <- order(tmp_val[1:idx])
    
    if( sum(order_val == idx:1) == idx){
      return(tmp[idx])
    }else{
      tmp_order <- order_val  - c(max(order_val), order_val[1: (length(order_val)-1)])
      if(sum(tmp_order > 0) > 0){
        tmp_idx <- min(which(tmp_order>0))-1
        return(tmp[tmp_idx])
      }else{
        return(tmp[1])
      }
    }
  }
  if( (type == "RR" & init.status == 1) ||(type == "LAD-M" & init.status == 1) ){
    
    ff4 <- function(a, r, h, c, s) return(mean(func( (r - a*h)/s,  c)))
    #upper_region = c(0.01,0.1,0.5,10,100,200,300,500,800,3000,5000)
    #upper_region = c(0.01,0.1,0.5,1,5)
    upper_region = c(0.01,0.1,0.5,5,10,100,200,300)
    
    tmp <- rep(NA, length(upper_region))
    tmp <-  tmp_val <- rep(NA, length(upper_region))
    for(i in 1:length(upper_region)){
      val = optimize(ff4, lower = -0.01, upper = upper_region[i], r = y.train - f.train, h = h.train, c = cc, s = ss)
      tmp[i] <- val$minimum
      tmp_val[i] <- val$objective
    }
    
    idx <- min(which(tmp_val == min(tmp_val)))
    order_val <- order(tmp_val[1:idx])
    
    if(sum(order_val == idx:1) == idx){ #continue going down
      return(tmp[idx])
    }else{
      tmp_order <- order_val  - c(max(order_val), order_val[1:(length(order_val)-1)])
      if(sum(tmp_order > 0) > 0){
        tmp_idx <- min(which(tmp_order>0))-1
        return(tmp[tmp_idx])
      }else{
        return(tmp[1])
      }
    }
  }
  
}


riemman <- function(uu, tt, rr = c(-1,1)){
  tmp <- sum((uu[-1] + uu[-length(uu)])/2 *(diff(tt)))
  if(rr[1] < tt[1]){
    tmp <- tmp + (tt[1] - rr[1])*uu[1]
  }
  if(rr[2]>tt[length(tt)]){
    tmp <- tmp + (rr[2] - tt[length(tt)])*uu[length(uu)]
  }
  return(tmp)
}
