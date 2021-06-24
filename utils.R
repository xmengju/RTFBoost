#' @import rpart stats lme4 splines graphics grid
mse <- function(a,b){
  return(mean((a-b)^2))
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
  res <- r^2
  res[abs(r) >cc] <- 2*cc*abs(r)[abs(r) >cc] - cc^2
  return (res)
}

func.huber.grad <- function(r, cc = 0.98) {
  res <- r
  res[abs(r) > cc] = sign(r)[abs(r) > cc]*cc
  return(res)
}

init.boost <- function(type){
  switch (type,
          "L2" = {
            func <- func.l2
            func.grad <- func.l2.grad
          },
          "LAD" = {
            func <- func.lad
            func.grad <- func.lad.grad
          },
          "RR" = {
            func <- func.tukey
            func.grad <- func.tukey.grad
          }
  )
  
  return (list(func = func, func.grad = func.grad))
}


cal.neggrad <- function(type, y_train,f_train_t, func.grad, init_status, ss, cc){
  
  if(type == "RR" & init_status == 0){
    tmp_numerator <- -func.grad((f_train_t - y_train)/ss, cc = cc)/ss
    tmp_denominator <- -sum(func.grad((f_train_t - y_train)/ss, cc = cc)*(y_train - f_train_t))/(ss^2)
    neg_grad <- tmp_numerator/(tmp_denominator)
  }
  if(type == "RR" & init_status == 1) {
    neg_grad <- - func.grad((f_train_t - y_train)/ss, cc = cc)/ss
  }
  
  if(type == "L2"){
   neg_grad <- -func.grad(f_train_t - y_train)
  }
  return(neg_grad)
}


cal.ss.rr <- function(f_train_t, y_train,  cc, bb) {
  
  ss <- RobStatTM::mscale(f_train_t - y_train,  tuning.chi=cc, delta = bb)
  return(ss)
}

mse <- function(x){
  
  return(mean(x^2))
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


cal.alpha <- function(f_train_t, h_train, y_train, func, type, init_status, ss,  cc){

    if(type == "L2"){
        return(1)
    }
    if(type == "LAD"){
      ff = function(a,r,h){
        return(mean(func(r - a*h)))
      }
      obj_val <- Inf
      min_val <- 0
      for(upper in c(1,5,10)){
        tmp <- optimize(ff, lower = -1, upper = upper, r = y_train - f_train_t, h = h_train)
        if(tmp$objective < obj_val){
          obj_val <- tmp$objective
          min_val <- tmp$minimum
        }
      }
      return(min_val)
    }
  
  if(type == "RR" & init_status == 0) {
    ff3 <- function(a, r, h) return(RobStatTM::mscale(r - a*h))
    upper_region = c(0.5,10,100,300)
    tmp <-  tmp_val <- rep(NA, length(upper_region))
    for(i in 1:length(upper_region)){
      val = optimize(ff3, lower = -1, upper = upper_region[i], r = y_train - f_train_t, h = h_train)
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
  if(type == "RR" & init_status == 1) {
    ff4 <- function(a, r, h, c, s) return(mean(func( (r - a*h)/s,  c)))
    upper_region = c(0.5,10,100,300)
    tmp <- rep(NA, length(upper_region))
    tmp <-  tmp_val <- rep(NA, length(upper_region))
    for(i in 1:length(upper_region)){
      val = optimize(ff4, lower = -1, upper = upper_region[i], r = y_train - f_train_t, h = h_train, c = cc, s = ss)
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

cal_error <- function(trim_prop, trim_c, err_type, f_t_test, y_test){
  
  if(err_type == "mse"){
    return(mse(f_t_test - y_test))
  }
  
  if(err_type == "tmse"){
    return(tmse(trim_prop, trim_c, f_t_test - y_test)$trmse)
  }
  
  if(err_type == "aad"){
    return(mean(abs(f_t_test - y_test)))
  }
  
}


# transform a basis matrix to an orthonormal basis matrix
compute.orthonormal <- function(B, grid, t_range){
  
  d <- ncol(B)
  Phi_i <- B
  Psi <- matrix(NA, nrow = nrow(B), ncol(B))
  Psi[,1] <-   B[,1]/ sqrt(riemman(B[,1]*B[,1], grid, t_range))
  
  for(i in 2:d){
    for(j in 1:(i-1)){
      Phi_i[,i] <-   Phi_i[,i]  -  riemman(Phi_i[,i]*Psi[,j], grid, t_range) * Psi[,j]
    }
    Psi[,i] <-   Phi_i[,i]/ sqrt(riemman(Phi_i[,i]*Phi_i[,i], grid, t_range))
  }
  return(Psi)
}

