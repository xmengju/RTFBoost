dat.generate.control <- function(n_train = 400,  n_val = 200, n_test = 1000, g_func_no = 1, SNR = 5){

  rr <- c(0,1)
  mattern_tmp <- readRDS("Code/Exp_A/mattern.rds")
  mu <-  mattern_tmp$mu
  phis <- mattern_tmp$phis
  K <- length(phis)
  
  switch(g_func_no,
         "1" = {
           g.func <- function(x,t){
             (sin((3*pi)*t/2) + sin(0.5*pi*t))*x
           }
           int_func <- function(x, t_eval){riemman(g.func(x, t_eval), t_eval, rr =rr)}
         },
         
         "2" = {
           g.func <- function(x,t){
             (x -mu(t)) *(phis[[1]](t)+ phis[[2]](t)) 
           }
           
           int_func <- function(x, t_eval){
             tmp <- riemman(g.func(x, t_eval), t_eval, rr = rr)
             sign(tmp)*(abs(tmp))^{1/3}
           }
         },
         "3" = {
           g.func <- function(x,t){
             x*log(abs(x))/2
           }
           int_func <- function(x, t_eval){ 5*exp(-abs(riemman(g.func(x, t_eval), t_eval, rr = rr)))}
         },
         "4" = {
           g.func <- function(x,t){
             x^2*sin(2*pi*t)
           }
           int_func <- function(x, t_eval){
             tmp <- riemman(g.func(x, t_eval), t_eval, rr = rr)
             return(5/(1+exp(-tmp)))
           }
         }, 
         "5" = {
           g.func.1 <-  function(x,t){
             cos(2*pi*t^2)*x
           }
           g.func.2 <-  function(x,t){
             sin(x)
           }
           threshold <- 0.5             
           rr_1 <- c(0,0.5)
           rr_2 <- c(0.5,1)
           
           int_func <- function(x, t_eval){5*(sqrt(abs(riemman(g.func.1(x[t_eval<threshold], t_eval[t_eval<threshold]), t_eval[t_eval<threshold], rr = rr_1))) +  
                                                + sqrt(abs(riemman(g.func.2(x[t_eval>=threshold], t_eval[t_eval>=threshold]), t_eval[t_eval>=threshold], rr = rr_2))))}
         }
  )
  
  return(list(n_train = n_train, n_val = n_val, n_test = n_test, g_func_no = g_func_no,  int_func = int_func,  SNR = SNR,
              mu = mu, phis = phis, K = K))
}


dat.generate <- function(seed, control = dat.generate.control(), type = "C0", direction = "symmetric"){
  
  n_train <- control$n_train;  n_val <- control$n_val;  n_test <- control$n_test; 
  n <- n_train + n_val + n_test
  
  g_func_no <- control$g_func_no
  int_func <- control$int_func
  SNR <- control$SNR
  
  # read values from saved files 
  sd_list <- readRDS("Code/Exp_A/sd_list.rds")
  sd_g <- sd_list[[control$g_func_no]]
  mattern_tmp <- readRDS("Code/Exp_A/mattern.rds")
  mu <-  mattern_tmp$mu
  phis <- mattern_tmp$phis
  K <- length(phis)
  
  set.seed(seed)
  
  t_eval <- seq(0,1, length.out = 1000)
  xi_matrix <- matrix(rnorm(n*4), ncol = 4)
  x_func <- generate_mattern(mu, phis, K)(xi_matrix, t_eval)

  train_idx <-1:n_train
  val_idx <-  (n_train+1): (n_train + n_val)
  test_idx <- (n_train + n_val+1):n
  
  
  tt <- seq(1, 1000, by = 10)
  tt_idx_t_eval <- (1:length(t_eval))[tt]
  
  if(type!= "C0"){
    outliers_train <- sample(1:n_train, round(0.1*n_train))
    outliers_val <-  sample(1:n_val, round(0.1*n_val))
  }
  
  if(type == "C1"){
    xi_matrix_out_train <- matrix(rnorm(length(outliers_train)*4), ncol = 4)
    xi_matrix_out_train[,2] <- rnorm(length(outliers_train), mean = 2, sd =0.25)
    
    new_mu <- function (tt) {
      1+ 2*cos(tt *2* pi) 
    }
    
    x_func_out_train <- generate_mattern(new_mu, control$phis, control$K)(xi_matrix_out_train, t_eval)
    x_func[outliers_train, ] <- x_func_out_train
    
    xi_matrix_out_val <- matrix(rnorm(length(outliers_val)*4), ncol = 4)
    xi_matrix_out_val[,2] <- rnorm(length(outliers_val), mean = 2, sd =0.25)
    
    x_func_out_val <- generate_mattern(control$mu, control$phis, control$K)(xi_matrix_out_val, t_eval)
    x_func[n_train+outliers_val,] <- x_func_out_val
  }
  if(type == "C2"){
    x_func[outliers_train, ] <- 2*x_func[outliers_train,] + 5
    x_func[n_train+outliers_val, ] <- 2*x_func[n_train+outliers_val,]  + 5
  }
  if(type == "C3"){

    outlier_num <- 10
    outlier_candidate <- 1:100
    
    t_out_train <- matrix(unlist(lapply(1:length(outliers_train), FUN= function(x){set.seed(x); sample(outlier_candidate, outlier_num)})), ncol = outlier_num, byrow = T)
    t_out_val <- matrix(unlist(lapply(1:length(outliers_val), FUN= function(x){set.seed(x); sample( outlier_candidate,  outlier_num)})), ncol = outlier_num, byrow = T)
    
    for(oo in 1:length(outliers_train)){

      tmp <- rbinom(outlier_num, 1,0.5)
      x_func[outliers_train[oo],  tt_idx_t_eval[t_out_train[oo,]]] <-  x_func[outliers_train[oo],   tt_idx_t_eval[t_out_train[oo,]]] + 
        tmp*rnorm(outlier_num, mean = 10, sd = 0.25) + (1-tmp)*rnorm( outlier_num, mean = -10, sd = 0.25)
    }
    
    for(oo in 1:length(outliers_val)){
      tmp <- rbinom(outlier_num, 1,0.5)
      x_func[n_train + outliers_val[oo],   tt_idx_t_eval[t_out_val[oo,]]] <- x_func[n_train+outliers_val[oo],   tt_idx_t_eval[t_out_val[oo,]]] + 
        tmp*rnorm( outlier_num, mean = 10, sd = 0.25) + (1-tmp)*rnorm( outlier_num, mean = -10, sd = 0.25)
    }
  }
  if(type == "C4"){
    
    interval_out_train <- unlist(lapply(1:length(outliers_train), FUN= function(x){set.seed(x); sample(seq(1, 100, by = 10), 1)}))
    interval_out_val <- unlist(lapply(1:length(outliers_val), FUN= function(x){set.seed(x); sample(seq(1, 100, by = 10), 1)}))
    
    for(oo in 1:length(interval_out_train)){
      x_func[outliers_train[oo],  tt_idx_t_eval[interval_out_train[oo]:(interval_out_train[oo]+9)]] <- 
        x_func[outliers_train[oo],  tt_idx_t_eval[interval_out_train[oo]:(interval_out_train[oo]+9)]] + rnorm(10, mean = 10, sd = 0.25) 
    }
    for(oo in 1:length(interval_out_val)){
      x_func[n_train + outliers_val[oo],  tt_idx_t_eval[interval_out_val[oo]:(interval_out_val[oo]+9)]] <- 
        x_func[n_train + outliers_val[oo],  tt_idx_t_eval[interval_out_val[oo]:(interval_out_val[oo]+9)]] + rnorm(10, mean = 10, sd = 0.25) 
    }
  }
  
  g <- unlist(lapply(1:n,  FUN = function(i) {int_func(x_func[i,], t_eval)}))
  S <- sd_g/sqrt(SNR)
  y <- g + S*rnorm(n)
  print(paste("constant is ", round(S, 3)))
  

  x_train <- x_func[train_idx,tt]
  x_val <- x_func[val_idx,tt]
  x_test <- x_func[test_idx,tt]
  
  y_train <- y[train_idx]
  y_val <- y[val_idx]
  y_test <- y[test_idx]
  
  dat <- list()
  dat$x <- list(x_train = x_train, x_val = x_val, x_test = x_test)
  dat$y <- list(y_train = y_train, y_val = y_val, y_test = y_test)
  dat$g <- list(g_train = g[train_idx], g_val = g[val_idx], g_test = g[test_idx])
  
  dat$t_eval <- t_eval
  dat$tt <- t_eval[tt]
  dat$S <- S
  
  dat2return <- dat
  
  if(type == "C0"){
  }else{

    dat2return$outliers_train <-  outliers_train 
    dat2return$outliers_val <-  outliers_val
    
    exm_val <- 30
    if(type %in% c("C1","C3","C4","C5", "C2","C6")){
      if(direction == "asymmetric"){
        dat2return$y$y_train[outliers_train] <- dat2return$g$g_train[outliers_train] + rnorm(length(outliers_train), -exm_val, 0.25)
        dat2return$y$y_val[outliers_val] <- dat2return$g$g_val[outliers_val] + rnorm(length(outliers_val), -exm_val, 0.25)
      }else{
        tmp_train <- rbinom(length(outliers_train), 1,0.5)
        dat2return$y$y_train[outliers_train] <- dat2return$g$g_train[outliers_train] + tmp_train*rnorm(length(outliers_train), exm_val, 0.25)+
          (1-tmp_train)*rnorm(length(outliers_train), - exm_val, 0.25)
        tmp_val <- rbinom(length(outliers_val), 1,0.5)
        dat2return$y$y_val[outliers_val] <- dat2return$g$g_val[outliers_val] + tmp_val*rnorm(length(outliers_val), exm_val, 0.25)+
          (1-tmp_val)*rnorm(length(outliers_val), -  exm_val, 0.25) 
      }
    }
  }
  
  dat2return$y$y_train <- dat2return$y$y_train
  dat2return$y$y_val <- dat2return$y$y_val
  
  return(dat2return)
}


cal_sd <- function(g_func_no = 1){
  
  set.seed(123)
  mattern_tmp <- readRDS("Code/Exp_A/mattern.rds")
  mu <-  mattern_tmp$mu
  phis <- mattern_tmp$phis
  K <- length(phis)
  control <- dat.generate.control(g_func_no = g_func_no)
  n <- 3000
  t_eval <- seq(0,1, length.out = 1000)
  xi_matrix <- matrix(rnorm(n*4), ncol = 4)
  x_func <- generate_mattern(mu, phis, K)(xi_matrix, t_eval)
  rr <- c(0,1)
  g <- unlist(lapply(1:n,  FUN = function(i) {control$int_func(x_func[i,], t_eval)}))
  return(sd(g))
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

generate_mattern <- function(mu, phis, K){
  
  # eigenvalues of the covariance function
  lambdas <- c(0.8, 0.3, 0.2, 0.1)
  # eigenfunctions
  
  fun <- function(xi_matrix,tt){
    xx <- matrix(1,nrow = nrow(xi_matrix), ncol = 1)%*%mu(tt)
    for(j in 1:K){
      xx <- xx + sqrt(lambdas[j]) *as.matrix(xi_matrix[,j])%*%(phis[[j]](tt))
    }
    return(xx)
  }
  return(fun)
}



#sd_list  <- list()
#for(g_func_no in c(1:5)){
#  sd_list[[g_func_no]] <- cal_sd(g_func_no = g_func_no)
#}
# 
#saveRDS(sd_list, file = "Code/Exp_A/sd_list.rds")










