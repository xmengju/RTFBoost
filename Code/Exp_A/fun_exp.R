set.control <- function(d, precision,  shrinkage, nknot, n_init, niter){
  
  control.l2 <-  RTFBoost.control(make_prediction = TRUE, 
                                  tree_control = TREE.control(tree_type  = "B", d = d, num_dir = 200),
                                  type = "L2", shrinkage  = shrinkage, precision = precision, 
                                  init_type = "mean", niter = niter, 
                                  nknot = nknot, save_f = FALSE, trace =  TRUE, save_tree = FALSE, error_type = c("mse"))
  
  control.lad <-  RTFBoost.control(make_prediction = TRUE, 
                                   tree_control = TREE.control(tree_type  = "B", d = d, num_dir = 200),
                                   type = "LAD", shrinkage  = shrinkage, precision = precision, 
                                   init_type = "median", niter = niter, 
                                   nknot = nknot, save_f = FALSE, trace = TRUE, save_tree = FALSE, error_type = c("mse"))
  
  control.rr <- RTFBoost.control(make_prediction = TRUE, eff_m= 0.95, bb = 0.5, 
                                 tree_control = TREE.control(tree_type  = "B", d = d, num_dir = 200),
                                 type = "RR", shrinkage  = shrinkage, precision = precision, 
                                 init_type = "LADTree", n_init = n_init, niter = niter, 
                                 nknot = nknot, save_f = FALSE, trace =  TRUE, save_tree = FALSE, error_type = c("mse"))
  
  return(list(control.l2 = control.l2, control.lad = control.lad, control.rr = control.rr))
  
}


exp.tree.methods<- function(x_train, y_train, x_val, y_val, x_test, y_test, niter, grid, t_range, control.tree.list){
  
  model.tree.list <- list()
  err_trains <- err_tests <- err_vals <-  matrix(NA, length(control.tree.list), niter)
  res_trains <- res_tests <- res_vals <- list()
  err_test <- rep(NA, length(control.tree.list))
  err_val <- rep(NA, length(control.tree.list))
  
  early_stop <- rep(NA, length(control.tree.list))
  time_vec <- rep(NA, length(control.tree.list))
  
  for(i in 1:length(control.tree.list)){
    
    print(paste(i, "out of", length(control.tree.list), "shrinkage", control.tree.list[[i]]$shrinkage, 
                "d",control.tree.list[[i]]$tree_control$d, "dd", control.tree.list[[i]]$tree_control$dd,
                "nknots", control.tree.list[[i]]$nknot))
    
    if(i <= 2){
       tt <- system.time(model.tree.list[[i]] <- RTFBoost(x_train = x_train, y_train = y_train, x_val = x_val,  y_val = y_val,
                                                       x_test = x_test, y_test = y_test, grid = grid, t_range  = t_range, 
                                                       control = control.tree.list[[i]]))
    }else{
       tt <- system.time(model.tree.list[[i]]  <- RTFBoost.validation(x_train = x_train, z_train = NULL, y_train = y_train,  x_val = x_val,  z_val = NULL, 
                                                    y_val = y_val, x_test = x_test, z_test = NULL, y_test = y_test, grid = grid,  t_range = t_range, 
                                               max_depth_init_set = c(1,2,3,4), min_leaf_size_init_set = c(10,20,30), control = control.tree.list[[i]]))
      when_init <- model.tree.list[[i]]$when_init
      err_cvs <-  model.tree.list[[i]]$err_cvs 
      params <-  model.tree.list[[i]]$params
    }
    
    time_vec[i] <- tt[3]
    err_trains[i,1:length(model.tree.list[[i]]$err_train)] <-model.tree.list[[i]]$err_train
    
    err_vals[i,1:length(model.tree.list[[i]]$err_val)] <-model.tree.list[[i]]$err_val
    err_val[i] <- model.tree.list[[i]]$err_val[model.tree.list[[i]]$early_stop]
    
    err_tests[i,1:model.tree.list[[i]]$early_stop] <-model.tree.list[[i]]$err_test[,1]
    err_test[i] <- model.tree.list[[i]]$err_test[model.tree.list[[i]]$early_stop,]

    early_stop[i] <- model.tree.list[[i]]$early_stop
    
    res_trains[[i]] <- model.tree.list[[i]]$f_train_t - y_train
    res_vals[[i]] <- model.tree.list[[i]]$f_val_t - y_val
    res_tests[[i]] <- model.tree.list[[i]]$f_test_t - y_test
    model.tree.list[[i]] <- NULL
    
  }
  return(list(params =   params, when_init = when_init,  err_cvs =  err_cvs , res_trains = res_trains, res_vals = res_vals, res_tests = res_tests, 
              time_vec  = time_vec, err_trains = err_trains, err_tests = err_tests, err_vals = err_vals, err_test = err_test, err_val = err_val, early_stop = early_stop))
}


do.exp <- function(seed, g_func_no, SNR, x_type, dd, control.tree.list, case_id){
  
  sd_g_list <- readRDS("Code/Exp_A/sd_g_list.rds")
  dat_gen_control <- dat.generate.control(x_type = x_type, SNR = SNR, g_func_no = g_func_no, n_train = 400, n_val = 200, n_test = 1000, sd_g_list = sd_g_list)
  dat <- dat.generate(seed = seed, control = dat_gen_control)
  
  x_train <- dat$x$x_train;  x_val <- dat$x$x_val;  x_test <- dat$x$x_test
  n_train <- nrow(x_train);  n_val <- nrow(x_val);  n_test <- nrow(x_test)
  
  g_train <- dat$g$g_train;  g_val <- dat$g$g_val;  g_test <- dat$g$g_test
  S <- dat$S
  e_train <- rnorm(n_train); e_val <- rnorm(n_val);  e_test <- rnorm(n_test)
  
  set.seed(seed)
  switch(as.character(dd), 
         "D0" = {
           out_train <- out_val <- NULL
         },
         "D1" = {
           out_train <- sample(n_train,size= 0.2*(n_train))
           out_val <- sample(n_val,size= 0.2*(n_val))
           
           tmp_train <- rbinom(length(out_train),1,0.5)
           tmp_val <- rbinom(length(out_val),1,0.5)
  
           e_train[out_train] <- tmp_train* rnorm(length(out_train), mean= 20, sd = 0.1) + (1-tmp_train)* rnorm(length(out_train), mean= -20, sd = 0.1)
           e_val[out_val] <-  tmp_val* rnorm(length(out_val), mean= 20, sd = 0.1) + (1-tmp_val)* rnorm(length(out_val), mean= -20, sd = 0.1)
         },
         "D2" = {
           out_train <- sample(n_train,size= 0.2*(n_train))
           out_val <- sample(n_val,size= 0.2*(n_val))
           e_train[out_train] <- rnorm(length(out_train), mean= 20, sd = 0.1) 
           e_val[out_val] <- rnorm(length(out_val), mean= 20, sd = 0.1) 
      }
    )
  
  y_train <-  g_train + S*e_train 
  y_val <-  g_val + S*e_val 
  y_test <-  g_test + S*e_test 
  
  grid <- dat$tt
  if(x_type == "ferraty"){
    t_range <- c(-1,1)
  }
  if(x_type == "mattern"){
    t_range <- c(0,1)
  }
  
  dat2return <- NULL
  if(case_id == 1){
    
    dat2return <- exp.tree.methods(x_train, y_train, x_val, y_val, x_test, y_test, niter,  grid, t_range, control.tree.list)
  }
  if(case_id == 0){
    u <- 1:length(y_train)
    range_beta <- 4:15
    range_eta <- 20
    norder  <- 4
    model_RobustFPLM_select <- FPLMBsplines(y = y_train, x = x_train, u = u, t = grid, range_freq = range_beta, range_spl = range_eta, 
                                       norder = norder, fLoss = "lmrob", trace = TRUE)
    
    range_beta <- 4:7
    model_RobustFPLM_fix <- FPLMBsplines(y = y_train, x = x_train, u = u, t = grid, range_freq = range_beta, range_spl = range_eta, 
                                            norder = norder, fLoss = "lmrob", trace = TRUE)
    
    dat2return <- list(select = FPLMBsplines.predict(model_RobustFPLM_select, newx = x_test, newy = y_test),
                       fix = FPLMBsplines.predict(model_RobustFPLM_fix, newx = x_test, newy = y_test))
  }

  dat2return$S <- dat$S
  dat2return$out_train <- out_train
  dat2return$out_val <- out_val
  return(dat2return)
}