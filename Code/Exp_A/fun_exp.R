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


exp.tree.methods<- function(x_train, y_train, x_val, y_val, x_test, y_test, niter, grid, t_range, control.tree.list, validation_tree){
  
  model.tree.list <- list()
  err_trains <- err_tests <- err_vals <-  matrix(NA, length(control.tree.list), niter)
  res_trains <- res_tests <- res_vals <- list()
  err_test <- rep(NA, length(control.tree.list))
  err_val <- rep(NA, length(control.tree.list))
  
  early_stop <- rep(NA, length(control.tree.list))
  time_vec <- rep(NA, length(control.tree.list))
  
  err_cvs <- NULL;  params <- NULL
  
  for(i in 1:length(control.tree.list)){
    
    print(paste(i, "out of", length(control.tree.list), "shrinkage", control.tree.list[[i]]$shrinkage, 
                "d",control.tree.list[[i]]$tree_control$d, "nknots", control.tree.list[[i]]$nknot))
    
    if(i <= 2){
       model.tree.list[[i]] <- RTFBoost(x_train = x_train, y_train = y_train, x_val = x_val,  y_val = y_val,
                                                       x_test = x_test, y_test = y_test, grid = grid, t_range  = t_range, 
                                                       control = control.tree.list[[i]])
    }else{
      if(validation_tree){
       model.tree.list[[i]]  <- RTFBoost.validation(x_train = x_train, z_train = NULL, y_train = y_train,  x_val = x_val,  z_val = NULL, 
                                                    y_val = y_val, x_test = x_test, z_test = NULL, y_test = y_test, grid = grid,  t_range = t_range, 
                                               max_depth_init_set = c(1,2,3,4), min_leaf_size_init_set = c(10,20,30), control = control.tree.list[[i]])
       err_cvs <-  model.tree.list[[i]]$err_cvs 
       params <-  model.tree.list[[i]]$params
       
      }else{
        tmp_control <-  control.tree.list[[i]]
        tmp_control$init_type <- "median"
        model.tree.list[[i]]  <- RTFBoost(x_train = x_train, y_train = y_train, x_val = x_val,  y_val = y_val,
                 x_test = x_test, y_test = y_test, grid = grid, t_range  = t_range, 
                 control =   tmp_control)
      }
      when_init <- model.tree.list[[i]]$when_init

    }
    
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
            err_trains = err_trains, err_tests = err_tests, err_vals = err_vals, err_test = err_test, err_val = err_val, early_stop = early_stop))
}


### choosing from c("FPPR", "FGAM", "MFLM", "RFSIR", "RFPLM", and "RTFBoost")

do.exp <- function(seed, g_func_no, SNR, d, control.tree.list, methods = c("RTFBoost"), type, validation_tree) {
  
  dat_gen_control <- dat.generate.control(SNR = SNR, g_func_no = g_func_no, n_train = 400, n_val = 200, n_test = 1000)
  dat <- dat.generate(seed = seed, control = dat_gen_control, type = type)
  
  x_train <- dat$x$x_train;  x_val <- dat$x$x_val;  x_test <- dat$x$x_test
  y_train <- dat$y$y_train;  y_val <- dat$y$y_val;  y_test <- dat$y$y_test
  
  n_train <- nrow(x_train);  n_val <- nrow(x_val);  n_test <- nrow(x_test)
  S <- dat$S
  grid <- dat$tt
  t_range <- c(0,1)
  
  set.seed(seed)
  
  
  dat2return <- list()
  
  if("RTFBoost" %in% methods){
    dat2return <- c(dat2return, list(RTFBoost = exp.tree.methods(x_train, y_train, x_val, y_val, x_test, y_test, niter,  grid, t_range, control.tree.list, validation_tree)))
  }
  
  if("RFPLM" %in% methods){
    u <- 1:length(y_train)
    range_beta <- 4:7
    range_eta <- 20
    norder  <- 4
    model_RobustFPLM <- FPLMBsplines(y = y_train, x = x_train, u = u, t = grid, range_freq = range_beta, range_spl = range_eta, 
                                     norder = norder, fLoss = "lmrob", trace = TRUE)
    
    dat2return <- c(dat2return, list(RFPLM = FPLMBsplines.predict(model_RobustFPLM, newx = x_test, newy = y_test)))
  }
  
  if("MFLM" %in% methods){
    u <- 1:length(y_train)
    range_beta <- 4:5
    range_eta <- 20
    norder  <- 4
    model_MFLM <- FPLMBsplines(y = y_train, x = x_train, u = u, t = grid, range_freq = range_beta, range_spl = range_eta, 
                               norder = norder, fLoss = "huang", trace = TRUE)
    
    dat2return <- c(dat2return, list(MFLM = FPLMBsplines.predict(model_MFLM, newx = x_test, newy = y_test))) 
  }
  
  if("FGAM" %in% methods){
    
    xx <- x_train
    yy <- y_train
    nbasises_FGAM <- 15
    model_FGAM <- fgam(yy ~ af(xx, splinepars=list(k=c(nbasises_FGAM,nbasises_FGAM),m=list(c(2,2),c(2,2)))), gamma = 1.2, method="REML")
    
    err_test_FGAM <-  mean((predict(model_FGAM ,newdata=list(xx =x_test),type='response')- y_test)^2)
    err_train_FGAM <-   mean((predict(model_FGAM ,newdata=list(xx =x_train),type='response')-y_train)^2)
    
    FGAM <- list (err_test  = err_test_FGAM, err_train =  err_train_FGAM)
    
    dat2return <- c(dat2return, list(FGAM =FGAM))
  }
  
  if("FPPR" %in% methods){
    niter_FPPR <- 15
    nknots_FPPR <- 3
    
    model_FPPR <- FPPR(x_train, y_train, x_val, y_val, x_test, y_test, grid = grid, t_range = t_range,  niter = niter_FPPR, nknots = nknots_FPPR)
    err_trains_FPPR <- model.FPPR$err_train
    err_tests_FPPR <- model.FPPR$err_test
    err_vals_FPPR <-  model.FPPR$err_val
    err_val_FPPR  <- model.FPPR$err_val[model.FPPR$early_stop]
    err_test_FPPR <- model.FPPR$err_test[model.FPPR$early_stop]
    Js_FPPR <- model.FPPR$J
    early_stop_FPPR <- model.FPPR$early_stop
    FPPR_obj <- list( err_trains_FPPR =  err_trains_FPPR,  err_tests_FPPR =  err_tests_FPPR,
                      err_vals_FPPR = err_vals_FPPR,  err_val_FPPR = err_val_FPPR, err_test_FPPR = err_test_FPPR,
                      Js_FPPR =  Js_FPPR ,  early_stop_FPPR =  early_stop_FPPR)
    dat2return <- c(dat2return, list(FPPR = FPPR_obj))
    
  }
  
  dat2return$dat <- dat
  
  
  return(dat2return)
}

