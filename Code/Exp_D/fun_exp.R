set.control <- function(d, precision,shrinkage, nknot, n_init, niter){
  
  control.l2 <-  RTFBoost.control(make_prediction = TRUE, 
                                  tree_control = TREE.control(tree_type  = "B", d = d, num_dir = 200),
                                  type = "L2", shrinkage  = shrinkage, precision = precision, 
                                  init_type = "mean", niter = niter, 
                                  nknot = nknot, save_f = FALSE, trace =  TRUE, save_tree = FALSE, error_type = c("mse"))
  
  control.lad <-  RTFBoost.control(make_prediction = TRUE, 
                                   tree_control = TREE.control(tree_type  = "B", d = d, num_dir = 200),
                                   type = "LAD", shrinkage  = shrinkage, precision = precision, 
                                   init_type = "median", niter = niter, 
                                   nknot = nknot, save_f = FALSE, trace =  TRUE, save_tree = FALSE, error_type = c("mse"))
  
  
  control.lad.m <-  RTFBoost.control(make_prediction = TRUE, 
                                   tree_control = TREE.control(tree_type  = "B", d = d, num_dir = 200),
                                   type = "LAD-M", shrinkage  = shrinkage, precision = precision, 
                                   init_type = "median",  n_init = n_init, niter = niter, 
                                   nknot = nknot, save_f = FALSE, trace = TRUE, save_tree = FALSE, error_type = c("mse"))
  
  control.rr <- RTFBoost.control(make_prediction = TRUE, eff_m= 0.95, bb = 0.2, 
                                 tree_control = TREE.control(tree_type  = "B", d = d, num_dir = 200),
                                 type = "RR", shrinkage  = shrinkage, precision = precision, 
                                 init_type = "LADTree", n_init = n_init, niter = niter, 
                                 nknot = nknot, save_f = FALSE, trace =  TRUE, save_tree = FALSE, error_type = c("mse"))
  
  return(list(control.l2 = control.l2,  control.lad = control.lad, control.lad.m = control.lad.m, control.rr = control.rr))
}


exp.tree.methods<- function(x_train, y_train, x_val, y_val, x_test, y_test,  grid, t_range, control.tree.list, validation_tree){
  
  model.tree.list <- list()
  err_trains <- err_tests <- err_vals <-  loss_vals <- matrix(NA, length(control.tree.list), control.tree.list[[1]]$niter)
  res_trains <- res_tests <- res_vals <- list()
  err_test <- rep(NA, length(control.tree.list))
  err_val <- rep(NA, length(control.tree.list))
  
  early_stop <- rep(NA, length(control.tree.list))
  time_vec <- rep(NA, length(control.tree.list))
  
  err_cvs <- NULL;  params <- NULL; errs_test <- NULL
  
  for(i in 1:length(control.tree.list)){
    
    print(paste(i, "out of", length(control.tree.list), "shrinkage", control.tree.list[[i]]$shrinkage, 
                "d",control.tree.list[[i]]$tree_control$d, "nknots", control.tree.list[[i]]$nknot))
    
    if(i <= 3){
       model.tree.list[[i]] <- RTFBoost(x_train = x_train, y_train = y_train, x_val = x_val,  y_val = y_val,
                                                       x_test = x_test, y_test = y_test, grid = grid, t_range  = t_range, 
                                                       control = control.tree.list[[i]])
       if(control.tree.list[[i]]$type == "LAD-M"){
         when_init_lad <- model.tree.list[[i]]$when_init
       }
    }else{
      if(validation_tree){
       model.tree.list[[i]]  <- RTFBoost.validation(x_train = x_train, z_train = NULL, y_train = y_train,  x_val = x_val,  z_val = NULL, 
                                                    y_val = y_val, x_test = x_test, z_test = NULL, y_test = y_test, grid = grid,  t_range = t_range, 
                                               max_depth_init_set = c(1,2,3), min_leaf_size_init_set = c(10,20,30), control = control.tree.list[[i]])
       err_cvs <-  model.tree.list[[i]]$err_cvs 
       params <-  model.tree.list[[i]]$params
       errs_test <-  model.tree.list[[i]]$errs_test
       
      }else{
        tmp_control <-  control.tree.list[[i]]
        tmp_control$init_type <- "median"
        model.tree.list[[i]]  <- RTFBoost(x_train = x_train, y_train = y_train, x_val = x_val,  y_val = y_val,
                 x_test = x_test, y_test = y_test, grid = grid, t_range  = t_range, 
                 control =   tmp_control)
      }
      when_init_rr <- model.tree.list[[i]]$when_init

    }
    
    err_trains[i,1:length(model.tree.list[[i]]$err_train)] <-model.tree.list[[i]]$err_train
    
    err_vals[i,1:length(model.tree.list[[i]]$err_val)] <-model.tree.list[[i]]$err_val
    err_val[i] <- model.tree.list[[i]]$err_val[model.tree.list[[i]]$early_stop]
    
    loss_vals[i,1:length(model.tree.list[[i]]$loss_val)] <- model.tree.list[[i]]$loss_val
    err_tests[i,1:model.tree.list[[i]]$early_stop] <-model.tree.list[[i]]$err_test[,1]
    err_test[i] <- model.tree.list[[i]]$err_test[model.tree.list[[i]]$early_stop,]

    early_stop[i] <- model.tree.list[[i]]$early_stop
    
    res_trains[[i]] <- model.tree.list[[i]]$f_train_t - y_train
    res_vals[[i]] <- model.tree.list[[i]]$f_val_t - y_val
    res_tests[[i]] <- model.tree.list[[i]]$f_test_t - y_test
    model.tree.list[[i]] <- NULL
    
  }
  
  
  return(list(errs_test = errs_test, params =   params, loss_vals = loss_vals, when_init_lad= when_init_lad, when_init_rr = when_init_rr,  err_cvs =  err_cvs , pred_train = res_trains, pred_val = res_vals, pred_test = res_tests, 
            err_trains = err_trains, err_tests = err_tests, err_vals = err_vals, err_test = err_test, err_val = err_val, early_stop = early_stop))
}


### choosing from c("FPPR", "FGAM", "MFLM", "RFSIR", "RFPLM", and "RTFBoost")
##  type: type of outliers 
##  validation_tree: using LADTree for initialization or not 
## 
do.exp <- function(seed, g_func_no, SNR, d, control.tree.list, methods = c("RTFBoost"), type,  validation_tree, competitor.control) {
  
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
    dat2return <- c(dat2return, list(RTFBoost = exp.tree.methods(x_train, y_train, x_val, y_val, x_test, y_test,  grid, t_range, control.tree.list, validation_tree)))
  }
  
  if("RFPLM" %in% methods){
    u <- 1:length(y_train)
    range_beta <- 4:7
    range_eta <- 20
    norder  <- 4
    model_RobustFPLM <- FPLMBsplines(y = y_train, x = x_train, u = u, t = grid, range_freq = range_beta, range_spl = range_eta, 
                                     norder = norder, fLoss = "lmrob", trace = TRUE)
    tmp_RobustFPLM <- FPLMBsplines.predict(model_RobustFPLM, newx = x_test, newy = y_test)
    dat2return <- c(dat2return, list(RFPLM = list(err_test =tmp_RobustFPLM$error, pred_test =  tmp_RobustFPLM$pred)))
  }
  
  if("MFLM" %in% methods){
    u <- 1:length(y_train)
    range_beta <- 4:7
    range_eta <- 20
    norder  <- 4
    
    model_MFLM <- FPLMBsplines(y = y_train, x = x_train, u = u, t = grid, range_freq = range_beta, range_spl = range_eta, 
                               norder = norder, fLoss = "huang", trace = TRUE)
    tmp_MFLM <- FPLMBsplines.predict(model_MFLM, newx = x_test, newy = y_test)
    dat2return <- c(dat2return, list(MFLM = list(err_test =tmp_MFLM$error, pred_test = tmp_MFLM$pred)))
  }
  
  if("FGAM" %in% methods){
    
    xx <- x_train
    yy <- y_train
    nbasises_FGAM <- competitor.control$nbasises_FGAM
    model_FGAM <- fgam(yy ~ af(xx, splinepars=list(k=c(nbasises_FGAM,nbasises_FGAM),m=list(c(2,2),c(2,2)))), gamma = 1.2, method="REML")
    
    pred_train <- predict(model_FGAM ,newdata=list(xx =x_train),type='response')
    pred_test <-  predict(model_FGAM ,newdata=list(xx =x_test),type='response')
    err_test <-  mean((pred_test- y_test)^2)
    err_train <-   mean((pred_train-y_train)^2)
    
    FGAM <- list(err_test  = err_test, err_train =  err_train, pred_train = pred_train, pred_test - pred_test)
    
    dat2return <- c(dat2return, list(FGAM =FGAM))
  }
  
  if("FPPR" %in% methods){
    niter_FPPR <- competitor.control$niter_FPPR
    nknots_FPPR <- competitor.control$nknots_FPPR
    
    model_FPPR <- FPPR(x_train, y_train, x_val, y_val, x_test, y_test, grid = grid, t_range = t_range,  niter = niter_FPPR, nknots = nknots_FPPR, err_type = "robust")
    err_trains_FPPR <- model_FPPR$err_train
    err_tests_FPPR <- model_FPPR$err_test
    err_vals_FPPR <-  model_FPPR$err_val
    err_val_FPPR  <- model_FPPR$err_val[model_FPPR$early_stop]
    err_test_FPPR <- model_FPPR$err_test[model_FPPR$early_stop]
    Js_FPPR <- model_FPPR$J
    early_stop_FPPR <- model_FPPR$early_stop
    FPPR_obj <- list( pred_train =  model_FPPR$f_train_t, pred_val =  model_FPPR$f_val_t, pred_test =  model_FPPR$f_test_t, 
                      err_trains =  err_trains_FPPR,  err_tests_FPPR =  err_tests_FPPR,
                      err_vals = err_vals_FPPR,  err_val = err_val_FPPR, err_test = err_test_FPPR,
                      J =  Js_FPPR ,  early_stop =  early_stop_FPPR)
    dat2return <- c(dat2return, list(FPPR = FPPR_obj))
    
  }
  
  if("RFSIR" %in% methods){
    
    qs <- competitor.control$qs
    k <- competitor.control$k
    model_RFSIR <- RFSIR(x_train, y_train, x_val, y_val, x_test, y_test, k= 10, qs = qs, make_prediction =  TRUE)
    dat2return <- c(dat2return, list(RFSIR = model_RFSIR))
    
  }
  dat2return$dat <- dat
  
  
  return(dat2return)
}

