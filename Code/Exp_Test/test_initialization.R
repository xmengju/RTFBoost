
rm(list = ls())


cal_ys <- function(x, t_range){
  
  deg <- 4; p <- ncol(x)
  grid0 <- seq(t_range[1],t_range[2], 1/(10*(p-1))) # in case of not evenly spaced
  knot <- quantile(grid0, (1:nknot)/(nknot+1) )
  delta <- sort(c(rep(range(grid0), deg), knot)) #exterior knots
  B<- spline.des(delta, grid, deg)$design
  B <- compute.orthonormal(B,grid, t_range)
  predictors <- t(apply(x, 1, function(xx){apply(B, 2, function(bb){riemman (xx*bb, grid, t_range)})}))

  
  ff <- function(predictors){
    
    aa = quantile(predictors[,1], c(0.2,0.75))
    return( (predictors[,1] < aa[1])*10  +  (predictors[,1] > aa[2])*(-10) )
  }
  
 return(ff(predictors))
}
try_parallel <- function(  g_func_no , seed,d, type){

  
  library(rpart)
  library(splines)
  library(MyFPLM)
  library(rpart)
  library(refund)
  library(fda)
  library(robustbase)
  library(np)
  library(RobStatTM)
  
  source("Code/RTFBoost.R")
  source("Code/TREE.init.R")
  source("Code/TREE.R")
  source("Code/user_rpart.R")
  source("Code/utils.R")
  source("Code/Exp_Test/fun_exp.R")
  source("Code/Exp_Test/dat_generate.R")
  source("Code/FSIM.R")
  source("Code/RFSIR.R")
  
  num_dir <- 200; type <- "D0"
  control <- dat.generate.control(n_train = 400,  n_val = 200, n_test = 1000, g_func_no = g_func_no, SNR = 5)
  dat <- dat.generate(seed = seed, control = control, type = type, direction = "symmetric")
  
  x_train <- dat$x$x_train;  x_val <- dat$x$x_val;  x_test <- dat$x$x_test
  n_train <- nrow(x_train);  n_val <- nrow(x_val);  n_test <- nrow(x_test)
  
  y_train <- 5*cal_ys(x_train, t_range) + rnorm(n_train)
  y_val <- 5*cal_ys(x_val, t_range) + rnorm(n_val)
  y_test <- 5*cal_ys(x_test, t_range) + rnorm(n_test)

  t_range <- c(0,1)
  grid <- dat$tt
  
  precision <- 6; niter <- 2000; nknot<- 3
  
  shrinkage <- 0.05
  d <- 1
  num_dir <- 200
  n_init <- 1000
  
  control.rr <-  RTFBoost.control(make_prediction = TRUE, eff_m = 0.95, bb = 0.2, 
                                   tree_control = TREE.control(tree_type  = "B", d = d, num_dir = num_dir),
                                   type = "RR", shrinkage  = shrinkage, precision = precision,
                                   init_type = "median",max_depth_init = 2, min_leaf_size_init = 10,  niter = niter, n_init = n_init, 
                                   nknot = nknot, save_f = FALSE, trace =  TRUE, save_tree = FALSE, error_type = c("mse"))
  
  model.rr <- RTFBoost(x_train = x_train, y_train = y_train, x_val = x_val,  y_val = y_val,
                        x_test = x_test, y_test = y_test, grid = grid, t_range  = t_range,
                        control = control.rr)
  shrinkage <- 0.05;
  n_init <- 500
  niter <- 1000
  
  control.l2 <-  RTFBoost.control(make_prediction = TRUE, 
                                    tree_control = TREE.control(tree_type  = "B", d = d, num_dir = num_dir),
                                    type = "L2", shrinkage  = shrinkage, precision = precision,
                                    init_type = "median",  n_init = n_init, niter = niter,
                                    nknot = nknot, save_f = FALSE, trace =  TRUE, save_tree = FALSE, error_type = c("mse"))
  
  model.l2 <- RTFBoost(x_train = x_train, y_train = y_train, x_val = x_val,  y_val = y_val,
                         x_test = x_test, y_test = y_test, grid = grid, t_range  = t_range,
                         control = control.l2)
  
  
  control.lad.pre <-  RTFBoost.control(make_prediction = TRUE, eff_m = 0.95, bb = 0.2, 
                                  tree_control = TREE.control(tree_type  = "B", d = d, num_dir = num_dir),
                                  type = "LAD", shrinkage  = 1, precision = precision,
                                  init_type = "median",  niter = 5, 
                                  nknot = nknot, save_f = FALSE, trace =  TRUE, save_tree = FALSE, error_type = c("mse"))
  
  model.lad.pre <- RTFBoost(x_train = x_train, y_train = y_train, x_val = x_val,  y_val = y_val,
                       x_test = x_test, y_test = y_test, grid = grid, t_range  = t_range,
                       control = control.lad.pre)
  
  niter <- 1000
  n_init  <- 500

  control.rr.after <-  RTFBoost.control(make_prediction = TRUE, eff_m = 0.95, bb = 0.5, 
                                  tree_control = TREE.control(tree_type  = "B", d = d, num_dir = num_dir),
                                  type = "RR", shrinkage  = shrinkage, precision = precision,
                                  init_type = "median",max_depth_init = 2, min_leaf_size_init = 10,  niter = niter, n_init = n_init, 
                                  nknot = nknot, save_f = FALSE, trace =  TRUE, save_tree = FALSE, error_type = c("mse"))
  
  model.rr.after <- RTFBoost(x_train = x_train, y_train = y_train - model.lad.pre$f_train_t, x_val = x_val,  y_val = y_val - model.lad.pre$f_val_t,
                       x_test = x_test, y_test = y_test - model.lad.pre$f_test_t, grid = grid, t_range  = t_range,
                       control = control.rr.after)
  
  
  plot(model.l2$err_train, ylim = c(0,100))
  lines(model.rr$err_train, col = "blue")
  lines(model.rr.after$err_train, col = "blue")
  
  plot(model.l2$err_test[,1], ylim = c(0,2000), xlim = c(0,2000))
  lines(model.rr$err_test[,1], col = "blue")
  lines(model.rr.after$err_test[,1], col = "blue")
  
  return(list(num_dir = num_dir,
              err_train_lad = model.lad$err_train, loss_train_lad = model.lad$loss_train,loss_val_lad = model.lad$loss_val,
              err_test_lad =model.lad$err_test, early_stop_lad = model.lad$early_stop,
              err_train_rr.2 = model.rr.2$err_train, loss_train_rr.2 = model.rr.2$loss_train,loss_val_rr.2 = model.rr.2$loss_val,
              err_test_rr.2 =model.rr.2$err_test, 
              err_train_rr.5 = model.rr.5$err_train, loss_train_rr.5 = model.rr.5$loss_train,loss_val_rr.5 = model.rr.5$loss_val,
              err_test_rr.5 =model.rr.5$err_test
  ))
}

