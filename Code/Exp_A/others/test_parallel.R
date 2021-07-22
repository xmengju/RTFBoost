## testing the performance of our method, without validaiton, and in parallel 


## testing the code on the server, and testing parallel on local ...... 

rm(list = ls())

try_parallel <- function(i, type =  "C2"){
  

  g_func_no <- 5
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
  source("Code/Exp_A/fun_exp.R")
  source("Code/Exp_A/dat_generate.R")
  source("Code/FSIM.R")
  source("Code/RFSIR.R")
  
  control <- dat.generate.control(n_train = 400,  n_val = 200, n_test = 1000, g_func_no =   g_func_no, SNR = 5)
  dat <- dat.generate(seed = i, control = control, type = type, direction = "symmetric")
  
  x_train <- dat$x$x_train;  x_val <- dat$x$x_val;  x_test <- dat$x$x_test
  n_train <- nrow(x_train);  n_val <- nrow(x_val);  n_test <- nrow(x_test)
  y_train <- dat$y$y_train;  y_val <- dat$y$y_val;  y_test <- dat$y$y_test
  
  t_range <- c(0,1)
  grid <- dat$tt
  
  xx <- x_train
  yy <- y_train
  nbasises_FGAM <- 15
  #model_FGAM <- fgam(yy ~ af(xx, splinepars=list(k=c(nbasises_FGAM,nbasises_FGAM),m=list(c(2,2),c(2,2)))), gamma = 1.2, method="REML")
  
  #pred_train <- predict(model_FGAM ,newdata=list(xx =x_train),type='response')
  #pred_test <-  predict(model_FGAM ,newdata=list(xx =x_test),type='response')
  #err_test_FGAM <-  mean((pred_test - y_test)^2)
  #err_train_FGAM <-   mean((pred_train - y_train)^2)

  d <- 1; shrinkage.l2 <- 0.05; precision <- 6; niter <- 1000; nknot<- 3
  #control.l2 <-  RTFBoost.control(make_prediction = TRUE,
  #                                tree_control = TREE.control(tree_type  = "B", d = d, num_dir = 200),
  #                                type = "L2", shrinkage  = shrinkage.l2, precision = precision,
  #                                init_type = "mean", niter = niter,
  #                                nknot = nknot, save_f = FALSE, trace =  TRUE, save_tree = FALSE, error_type = c("mse"))

 # model.l2 <- RTFBoost(x_train = x_train, y_train = y_train, x_val = x_val,  y_val = y_val,
 #                      x_test = x_test, y_test = y_test, grid = grid, t_range  = t_range,
#                       control = control.l2)

  shrinkage.lad <- 0.05
  control.lad <-  RTFBoost.control(make_prediction = TRUE,
                                  tree_control = TREE.control(tree_type  = "B", d = d, num_dir = 200),
                                  type = "LAD", shrinkage  = shrinkage.lad, precision = precision,
                                  init_type = "mean", niter = niter,
                                  nknot = nknot, save_f = FALSE, trace =  TRUE, save_tree = FALSE, error_type = c("mse"))
  
  model.lad <- RTFBoost(x_train = x_train, y_train = y_train, x_val = x_val,  y_val = y_val,
                       x_test = x_test, y_test = y_test, grid = grid, t_range  = t_range,
                       control = control.lad)
  
  
  
  shrinkage.rr <- 0.05; n_init <- 500
  niter <- 1000
  control.rr <- RTFBoost.control(make_prediction = TRUE, eff_m= 0.95, bb = 0.5,
                                 tree_control = TREE.control(tree_type  = "B", d = d, num_dir = 200),
                                 type = "RR", shrinkage  = shrinkage.rr, precision = precision,
                                 init_type = "median", n_init = n_init, niter = niter,
                                 nknot = nknot, save_f = FALSE, trace =  TRUE, save_tree = FALSE, error_type = c("mse"))

  model.rr <- RTFBoost(x_train = x_train, y_train = y_train, x_val = x_val,  y_val = y_val,
                       x_test = x_test, y_test = y_test, grid = grid, t_range  = t_range,
                       control = control.rr)

  return(list(when_init = model.rr$when_init, err_test_rr =model.rr$err_test,
              err_test_lad =model.lad$err_test, early_stop_lad = model.lad$early_stop, early_stop_rr = model.rr$early_stop))
  
}

save_data <-try_parallel(6)

library(doParallel)

save_data_all <- list()
for(type in c("C0","C1","C2","C3","C4","C5")){
  cl <- makeCluster(4)
  clusterExport(cl=cl, varlist=c("type", "try_parallel"))
  system.time(save_data <- parLapply(cl, 1:4, function(x) try_parallel(x,type)))
  stopCluster(cl)
  save_data_all[[type]] <- save_data
}




tmp <- save_data_all[["C4"]] 
tmp <- save_data
res <- NULL

for(i in 1:4){
  res =  rbind(res, c(tmp[[i]]$err_test_FGAM,tail(tmp[[i]]$err_test_l2,1), tail(tmp[[i]]$err_test_lad,1), tail(tmp[[i]]$err_test_rr,1)))
}
res

tmp[[i]]$early_stop_rr



save_data[[1]]$when_init

  