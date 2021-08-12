## testing the performance of our method, without validaiton, and in parallel 


## testing the code on the server, and testing parallel on local ...... 

rm(list = ls())

try_parallel <- function(d, seed, g_func_no,  type =  "C0"){
  
  num_dir <- 200
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
  source("Code/Exp_C/fun_exp.R")
  source("Code/Exp_C/dat_generate.R")
  source("Code/FSIM.R")
  source("Code/RFSIR.R")
  

  control <- dat.generate.control(n_train = 400,  n_val = 200, n_test = 1000, g_func_no =   g_func_no, SNR = 5)
  dat <- dat.generate(seed = seed, control = control, type = type, direction = "symmetric")
  
  x_train <- dat$x$x_train;  x_val <- dat$x$x_val;  x_test <- dat$x$x_test
  n_train <- nrow(x_train);  n_val <- nrow(x_val);  n_test <- nrow(x_test)
  y_train <- dat$y$y_train;  y_val <- dat$y$y_val;  y_test <- dat$y$y_test
  
  t_range <- c(0,1)
  grid <- dat$tt


  precision <- 6; niter <- 1000; nknot<- 3

  
  u <- 1:length(y_train)
  range_beta <- 4:7
  range_eta <- 20
  norder  <- 4
  model_RobustFPLM <- FPLMBsplines(y = y_train, x = x_train, u = u, t = grid, range_freq = range_beta, range_spl = range_eta, 
                                   norder = norder, fLoss = "huang", trace = TRUE)
  tmp_RobustFPLM <- FPLMBsplines.predict(model_RobustFPLM, newx = x_test, newy = y_test)
  return(err_test =tmp_RobustFPLM$error)
}
  
  #shrinkage.lad <- 0.05
  #control.lad <-  RTFBoost.control(make_prediction = TRUE,
   #                               tree_control = TREE.control(tree_type  = "B", d = d, num_dir = num_dir),
  #                                type = "LAD", shrinkage  = shrinkage.lad, precision = precision,
  #                                init_type = "mean", niter = niter,
  #                                nknot = nknot, save_f = FALSE, trace =  TRUE, save_tree = FALSE, error_type = c("mse"))
  
 # model.lad <- RTFBoost(x_train = x_train, y_train = y_train, x_val = x_val,  y_val = y_val,
#                       x_test = x_test, y_test = y_test, grid = grid, t_range  = t_range,
#                       control = control.lad)
 
  
#  shrinkage.rr <- 0.05;
 # n_init <- 500
#  niter <- 1000
 # control.rr <- RTFBoost.control(make_prediction = TRUE, eff_m= 0.95, bb = 0.2,
  #                               tree_control = TREE.control(tree_type  = "B", d = d, num_dir = num_dir),
  #                               type = "RR", shrinkage  = shrinkage.rr, precision = precision,
  #                               init_type = "median", n_init = n_init, niter = niter,
  #                               nknot = nknot, save_f = FALSE, trace =  TRUE, save_tree = FALSE, error_type = c("mse"))
#  
 # model.rr <- RTFBoost(x_train = x_train, y_train = y_train, x_val = x_val,  y_val = y_val,
  #                     x_test = x_test, y_test = y_test, grid = grid, t_range  = t_range,
  #                     control = control.rr)
 # 
#return(list(num_dir = num_dir, 
#            err_train_lad = model.lad$err_train, loss_train_lad = model.lad$loss_train,loss_val_lad = model.lad$loss_val,
#              err_test_lad =model.lad$err_test, early_stop_lad = model.lad$early_stop,
#              err_train_rr = model.rr$err_train, loss_train_rr = model.rr$loss_train,loss_val_rr = model.rr$loss_val,
#              err_test_rr =model.rr$err_test, early_stop_lad = model.rr$early_stop
#              ))

#save_data <-try_parallel(6)

library(doParallel)

try_parallel(d, seed = 5, g_func_no,  type =  type)

seed <- 2
save_data_all <- list()
g_func_no <- 3
for(d in 1){
  save_data_all[[d]] <- list()
  for(type in c("C0","C1","C2","C3","C4","C5")){
    print(c(d, type))
    cl <- makeCluster(5)
    clusterExport(cl=cl, varlist=c("try_parallel","g_func_no",'d',"type"))
    system.time(save_data <- parLapply(cl, 1:5, function(x)try_parallel(d, seed = x, g_func_no,  type =  type)))
    stopCluster(cl)
    
    save_data_all[[d]][[type]] <- save_data
   }
}

res <- NULL
d <- 1
type <- "C2"
for(i in 1:5){
  res <- rbind(res, c(tail(save_data_all[[d]][["C0"]][[i]],1),tail(save_data_all[[d]][["C1"]][[i]],1),tail(save_data_all[[d]][["C2"]][[i]],1),
                      tail(save_data_all[[d]][["C3"]][[i]],1),tail(save_data_all[[d]][["C4"]][[i]],1), tail(save_data_all[[d]][["C5"]][[i]],1)))
}

save_data_all_3[[4]]$num_dir
save_data <- save_data_all[[4]]

par(mfrow = c(1,1))
plot(save_data[[1]]$err_test_l2[,1],type = "l", lwd= 3,ylim = c(0,2))
lines(save_data[[2]]$err_test_l2[,1], lwd = 3, col = "red")
lines(save_data[[3]]$err_test_l2[,1],  lwd = 3,col = "green")
lines(save_data[[4]]$err_test_l2[,1],  lwd = 3,col = "blue")


plot(save_data[[1]]$err_test_lad[,1],type = "l", lwd= 3,ylim = c(0,2))
lines(save_data[[2]]$err_test_lad[,1], col = "red")
lines(save_data[[3]]$err_test_lad[,1], col = "green")
lines(save_data[[4]]$err_test_lad[,1], col = "blue")

lines(save_data[[1]]$err_test_rr[,1],type = "l", lty = 3, lwd= 3,ylim = c(0,2))
lines(save_data[[2]]$err_test_rr[,1],lwd= 3,  lty = 3,col = "red")
lines(save_data[[3]]$err_test_rr[,1],lwd= 3,  lty = 3,col = "green")
lines(save_data[[4]]$err_test_rr[,1], lwd= 3,  lty = 3,col = "blue")


plot(save_data[[1]]$loss_val_rr,type = "l", lty = 3, lwd= 3)
lines(save_data[[2]]$loss_val_rr,lwd= 3,  lty = 3,col = "red")
lines(save_data[[3]]$loss_val_rr,lwd= 3,  lty = 3,col = "green")
lines(save_data[[4]]$loss_val_rr, lwd= 3,  lty = 3,col = "blue")

plot(save_data[[1]]$loss_train_rr,type = "l", lty = 3, lwd= 3)
lines(save_data[[2]]$loss_train_rr,lwd= 3,  lty = 3,col = "red")
lines(save_data[[3]]$loss_train_rr,lwd= 3,  lty = 3,col = "green")
lines(save_data[[4]]$loss_train_rr, lwd= 3,  lty = 3,col = "blue")





# 
# tmp <- save_data_all[["C4"]] 
# res <- NULL
# 
# for(i in 1:10){
#   print(c(tmp[[i]]$when_init,tmp[[i]]$early_stop_lad, tmp[[i]]$early_stop_rr))
#   res =  rbind(res, c(tmp[[i]]$err_test_FGAM,tail(tmp[[i]]$err_test_l2,1), tail(tmp[[i]]$err_test_lad,1), tail(tmp[[i]]$err_test_rr,1)))
# }
# boxplot(unlist(res[,1]), unlist(res[,2]))
# 
# tmp[[i]]$early_stop_rr
# 
# 
# 
# save_data[[1]]$when_init
# 
#   