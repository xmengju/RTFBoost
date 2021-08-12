
rm(list = ls())

try_parallel <- function(g_func_no , seed, type){

  
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
  
  control <- dat.generate.control(n_train = 400,  n_val = 200, n_test = 1000, g_func_no = g_func_no, SNR = 5)
  dat <- dat.generate(seed = seed, control = control, type = type, direction = "asymmetric")
  
  x_train <- dat$x$x_train;  x_val <- dat$x$x_val;  x_test <- dat$x$x_test
  n_train <- nrow(x_train);  n_val <- nrow(x_val);  n_test <- nrow(x_test)
  y_train <- dat$y$y_train;  y_val <- dat$y$y_val;  y_test <- dat$y$y_test
  
  t_range <- c(0,1)
  grid <- dat$tt
  

  u <- 1:length(y_train)
  range_beta <- 4:7
  range_eta <- 20
  norder  <- 4
  model_RobustFPLM <- FPLMBsplines(y = y_train, x = x_train, u = u, t = grid, range_freq = range_beta, range_spl = range_eta, 
                                   norder = norder, fLoss = "lmrob", trace = TRUE)
  tmp_RobustFPLM <- FPLMBsplines.predict(model_RobustFPLM, newx = x_test, newy = y_test)

  return(list(err_test =tmp_RobustFPLM$error))
}


res <- try_parallel(g_func_no = 1 , seed = 1, type = "C0")

  

save_data_all <- list()

g_func_no <- 5
for(type in paste("C", 0:5, sep = "")){
  cl <- makeCluster(5)
  clusterExport(cl=cl, varlist=c("try_parallel","type", "g_func_no"))
  system.time(save_data <- parLapply(cl,1:20, function(x) try_parallel(g_func_no,x,type)))
  stopCluster(cl)
  save_data_all[[type]] <- save_data
}
res <- NULL
for(type in paste("C", 0:5, sep = "")){
  res <- rbind(res, unlist(save_data_all[[type]]))
}

boxplot(t(res))


