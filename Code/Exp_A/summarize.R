## summary code
robust_mse <- function(x){
  return( median(x)^2 + mscale(x))
}

summarize0 <- function(g_func_no, type, SNR = 5, nknot = 3, seeds = 1:20, summarize_type = "pred"){
  
  dir <- paste("Results/Results_0_type_", type, "_g_", g_func_no, "_SNR_", SNR,"_nknot", nknot, sep = "")
  file_list <- list.files(dir)
  
  res <- NULL
  
  for(seed in seeds){
    
    file_name <-  paste("g_func_no_", g_func_no, "_case_id_0_d_1",  "_type_",type, "_seed_", seed,".RData", sep = "")
    if(file_name %in% file_list){
      load(paste(dir,"/",file_name,sep = ""))
      res <- rbind(res, c(dat2save$FPPR$err_test, dat2save$FGAM$err_test,dat2save$RFPLM$err_test, dat2save$MFLM$err_test, 
                          dat2save$RFSIR$err_test))
    }
  }
  
  colnames(res) <- c("FPPR","FGAM","RFPLM", "MFLM","RFSIR")
  return(res)
}

summarize1 <- function(g_func_no, type, SNR = 5, nknot = 3, seeds = 1:20, summarize_type = "pred"){
  
  dir <- paste("Results/Results_1_type_", type, "_g_", g_func_no, "_SNR_", SNR,"_nknot", nknot, sep = "")
  file_list <- list.files(dir)
  
  res <- NULL
  res_early <- NULL
  for(seed in seeds){
    tmp_val <- NULL
    tmp_test <- NULL
    tmp_early <- NULL
    for(d in 1:4){
      file_name <-  paste("g_func_no_", g_func_no, "_case_id_1_d_", d, "_type_", type,"_seed_", seed,".RData",sep = "")
      if(file_name %in% file_list){
        load(paste(dir,"/",file_name,sep = ""))
        tmp_val <- rbind(tmp_val, c(robust_mse(dat2save$RTFBoost$pred_val[[1]]),robust_mse(dat2save$RTFBoost$pred_val[[2]]),
                                    robust_mse(dat2save$RTFBoost$pred_val[[3]])))
        tmp_test <- rbind(tmp_test,  dat2save$RTFBoost$err_test)
        tmp_early <- rbind(tmp_early, c(dat2save$RTFBoost$when_init,dat2save$RTFBoost$early_stop))
      }
    }
    
    idx <- apply(tmp_val, 2, which.min)
    tmp <- rep(NA, 3)
    for(j in 1:3){
      tmp[j] <-tmp_test[idx[j],j]
    }
    res <- rbind(res, tmp)
    res_early  <- rbind(res_early,  c( tmp_early[idx[3],1], tmp_early[idx[1],2], tmp_early[idx[2],3], tmp_early[idx[3],4]))
  }
  
  dat2return <- list(res = res, res_early = res_early)
  colnames(res) <- c("RTFBoost(L2)", "RTFBoost(LAD)", "RTFBoost(RR)")
  if(summarize_type == "pred"){
      return(res)
  }else{
    return(dat2return)
  }
}

summarize1_same_d_median <- function(g_func_no, type, d, SNR = 5, nknot = 3,  seeds = 1:20){
  
  dir <- paste("Results/Results_1_type_", type, "_g_", g_func_no, "_SNR_", SNR,"_nknot", nknot, sep = "")
  file_list <- list.files(dir)
  
  res_val <- res_test <- NULL
  res_early <- NULL
  res_rr <- NULL
  for(seed in seeds){
  
    file_name <-  paste("g_func_no_", g_func_no, "_case_id_1_d_", d, "_type_", type,"_seed_", seed,".RData",sep = "")
    if(file_name %in% file_list){
        load(paste(dir,"/",file_name,sep = ""))
        res_val <- rbind(res_val, c(robust_mse(dat2save$RTFBoost$pred_val[[1]]),robust_mse(dat2save$RTFBoost$pred_val[[2]]),
                                    robust_mse(dat2save$RTFBoost$pred_val[[3]])))
        res_test <- rbind(res_test,  dat2save$RTFBoost$err_test)
        res_rr <- rbind(res_rr,  matrix(dat2save$RTFBoost$errs_test, nrow = 1))
        res_early <- rbind(res_early, c(dat2save$RTFBoost$when_init,dat2save$RTFBoost$early_stop))
      }
    }
    
  dat2return <- list(res_test = res_test, res_val = res_val, res_rr = res_rr, res_early = res_early)
}
  


d = 1

res <- summarize1_same_d_median(g_func_no = 4, type = "C3", d = d, SNR = 5, nknot = 3,  seeds = 1:20)
boxplot(res$res_test)

res$res_early

boxplot(res$res_test) 

boxplot(res$res_rr) 

g_func_no = 4; type = "C5"; d <- 1; SNR = 5; nknot <- 3; seed <- 6
dir <- paste("Results/Results_1_type_", type, "_g_", g_func_no, "_SNR_", SNR,"_nknot", nknot, sep = "")
file_name <-  paste("g_func_no_", g_func_no, "_case_id_1_d_", d, "_type_", type,"_seed_", seed,".RData",sep = "")
load(paste(dir,"/",file_name,sep = ""))
plot(dat2save$RTFBoost$err_vals[3,])
plot(dat2save$RTFBoost$err_vals[2,])


#rr <- summarize1(g_func_no = 2, type= "C3", SNR = 5, nknot = 3, seeds = 1:30, summarize_type = "others")
 # rr$res_early
## pre analysis 

save_data[[1]]
