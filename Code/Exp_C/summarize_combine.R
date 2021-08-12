### analyze partial results 
SNR <- 5; nknot <- 3
robust_mse <- function(x){
  return( median(x)^2 + mscale(x)^2)
}


summarize_1_combine <- function(g_func_no, type, SNR = 5, nknot = 3, seeds = 1:20, summarize_type = "median"){
  
  dir_1 <- paste("Results_cedar6/Results_1_type_", type, "_g_", g_func_no, "_SNR_", SNR,"_nknot", nknot, sep = "")
  dir_2 <- paste("Results_LAD/Results_1_type_", type, "_g_", g_func_no, "_SNR_", SNR,"_nknot", nknot, sep = "")
  dir_3 <- paste("Results_cedar9/Results_1_type_", type, "_g_", g_func_no, "_SNR_", SNR,"_nknot", nknot, sep = "")
  
  file_list_1 <- list.files(dir_1)
  file_list_2 <- list.files(dir_2)
  file_list_3 <- list.files(dir_3)
  
  res <- NULL
  res_early <- NULL
  res_d <- NULL
  for(seed in seeds){
    tmp_val <- NULL
    tmp_test <- NULL
    tmp_early <- NULL
    for(d in 1:4){
      file_name <-  paste("g_func_no_", g_func_no, "_case_id_1_d_", d, "_type_", type,"_seed_", seed,".RData",sep = "")
      if( (file_name %in% file_list_1) & (file_name %in% file_list_2) & (file_name %in% file_list_3)){
        load(paste(dir_1,"/",file_name,sep = ""))
        dat2save1 <- dat2save
        load(paste(dir_2,"/",file_name,sep = ""))
        dat2save2 <- dat2save
        load(paste(dir_3,"/",file_name,sep = ""))
        dat2save3 <- dat2save
        
        # L2, LAD, LAD-M, RR.5, RR.2
        tmp_val <- rbind(tmp_val, c(robust_mse(dat2save1$RTFBoost$pred_val[[1]]),robust_mse(dat2save2$RTFBoost$pred_val[[1]]), 
                                    robust_mse(dat2save3$RTFBoost$pred_val[[2]]),  robust_mse(dat2save3$RTFBoost$pred_val[[3]]),
                                    robust_mse(dat2save1$RTFBoost$pred_val[[3]])))
        if(summarize_type == "median"){
          tmp_test <- rbind(tmp_test,  c(dat2save1$RTFBoost$err_test[1], dat2save2$RTFBoost$err_test[1],
                                         dat2save3$RTFBoost$err_test[2], dat2save3$RTFBoost$errs_test[1],
                                         dat2save1$RTFBoost$errs_test[1]))
        }else{
          tmp_test <- rbind(tmp_test,  c(dat2save1$RTFBoost$err_test[1], dat2save2$RTFBoost$err_test[1],
                                         dat2save3$RTFBoost$err_test[2], dat2save3$RTFBoost$err_test[3], 
                                         dat2save1$RTFBoost$err_test[3]))
        }
        #when_init, LAD-M, RR.5, RR.2
        #early_stop
        tmp_early <- rbind(tmp_early, c(dat2save3$RTFBoost$when_init_lad, dat2save3$RTFBoost$when_init_rr,
                                        dat2save1$RTFBoost$when_init_rr,
                                        dat2save1$RTFBoost$early_stop[1], dat2save2$RTFBoost$early_stop[1],
                                        dat2save3$RTFBoost$early_stop[2],  dat2save3$RTFBoost$early_stop[3],
                                        dat2save1$RTFBoost$early_stop[3]))
      
        }
    }
    
    idx <- apply(tmp_val, 2, which.min)
    tmp <- rep(NA, 5)
    for(j in 1:5){
      tmp[j] <-tmp_test[idx[j],j]
    }

    
    res_d <- rbind(res_d, idx)
    res <- rbind(res, tmp)
    res_early  <- rbind(res_early,  c(tmp_early[idx[3],1], tmp_early[idx[4],2],tmp_early[idx[5],3],
                                      tmp_early[idx[1],4], tmp_early[idx[2],5],  tmp_early[idx[3],6],  tmp_early[idx[4],7],
                                      tmp_early[idx[5],8]))
  }
  
  colnames(res) <- c("RTFBoost(L2)", "RTFBoost(LAD)", "RTFBoost(LAD-M)","RTFBoost(RR.5)", "RTFBoost(RR.2)")
  dat2return <- list(res = res, res_early = res_early,   res_d  =   res_d )
  
}
  

summarize0 <- function(g_func_no, type, SNR = 5, nknot = 3, seeds = 1:20, summarize_type = "pred"){
  
  dir <- paste("Results_cedar6/Results_0_type_", type, "_g_", g_func_no, "_SNR_", SNR,"_nknot", nknot, sep = "")
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



summarize1_same_d_median <- function(g_func_no, type, d, SNR = 5, nknot = 3,  seeds = 1:20, method = "LAD"){
  
  if(method == "LAD"){
    dir <- paste("Results_LAD/Results_1_type_", type, "_g_", g_func_no, "_SNR_", SNR,"_nknot", nknot, sep = "")
  }else{
    dir <- paste("Results/Results_1_type_", type, "_g_", g_func_no, "_SNR_", SNR,"_nknot", nknot, sep = "")
  }
  file_list <- list.files(dir)
  
  res_val <- res_test <- NULL
  res_early <- NULL
  res_rr <- NULL
  for(seed in seeds){
    
    file_name <-  paste("g_func_no_", g_func_no, "_case_id_1_d_", d, "_type_", type,"_seed_", seed,".RData",sep = "")
    if(file_name %in% file_list){
      load(paste(dir,"/",file_name,sep = ""))
      
      if(method == "LAD"){
        res_val <- rbind(res_val, c(robust_mse(dat2save$RTFBoost$pred_val[[1]])))
        }else{
          res_val <- rbind(res_val, c(robust_mse(dat2save$RTFBoost$pred_val[[1]]),robust_mse(dat2save$RTFBoost$pred_val[[2]]),
                                      robust_mse(dat2save$RTFBoost$pred_val[[3]])))
        }
      res_test <- rbind(res_test,  dat2save$RTFBoost$err_test)

      if(method == "LAD"){
        res_early <- rbind(res_early, c(dat2save$RTFBoost$early_stop))
        
      }else{
        res_early <- rbind(res_early, c(dat2save$RTFBoost$when_init_rr,dat2save$RTFBoost$early_stop))
      }
    }
  }
  
  dat2return <- list(res_test = res_test, res_val = res_val,res_early = res_early)
}



# y_lim <- list(c(0.1,0.4),  c(0.15,0.6),  c(0.25,0.8), 
#               c(0.25,1.1),  c(0.45,1.5))
# 
# g_func_no = 2
# type = "C1"
# pdf("test.pdf", width = 18, height = 10)
# for(g_func_no in 1:5){
#   par(mfrow = c(2,1))
#   for(type in paste("C",0:5, sep = "")){
#     res1 <- summarize1_same_d_median(g_func_no, type, d = 1, SNR = 5, nknot = 3,  seeds = 1:20, method = "LAD")
#     res2 <- summarize1_same_d_median(g_func_no, type, d= 2, SNR = 5, nknot = 3,  seeds = 1:20, method = "LAD")
#     res3 <- summarize1_same_d_median(g_func_no, type, d= 3, SNR = 5, nknot = 3,  seeds = 1:20, method = "LAD")
#     res4 <- summarize1_same_d_median(g_func_no, type, d =4, SNR = 5, nknot = 3,  seeds = 1:20, method = "LAD")
#     
#     res11 <- summarize1_same_d_median(g_func_no, type, d = 1, SNR = 5, nknot = 3,  seeds = 1:20, method = "others")
#     res22 <- summarize1_same_d_median(g_func_no, type, d = 2, SNR = 5, nknot = 3,  seeds = 1:20, method = "others")
#     res33 <- summarize1_same_d_median(g_func_no, type, d = 3, SNR = 5, nknot = 3,  seeds = 1:20, method = "others")
#     res44 <- summarize1_same_d_median(g_func_no, type, d = 4, SNR = 5, nknot = 3,  seeds = 1:20, method = "others")
#     
#     boxplot(c(res1$res_test),c(res2$res_test),c(res3$res_test),c(res4$res_test),
#             c(res11$res_test[,2]),c(res22$res_test[,2]),c(res33$res_test[,2]),c(res44$res_test[,2]),
#             c(res11$res_test[,3]),c(res22$res_test[,3]),c(res33$res_test[,3]),c(res44$res_test[,3]),
#             main = paste( "r", g_func_no, ": ", type, sep = ""), ylim =y_lim[[g_func_no]],names = c("LAD.1","LAD.2","LAD.3","LAD.4",
#                                                                                    "LAD-M.1","LAD-M.2","LAD-M.3","LAD-M.4",   "RR(.2).1","RR(.2).2","RR(.2).3","RR(.2).4")  )
#   }
# }
# dev.off()



