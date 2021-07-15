### summary code 

robust_mse <- function(x){

  return( median(x)^2 + mscale(x))
}

mse <- function(x){
  mean(x^2)
}
summarize_0 <- function(g_func_no, SNR = 5, dd = "D0", nknot = 3, seeds = 1:30){
  dir <- paste("Results/Results_0_x_mattern_g_", g_func_no, "_SNR_", SNR,"_nknot", nknot, sep = "")
  file_list <- list.files(dir)
  
  res <- NULL
  for(seed in seeds){

    file_name <-  paste("g_func_no_", g_func_no, "_case_id_0_d_1",  "_dd_", dd,"_seed_", seed,".RData", sep = "")
    if(file_name %in% file_list){
      load(paste(dir,"/",file_name,sep = ""))
      res <- rbind(res, dat2save$fix$error)
    }
  }
  
  return(res)
  
}
g_func_no <- 1
SNR = 5
dd = "D0"
nknot = 3  
summarize_1 <- function(g_func_no, SNR = 5, dd = "D0", nknot = 3, seeds = 1:30){
  
  dir <- paste("Results/Results_1_x_mattern_g_", g_func_no, "_SNR_", SNR,"_nknot", nknot, sep = "")
  file_list <- list.files(dir)
  
  res <- NULL
  for(seed in seeds){
      tmp_val <- NULL
      tmp_test <- NULL
    for(d in 1:4){
        file_name <-  paste("g_func_no_", g_func_no, "_case_id_1_d_", d, "_dd_", dd,"_seed_", seed,".RData",sep = "")
        if(file_name %in% file_list){
          load(paste(dir,"/",file_name,sep = ""))
        if(dd == "D0"){
          tmp_val <- rbind(tmp_val, c(mse(dat2save$res_vals[[1]]),mse(dat2save$res_vals[[2]]), mse(dat2save$res_vals[[3]])))
        }else{
          tmp_val <- rbind(tmp_val, c(robust_mse(dat2save$res_vals[[1]]),robust_mse(dat2save$res_vals[[2]]), robust_mse(dat2save$res_vals[[3]])))
        }
        tmp_test <- rbind(tmp_test,  dat2save$err_test)
        }
    }
        
    idx <- apply(tmp_val, 2, which.min)
    tmp <- rep(NA, 3)
    for(j in 1:3){
      tmp[j] <-tmp_test[idx[j],j]
    }
    res <- rbind(res, tmp)
  }
  
  return(res)
}


y_lim <- list()
y_lim[[1]] <- c(0,0.7)
y_lim[[7]] <- c(0.2,1.5)
y_lim[[8]] <- c(0.2,1.5)
y_lim[[9]] <- c(0.4,2)

pdf("sim.pdf", width = 12, height = 3)
for(g_func_no  in c(1,7,8,9)){
  
  par(mfrow = c(1,3))
  for(dd in c("D0","D1","D2")){
    res_0 <- summarize_0(g_func_no, SNR = 5, dd = dd, nknot = 3, seeds = 1:30)
    res_1 <- summarize_1(g_func_no, SNR = 5, dd = dd, nknot = 3, seeds = 1:30)
    switch(as.character(g_func_no), 
           "1" = {
             r_no <- 1
            },
           "7" = {
             r_no <- 2
            },
           "8" = {
             r_no <- 3
            },
           "9" = {
             r_no <- 4
             } 
          )
    main_title = paste("(r", r_no, ", ", dd,")", sep = "")
    boxplot(res_1[,1], res_1[,2], res_1[,3], res_0, ylim = c(y_lim[[g_func_no]][1],y_lim[[g_func_no]][2]), main = main_title, names=c("TFBoost(L2)","TFBoost(LAD)","TFBoost(RR)","RobustFPLM"))
  }
}
dev.off()

