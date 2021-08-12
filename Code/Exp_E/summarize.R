## summary code
robust_mse <- function(x){
  return( median(x)^2 + mscale(x)^2)
}

summarize0 <- function(g_func_no, type, SNR = 5, nknot = 3, seeds = 1:20, summarize_type = "pred"){
  
  dir <- paste("Results_cedar11/Results_0_type_", type, "_g_", g_func_no, "_SNR_", SNR,"_nknot", nknot, sep = "")
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

summarize1 <- function(g_func_no, type, SNR = 5, nknot = 3, seeds =c(1:3, 5:20), summarize_type = "pred"){
  
  dir <- paste("Results_cedar11/Results_1_type_", type, "_g_", g_func_no, "_SNR_", SNR,"_nknot", nknot, sep = "")
  file_list <- list.files(dir)
  
  res <- NULL
  res_early <- NULL
  res_d <- NULL
  for(seed in seeds){
    tmp_val <- NULL
    tmp_test <- NULL
    tmp_early <- NULL
    for(d in 1:4){
      file_name <-  paste("g_func_no_", g_func_no, "_case_id_1_d_", d, "_type_", type,"_seed_", seed,".RData",sep = "")
      if(file_name %in% file_list){
        load(paste(dir,"/",file_name,sep = ""))
        tmp_val <- rbind(tmp_val, c(robust_mse(dat2save$RTFBoost$pred_val[[1]]),robust_mse(dat2save$RTFBoost$pred_val[[2]]),
                                    robust_mse(dat2save$RTFBoost$pred_val[[3]]),  robust_mse(dat2save$RTFBoost$pred_val[[4]]),
                                    robust_mse(dat2save$RTFBoost$pred_val[[5]])))
        tmp_test <- rbind(tmp_test,  c(dat2save$RTFBoost$err_test))
        tmp_early <- rbind(tmp_early, dat2save$RTFBoost$early_stop)
      }
    }
    
    idx <- apply(tmp_val, 2, which.min)
    tmp <- rep(NA, 5)
    for(j in 1:5){
      tmp[j] <-tmp_test[idx[j],j]
    }
    res_d <- rbind(res_d, idx)
    res <- rbind(res, tmp)
    res_early  <- rbind(res_early,  c( tmp_early[idx[1],1], tmp_early[idx[2],2], tmp_early[idx[3],3], tmp_early[idx[4],4], tmp_early[idx[5],5]))
  }
  
  colnames(res) <- c("RTFBoost(L2)", "RTFBoost(LAD)", "RTFBoost(LAD-M)", "RTFBoost(RR.2)", "RTFBoost(RR.5)")
  dat2return <- list(res = res, res_early = res_early,   res_d  =   res_d )
  
  if(summarize_type == "pred"){
      return(res)
  }else{
    return(dat2return)
  }
}



make_figure <- function(g_func_no, method_type = "RTFBoost"){
  
  y_lim <- list(c(0.1,0.4),  c(0.15,0.6),  c(0.25,0.8), 
                c(0.4,1.1),  c(0.45,1.5))
  tmp <- NULL
  for(type in c("C0","C1","C2","C3","C4","C5")){
    res1 <- summarize1(g_func_no = g_func_no, type = type, summarize_type = "all")
    res1 <- data.frame(res1$res, rep(type, nrow(res1$res)))
    
    colnames(res1) <- c("TFBoost(L2)", "TFBoost(LAD)", "TFBoost(LAD-M)", "TFBoost(RR.2)", "TFBoost(RR.5)","type")
    tmp  <- rbind(tmp, res1)
  }
  dat2plot <- melt(tmp, id.vars = "type", variable.name = "method")
  
  tmp <- NULL
  for(type in c("C0","C1","C2","C3","C4","C5")){
    res0 <- summarize0(g_func_no = g_func_no, type = type, summarize_type = "all")
    res0 <- data.frame(res0, rep(type, nrow(res0)))
    colnames(res0) <- c( "FPPR","FGAM","RFPLM","MFLM", "RFSIR", "type")
    tmp  <- rbind(tmp, res0)
  }
  
  dat2plot <- rbind(dat2plot, melt(tmp, id.vars = "type", variable.name = "method"))
  dat2plot$value <- as.numeric(dat2plot$value)
  p <- ggplot(dat2plot, aes(fill=method, y=value)) + 
    geom_boxplot() + ylab("MSPE") + xlab(" ")  + ggtitle(paste("r", g_func_no, sep = ""))+
    facet_wrap(vars(type), ncol = 3) + theme(axis.ticks.x=element_blank(),axis.text.x=element_blank(), plot.title = element_text(hjust = 0.5), legend.position="bottom") +  scale_fill_brewer(palette="Spectral")
  
  return(p +coord_cartesian(ylim= c(y_lim[[g_func_no]][1],y_lim[[g_func_no]][2])))
  
}

library(reshape2)
library(ggplot2)
pp <- list()
for(g_func_no in 1:5){
  pp[[g_func_no]] <- make_figure(g_func_no = g_func_no, method_type = "RTFBoost")
}

pdf("make_figures_11.pdf", width = 12, height = 8)
for(i in 1:5){
  print(pp[[i]])
}
dev.off()



  
