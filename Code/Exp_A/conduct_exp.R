rm(list = ls())
library(rpart)

source("Code/RTFBoost.R")
source("Code/Tree.init.R")
source("Code/Tree.R")
source("Code/user_rpart.R")
source("Code/utils.R")
source("Code/Exp_A/fun_exp.R")
source("Code/Exp_A/dat_gen.R")


args <- commandArgs(trailingOnly = TRUE) 
x_type <- c("mattern") 
case_id <- 1
nknot <-3
g_func_nos <- c(1,7,8,9)
ds <- c(1,2,3,4)
SNRs <- c(5)
shrinkage <- 0.05
niter <- 1000
dds <- c("D0","D1","D2")

if(case_id == 1){
  niter <- 1000  # doesn't matter
  ncol_num <- 2
}

if(length(args)== 0){
  g_func_nos <- c(1,7)
  SNR <- 5
  ds <- c(1,3)
  test <- TRUE
  RNGkind(sample.kind = "default")
  exp_id <- 1  # row index in the conduct_sheet 
  case_id <- 1 # shrinkage id  0 means only run FPPR and FAM
  niter <- 2
  if(case_id == 1){
    n_init <- 500
    niter <- 1000  
    ncol_num <- 2
  }
  if(case_id == 0){
    n_init <- 500 # doesn't matter
    niter <- 1000  # doesn't matter
    ncol_num <- 2
  }
}else{
  exp_id <- as.numeric(args[1]) 
  case_id <- as.numeric(args[2]) 
  if(case_id == 1){
    n_init <- 500
    niter <- 1000  
    ncol_num <- 2
  }
  if(case_id == 0){
    n_init <- 500 # doesn't matter
    niter <- 1000  # doesn't matter
    ncol_num <- 2
  }
}


seeds <- 1:10

if(case_id == 0){
  ds <- 1
}

all_exps <- expand.grid(g_func_nos, ds, dds, seeds, SNRs)
colnames(all_exps) <- c("g_func_no", "d","dd","seed","SNR" )

if(nrow(all_exps)%%ncol_num!=0){
  
  if(nrow(all_exps)%%ncol_num!=0){
    conduct_sheet <- matrix(c(1:nrow(all_exps), rep(NA, ncol_num - nrow(all_exps)%%ncol_num)), ncol = ncol_num)
  }else{
    conduct_sheet <- matrix(c(1:nrow(all_exps)), ncol = ncol_num)
  }
}else{
  conduct_sheet <- matrix(1:nrow(all_exps), ncol = ncol_num)
}

print(dim(conduct_sheet))

conduct.exp <- function(exp_id = 1, conduct_sheet){
  
  for(i in 1:ncol(conduct_sheet)) {
    
    print(c(i,"th trial"))
    if(!is.na(conduct_sheet[exp_id,i])){
      
      g_func_no <- all_exps[conduct_sheet[exp_id, i],]$g_func_no
      seed <- all_exps[conduct_sheet[exp_id, i],]$seed
      d <-  all_exps[conduct_sheet[exp_id, i],]$d
      dd <-  all_exps[conduct_sheet[exp_id, i],]$dd
      SNR <- all_exps[conduct_sheet[exp_id, i],]$SNR
      
      if(x_type == "ferraty" && g_func_no == 6){
        precision <- 5
      }else{
        if(x_type == "mattern" && (g_func_no %in% c(4,6))){
          precision <- 7
        }else{
          precision <- 6
        }
      }
      
      if(test){
        dir_name <- paste("Results_tmp/Results_", case_id, "_x_", x_type, "_g_",g_func_no,  "_SNR_", SNR, "_nknot", nknot,  sep = "")
      }else{
        dir_name <- paste("Results/Results_", case_id,"_x_", x_type, "_g_",g_func_no, "_SNR_", SNR, "_nknot", nknot,  sep = "")
      }
      
      if(dir.exists(dir_name) == FALSE){
        dir.create(dir_name)
      }
      
      start_time <- Sys.time()
      control.tree.list <- set.control(d, precision,  shrinkage, nknot, n_init, niter)
      dat2save <- tryCatch(do.exp(seed = seed, g_func_no = g_func_no, SNR = SNR, x_type = x_type, dd = dd, 
                                  control.tree.list = control.tree.list, case_id = case_id))
      end_time <- Sys.time()
      
      dat2save$time <- end_time - start_time
      print(dat2save$time)
      save(file = paste(dir_name, "/g_func_no_", g_func_no, "_case_id_",case_id, "_d_", d,  "_seed_", seed, ".RData", sep = ""), dat2save)
      
    }
  }
}


conduct.exp(exp_id = exp_id, conduct_sheet)
